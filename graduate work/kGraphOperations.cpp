#include "kGraphOperations.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>
#include <queue>
#include <map>
#include <optional>


std::ofstream output("C:/Users/User/source/repos/graduate work/graduate work/output.txt");
std::ofstream output1("C:/Users/User/source/repos/graduate work/graduate work/output1.txt", std::ios_base::app);
int NumberOfRec = 0;
std::vector<double> sumReliab;
double globsumReliab = 0;
std::vector<std::vector<double>> BlockReliab;


std::optional<int> kGraph::FindLastUnreliableEdge() const {
    for (int i = static_cast<int>(this->PArray.size()) - 1; i >= 0; --i)
        if (this->PArray[i] < 1.0)
            return i;
    return std::nullopt;
}

std::pair<int, int> kGraph::FindReverseEdgeIndices(int edgeIndex) const {
	int from = 1;
	while (from < this->KAO.size() && this->KAO[from] <= edgeIndex)
		++from;
	int to = this->FO[edgeIndex];
	int reverseStart = this->KAO[to - 1];
	int reverseEnd = this->KAO[to];
	for (int k = reverseStart; k < reverseEnd; ++k)
		if (this->FO[k] == from)
			return { from, k };
	return { from, -1 }; // Если обратное ребро не найдено
}

void Factoring(kGraph G, const int variant, const int d, double Reliab) {
	//��������� ����� � sumReliab[0]
	//���������, variant=0 - ����� ��������, variant=1 - ����� ������������� �����
	NumberOfRec++;

	if ((!variant && !G.KComponent()) || (variant && !G.CheckDistanceFloyd(d)))
		return;// add 0 to sum

	auto optIndex = G.FindLastUnreliableEdge();
	if (!optIndex) {
		sumReliab[0] += Reliab;
		return;
	}

	int edgeIdx = *optIndex;
	double p = G.PArray[edgeIdx];
	G.PArray[edgeIdx] = 1.0;

	auto [from, reverseIdx] = G.FindReverseEdgeIndices(edgeIdx);
	if (reverseIdx != -1) G.PArray[reverseIdx] = 1.0;

	Factoring(G, 1, d, Reliab * p);

	if (reverseIdx != -1) G.PArray[reverseIdx] = p;
	G.PArray[edgeIdx] = p;

	kGraph T = G.DeleteEdgeK(G.FO[edgeIdx], from);
	Factoring(T, 0, d, Reliab * (1 - p));
}

void Factoring2Vert(kGraph G, int x, int y, int variant, int d, double Reliab) {
	//результат будет в globsumReliab
	NumberOfRec++;

	if (!variant && G.DistanceDijkstra(x, y) > d)
		return;

	auto optIndex = G.FindLastUnreliableEdge();
	if (!optIndex) {
		globsumReliab += Reliab;
		return;
	}

	int edgeIdx = *optIndex;
	double p = G.PArray[edgeIdx];
	G.PArray[edgeIdx] = 1.0;

	auto [from, reverseIdx] = G.FindReverseEdgeIndices(edgeIdx);
	if (reverseIdx != -1) G.PArray[reverseIdx] = 1.0;

	Factoring2Vert(G, x, y, 1, d, Reliab * p);

	if (reverseIdx != -1) G.PArray[reverseIdx] = p;
	G.PArray[edgeIdx] = p;

	kGraph T = G.DeleteEdgeK(G.FO[edgeIdx], from);
	Factoring2Vert(T, x, y, 0, d, Reliab * (1 - p));
}

void Factoring2VertM(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound) {
	//��������� ������ ����������� � ��������������� ������������ �� �������
	//���������, variant=0 - ����� ��������, variant=1 - ����� ������������� �����
	NumberOfRec++;

	int dist = G.DistanceDijkstra(x, y);
	if (!variant && dist > d) {
		if (dist <= UpperBound) {
			NumberOfRec--;
			Factoring2VertM(G, x, y, 0, d + 1, Reliab, LowerBound, UpperBound);
		}
		return;// add 0 to sum if distance from x to y >d
	}

	auto optIndex = G.FindLastUnreliableEdge();
	if (!optIndex) {
		sumReliab[d - LowerBound] += Reliab;
		return;
	}
	
	int edgeIdx = *optIndex;
	double p = G.PArray[edgeIdx];
	G.PArray[edgeIdx] = 1.0;

	int j;
	for (j = 1; j < G.KAO.size(); j++)
		if (G.KAO[j] > edgeIdx)
			break;
	int k;
	for (k = G.KAO[G.FO[edgeIdx] - 1]; k < G.KAO[G.FO[edgeIdx]]; k++)
		if (G.FO[k] == j)
			break;
	
	G.PArray[k] = 1;
	Factoring2VertM(G, x, y, 1, d, Reliab * p, LowerBound, UpperBound);
	G.PArray[edgeIdx] = p;
	G.PArray[k] = p;
	kGraph T = G.DeleteEdgeK(G.FO[edgeIdx], j);
	Factoring2VertM(T, x, y, 0, d, Reliab * (1 - p), LowerBound, UpperBound);
}

void Factoring2VertMParallel(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound) {
	//��������� ������ ����������� � ��������������� ������������ �� �������
	//���������, variant=0 - ����� ��������, variant=1 - ����� ������������� �����
	NumberOfRec++;
	int dist = G.DistanceDijkstra(x, y);
	if (!variant && dist > d) {
		if (dist <= UpperBound) {
			NumberOfRec--;
			Factoring2VertMParallel(G, x, y, 0, d + 1, Reliab, LowerBound, UpperBound);
		}
		return;// add 0 to sum if distance from x to y >d
	}
	int i = G.PArray.size();
	for (int j = G.PArray.size() - 1; j >= 0; j--)
		if (G.PArray[j] < 1) {
			i = j;
			break;
		}
	if (i == G.FO.size()) {
		BlockReliab[omp_get_thread_num()][d - LowerBound] += Reliab;
		return; // add Reliab to sum
	}
	double p = G.PArray[i];
	G.PArray[i] = 1;
	int j;
	for (j = 1; j < G.KAO.size(); j++)
		if (G.KAO[j] > i)
			break;
	int k;
	for (k = G.KAO[G.FO[i] - 1]; k < G.KAO[G.FO[i]]; k++)
		if (G.FO[k] == j)
			break;
	G.PArray[k] = 1;
	Factoring2VertMParallel(G, x, y, 1, d, Reliab * p, LowerBound, UpperBound);
	G.PArray[i] = p;
	G.PArray[k] = p;
	kGraph T = G.DeleteEdgeK(G.FO[i], j);
	Factoring2VertMParallel(T, x, y, 0, d, Reliab * (1 - p), LowerBound, UpperBound);
}

void GraphMerging(int k) {
	std::ifstream fin("C:/Users/User/source/repos/graduate work/graduate work/GraphsToMerge.txt");
	if (!fin) {
		std::cout << "Error!\n";
		throw std::runtime_error("OPEN_ERROR");
	}
	kGraph G1 = kGraphFileInput(fin);
	kGraph G2 = kGraphFileInput(fin);
	fin.close();
	kGraph G = UnionGraphs(G1, G2, k);
	for (int i = 1; i < G1.KAO.size() / 2; i++)
		G.ChangeVertex(i, G1.KAO.size() - 1);
	std::ofstream fout("C:/Users/User/source/repos/graduate work/graduate work/MergedGraphs.txt", std::ofstream::trunc);
	if (!fout) {
		std::cout << "Error!\n";
		throw std::runtime_error("OPEN_ERROR");
	}
	G.kGraphFileOutput(fout);
}

kGraph UnionGraphs(kGraph G1, kGraph G2, int k) {
	//{��������� ����, ���������� ������������ G1 � G2 �� �� ������ k ��������}
	//{NN - ������ � �������� �-� G2 � G; ������ �-� G1 � G ��������� � �� �������� � G1}
	int N1 = G1.KAO.size() - 1, N2 = G2.KAO.size() - 1;
	std::vector<int> NN(N2 + 1);
	for (int i = 1; i < k + 1; i++)
		NN[i] = i;
	for (int i = k + 1; i < N2 + 1; i++)
		NN[i] = N1 + i - k;

	std::vector<int> KAO(N1 + N2 - k + 1), FO(0);
	int l = 0;
	bool boolean;
	for (int i = 1; i < k + 1; i++) {
		KAO[i] = KAO[i - 1];
		for (int j = G1.KAO[i - 1]; j < G1.KAO[i]; j++) {
			KAO[i]++;
			FO.resize(++l);
			FO[l - 1] = G1.FO[j];
		}
		for (int j = G2.KAO[i - 1]; j < G2.KAO[i]; j++) {
			boolean = true;
			for (int s = KAO[i - 1]; s < KAO[i]; s++)
				if (FO[s] == NN[G2.FO[j]])
					boolean = false;
			if (boolean) {
				KAO[i]++;
				FO.resize(++l);
				FO[l - 1] = NN[G2.FO[j]];
			}
		}
	}
	for (int i = k + 1; i < N1 + 1; i++) {
		KAO[i] = KAO[i - 1];
		for (int j = G1.KAO[i - 1]; j < G1.KAO[i]; j++) {
			KAO[i]++;
			FO.resize(++l);
			FO[l - 1] = G1.FO[j];
		}
	}
	for (int i = N1 + 1; i < N1 + N2 - k + 1; i++) {
		KAO[i] = KAO[i - 1];
		for (int j = G2.KAO[i - N1 + k - 1]; j < G2.KAO[i - N1 + k]; j++) {
			KAO[i]++;
			FO.resize(++l);
			FO[l - 1] = NN[G2.FO[j]];
		}
	}
	std::vector<double> PArray(FO.size());
	for (int i = 0; i < PArray.size(); i++)
		PArray[i] = 0.9; // �����, ���� ����� �����������, ��������������� ������
	std::vector<int> Targets(KAO.size());
	//Targets ���� ���� ��������� ��������������� �������
	kGraph Result(KAO, FO, PArray, Targets);

	return Result;
}

void kGraph::ChangeVertex(int u, int v) {
	//������ � ����� ������� u � v ������� (����������������)
	if (u > this->KAO.size() - 1 || v > this->KAO.size() - 1 || u == v)
		return;
	if (u > v)
		std::swap(u, v);
	std::vector<int> KAO(this->KAO.size()), FO(this->FO.size()), Targets(this->Targets.size());
	std::vector<double> PArray(this->PArray.size());
	for (int i = 0; i < u; i++) {
		KAO[i] = this->KAO[i];
		Targets[i] = this->Targets[i];
	}
	KAO[u] = this->KAO[u - 1] + this->KAO[v] - this->KAO[v - 1];
	Targets[u] = this->Targets[v];
	for (int i = u + 1; i < v; i++) {
		KAO[i] = KAO[i - 1] + this->KAO[i] - this->KAO[i - 1];
		Targets[i] = this->Targets[i];
	}
	KAO[v] = KAO[v - 1] + this->KAO[u] - this->KAO[u - 1];
	Targets[v] = this->Targets[u];
	for (int i = v + 1; i < this->KAO.size(); i++) {
		KAO[i] = this->KAO[i];
		Targets[i] = this->Targets[i];
	}
		
	int j = 0;
	for (int i = this->KAO[v - 1]; i < this->KAO[v]; i++) {
		FO[KAO[u - 1] + j] = this->FO[i];
		PArray[KAO[u - 1] + j] = this->PArray[i];
		j++;
	}
	j = 0;
	for (int i = this->KAO[u - 1]; i < this->KAO[u]; i++) {
		FO[KAO[v - 1] + j] = this->FO[i];
		PArray[KAO[v - 1] + j] = this->PArray[i];
		j++;
	}
	for (int ver = u + 1; ver < v; ver++) {
		j = 0;
		for (int i = this->KAO[ver - 1]; i < this->KAO[ver]; i++) {
			FO[KAO[ver - 1] + j] = this->FO[i];
			PArray[KAO[ver - 1] + j] = this->PArray[i];
			j++;
		}
	}
	for (int ver = 1; ver < u; ver++)
		for (int i = this->KAO[ver - 1]; i < this->KAO[ver]; i++) {
			FO[i] = this->FO[i];
			PArray[i] = this->PArray[i];
		}
	for (int ver = v + 1; ver < this->KAO.size(); ver++)
		for (int i = this->KAO[ver - 1]; i < this->KAO[ver]; i++) {
			FO[i] = this->FO[i];
			PArray[i] = this->PArray[i];
		}
	for (int i = 0; i < this->FO.size(); i++)
		if (FO[i] == u)
			FO[i] = v;
		else if (FO[i] == v)
			FO[i] = u;
	kGraph Result(KAO, FO, PArray, Targets);
	*this = Result;
}

// НОВАЯ вспомогательная функция для поиска пути в графе блоков (BFS)
std::vector<int> find_path_in_block_graph(
	int start_block_node_idx, // Индекс стартового узла в block_nodes
	int end_block_node_idx,   // Индекс целевого узла в block_nodes
	int num_block_nodes,
	const std::vector<BlockGraphEdge>& block_edges,
	std::vector<int>& out_path_articulation_points_original_ids // ID точек сочленения на пути
) {
	std::vector<std::vector<std::pair<int, int>>> adj(num_block_nodes); // смежный узел, ID точки сочленения
	for (const auto& edge : block_edges) {
		adj[edge.from_block_node_idx].push_back({ edge.to_block_node_idx, edge.articulation_point_original_id });
		adj[edge.to_block_node_idx].push_back({ edge.from_block_node_idx, edge.articulation_point_original_id }); // Для неориентированного графа блоков
	}

	std::vector<int> path_nodes;
	out_path_articulation_points_original_ids.clear();

	if (start_block_node_idx == end_block_node_idx) {
		path_nodes.push_back(start_block_node_idx);
		return path_nodes;
	}

	std::queue<std::vector<std::pair<int, int>>> q; // {node_idx, articulation_point_id_to_reach_this_node}
	q.push({ {start_block_node_idx, 0} }); // Начальная точка сочленения 0 (не используется)

	std::vector<std::vector<std::pair<int, int>>> paths_to_node(num_block_nodes); // Хранит путь к каждому узлу
	paths_to_node[start_block_node_idx] = { {start_block_node_idx, 0} };
	std::vector<bool> visited(num_block_nodes, false);
	visited[start_block_node_idx] = true;

	while (!q.empty()) {
		std::vector<std::pair<int, int>> current_path_info = q.front();
		q.pop();
		int u_node_idx = current_path_info.back().first;

		if (u_node_idx == end_block_node_idx) {
			for (const auto& p_info : current_path_info) {
				path_nodes.push_back(p_info.first);
				if (p_info.second != 0) { // Пропускаем фиктивную точку сочленения для старта
					out_path_articulation_points_original_ids.push_back(p_info.second);
				}
			}
			return path_nodes;
		}

		for (auto& edge_info : adj[u_node_idx]) {
			int v_node_idx = edge_info.first;
			int ap_id = edge_info.second;
			if (!visited[v_node_idx]) {
				visited[v_node_idx] = true;
				std::vector<std::pair<int, int>> new_path_info = current_path_info;
				new_path_info.push_back({ v_node_idx, ap_id });
				q.push(new_path_info);
				paths_to_node[v_node_idx] = new_path_info; // Сохраняем путь
			}
		}
	}
	return {}; // Путь не найден
}

void kGraph::ReliabilityDiamConstr(int d) {
	//{������ ����������� ��������� � ������������ �� ������� d ������� ���������.
	// ������������ ��������� ����-�� ��������� � �������� �-��, ���������� ������,
	// ���������� �-� �������� �� ����������}
	Nconst = this->KAO.size() * this->KAO.size();
	//double t = Factoring(G, 0, d);
	//return t;
	Factoring(*this, 0, d, 1);
}

GraphMethodResult kGraph::ReliabilityDiamConstr2Vert(int x, int y, int d) {
    Nconst = this->KAO.size() * this->KAO.size();
    clock_t start_time = clock();
    globsumReliab = 0; // Сбрасываем глобальную переменную
    NumberOfRec = 0;
    Factoring2Vert(*this, x, y, 0, d, 1.0);
    double time_sec = (clock() - start_time) / 1000.0000;
    return GraphMethodResult(globsumReliab, NumberOfRec, time_sec);
}

GraphMethodResult kGraph::ReliabilityDiamConstr2VertDecomposeSimpleFacto(int x, int y, int UpperBound)
{
	Nconst = this->KAO.size() * this->KAO.size();
	clock_t start_time = clock();
	NumberOfRec = 0;
	sumReliab.clear();
	BlockReliab.clear();

	std::vector<int> spisok = this->DecomposeOnBlocksK();
	int BlockNum = spisok.back();

	if (BlockNum == 1) {
		int BlockDiam = this->DistanceDijkstra(x, y);
		
		if (BlockDiam > UpperBound) {
			double time_sec = (clock() - start_time) / 1000.0000;
			std::cout << "Error: BlockDiam > UpperBound" << std::endl;
			return GraphMethodResult(0.0, NumberOfRec, time_sec);
		}

		Factoring2Vert(*this, x, y, 0, UpperBound, 1);
		double reliability = globsumReliab;
		double time_sec = (clock() - start_time) / 1000.0000;
		return GraphMethodResult(reliability, NumberOfRec, time_sec);
	}

	std::vector<int> BlockDiam(BlockNum);
	int diamsum = 0;

	for (int i = 0; i < BlockNum; i++) {
		kGraph Restored = this->RestoreBlockK(i + 1, spisok);
		int local_x = 0, local_y = 0;
		for (int j = 1; j < Restored.Targets.size(); j++) {
			if (Restored.Targets[j]) {
				if (local_x == 0) local_x = j;
				else if (local_y == 0) local_y = j;
			}
		}
		BlockDiam[i] = Restored.DistanceDijkstra(local_x, local_y);
		diamsum += BlockDiam[i];
	}

	if (diamsum > UpperBound) {
		double time_sec = (clock() - start_time) / 1000.0000;
		std::cout << "Error: diamsum > UpperBound" << std::endl;
		return GraphMethodResult(0.0, NumberOfRec, time_sec);
	}

	int gap = UpperBound - diamsum;
	BlockReliab.resize(BlockNum, std::vector<double>(gap + 1, 0.0));

	for (int i = 0; i < BlockNum; i++) {
		kGraph Restored = this->RestoreBlockK(i + 1, spisok);
		int local_x = 0, local_y = 0;
		for (int j = 1; j < Restored.Targets.size(); j++) {
			if (Restored.Targets[j]) {
				if (local_x == 0) local_x = j;
				else if (local_y == 0) local_y = j;
			}
		}
		for (int j = 0; j < gap + 1; j++) {
			globsumReliab = 0;
			Factoring2Vert(Restored, local_x, local_y, 0, BlockDiam[i] + j, 1);
			BlockReliab[i][j] = globsumReliab;
		}
	}

	std::vector<double> curBlockRel(gap + 1, 0.0);
	//for (int j = 1; j < BlockReliab[BlockNum - 1].size(); j++) {
	//	BlockReliab[BlockNum - 1][j] += BlockReliab[BlockNum - 1][j - 1];
	//}
	//получили полные надежности вместо точных, которые не применимы в формуле
	for (int i = 0; i < BlockNum - 1; i++) {
		for (int j = gap; j > 0; j--) {
			BlockReliab[i][j] -= BlockReliab[i][j - 1];
		}
	}

	for (int i = BlockNum - 1; i > 0; i--) {
		for (int diam = 0; diam < gap + 1; diam++) {
			for (int j = 0; j < diam + 1; j++) {
				curBlockRel[diam] += BlockReliab[i - 1][j] * BlockReliab[i][diam - j];
			}
		}
		for (int j = 0; j < gap + 1; j++) {
			BlockReliab[i - 1][j] = curBlockRel[j];
			curBlockRel[j] = 0.0;
		}
	}

	double time_sec = (clock() - start_time) / 1000.0000;
	
	double reliability = BlockReliab[0].back();
	
	return GraphMethodResult(reliability, NumberOfRec, time_sec);
}

void kGraph::ReliabilityDiamConstr2VertM(int x, int y, const int& UpperBound){
//    {Расчет в-ти св-ти двух вершин x,y с огр-м на диаметр d методом модифицированной факторизации.
//     Используется разделение ветвей, встроенная ф-я проверки на расстояния;
//     не используется выделение компонентв с двумя в-ми (так быстрее)}
	Nconst = this->KAO.size() * this->KAO.size();
	clock_t start_time = clock();
	int TerminalShortestPathLen = this->DistanceDijkstra(x, y);
	if (TerminalShortestPathLen > UpperBound) {
		output << "Error 'TerminalShortestPath > UpperBound'";
		return;
	}
	int gap = UpperBound - TerminalShortestPathLen;
	
	sumReliab.resize(gap + 1);
	Factoring2VertM(*this, x, y, 0, UpperBound, 1.0, TerminalShortestPathLen, UpperBound);
	output << std::setprecision(16) << sumReliab[gap] << std::endl;
	output << "Recursions: " << NumberOfRec << std::endl;
	output1 << (clock() - start_time) / 1000.0000 << std::endl;
}

GraphMethodResult kGraph::ReliabilityDiamConstr2VertMDecompose(int x, int y, const int& UpperBound) {
	Nconst = this->KAO.size() * this->KAO.size();
	clock_t start_time = clock();
	NumberOfRec = 0;
	sumReliab.clear();
	BlockReliab.clear();

	std::vector<int> spisok = this->DecomposeOnBlocksK();
	int BlockNum = spisok.back();

	if (BlockNum == 1) {
		int BlockDiam = this->DistanceDijkstra(x, y);
		sumReliab.resize(UpperBound - BlockDiam + 1);
		Factoring2VertM(*this, x, y, 0, UpperBound, 1, BlockDiam, UpperBound);
		for (int i = 1; i < sumReliab.size(); i++) {
			sumReliab[i] += sumReliab[i - 1];
		}
		double time_sec = (clock() - start_time) / 1000.0000;
		double reliability = sumReliab.back();
		return GraphMethodResult(reliability, NumberOfRec, time_sec);
	}

	std::vector<int> BlockDiam(BlockNum);
	int diamsum = 0;

	for (int i = 0; i < BlockNum; i++) {
		kGraph Restored = this->RestoreBlockK(i + 1, spisok);
		int local_x = 0, local_y = 0;
		for (int j = 1; j < Restored.Targets.size(); j++) {
			if (Restored.Targets[j]) {
				if (local_x == 0) local_x = j;
				else if (local_y == 0) local_y = j;
			}
		}
		BlockDiam[i] = Restored.DistanceDijkstra(local_x, local_y);
		diamsum += BlockDiam[i];
	}

	if (diamsum > UpperBound) {
		std::cout << "Diameter sum overflow (diamsum = " + std::to_string(diamsum) + ", UpperBound = " + std::to_string(UpperBound) + ")\n";
		double time_sec = (clock() - start_time) / 1000.0000;
		return GraphMethodResult(0.0, NumberOfRec, time_sec);
	}

	int gap = UpperBound - diamsum;
	BlockReliab.resize(BlockNum);

	for (int i = 0; i < BlockNum; i++) {
		kGraph Restored = this->RestoreBlockK(i + 1, spisok);
		int local_x = 0, local_y = 0;
		for (int j = 1; j < Restored.Targets.size(); j++) {
			if (Restored.Targets[j]) {
				if (local_x == 0) local_x = j;
				else if (local_y == 0) local_y = j;
			}
		}
		sumReliab.clear();
		sumReliab.resize(gap + 1, 0.0);
		Factoring2VertM(Restored, local_x, local_y, 0, BlockDiam[i], 1, BlockDiam[i], BlockDiam[i] + gap);
		BlockReliab[i] = sumReliab;
	}

	std::vector<double> curBlockRel(gap + 1, 0.0);
	for (int j = 1; j < BlockReliab[BlockNum - 1].size(); j++) {
		BlockReliab[BlockNum - 1][j] += BlockReliab[BlockNum - 1][j - 1];
	}

	for (int i = BlockNum - 1; i > 0; i--) {
		for (int diam = 0; diam < gap + 1; diam++) {
			for (int j = 0; j < diam + 1; j++) {
				curBlockRel[diam] += BlockReliab[i - 1][j] * BlockReliab[i][diam - j];
			}
		}
		for (int j = 0; j < gap + 1; j++) {
			BlockReliab[i - 1][j] = curBlockRel[j];
			curBlockRel[j] = 0.0;
		}
	}

	double time_sec = (clock() - start_time) / 1000.0000;
	
	double reliability = BlockReliab[0].back();

	return GraphMethodResult(reliability, NumberOfRec, time_sec);
}

void kGraph::ReliabilityDiamConstr2VertMDecomposeParallel(int x, int y, const int& UpperBound) {
	//��� ������ Release � VS ����������������� �� ��������
	Nconst = this->KAO.size() * this->KAO.size();
	clock_t start_time = clock();
	std::vector<int> spisok = this->DecomposeOnBlocksK();
	int BlockNum = spisok.back();
	//���� 1 ����
	if (BlockNum == 1) {
		int BlockDiam = this->DistanceDijkstra(x, y);
		sumReliab.resize(UpperBound - BlockDiam + 1);
		Factoring2VertM(*this, x, y, 0, UpperBound, 1, BlockDiam, UpperBound);
		for (int i = 1; i < sumReliab.size(); i++)
			sumReliab[i] += sumReliab[i - 1];
		output << std::setprecision(16) << sumReliab.back() << std::endl;
		output << "Recursions: " << NumberOfRec << std::endl;
		output << "Decompose Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
		return;
	}

	std::vector<int> BlockDiam(BlockNum);
	int diamsum = 0;  //����� ���� ���������

	for (int i = 0; i < BlockNum; i++) {
		kGraph Restored = this->RestoreBlockK(i + 1, spisok);
		x = 0; y = 0;
		for (int j = 1; j < Restored.Targets.size(); j++)
			if (Restored.Targets[j]) {
				x = j;
				break;
			}
		for (int j = x + 1; j < Restored.Targets.size(); j++)
			if (Restored.Targets[j]) {
				y = j;
				break;
			}
		BlockDiam[i] = Restored.DistanceDijkstra(x, y);
		diamsum += BlockDiam[i];
	}

	output << "Diameter calculations(ms): " << (clock() - start_time) << std::endl;

	//diamsum - ����������� ���������� ����� �������� ���������
	if (diamsum > UpperBound) {
		output << "Diameter sum overflow (diamsum = " + std::to_string(diamsum) + ", UpperBound = " + std::to_string(UpperBound) + ")\n";
		return;
	}
	int gap = UpperBound - diamsum;
	start_time = clock();

	//omp_set_num_threads(BlockNum);
	BlockReliab.resize(BlockNum,std::vector<double>(gap + 1));
	
	#pragma omp parallel for num_threads(BlockNum)
			for (int i = 0; i < BlockNum; i++) {
				//std::cout << omp_get_num_threads() << std::endl;
					kGraph Restored = this->RestoreBlockK(i + 1, spisok);
					for (int j = 1; j < Restored.Targets.size(); j++)
						if (Restored.Targets[j]) {
							x = j;
							break;
						}
					for (int j = x + 1; j < Restored.Targets.size(); j++)
						if (Restored.Targets[j]) {
							y = j;
							break;
						}
					Factoring2VertMParallel(Restored, x, y, 0, BlockDiam[i], 1, BlockDiam[i], gap + BlockDiam[i]);
					//std::cout << omp_get_thread_num() << "finished" << std::endl;
			}

	clock_t reliabcalc = clock();
	//����� n ������� n - 1 ����� ����������
	std::vector<double> curBlockRel(gap + 1);
	for (int j = 1; j < BlockReliab[BlockNum - 1].size(); j++)
		BlockReliab[BlockNum - 1][j] += BlockReliab[BlockNum - 1][j - 1];

	for (int i = spisok.back() - 1; i > 0; i--) {
		//int d1 = BlockDiam[i - 1], d2 = BlockDiam[i];
		for (int diam = 0; diam < gap + 1; diam++) {
			for (int j = 0; j < diam + 1; j++)
				curBlockRel[diam] += BlockReliab[i - 1][j] * BlockReliab[i][diam - j];
		}
		for (int j = 0; j < BlockReliab[i - 1].size(); j++) {
			BlockReliab[i - 1][j] = curBlockRel[j];
			curBlockRel[j] = 0;
		}
	}
	//for (int i = 1; i < BlockReliab[0].size(); i++)
	//	BlockReliab[0][i] += BlockReliab[0][i - 1];
	output << "reliabcalc Time(ms): " << (clock() - reliabcalc) << std::endl;
	output << std::setprecision(16) << BlockReliab[0].back() << std::endl;
	output << "Recursions: " << NumberOfRec << std::endl;
	output << "Decompose Parallel Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
	output1 << (clock() - start_time) / 1000.0000 << std::endl;
}

GraphMethodResult kGraph::ReliabilityDiamConstr2VertRecDecompositionSF(int s, int t, int UpperBound) {
	Nconst = (this->KAO.size() > 1) ? (this->KAO.size() * this->KAO.size()) : 100000;
	clock_t start_time = clock();
	
	int terminaldiam = this->DistanceDijkstra(s, t);

	if (terminaldiam > UpperBound) {
		double time_sec = (clock() - start_time) / 1000.0000;
		std::cout << "Error: terminaldiam > UpperBound" << std::endl;
		return GraphMethodResult(0.0, NumberOfRec, time_sec);
	}
	//!TODO убедиться в необходимости реализации отдельной функции
}


GraphMethodResult kGraph::ReliabilityDiamConstr2VertRecursiveDecomposition(int s_node, int t_node, int UpperBound) {
	Nconst = (this->KAO.size() > 1) ? (this->KAO.size() * this->KAO.size()) : 100000;
	clock_t start_time = clock();
	globsumReliab = 0;

	// 1. Декомпозиция графа на блоки
	std::vector<int> spisok_decomp_original = this->DecomposeOnBlocksK();

	int BlockNum = spisok_decomp_original.back();

	int diamsum = 0;
	if (spisok_decomp_original.empty()) {
		std::cout << "Can't decompose a graph";
		return GraphMethodResult(0, 0, 0);
	}
	if (BlockNum == 1) {
		// Прямой расчёт для одного блока
		int BlockDiam = this->DistanceDijkstra(s_node, t_node);
		if (BlockDiam > UpperBound) {
			double time_sec = (clock() - start_time) / 1000.0000;
			std::cout << "Error: BlockDiam > UpperBound" << std::endl;
			return GraphMethodResult(0.0, NumberOfRec, time_sec);
		}

		globsumReliab = 0.0;
		Factoring2Vert(*this, s_node, t_node, 0, UpperBound, 1.0);

		double reliability = globsumReliab;
		double time_sec = (clock() - start_time) / 1000.0000;
		return GraphMethodResult(reliability, NumberOfRec, time_sec);
	}

	std::vector<int> BlockDiam(BlockNum);

	for (int i = 0; i < BlockNum; i++) {
		kGraph Restored = this->RestoreBlockK(i + 1, spisok_decomp_original);
		int local_x = 0, local_y = 0;
		for (int j = 1; j < Restored.Targets.size(); j++) {
			if (Restored.Targets[j]) {
				if (local_x == 0) local_x = j;
				else if (local_y == 0) local_y = j;
			}
		}
		BlockDiam[i] = Restored.DistanceDijkstra(local_x, local_y);
		diamsum += BlockDiam[i];
	}

	if (diamsum > UpperBound) {
		std::cout << "Diameter sum overflow (diamsum = " + std::to_string(diamsum) + ", UpperBound = " + std::to_string(UpperBound) + ")\n";
		double time_sec = (clock() - start_time) / 1000.0000;
		return GraphMethodResult(0.0, NumberOfRec, time_sec);
	}

	// 2. Построение графа блоков
	std::vector<BlockGraphNode> block_nodes(BlockNum);
	std::map<int, int> original_block_id_to_node_idx;
	for (int i = 0; i < BlockNum; ++i) {
		block_nodes[i].original_block_id = i + 1;
		original_block_id_to_node_idx[i + 1] = i;
	}

	std::vector<BlockGraphEdge> block_edges;
	std::vector<bool> articulation_point_processed(this->KAO.size(), false);

	for (int v_orig = 1; v_orig < this->KAO.size(); ++v_orig) {
		std::vector<int> blocks_for_v = get_blocks_containing_vertex(v_orig, spisok_decomp_original);
		if (blocks_for_v.size() > 1) {
			if (articulation_point_processed[v_orig]) continue;
			articulation_point_processed[v_orig] = true;
			for (size_t i = 0; i < blocks_for_v.size(); ++i) {
				for (size_t j = i + 1; j < blocks_for_v.size(); ++j) {
					int block1_orig_id = blocks_for_v[i];
					int block2_orig_id = blocks_for_v[j];
					if (original_block_id_to_node_idx.count(block1_orig_id) && original_block_id_to_node_idx.count(block2_orig_id)) {
						block_edges.push_back({ original_block_id_to_node_idx[block1_orig_id],
											   original_block_id_to_node_idx[block2_orig_id],
											   v_orig });
					}
				}
			}
		}
	}

	// 3. Определение стартового и конечного блоков
	std::vector<int> s_containing_orig_blocks = get_blocks_containing_vertex(s_node, spisok_decomp_original);
	std::vector<int> t_containing_orig_blocks = get_blocks_containing_vertex(t_node, spisok_decomp_original);

	if (s_containing_orig_blocks.empty() || t_containing_orig_blocks.empty()) {
		double time_sec = (clock() - start_time) / 1000.0000;
		return GraphMethodResult(0.0, NumberOfRec, time_sec);
	}

	int start_block_node_idx = original_block_id_to_node_idx[s_containing_orig_blocks[0]];
	int end_block_node_idx = original_block_id_to_node_idx[t_containing_orig_blocks[0]];

	// 4. Найти путь блоков
	std::vector<int> path_articulation_points_ids;
	std::vector<int> block_path_node_indices = find_path_in_block_graph(
		start_block_node_idx, end_block_node_idx, block_nodes.size(), block_edges, path_articulation_points_ids);

	if (block_path_node_indices.empty()) {
		double time_sec = (clock() - start_time) / 1000.0000;
		return GraphMethodResult(0.0, NumberOfRec, time_sec);
	}

	// 5. Подготовка упорядоченной цепи блоков
	std::vector<int> ordered_original_block_ids_on_path;
	for (int node_idx : block_path_node_indices) {
		ordered_original_block_ids_on_path.push_back(block_nodes[node_idx].original_block_id);
	}

	std::vector<int> final_aps_for_solver(ordered_original_block_ids_on_path.size() + 1, 0);
	for (size_t i = 0; i < path_articulation_points_ids.size(); ++i) {
		final_aps_for_solver[i + 1] = path_articulation_points_ids[i];
	}

	// 6. Вызов рекурсивной функции
	std::vector<double> result_reliabilities = this->solve_recursive_for_block_chain_ordered(
		0, s_node, t_node, UpperBound, spisok_decomp_original,
		ordered_original_block_ids_on_path, final_aps_for_solver);

	double time_sec = (clock() - start_time) / 1000.0000;
	double reliability = result_reliabilities.empty() ? 0.0 : result_reliabilities[std::min((int)result_reliabilities.size() - 1, UpperBound)];
	return GraphMethodResult(reliability, NumberOfRec, time_sec);
}

bool kGraph::CheckDistanceFloyd(const int d) {
	//������� ������ ������ ������� ���������� � ��������� ������
	int N = this->KAO.size();
	std::vector<std::vector<int> > M(N + 1);

	for (int i = 0; i < N; i++)
		M[i].resize(N + 1);
	for (int i = 1; i < N + 1; i++) {
		for (int j = i + 1; j < N + 1; j++) {
			if (this->CheckEdge(i, j))
				M[i][j] = 1;
			else M[i][j] = Nconst;
		}
		for (int k = 1; k < N + 1; k++) { 
			for (i = 1; i < k; i++) {
				for (int j = i + 1; j < k; j++)
					if (M[i][j] > M[i][k] + M[j][k])
						M[i][j] = M[i][k] + M[j][k];
				for (int j = k + 1; j < N + 1; j++) 
					if (M[i][j] > M[i][k] + M[k][j])
						M[i][j] = M[i][k] + M[k][j];
			}
			for (i = k + 1; i < N + 1; i++) { 
				for (int j = i + 1; j < N + 1; j++) 
					if (M[i][j] > M[k][i] + M[k][j])
						M[i][j] = M[k][i] + M[k][j];
			}
		}
	}
	for (int i = 1; i < N; i++)
		if (this->Targets[i])
			for (int j = i + 1; j < N; j++)
				if (this->Targets[j] && (M[i][j] > d))
					return false;
	return true;
}

void kGraph::SearchVertice(std::vector<int>& DfNumber, std::vector<int>& TreeEdges, std::vector<int>& Low, std::vector<int>& Stek, int& r, int& l, int v, std::vector<int>& DOB) {
	int last = 0;
	DfNumber[v] = l;
	Low[v] = l++;
	for (int i = this->KAO[v - 1]; i < this->KAO[v]; i++) {
		if ((!DfNumber[this->FO[i]]) || ((DfNumber[this->FO[i]] < DfNumber[v]) && (!TreeEdges[i]))) {
			last = Stek.size();
			Stek.push_back(i);
			Stek.push_back(this->SearchEdge(this->FO[i], v));
		}
		if (!DfNumber[this->FO[i]]) {
			TreeEdges[i] = 1;
			TreeEdges[this->SearchEdge(this->FO[i], v)] = 1;
			this->SearchVertice(DfNumber, TreeEdges, Low, Stek, r, l, this->FO[i], DOB);

			if (Low[this->FO[i]] >= DfNumber[v]) {
				for (int j = last; j < Stek.size(); j++)
					DOB[Stek[j]] = r;
				r++;
				Stek.resize(last);
			}
			Low[v] = std::min(Low[v], Low[this->FO[i]]);
			//if (Low[v] > Low[this->FO[i]]) 
			//	Low[v] = Low[this->FO[i]];
		}
		else if ((!TreeEdges[i]) && (DfNumber[this->FO[i]] < DfNumber[v]) && (Low[v] > DfNumber[this->FO[i]]))
			Low[v] = DfNumber[this->FO[i]];
	}
}

std::vector<int> kGraph::DecomposeOnBlocks() {
	//Result[i] - ����� ����� � ������ ������ ����� � ������� i(����� ����� �� ������� FO)
	//Result[Result.size() - 1] - ��� - �� ������
	std::vector<int> DfNumber(this->KAO.size()), TreeEdges(this->FO.size()), Low(this->KAO.size()), Stek;
	int r = 1, l = 1;
	std::vector<int> Result(this->FO.size() + 1);
	this->SearchVertice(DfNumber, TreeEdges, Low, Stek, r, l, 1, Result);
	Result[this->FO.size()] = r - 1;
	return Result;
}

std::vector<int> kGraph::DecomposeOnBlocksK() {
	//{���������� K - ����� �� �����.
	//Result[i] - ����� ����� � ������ ������ ����� � ������� i(����� ����� �� ������� FO)
	//Result[length(Result) - 1] - ��� - �� ������, ������� �������� ���.� - �� ��� �������� ����������
	//�����, �������� � �����, ������� �� �������� ���.� - � � �� �������� ���������� ������ ���� � ������� 0
	//����� ���������� � ��������� ������ ��������� � ������� ������� ����� G
	//���� � ����� ���� ���.� - ��, �� Result[i] = 0 ��� ���� i
	//���� � ����� ��� �����, �� Result ������� �� ������ ��������, Result[0] = 0
	//���������� ������� ����������� � ��������� ����� � ����� ��� ������� ������}
	std::vector<int> S = this->DecomposeOnBlocks();
	if (S.size() == 1)
		return S;
	std::vector<int> TargetBlocks(S.back() + 1);
	for (int i = 1; i < this->Targets.back() + 1; i++) //����� �������� �������
		if (this->Targets[i] == 1) {
			bool boolean = true;
			int k = S[this->KAO[i - 1]];
			for (int j = this->KAO[i - 1] + 1; j < this->KAO[i]; j++)
				if (S[j] != k)
					boolean = false;
			if (boolean)
				TargetBlocks[k] = 1;
		}

	for (int i = 1; i < TargetBlocks.size(); i++)
		if (TargetBlocks[i] == 0 && ConnectivityWithoutBlock(i, S) == false)
			TargetBlocks[i] = 1;

	// ����� ��� ���������� targetsOrder
	// std::vector<int> TargetOrder(this->Targets.size());
	// TargetOrder[!y!] = *amount of targetBlocks*;
	for (int i = 1; i < this->KAO.size(); i++)
		if (this->Targets[i] == 0) {
			bool boolean = true;
			int k = 0;
			for (int j = this->KAO[i - 1]; j < this->KAO[i]; j++)
				if (TargetBlocks[S[j]] == 1)
					if (k == 0)
						k = S[j];
					else if (S[j] != k)
						boolean = false;
			if (boolean == false)
				this->Targets[i] = 1;
		}
	std::vector<int> NewNumbers(S.back() + 1);
	for (int i = 0; i < NewNumbers.size(); i++)
		NewNumbers[i] = i;
	int k = 0;
	for (int i = 1; i < TargetBlocks.size(); i++)
		if (TargetBlocks[i] == 0) {
			k++;
			NewNumbers[i] = 0;
			for (int j = i + 1; j < NewNumbers.size(); j++)
				NewNumbers[j]--;
		}
	//S[] - ������ � ������� ������� ������ ����� ����������� ������ �����
	//NewNumbers[] - �������� ��� ������� ����� ����������. ���� 0, ���� ���� ����������, ���� �����- ������� ����� � ����
	// �� ���� � ����� ������� ����� ������������� �����
	//!TODO ��������� �� ���������� ���������
	//���� 
	// NewNumbers[17] = 1;
	// NewNumbers[20] = 2;
	// NewNumbers[24] = 3;
	// NewNumbers[50] = 4;
	//�����
	//NewNumbers[24] = 1;
	//NewNumbers[50] = 2;
	//NewNumbers[20] = 3;
	//NewNumbers[17] = 4;
	for (int i = 0; i < S.size() - 1; i++)
		S[i] = NewNumbers[S[i]];
	S.back() -= k;
	return S;
}

kGraph kGraph::DeleteEdgeK(const int u, const int v) {
	int e = this->SearchEdge(u, v);
	int f = this->SearchEdge(v, u);
	kGraph Result;

	if (e >= this->FO.size() || f >= this->FO.size())
		return *this; // если одно из рёбер не найдено, вернём копию

	// Удалим оба ребра из FO и PArray
	Result.FO.reserve(this->FO.size() - 2);
	Result.PArray.reserve(this->PArray.size() - 2);
	for (int i = 0; i < this->FO.size(); ++i) {
		if (i != e && i != f) {
			Result.FO.push_back(this->FO[i]);
			Result.PArray.push_back(this->PArray[i]);
		}
	}

	// Пересчитаем KAO
	Result.KAO.resize(this->KAO.size());
	for (int i = 1; i < (int)this->KAO.size(); ++i) {
		int count = 0;
		for (int j = this->KAO[i - 1]; j < this->KAO[i]; ++j) {
			if (j != e && j != f)
				++count;
		}
		Result.KAO[i] = Result.KAO[i - 1] + count;
	}

	// Копируем Targets (если это необходимо)
	Result.Targets = this->Targets;

	return Result;
}

kGraph kGraph::InducedKGraph(const std::vector<int> spisok, int Number) const {
    int oldN = (int)spisok.size() - 1;

    // 1) Построим отображение старых ID -> новых [1..k]
    std::vector<int> S(spisok.size(), 0);
    int k = 0;
    for (int i = 1; i <= oldN; ++i) {
        if (spisok[i] == Number || spisok[i] == 0) {
            S[i] = ++k;
        }
    }

    // 2) Временные списки смежности
    // newAdj[v] = вектор пар (u, prob) для новой вершины v (1-based)
    std::vector<std::vector<std::pair<int,double>>> newAdj(k+1);
    for (int oldV = 1; oldV <= oldN; ++oldV) {
        int newV = S[oldV];
        if (newV == 0) continue;

        // сканируем все соседи старой вершины oldV
        for (int idx = KAO[oldV-1]; idx < KAO[oldV]; ++idx) {
            int oldU = FO[idx];
            int newU = S[oldU];
            // включаем ребро, только если новый U существует и Boolka разрешает
            if (newU > 0 && Boolka(spisok, Number, oldV, idx)) {
                newAdj[newV].emplace_back(newU, PArray[idx]);
            }
        }
    }

    // 3) Переводим newAdj в CSR (KAO2, FO2, P2)
    std::vector<int> KAO2(k+1), FO2;
    std::vector<double> P2;
    KAO2[0] = 0;
    for (int v = 1; v <= k; ++v) {
        KAO2[v] = KAO2[v-1] + (int)newAdj[v].size();
    }
    int total = KAO2[k];
    FO2.reserve(total);
    P2 .reserve(total);
    for (int v = 1; v <= k; ++v) {
        for (auto &pr : newAdj[v]) {
            FO2.push_back(pr.first);
            P2.push_back(pr.second);
        }
    }

    // 4) Targets2
    std::vector<int> T2(k+1, 0);
    for (int oldV = 1; oldV <= oldN; ++oldV) {
        int newV = S[oldV];
        if (newV > 0) {
            T2[newV] = this->Targets[oldV];
        }
    }

    return kGraph(std::move(KAO2), std::move(FO2), std::move(P2), std::move(T2));
}

bool kGraph::KComponent() {
	// Проверка: если целевая вершина всего одна, делаем граф тривиальным
	int target_count = 0;
	for (int i = 1; i < this->Targets.size(); ++i) {
		if (this->Targets[i] == 1)
			++target_count;
	}

	if (target_count == 1) {
		this->KAO = { 0, 0 };
		this->FO.clear();
		this->PArray.clear();
		this->Targets = { 0, 0 };
		return true;
	}

	// BFS для обхода компоненты связности по целевым вершинам
	std::vector<int> spot(this->KAO.size(), 2); // 2 — не посещено
	std::queue<int> q;

	// Найдём первую целевую вершину и начнем обход
	int start = -1;
	for (int i = 1; i < this->Targets.size(); ++i) {
		if (this->Targets[i] == 1) {
			start = i;
			break;
		}
	}
	if (start == -1) return false; // на всякий случай

	q.push(start);
	spot[start] = 1;

	int reachable_targets = 1;
	int reachable_total = 1;

	while (!q.empty()) {
		int u = q.front(); q.pop();
		for (int i = this->KAO[u - 1]; i < this->KAO[u]; ++i) {
			int v = this->FO[i];
			if (spot[v] == 2) {
				spot[v] = 1;
				q.push(v);
				if (this->Targets[v] == 1)
					++reachable_targets;
				++reachable_total;
			}
		}
	}

	bool all_targets_reached = (reachable_targets == target_count);

	// Если граф связный, но содержит лишние вершины — обрезаем его
	if (all_targets_reached && reachable_total < (int)this->KAO.size() - 1) {
		*this = this->InducedKGraph(spot, 1); // оставить только достигнутые вершины
	}

	return all_targets_reached;
}

bool kGraph::Boolka(std::vector<int> spisok, int Number, int i, int j) const{
	if (Number == 1)
		return true;
	return spisok[i] || spisok[this->FO[j]];
}

bool kGraph::ConnectivityWithoutBlock(int Block, std::vector<int> S) {
	std::vector<int> A(1);
	for (int i = 1; i < this->Targets.size(); i++)
		if (this->Targets[i] == 1) {
			A[0] = i;
			break;
		}
	int l = A.size(), sum = 1;
	std::vector<int> Spot(this->KAO.size()), B;
	Spot[A[0]] = 1;
	while (l > 0) {
		B.resize(0);
		for (int i = 0; i < A.size(); i++)
			for (int j = this->KAO[A[i] - 1]; j < this->KAO[A[i]]; j++)
				if ((!Spot[this->FO[j]]) && (S[j] != Block)) {
					B.resize(B.size() + 1, this->FO[j]);
					Spot[this->FO[j]] = 1;
					sum += this->Targets[this->FO[j]] == 1;
				}
		A.resize(0);
		l = B.size();
		if (l > 0)
			for (int i = 0; i < l; i++)
				A.push_back(B[i]);
	}

	int sum1 = 0;
	for (int i = 1; i < this->Targets.size(); i++)
		if (this->Targets[i] == 1)
			sum1++;
	if (this->KAO.size() == 1)
		return true;
	return sum == sum1;
}

kGraph kGraph::RestoreBlockK(int Number, const std::vector<int>& spisok) {
	//{��������������� ���� K - ����� G � ������� Number,
	// ��������� spisok �������������� �������������� ���������� DecopmoseOnBlocksK
	// ��� �� ��������� ��� ��������� ������ ���.� - � ������������ ������� ����������}
	std::vector<int> S(this->KAO.size()), BackS(this->KAO.size());
	int j = 0;
	for (int i = 0; i < spisok.size() - 1; i++)
		if ((spisok[i] == Number) && !S[this->FO[i]]) {
			j++;
			S[this->FO[i]] = j;
			BackS[j] = this->FO[i];
		}
	std::vector<int> FO(this->FO.size()), KAO(j + 1), Targets(j + 1);
	std::vector<double> Parray(this->FO.size());
	BackS.resize(j + 1);
	for (int i = 1; i < Targets.size(); i++)
		Targets[i] = this->Targets[BackS[i]];
	int l = 0;
	for (int i = 1; i < BackS.size(); i++)
		for (j = this->KAO[BackS[i] - 1]; j < this->KAO[BackS[i]]; j++)
			if (spisok[j] == Number) {
				FO[l] = S[this->FO[j]];
				Parray[l] = this->PArray[j];
				l++;
				KAO[i] = l;
			}
	FO.resize(l);
	Parray.resize(l);
	kGraph Result(KAO, FO, Parray, Targets);
	return Result;
}

kGraph kGraph::RestoreBlock(int Number, const std::vector<int>& spisok) {
	std::vector<int> S(this->KAO.size()), FO;
	std::vector<double> Parray(0);
	int j = 0;
	int minFO = this->KAO.size();
	for (int i = 0; i < spisok.size() - 1; i++)
		if (spisok[i] == Number) {
			FO.push_back(this->FO[i]);
			Parray.push_back(this->PArray[i]);
			if (!S[this->FO[i]]) {
				j++;
				minFO = std::min(minFO, this->FO[i]);
			}
			S[this->FO[i]]++;
			//BackS[j] = this->FO[i];
		}
	//std::vector<int> KAO(j + 1), Targets(j + 1);
	std::vector<int> KAO(1), Targets(1);


	for (int i = 1; i < S.size(); i++) {
		if (S[i]) {
			KAO.push_back(KAO.back() + S[i]);
			Targets.push_back(this->Targets[i]);
		}
	}
	//for (int i = 1; i < KAO.size(); i++) {
	//    KAO[i] += S[i];
	//    if (S[i] && this->Targets[i])
	//        Targets[i] = 1;
	//}

	for (int i = 0; i < FO.size(); i++)
		FO[i] -= minFO - 1;
	kGraph Result(KAO, FO, Parray, Targets);
	return Result;
}

kGraph kGraphFileInput(std::ifstream& fin) {

	if (!fin) {
		std::cout << "Error!\n";
		throw std::runtime_error("OPEN_ERROR");
	}
	std::string line;
	std::getline(fin, line, '\n');
	std::istringstream skao(line);
	std::string temp;
	std::vector<int> KAO;
	while (std::getline(skao, temp, ','))
		KAO.push_back(std::stoi(temp));

	std::getline(fin, line, '\n');
	std::istringstream sfo(line);
	std::vector<int> FO;
	while (std::getline(sfo, temp, ','))
		FO.push_back(std::stoi(temp));

	std::getline(fin, line, '\n');
	std::istringstream stargets(line);
	std::vector<int> Targets;
	while (std::getline(stargets, temp, ','))
		Targets.push_back(std::stoi(temp));
	std::getline(fin, line);

	size_t found = line.find(',');
	if (found != std::string::npos)
		line[found] = '.';
	double p = std::stod(line);
	std::vector<double> Parray(FO.size(), p);
	kGraph G(KAO, FO, Parray, Targets);
	return G;
}

void kGraph::kGraphFileOutput(std::ofstream& fout) {
	for (int i = 0; i < this->KAO.size() - 1; i++) {
		fout << KAO[i] << ",";
	}
	fout << this->KAO.back() << std::endl;
	for (int i = 0; i < this->FO.size() - 1; i++) {
		fout << FO[i] << ",";
	}
	fout << this->FO.back() << std::endl;
	for (int i = 0; i < this->Targets.size() - 1; i++) {
		fout << Targets[i] << ",";
	}
	fout << this->Targets.back() << std::endl;
	fout << this->PArray[0];
}

// Вспомогательная функция для определения, каким блокам принадлежит вершина (на основе ребер)
std::vector<int> kGraph::get_blocks_containing_vertex(
	int vertex_original_id,
	const std::vector<int>& spisok_decomposition
) const {
	std::vector<int> blocks;
	std::vector<bool> found_block(spisok_decomposition.back() + 1, false); // spisok_decomposition.back() - общее число блоков

	if (vertex_original_id == 0 || vertex_original_id >= this->KAO.size()) return blocks; // Неверный ID

	for (int edge_idx = this->KAO[vertex_original_id - 1]; edge_idx < this->KAO[vertex_original_id]; ++edge_idx) {
		if (edge_idx < spisok_decomposition.size() - 1) { // Последний элемент spisok_decomposition - общее число блоков
			int block_num = spisok_decomposition[edge_idx];
			if (block_num > 0 && !found_block[block_num]) {
				blocks.push_back(block_num);
				found_block[block_num] = true;
			}
		}
	}
	return blocks;
}



// НОВАЯ версия solve_recursive_for_block_chain, работающая с упорядоченной цепью
std::vector<double> kGraph::solve_recursive_for_block_chain_ordered(
    int current_block_path_idx,         // Индекс текущего блока в ordered_original_block_ids_on_path
    int entry_node_original_id,         // ID узла входа в этот блок (в нумерации G_original)
    int target_node_original_id,        // Конечный узел t (в нумерации G_original)
    int max_len_budget,                 // Максимально допустимая длина пути для текущего вызова
    const std::vector<int>& spisok_decomposition_G_original, // Для get_block_graph_and_map_ids
    const std::vector<int>& ordered_original_block_ids_on_path, // Упорядоченная цепь ID блоков
    const std::vector<int>& ordered_articulation_points_orig_ids // ordered_aps[i] - точка между блоком i-1 и i в пути
) const { // Сделал const

    if (max_len_budget < 0 || current_block_path_idx >= ordered_original_block_ids_on_path.size()) return {};

    int original_block_id_to_restore = ordered_original_block_ids_on_path[current_block_path_idx];

    kGraph current_block_graph;
    std::vector<int> map_new_to_orig, map_orig_to_new;
    // Важно: get_block_graph_and_map_ids вызывается от G_original
    this->get_block_graph_and_map_ids(original_block_id_to_restore, spisok_decomposition_G_original,
                                current_block_graph, map_new_to_orig, map_orig_to_new);

    if (current_block_graph.KAO.size() <= 1 || map_orig_to_new[entry_node_original_id] == 0) {
        return std::vector<double>(max_len_budget + 1, 0.0);
    }
    int s_in_block = map_orig_to_new[entry_node_original_id];

    // БАЗОВЫЙ СЛУЧАЙ: текущий блок является последним в пути И содержит target_node_original_id
    bool is_last_block_in_path = (current_block_path_idx == ordered_original_block_ids_on_path.size() - 1);
    
    std::vector<int> blocks_of_target = this->get_blocks_containing_vertex(target_node_original_id, spisok_decomposition_G_original);
    bool target_in_current_original_block = false;
    for(int b_id : blocks_of_target) {
        if (b_id == original_block_id_to_restore) {
            target_in_current_original_block = true;
            break;
        }
    }

    if (is_last_block_in_path && target_in_current_original_block) {
        if (map_orig_to_new[target_node_original_id] == 0) {
            return std::vector<double>(max_len_budget + 1, 0.0);
        }
        int t_in_block = map_orig_to_new[target_node_original_id];
        std::vector<double> R_cumulative(max_len_budget + 1, 0.0);
        double global_globsumReliab_backup = globsumReliab;

        int d_min_this_block_path = current_block_graph.DistanceDijkstra(s_in_block, t_in_block);
        for (int d = 0; d <= max_len_budget; ++d) {
            if (d < d_min_this_block_path && d_min_this_block_path != Nconst) {
                R_cumulative[d] = 0.0;
            } else {
                globsumReliab = 0.0;
                Factoring2Vert(current_block_graph, s_in_block, t_in_block, 0, d, 1.0);
                R_cumulative[d] = globsumReliab;
            }
        }
        for (int d = 1; d <= max_len_budget; ++d) {
            if (R_cumulative[d] < R_cumulative[d - 1]) R_cumulative[d] = R_cumulative[d-1];
        }
        globsumReliab = global_globsumReliab_backup;
        return R_cumulative;
    }

    // РЕКУРСИВНЫЙ ШАГ
    if (is_last_block_in_path && !target_in_current_original_block) { // Последний блок, но t не в нем - ошибка пути
         return std::vector<double>(max_len_budget + 1, 0.0);
    }
    // Если не последний блок, нужна точка выхода
    // ordered_articulation_points_orig_ids[current_block_path_idx + 1] это точка между current и next
    if (current_block_path_idx + 1 >= ordered_articulation_points_orig_ids.size() || ordered_articulation_points_orig_ids[current_block_path_idx + 1] == 0) {
        return std::vector<double>(max_len_budget + 1, 0.0); // Нет точки выхода
    }
    int exit_node_original_id = ordered_articulation_points_orig_ids[current_block_path_idx + 1];
    if (map_orig_to_new.empty() || exit_node_original_id >= map_orig_to_new.size() || map_orig_to_new[exit_node_original_id] == 0) {
        return std::vector<double>(max_len_budget + 1, 0.0); // Точка выхода не в текущем восстановленном блоке
    }
    int x_in_block = map_orig_to_new[exit_node_original_id];


    // 1. Рассчитать R_cumul_s_x (s_in_block -> x_in_block) в current_block_graph
    std::vector<double> R_cumul_s_x(max_len_budget + 1, 0.0);
    double global_globsumReliab_backup_sx = globsumReliab;
    int d_min_s_x = current_block_graph.DistanceDijkstra(s_in_block, x_in_block);
    for (int d = 0; d <= max_len_budget; ++d) {
        if (d < d_min_s_x && d_min_s_x != Nconst) {
             R_cumul_s_x[d] = 0.0;
        } else {
            globsumReliab = 0.0;
            Factoring2Vert(current_block_graph, s_in_block, x_in_block, 0, d, 1.0); // Здесь происходит пересчет для одних и тех же значений диаметра??
            R_cumul_s_x[d] = globsumReliab;
        }
    }
    for (int d = 1; d <= max_len_budget; ++d) {
        if (R_cumul_s_x[d] < R_cumul_s_x[d-1]) R_cumul_s_x[d] = R_cumul_s_x[d-1];
    }
    globsumReliab = global_globsumReliab_backup_sx;

    std::vector<double> R_bar_s_x(max_len_budget + 1);
    R_bar_s_x[0] = R_cumul_s_x[0];
    for (int d = 1; d <= max_len_budget; ++d) {
        R_bar_s_x[d] = R_cumul_s_x[d] - R_cumul_s_x[d - 1];
    }

    // 2. Рекурсивный вызов для остальной части цепи
    std::vector<double> R_cumul_x_t_rest = this->solve_recursive_for_block_chain_ordered(
        current_block_path_idx + 1, exit_node_original_id, target_node_original_id,
        max_len_budget, spisok_decomposition_G_original,
        ordered_original_block_ids_on_path, ordered_articulation_points_orig_ids);

    // 3. Комбинирование
    std::vector<double> final_R_for_this_call(max_len_budget + 1, 0.0);
    for (int d_target = 0; d_target <= max_len_budget; ++d_target) {
        double sum_for_this_d_target = 0.0;
        for (int len1 = 0; len1 <= d_target; ++len1) {
            if (R_bar_s_x[len1] == 0.0) continue;
            int len_for_rest_budget = d_target - len1;
            if (len_for_rest_budget < 0) continue;
            double prob_rest_of_chain = 0.0;
            if (len_for_rest_budget < R_cumul_x_t_rest.size()) {
                prob_rest_of_chain = R_cumul_x_t_rest[len_for_rest_budget];
            }
            sum_for_this_d_target += R_bar_s_x[len1] * prob_rest_of_chain;
        }
        final_R_for_this_call[d_target] = sum_for_this_d_target;
    }
    for (int d = 1; d <= max_len_budget; ++d) {
        if (final_R_for_this_call[d] < final_R_for_this_call[d-1]) final_R_for_this_call[d] = final_R_for_this_call[d-1];
    }
    return final_R_for_this_call;
}

// Вспомогательная функция для восстановления графа блока и карт вершин
// Принимает *this (исходный граф) как G_original неявно
void kGraph::get_block_graph_and_map_ids(
    int block_idx_to_restore,
    const std::vector<int>& spisok_decomposition,
    kGraph& out_block_graph,
    std::vector<int>& out_map_new_id_to_original_id,
    std::vector<int>& out_map_original_id_to_new_id
) const { // Добавлено const
    out_map_original_id_to_new_id.assign(this->KAO.size(), 0);
    out_map_new_id_to_original_id.assign(1, 0); // 0-й элемент не используется (для 1-based new_id)

    int new_vertex_id_counter = 0;
    std::vector<bool> is_vertex_in_block_definitively(this->KAO.size(), false);

    // Определяем вершины, которые точно принадлежат этому блоку
    // (т.е. не являются точками сочленения с другими блоками, которые мы сейчас не рассматриваем,
    // или являются точками сочленения, но мы их включаем как часть этого блока)
    for (int v_orig = 1; v_orig < this->KAO.size(); ++v_orig) {
        bool part_of_block = false;
        if (this->KAO[v_orig-1] < this->KAO[v_orig]) { // Если вершина не изолирована
             for (int edge_idx = this->KAO[v_orig - 1]; edge_idx < this->KAO[v_orig]; ++edge_idx) {
                if (edge_idx < spisok_decomposition.size() -1 && spisok_decomposition[edge_idx] == block_idx_to_restore) {
                    part_of_block = true;
                    break;
                }
            }
        } else if (this->Targets[v_orig] && block_idx_to_restore == 1) { // Особый случай для изолированных целевых вершин в первом блоке
            // Это эвристика, может потребовать уточнения.
            // Если s или t - изолированные вершины, DecomposeOnBlocksK может их не включить в реберные блоки.
            // Здесь мы предполагаем, что если Targets[v_orig] = 1 и мы восстанавливаем блок 1, то это может быть s или t.
            // Лучше, если DecomposeOnBlocksK корректно обрабатывает такие случаи или если входные данные не имеют изолированных s/t.
        }


        if (part_of_block) {
            is_vertex_in_block_definitively[v_orig] = true;
            new_vertex_id_counter++;
            out_map_original_id_to_new_id[v_orig] = new_vertex_id_counter;
            out_map_new_id_to_original_id.push_back(v_orig);
        }
    }
    
    // Если блок пуст (например, неверный block_idx_to_restore или пустая декомпозиция)
    if (new_vertex_id_counter == 0) {
        out_block_graph.KAO.assign(1,0);
        out_block_graph.FO.clear();
        out_block_graph.PArray.clear();
        out_block_graph.Targets.assign(1,0);
        // Карты уже инициализированы правильно (пустые или с 0)
        return;
    }


    std::vector<int> KAO_new(new_vertex_id_counter + 1);
    std::vector<int> FO_new_edges;
    std::vector<double> Parray_new_edges;
    std::vector<int> Targets_new(new_vertex_id_counter + 1, 0);

    KAO_new[0] = 0;
    for (int new_v_idx = 1; new_v_idx <= new_vertex_id_counter; ++new_v_idx) {
        KAO_new[new_v_idx] = KAO_new[new_v_idx - 1];
        int original_v_id_current = out_map_new_id_to_original_id[new_v_idx];
        
        if (original_v_id_current < this->Targets.size()) {
             Targets_new[new_v_idx] = this->Targets[original_v_id_current];
        } else {
            // Это может произойти, если KAO.size() больше Targets.size(), что странно.
            // Обычно они должны быть одинакового размера (N+1).
        }

        for (int edge_orig_idx = this->KAO[original_v_id_current - 1]; edge_orig_idx < this->KAO[original_v_id_current]; ++edge_orig_idx) {
            if (edge_orig_idx < spisok_decomposition.size() -1 && spisok_decomposition[edge_orig_idx] == block_idx_to_restore) {
                int original_adj_v_id = this->FO[edge_orig_idx];
                if (out_map_original_id_to_new_id[original_adj_v_id] != 0) { // Если смежная вершина также в этом блоке
                    FO_new_edges.push_back(out_map_original_id_to_new_id[original_adj_v_id]);
                    Parray_new_edges.push_back(this->PArray[edge_orig_idx]);
                    KAO_new[new_v_idx]++;
                }
            }
        }
    }
    out_block_graph.KAO = KAO_new;
    out_block_graph.FO = FO_new_edges;
    out_block_graph.PArray = Parray_new_edges;
    out_block_graph.Targets = Targets_new;
}


void kGraph::convertEdgeListToKAOFO(const std::string& inputPath,
	const std::string& outputPath,
	double reliability) {
	// Базовая конверсия
	Graph::convertEdgeListToKAOFO(inputPath, outputPath, reliability);

	// Обработка Targets
	int n = static_cast<int>(KAO.size()) - 2;
	Targets.assign(n, 0);
	std::ofstream outFile(outputPath, std::ios::app);
	for (int i = 0; i < n; ++i)
		outFile << Targets[i] << (i + 1 < n ? ',' : '\n');
}

void kGraph::convertKAOFOToEdgeList(const std::string& inputPath,
	const std::string& outputPath) {
	// Базовый метод восстанавливает список рёбер
	Graph::convertKAOFOToEdgeList(inputPath, outputPath);
	// TODO: чтение Targets, если потребуется
}