#include "kGraphOperations.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>

std::ofstream output("C:/Users/User/source/repos/graduate work/graduate work/output.txt");
std::ofstream output1("C:/Users/User/source/repos/graduate work/graduate work/output1.txt", std::ios_base::app);
int NumberOfRec = 0;
std::vector<double> sumReliab;
double globsumReliab = 0;
std::vector<std::vector<double>> BlockReliab;


void Factoring(kGraph G, const int variant, const int d, double Reliab) {
	//результат лежит в sumReliab[0]
	//Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра
	int i, j, k;
	NumberOfRec++;
	if (!variant) {
		if (!G.KComponent())
			return; // add 0 to sum
	}
	else if (!G.CheckDistanceFloyd(d))
		return; // add 0 to sum

	i = G.PArray.size();  //можно ли заменить PArray на FO? 
	for (j = G.PArray.size() - 1; j >= 0; j--)
		if (G.PArray[j] < 1) {
			i = j;
			break;
		}
	if (i == G.FO.size()) {
		sumReliab[0] += Reliab;
		return; // add Reliab to sum
	}

	double p = G.PArray[i];
	G.PArray[i] = 1;
	for (j = 1; j < G.KAO.size(); j++)
		if (G.KAO[j] > i)
			break;
	for (k = G.KAO[G.FO[i] - 1]; k < G.KAO[G.FO[i]]; k++)
		if (G.FO[k] == j)
			break;
	G.PArray[k] = 1;
	Factoring(G, 1, d, Reliab * p);

	kGraph T = G.DeleteEdgeK(G.FO[i], j);
	Factoring(T, 0, d, Reliab * (1 - p));
}

void Factoring2Vert(kGraph G, const int x, const int y, const int variant, const int d, double Reliab) {
	//результат лежит в globsumReliab
	//Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра
	NumberOfRec++;
	if (!variant && G.DistanceDijkstra(x, y) > d)
		return;// add 0 to sum if distance from x to y > d

	int i = G.PArray.size();//можно ли заменить PArray на FO?
	for (int j = G.PArray.size() - 1; j >= 0; j--)
		if (G.PArray[j] < 1) {
			i = j;
			break;
		}
	if (i == G.FO.size()) {
		globsumReliab += Reliab;
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
	Factoring2Vert(G, x, y, 1, d, Reliab * p);
	G.PArray[i] = p;
	G.PArray[k] = p;
	kGraph T = G.DeleteEdgeK(G.FO[i], j);
	Factoring2Vert(T, x, y, 0, d, Reliab * (1 - p));
}

void Factoring2VertM(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound) {
	//Заполняет вектор надежностей с соответствующим ограничением на диаметр
	//Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра
	NumberOfRec++;
	int dist = G.DistanceDijkstra(x, y);
	if (!variant && dist > d) {
		if (dist <= UpperBound) {
			NumberOfRec--;
			Factoring2VertM(G, x, y, 0, d + 1, Reliab, LowerBound, UpperBound);
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
		sumReliab[d - LowerBound] += Reliab;
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
	Factoring2VertM(G, x, y, 1, d, Reliab * p, LowerBound, UpperBound);
	G.PArray[i] = p;
	G.PArray[k] = p;
	kGraph T = G.DeleteEdgeK(G.FO[i], j);
	Factoring2VertM(T, x, y, 0, d, Reliab * (1 - p), LowerBound, UpperBound);
}

void Factoring2VertMParallel(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound) {
	//Заполняет вектор надежностей с соответствующим ограничением на диаметр
	//Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра
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
	//{Формирует граф, полученный объединением G1 и G2 по их первым k вершинам}
	//{NN - массив с номерами в-н G2 в G; номера в-н G1 в G совпадает с их номерами в G1}
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
		PArray[i] = 0.9; // плохо, надо брать вероятности, соответствующие ребрам
	std::vector<int> Targets(KAO.size());
	//Targets тоже надо заполнять соответствующим образом
	kGraph Result(KAO, FO, PArray, Targets);

	return Result;
}

void kGraph::ChangeVertex(int u, int v) {
	//Меняет в графе вершины u и v местами (перенумеровывает)
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

void kGraph::ReliabilityDiamConstr(int d) {
	//{Расчет вероятности связности с ограничением на диаметр d методом ветвления.
	// Используется выделение комп-ты связности с целевыми в-ми, разделение ветвей,
	// встроенная ф-я проверки на расстояния}
	Nconst = this->KAO.size() * this->KAO.size();
	//double t = Factoring(G, 0, d);
	//return t;
	Factoring(*this, 0, d, 1);
}

void kGraph::ReliabilityDiamConstr2Vert(int x, int y, int d) {
	//    {Расчет в-ти св-ти двух вершин x,y с огр-м на диаметр d методом ветвления.
	//     Используется разделение ветвей, встроенная ф-я проверки на расстояния;
	//     не используется выделение компонентв с двумя в-ми (так быстрее)}
	Nconst = this->KAO.size() * this->KAO.size();
	clock_t start_time = clock();
	Factoring2Vert(*this, x, y, 0, d, 1.0);
	output << std::setprecision(16) << globsumReliab << std::endl;
	output << "Recursions: " << NumberOfRec << std::endl;
	output1 << (clock() - start_time) / 1000.0000 << std::endl;
}

void kGraph::ReliabilityDiamConstr2VertDecomposeSimpleFacto(int x, int y, const int& UpperBound) {
	Nconst = this->KAO.size() * this->KAO.size();
	clock_t start_time = clock();
	std::vector<int> spisok = this->DecomposeOnBlocksK();
	int BlockNum = spisok[spisok.size() - 1];
	//если 1 блок
	if (BlockNum == 1) {
		Factoring2Vert(*this, x, y, 0, UpperBound, 1);
		output << std::setprecision(16) << globsumReliab << std::endl;
		output << "Recursions: " << NumberOfRec << std::endl;
		output << "Decompose SimpleFacto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
		return;
	}

	//если >1
	//spisok = SpisokSort;
	//std::vector<int> newSpisok(spisok);
	//int blocknum;
	//for (int i = 1; i < spisok[spisok.size() - 1] + 1; i++) {
	//	blocknum = spisok[KAO[x] - 1];
	//}
	std::vector<int> BlockDiam(BlockNum);
	int diamsum = 0;  //сумма всех диаметров

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

	//diamsum - минимальное расстояние между целевыми вершинами
	if (diamsum > UpperBound) {
		output << "Error: diamsum > UpperBound" << std::endl;
		return;
	}
	int gap = UpperBound - diamsum;
	BlockReliab.resize(BlockNum, std::vector<double>(gap + 1));
	start_time = clock();
	for (int i = 0; i < BlockNum; i++) {
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
		for (int j = 0; j < gap + 1; j++) {
			Factoring2Vert(Restored, x, y, 0, BlockDiam[i] + j, 1);
			BlockReliab[i][j] = globsumReliab;
			globsumReliab = 0;
		}
	}
	clock_t reliabcalc = clock();
	//между n блоками n - 1 точек сочленения
	std::vector<double> curBlockRel(gap + 1);
	//for (int j = 1; j < BlockReliab[BlockNum - 1].size(); j++)
	//	BlockReliab[BlockNum - 1][j] += BlockReliab[BlockNum - 1][j - 1];

	for (int i = BlockNum - 1; i > 0; i--) {
		// d1 = BlockDiam[i - 1], d2 = BlockDiam[i];
		for (int diam = 0; diam < gap + 1; diam++) {
			for (int j = 0; j < diam + 1; j++)
				curBlockRel[diam] = BlockReliab[i - 1][j] * BlockReliab[i][diam - j];
		}
		for (int j = 0; j < BlockReliab[i - 1].size(); j++) {
			BlockReliab[i - 1][j] = curBlockRel[j];
			curBlockRel[j] = 0;
		}
	}
	//for (int i = 1; i < BlockReliab[0].size(); i++)
	//	BlockReliab[0][i] += BlockReliab[0][i - 1];
	output << "reliabcalc Time(ms): " << (clock() - reliabcalc) << std::endl;
	output << std::setprecision(16) << BlockReliab[0][BlockReliab[0].size() - 1] << std::endl;
	output << "Recursions: " << NumberOfRec << std::endl;
	output << "Decompose Simple Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
	output1 << (clock() - start_time) / 1000.0000 << std::endl;
}

void kGraph::ReliabilityDiamConstr2VertMDecompose(int x, int y, const int& UpperBound) {
	//нет правильной нумерации блоков. вследствие чего для GEANTа блоки нумеровались вручную
	Nconst = this->KAO.size() * this->KAO.size();
	clock_t start_time = clock();
	std::vector<int> spisok = this->DecomposeOnBlocksK();
	int BlockNum = spisok[spisok.size() - 1];
	//если 1 блок
	if (BlockNum == 1) {
		int BlockDiam = this->DistanceDijkstra(x, y);
		sumReliab.resize(UpperBound - BlockDiam + 1);
		Factoring2VertM(*this, x, y, 0, UpperBound, 1, BlockDiam, UpperBound);
		for (int i = 1; i < sumReliab.size(); i++)
			sumReliab[i] += sumReliab[i - 1];
		output << std::setprecision(16) << sumReliab[sumReliab.size() - 1] << std::endl;
		output << "Recursions: " << NumberOfRec << std::endl;
		output << "Decompose Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
		return;
	}

	//если >1
	//spisok = SpisokSort;
	//std::vector<int> newSpisok(spisok);
	//int blocknum;
	//for (int i = 1; i < spisok[spisok.size() - 1] + 1; i++) {
	//	blocknum = spisok[KAO[x] - 1];
	//}
	std::vector<int> BlockDiam(BlockNum);
	int diamsum = 0;  //сумма всех диаметров

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
	
	//diamsum - сумма диаметров всех блоков
	if (diamsum > UpperBound) {
		output << "Diameter sum overflow (diamsum = " + std::to_string(diamsum) + ", UpperBound = " + std::to_string(UpperBound) + ")\n";
		return;
	}

	int gap = UpperBound - diamsum;
	start_time = clock();
	for (int i = 0; i < BlockNum; i++) {
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
		sumReliab.resize(0);
		sumReliab.resize(gap + 1);
		Factoring2VertM(Restored, x, y, 0, BlockDiam[i], 1, BlockDiam[i], gap + BlockDiam[i]);
		BlockReliab.push_back(sumReliab);//получили sumreliab для блока, сохраняем его в BlockReliab
	}
	clock_t reliabcalc = clock();
	//между n блоками n - 1 точек сочленения
	std::vector<double> curBlockRel(gap + 1);
	for (int j = 1; j < BlockReliab[BlockNum - 1].size(); j++)
		BlockReliab[BlockNum - 1][j] += BlockReliab[BlockNum - 1][j - 1];

	for (int i = BlockNum - 1; i > 0; i--) {
		//int d1 = BlockDiam[i - 1], d2 = BlockDiam[i];
		for (int diam = 0; diam < gap + 1; diam++) {
			for (int j = 0; j < diam + 1; j++)
				curBlockRel[diam] += BlockReliab[i - 1][j] * BlockReliab[i][diam - j];
		}
		for (int j = 0; j < gap + 1; j++) {
			BlockReliab[i - 1][j] = curBlockRel[j];
			curBlockRel[j] = 0;
		}
	}
	//for (int i = 1; i < BlockReliab[0].size(); i++)
	//	BlockReliab[0][i] += BlockReliab[0][i - 1];
	output << "reliabcalc Time(ms): " << (clock() - reliabcalc) << std::endl;
	output << std::setprecision(16) << BlockReliab[0][BlockReliab[0].size() - 1] << std::endl;
	output << "Recursions: " << NumberOfRec << std::endl;
	output << "Decompose Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
	//output1 << std::setprecision(16) << BlockReliab[0][BlockReliab[0].size() - 1] << std::endl;
	output1 << (clock() - start_time) / 1000.0000 << std::endl;
}

void kGraph::ReliabilityDiamConstr2VertMDecomposeParallel(int x, int y, const int& UpperBound) {
	//для сборки Release в VS распараллеливание не работает
	Nconst = this->KAO.size() * this->KAO.size();
	clock_t start_time = clock();
	std::vector<int> spisok = this->DecomposeOnBlocksK();
	int BlockNum = spisok[spisok.size() - 1];
	//если 1 блок
	if (BlockNum == 1) {
		int BlockDiam = this->DistanceDijkstra(x, y);
		sumReliab.resize(UpperBound - BlockDiam + 1);
		Factoring2VertM(*this, x, y, 0, UpperBound, 1, BlockDiam, UpperBound);
		for (int i = 1; i < sumReliab.size(); i++)
			sumReliab[i] += sumReliab[i - 1];
		output << std::setprecision(16) << sumReliab[sumReliab.size() - 1] << std::endl;
		output << "Recursions: " << NumberOfRec << std::endl;
		output << "Decompose Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
		return;
	}

	std::vector<int> BlockDiam(BlockNum);
	int diamsum = 0;  //сумма всех диаметров

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

	//diamsum - минимальное расстояние между целевыми вершинами
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
	//между n блоками n - 1 точек сочленения
	std::vector<double> curBlockRel(gap + 1);
	for (int j = 1; j < BlockReliab[BlockNum - 1].size(); j++)
		BlockReliab[BlockNum - 1][j] += BlockReliab[BlockNum - 1][j - 1];

	for (int i = spisok[spisok.size() - 1] - 1; i > 0; i--) {
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
	output << std::setprecision(16) << BlockReliab[0][BlockReliab[0].size() - 1] << std::endl;
	output << "Recursions: " << NumberOfRec << std::endl;
	output << "Decompose Parallel Facto Time(sec): " << (clock() - start_time) / 1000.0000 << std::endl;
	output1 << (clock() - start_time) / 1000.0000 << std::endl;
}

bool kGraph::CheckDistanceFloyd(const int d) {
	//Методом Флойда строит матрицу расстояний и проверяет нужные
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
	//Result[i] - номер блока в котрый попало ребро с номером i(номер ребра по массиву FO)
	//Result[Result.size() - 1] - кол - во блоков
	std::vector<int> DfNumber(this->KAO.size()), TreeEdges(this->FO.size()), Low(this->KAO.size()), Stek;
	int r = 1, l = 1;
	std::vector<int> Result(this->FO.size() + 1);
	this->SearchVertice(DfNumber, TreeEdges, Low, Stek, r, l, 1, Result);
	Result[this->FO.size()] = r - 1;
	return Result;
}

std::vector<int> kGraph::DecomposeOnBlocksK() {
	//{Разложение K - графа на блоки.
	//Result[i] - номер блока в котрый попало ребро с номером i(номер ребра по массиву FO)
	//Result[length(Result) - 1] - кол - во блоков, которые содержат цел.в - ны или являются связующими
	//Ребра, входящие в блоки, которые не содержат цел.в - н и не являются связующими входят блок с номером 0
	//Точки сочленения в связующих блоках заносятся в целевые вершины графа G
	//Если в графе одна цел.в - на, то Result[i] = 0 для всех i
	//Если в графе нет ребер, то Result состоит из одного элемента, Result[0] = 0
	//Применение функции некорректно к несвязому графу и графу без целевых вершин}
	std::vector<int> S = this->DecomposeOnBlocks();
	if (S.size() == 1)
		return S;
	std::vector<int> TargetBlocks(S[S.size() - 1] + 1);
	for (int i = 1; i < this->Targets[this->Targets.size() - 1] + 1; i++) //очень странная строчка
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

	// место для объявления targetsOrder
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
	std::vector<int> NewNumbers(S[S.size() - 1] + 1);
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
	//S[] - массив в котором указано какому блоку принадлежит каждое ребро
	//NewNumbers[] - содержит для каждого блока информацию. Либо 0, если блок отсекается, либо число- позиция блока в цепи
	// то есть в каком порядке будут обсчитываться блоки
	//!TODO исправить на нормальную нумерацию
	//было 
	// NewNumbers[17] = 1;
	// NewNumbers[20] = 2;
	// NewNumbers[24] = 3;
	// NewNumbers[50] = 4;
	//стало
	//NewNumbers[24] = 1;
	//NewNumbers[50] = 2;
	//NewNumbers[20] = 3;
	//NewNumbers[17] = 4;
	for (int i = 0; i < S.size() - 1; i++)
		S[i] = NewNumbers[S[i]];
	S[S.size() - 1] -= k;
	return S;
}

kGraph kGraph::DeleteEdgeK(const int u, const int v) {
	//    Удаляет из графа ребро(u, v)
	int e, f;
	e = this->SearchEdge(u, v);
	f = this->SearchEdge(v, u);
	kGraph Result;
	Result.Targets.resize(this->Targets.size());
	if ((e < this->FO.size()) && (f < this->FO.size())) {
		Result.KAO.resize(this->KAO.size());
		Result.FO.resize(this->FO.size() - 2);
		Result.PArray.resize(this->PArray.size() - 2);
		if (e < f) {
			for (int i = 0; i < u; i++)
				Result.KAO[i] = this->KAO[i];
			for (int i = u; i < v; i++)
				Result.KAO[i] = this->KAO[i] - 1;
			for (int i = v; i < this->KAO.size(); i++)
				Result.KAO[i] = this->KAO[i] - 2;
			if (e > 0)
				for (int i = 0; i < e; i++) {
					Result.FO[i] = this->FO[i];
					Result.PArray[i] = this->PArray[i];
				}
			if (f > e + 1)
				for (int i = e; i < f - 1; i++) {
					Result.FO[i] = this->FO[i + 1];
					Result.PArray[i] = this->PArray[i + 1];
				}
			if (this->FO.size() > f + 1)
				for (int i = f - 1; i < this->FO.size() - 2; i++) {
					Result.FO[i] = this->FO[i + 2];
					Result.PArray[i] = this->PArray[i + 2];
				}
		}
		else {
			for (int i = 0; i < v; i++)
				Result.KAO[i] = this->KAO[i];
			for (int i = v; i < u; i++)
				Result.KAO[i] = this->KAO[i] - 1;
			for (int i = u; i < this->KAO.size(); i++)
				Result.KAO[i] = this->KAO[i] - 2;
			if (f > 0)
				for (int i = 0; i < f; i++) {
					Result.FO[i] = this->FO[i];
					Result.PArray[i] = this->PArray[i];
				}
			if (e > f + 1)
				for (int i = f; i < e - 1; i++) {
					Result.FO[i] = this->FO[i + 1];
					Result.PArray[i] = this->PArray[i + 1];
				}
			if (this->FO.size() > e + 1)
				for (int i = e - 1; i < this->FO.size() - 2; i++) {
					Result.FO[i] = this->FO[i + 2];
					Result.PArray[i] = this->PArray[i + 2];
				}
		}
	}
	Result.Targets.resize(this->Targets.size());
	for (int i = 0; i < Result.Targets.size(); i++)
		Result.Targets[i] = this->Targets[i];
	return Result;
}

kGraph kGraph::InducedKGraph(std::vector<int> spisok, int Number) {
	//{Восстанавливает подграф K-графа G, в который входят  ребра с номером Number в списке spisok}
	std::vector<int> S(spisok.size());
	int k = 0;
	for (int i = 1; i < spisok.size(); i++)
		if ((spisok[i] == Number) || (!spisok[i])) {
			k++;
			S[i] = k;
		}
		else
			S[i] = 0;
	std::vector<int> KAO(k + 1), FO(this->FO.size());
	std::vector<double>PArray(this->FO.size());
	int l = 0;
	for (int i = 1; i < spisok.size(); i++)
		if (S[i] != 0) {
			KAO[S[i]] = KAO[S[i] - 1];
			for (int j = this->KAO[i - 1]; j < this->KAO[i]; j++)
				if ((S[this->FO[j]]) && this->Boolka(spisok, Number, i, j)) {
					KAO[S[i]]++;
					FO[l] = S[this->FO[j]];
					PArray[l] = this->PArray[j];
					l++;
				}
		}
	FO.resize(l);
	PArray.resize(l);
	std::vector<int> Targets(KAO.size());
	kGraph Result(KAO, FO, PArray, Targets);
	for (int i = 1; i < this->Targets.size(); i++)
		if (S[i] > 0)
			Result.Targets[S[i]] = this->Targets[i];

	return Result;
}

bool kGraph::KComponent() {
	//{Выделяет компоненту связности графа G, содержащую все целевые вершины.
	//Если это возможно, то компонента сохраняется как G, а результат - истина.
	//В противном случае граф не изменяется, а результат - ложь.
	//Переформирование графа G происходит лишь в том случае, когда он несвязен.
	//Если граф не содержит целевых вершин, то результат - ложь.
	//Применение функции к пустому графу некорректно}
	int sum1 = 0;
	for (int i = 1; i < this->Targets.size(); i++)
		if (this->Targets[i] == 1)
			sum1++;
	if (sum1 == 1) {
		// как ищется компонента связности?
		// идем по массиву Таргетс со 2го до последнего элемента
		// и если там есть всего 1 целевая вершина, то возвращаем безвершинный граф
		// а если первый элемент таргетс 1
		this->KAO.resize(2);
		this->KAO[0] = 0;
		this->KAO[1] = 0;
		this->FO.resize(0);
		this->PArray.resize(0);
		this->Targets.resize(2);
		this->Targets[0] = 0;
		this->Targets[1] = 0;
		return true;
	}
	else {
		std::vector<int> A, B, Spot;
		int sum = 1;
		int SumAll = 1;
		A.resize(1);
		for (int i = 1; i < this->Targets.size(); i++)
			if (this->Targets[i] == 1) {
				A[0] = i;
				break;
			}
		int l = A.size();
		Spot.resize(this->KAO.size());
		for (int i = 1; i < this->KAO.size(); i++)
			Spot[i] = 2;
		Spot[A[0]] = 1;
		while (l > 0) {
			B.resize(0);
			for (int i : A)
				for (int j = this->KAO[i - 1]; j < this->KAO[i]; j++)
					if (Spot[this->FO[j]] == 2) {
						B.resize(B.size() + 1);
						B[B.size() - 1] = this->FO[j];
						Spot[this->FO[j]] = 1;
						if (this->Targets[this->FO[j]] == 1)
							sum++;
						SumAll++;
					}

			A.resize(B.size());
			l = A.size();
			if (l > 0)
				for (int i = 0; i < l; i++)
					A[i] = B[i];
		}
		bool Result = (sum == sum1);
		if (Result && (SumAll < (this->KAO.size() - 1)))
			*this = this->InducedKGraph(Spot, 1); //requires attention!
		return Result;
	}
}

bool kGraph::Boolka(std::vector<int> spisok, int Number, int i, int j) {
	if (Number == 1)
		return true;
	//    if (!(spisok[i] || spisok[this->FO[j]]))
	//        return false;
	//    return true;
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
	//{Восстанавливает блок K - графа G с номером Number,
	// используя spisok предварительно сформированный процедурой DecopmoseOnBlocksK
	// Эта же процедура уже пополнила спсиок цел.в - н необходимыми точками сочленения}
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
			KAO.push_back(KAO[KAO.size() - 1] + S[i]);
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
	fout << this->KAO[this->KAO.size() - 1] << std::endl;
	for (int i = 0; i < this->FO.size() - 1; i++) {
		fout << FO[i] << ",";
	}
	fout << this->FO[this->FO.size() - 1] << std::endl;
	for (int i = 0; i < this->Targets.size() - 1; i++) {
		fout << Targets[i] << ",";
	}
	fout << this->Targets[this->Targets.size() - 1] << std::endl;
	fout << this->PArray[0];
}