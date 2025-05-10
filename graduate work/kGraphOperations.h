#pragma once
#include <fstream>
#include "GraphOperations.h"

extern std::ofstream output, output1;
extern int NumberOfRec;
extern std::vector<double> sumReliab;
extern double globsumReliab;


struct kGraph : public Graph {
	std::vector<int> Targets;

	kGraph() = default;

	kGraph(std::vector<int> KAO, std::vector<int> FO, std::vector<double> PArray, std::vector<int> Targets) : Graph() {
		this->KAO = std::move(KAO);
		this->FO = std::move(FO);
		this->PArray = std::move(PArray);
		this->Targets = std::move(Targets);
	}
	//    ~kGraph()= default;

	void ReliabilityDiamConstr(int d);
	//{������ ����������� ��������� � ������������ �� ������� d ������� ���������.
	// ������������ ��������� ����-�� ��������� � �������� �-��, ���������� ������,
	// ���������� �-� �������� �� ����������}

	void ReliabilityDiamConstr2Vert(int x, int y, int d);
	//{������ �-�� ��-�� ���� ������ x,y � ���-� �� ������� d ������� ���������.
	// ������������ ���������� ������, ���������� �-� �������� �� ����������;
	// �� ������������ ��������� ���������� � ����� �-�� (��� �������)}
	void ReliabilityDiamConstr2VertDecomposeSimpleFacto(int x, int y, const int& UpperBound);

	void ReliabilityDiamConstr2VertM(int x, int y,/* const int& LowerBound,*/ const int& UpperBound);

	void ReliabilityDiamConstr2Vert2Blocks(int x, int y, const int& LowerBound, const int& UpperBound);

	void ReliabilityDiamConstr2VertMDecompose(int x, int y, const int& UpperBound);
	
	void ReliabilityDiamConstr2VertMDecomposeParallel(int x, int y, const int& UpperBound);

	bool CheckDistanceFloyd(const int d);

	void SearchVertice(std::vector<int>& DfNumber, std::vector<int>& TreeEdges, std::vector<int>& Low, std::vector<int>& Stek, int& r, int& l, int v, std::vector<int>& DOB);

	std::vector<int> DecomposeOnBlocks();
		//Result[i] - ����� ����� � ������ ������ ����� � ������� i(����� ����� �� ������� FO)
		//Result[Result.size() - 1] - ��� - �� ������

	std::vector<int> DecomposeOnBlocksK();
	//{���������� K - ����� �� �����.
	//Result[i] - ����� ����� � ������ ������ ����� � ������� i(����� ����� �� ������� FO)
	//Result[length(Result) - 1] - ��� - �� ������, ������� �������� ���.� - �� ��� �������� ����������
	//�����, �������� � �����, ������� �� �������� ���.� - � � �� �������� ���������� ������ ���� � ������� 0
	//����� ���������� � ��������� ������ ��������� � ������ ������� ����� G
	//���� � ����� ���� ���.� - ��, �� Result[i] = 0 ��� ���� i
	//���� � ����� ��� �����, �� Result ������� �� ������ ��������, Result[0] = 0
	//���������� ������� ����������� � ��������� ����� � ����� ��� ������� ������}

	kGraph RestoreBlockK(int Number, const std::vector<int>& spisok);
	//{��������������� ���� K - ����� G � ������� Number,
	// ��������� spisok �������������� �������������� ���������� DecopmoseOnBlocksK
	// ��� �� ��������� ��� ��������� ������ ���.� - � ������������ ������� ����������}

	kGraph RestoreBlock(int Number, const std::vector<int>& spisok);

	bool ConnectivityWithoutBlock(int Block, std::vector<int> S);

	kGraph DeleteEdgeK(const int u, const int v);
	//    ������� �� ����� �����(u, v)

	kGraph InducedKGraph(std::vector<int> spisok, int Number);
	//{��������������� ������� K-����� G, � ������� ������  ����� � ������� Number � ������ spisok}

	bool KComponent();
	//{�������� ���������� ��������� ����� G, ���������� ��� ������� �������.
	//���� ��� ��������, �� ���������� ����������� ��� G, � ��������� - ������.
	//� ��������� ������ ���� �� ����������, � ��������� - ����.
	//���������������� ����� G ���������� ���� � ��� ������, ����� �� ��������.
	//���� ���� �� �������� ������� ������, �� ��������� - ����.
	//���������� ������� � ������� ����� �����������}

	bool Boolka(std::vector<int> spisok, int Number, int i, int j);

	void kGraphFileOutput(std::ofstream& fout);

	void ChangeVertex(int u, int v);
};

// Вспомогательная структура для представления графа блоков
struct BlockGraphNode {
    int original_block_id; // ID из DecomposeOnBlocksK
    std::vector<int> vertices_original_ids; // Вершины исходного графа, принадлежащие этому блоку
    // Можно добавить другую мета-информацию о блоке, если нужно
};

struct BlockGraphEdge {
    int from_block_node_idx; // Индекс узла в векторе узлов графа блоков
    int to_block_node_idx;   // Индекс узла в векторе узлов графа блоков
    int articulation_point_original_id; // ID точки сочленения
};

void Factoring(kGraph G, const int variant, const int d, double Reliab);
//��������� ����� � sumReliab[0]
//���������, variant=0 - ����� ��������, variant=1 - ����� ������������� �����

void Factoring2Vert(kGraph G, const int x, const int y, const int variant, const int d, double Reliab);
//��������� ����� � globsumReliab
//���������, variant=0 - ����� ��������, variant=1 - ����� ������������� �����

void Factoring2VertM(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound);
//��������� ������ ����������� � ��������������� ������������ �� �������
//���������, variant=0 - ����� ��������, variant=1 - ����� ������������� �����

void Factoring2VertMParallel(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound);

void GraphMerging(int k);
//��������� �� ������ � ���������� ��� ����� �� ������ k �������� � ������� ������� UnionGraphs

kGraph UnionGraphs(kGraph G1, kGraph G2, int k);
//{��������� ����, ���������� ������������ G1 � G2 �� �� ������ k ��������}

kGraph kGraphFileInput(std::ifstream& fin);