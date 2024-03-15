#pragma once

#include "GraphOperations.h"

extern std::ofstream output("C:/Users/User/source/repos/graduate work/graduate work/output.txt", std::ofstream::trunc);
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

	void ReliabilityDiamConstr2VertM(int x, int y, int d);

	void ReliabilityDiamConstr2VertDecompose(int x, int y, int d);

	
	bool CheckDistanceFloyd(const int d);

	void SearchVertice(std::vector<int>& DfNumber, std::vector<int>& TreeEdges, std::vector<int>& Low, std::vector<int>& Stek, int& r, int& l, int v, std::vector<int>& DOB);

	std::vector<int> DecomposeOnBlocks();
		//Result[i] - ����� ����� � ������ ������ ����� � ������� i(����� ����� �� ������� FO)
		//Result[Result.size() - 1] - ��� - �� ������

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
};

void Factoring(kGraph G, const int variant, const int d, double Reliab);
//��������� ����� � sumReliab[0]
//���������, variant=0 - ����� ��������, variant=1 - ����� ������������� �����

void Factoring2Vert(kGraph G, const int x, const int y, const int variant, const int d, double Reliab);
//��������� ����� � sumReliab[0]
//���������, variant=0 - ����� ��������, variant=1 - ����� ������������� �����

void Factoring2VertM(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound);
//��������� ������ ����������� � ��������������� ������������ �� �������
//���������, variant=0 - ����� ��������, variant=1 - ����� ������������� �����
