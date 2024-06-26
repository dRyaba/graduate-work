#pragma once
#include <vector>

extern int Nconst;
extern std::vector<int> cutPoints;

struct Graph {
	std::vector<int> KAO, FO;
	std::vector<double> PArray;

	Graph() = default;

	Graph(std::vector<int> KAO, std::vector<int> FO, std::vector<double> PArray) {
		this->KAO = std::move(KAO);
		this->FO = std::move(FO);
		this->PArray = std::move(PArray);
	}
	//    ~Graph()= default;

	void CutPointsSearch(int v, int p);
	// ������� ����� ���������� � ���������� �� � cutPoints

	int SearchEdge(const int i, const int j);
	//��������� ����� ����� �� i � j � ������� FO

	int DistanceDijkstra(const int x, const int y);

	bool CheckEdge(int i, int j);

	Graph DeleteEdge(int u, int v);

	std::vector<int> CutDecompose(const std::vector<int> V);

	//�� ���� ���������� CutSearch �������� ������ � ������� ��� ������ ������� 
	// ������� � ����� ���������� ��������� ��� ���������
	//��� ��� ����������� ��� ���������� ����������
	//���� ��� �������  reliability2vert � ������� ������� �� ����� ���������� � �� ����� ���������� �� ������
	//� � ������ �� ���� ���������� ����� ������������� ������� ������ �� ��������������� ���������� ���������
	//��� ��� �������?
	//��� ������ ���� ����� ������������� ��������� ����� �� ��� � ������� ���������� ���������
	//������-�� �������� � �� �����, �� ����� ������������ ��������� �� ����� ����, ��� ��� ����� ��������� ����� ����� ����������,
	//������� ������� ���������� � ��� ��� ��������  
	std::vector<int> CutSearch();

	std::vector<int> CutDecomposeOnTwo();

	Graph ChangVertex(int u, int v);
	//������ � ����� ������� u � v �������(����������������)
};