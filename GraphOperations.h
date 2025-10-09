#pragma once
#include <vector>
#include <fstream>

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

	int SearchEdge(const int i, const int j) const;
	//��������� ����� ����� �� i � j � ������� FO

	int DistanceDijkstra(const int x, const int y) const;

	bool CheckEdge(int i, int j) const;

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

	// ������ ������ ���� � �������������� � ������ KAO/FO
	virtual void convertEdgeListToKAOFO(const std::string& inputPath,
		const std::string& outputPath,
		double reliability);

	// �������� �������������� KAO/FO -> ������ ����
	virtual void convertKAOFOToEdgeList(const std::string& inputPath,
		const std::string& outputPath);
};
