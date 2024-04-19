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
	// Ќаходит точки сочленени€ и записывает их в cutPoints

	int SearchEdge(const int i, const int j);
	//вычисл€ет номер ребра из i в j в массиве FO

	int DistanceDijkstra(const int x, const int y);

	bool CheckEdge(int i, int j);

	Graph DeleteEdge(int u, int v);

	std::vector<int> CutDecompose(const std::vector<int> V);

	//то есть используем CutSearch получаем вектор в котором про каждую вершину 
	// сказано в какой компоненте св€зности она находитс€
	//как это исползовать при вычислении надежности
	//надо два запуска  reliability2vert с искомой вершины до точки сочленени€ и от точки сочленени€ до другой
	//и в каждом из этих вычислений нужно ииспользовать вершины только из соответствующей компоненты св€зности
	//как это сделать?
	//при выборе рЄбер нужно дополнительно провер€ть лежит ли оно в текущей компоненте св€зности
	//вообще-то проверка и не нужна, тк между компонентами св€зности не будет рЄбер, они все будут проходить через точку сочленени€,
	//поэтому решение оптимально и без доп проверок  
	std::vector<int> CutSearch();

	std::vector<int> CutDecomposeOnTwo();

	Graph ChangVertex(int u, int v);
	//ћен€ет в графе вершины u и v местами(перенумеровывает)
};