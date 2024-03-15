#pragma once

static int Nconst;

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

	int SearchEdge(const int i, const int j);
	//вычисляет номер ребра из i в j в массиве FO

	int DistanceDijkstra(const int x, const int y);

	bool CheckEdge(int i, int j);

	Graph DeleteEdge(int u, int v);

	std::vector<int> CutDecompose(const std::vector<int> V);

	//то есть используем CutSearch получаем вектор в котором про каждую вершину 
	// сказано в какой компоненте связности она находится
	//как это исползовать при вычислении надежности
	//надо два запуска  reliability2vert с искомой вершины до точки сочленения и от точки сочленения до другой
	//и в каждом из этих вычислений нужно ииспользовать вершины только из соответствующей компоненты связности
	//как это сделать?
	//при выборе рёбер нужно дополнительно проверять лежит ли оно в текущей компоненте связности
	//вообще-то проверка и не нужна, тк между компонентами связности не будет рёбер, они все будут проходить через точку сочленения,
	//поэтому решение оптимально и без доп проверок  
	std::vector<int> CutSearch();

	std::vector<int> CutDecomposeOnTwo();
};