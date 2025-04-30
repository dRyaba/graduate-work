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
	//{Расчет вероятности связности с ограничением на диаметр d методом ветвления.
	// Используется выделение комп-ты связности с целевыми в-ми, разделение ветвей,
	// встроенная ф-я проверки на расстояния}

	void ReliabilityDiamConstr2Vert(int x, int y, int d);
	//{Расчет в-ти св-ти двух вершин x,y с огр-м на диаметр d методом ветвления.
	// Используется разделение ветвей, встроенная ф-я проверки на расстояния;
	// не используется выделение компонентв с двумя в-ми (так быстрее)}
	void ReliabilityDiamConstr2VertDecomposeSimpleFacto(int x, int y, const int& UpperBound);

	void ReliabilityDiamConstr2VertM(int x, int y,/* const int& LowerBound,*/ const int& UpperBound);

	void ReliabilityDiamConstr2Vert2Blocks(int x, int y, const int& LowerBound, const int& UpperBound);

	void ReliabilityDiamConstr2VertMDecompose(int x, int y, const int& UpperBound);
	
	void ReliabilityDiamConstr2VertMDecomposeParallel(int x, int y, const int& UpperBound);

	bool CheckDistanceFloyd(const int d);

	void SearchVertice(std::vector<int>& DfNumber, std::vector<int>& TreeEdges, std::vector<int>& Low, std::vector<int>& Stek, int& r, int& l, int v, std::vector<int>& DOB);

	std::vector<int> DecomposeOnBlocks();
		//Result[i] - номер блока в котрый попало ребро с номером i(номер ребра по массиву FO)
		//Result[Result.size() - 1] - кол - во блоков

	std::vector<int> DecomposeOnBlocksK();
	//{Разложение K - графа на блоки.
	//Result[i] - номер блока в котрый попало ребро с номером i(номер ребра по массиву FO)
	//Result[length(Result) - 1] - кол - во блоков, которые содержат цел.в - ны или являются связующими
	//Ребра, входящие в блоки, которые не содержат цел.в - н и не являются связующими входят блок с номером 0
	//Точки сочленения в связующих блоках заносятся в цлевые вершины графа G
	//Если в графе одна цел.в - на, то Result[i] = 0 для всех i
	//Если в графе нет ребер, то Result состоит из одного элемента, Result[0] = 0
	//Применение функции некорректно к несвязому графу и графу без целевых вершин}

	kGraph RestoreBlockK(int Number, const std::vector<int>& spisok);
	//{Восстанавливает блок K - графа G с номером Number,
	// используя spisok предварительно сформированный процедурой DecopmoseOnBlocksK
	// Эта же процедура уже пополнила спсиок цел.в - н необходимыми точками сочленения}

	kGraph RestoreBlock(int Number, const std::vector<int>& spisok);

	bool ConnectivityWithoutBlock(int Block, std::vector<int> S);

	kGraph DeleteEdgeK(const int u, const int v);
	//    Удаляет из графа ребро(u, v)

	kGraph InducedKGraph(std::vector<int> spisok, int Number);
	//{Восстанавливает подграф K-графа G, в который входят  ребра с номером Number в списке spisok}

	bool KComponent();
	//{Выделяет компоненту связности графа G, содержащую все целевые вершины.
	//Если это возможно, то компонента сохраняется как G, а результат - истина.
	//В противном случае граф не изменяется, а результат - ложь.
	//Переформирование графа G происходит лишь в том случае, когда он несвязен.
	//Если граф не содержит целевых вершин, то результат - ложь.
	//Применение функции к пустому графу некорректно}

	bool Boolka(std::vector<int> spisok, int Number, int i, int j);

	void kGraphFileOutput(std::ofstream& fout);

	void ChangeVertex(int u, int v);
};

void Factoring(kGraph G, const int variant, const int d, double Reliab);
//результат лежит в sumReliab[0]
//Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра

void Factoring2Vert(kGraph G, const int x, const int y, const int variant, const int d, double Reliab);
//результат лежит в globsumReliab
//Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра

void Factoring2VertM(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound);
//Заполняет вектор надежностей с соответствующим ограничением на диаметр
//Ветвление, variant=0 - после удаления, variant=1 - после обнадеживания ребра

void Factoring2VertMParallel(kGraph G, const int x, const int y, const int variant, const int d, double Reliab, int LowerBound, int UpperBound);

void GraphMerging(int k);
//Считывает из файлов и объединяет два графа по первым k вершинам с помощью функции UnionGraphs

kGraph UnionGraphs(kGraph G1, kGraph G2, int k);
//{Формирует граф, полученный объединением G1 и G2 по их первым k вершинам}

kGraph kGraphFileInput(std::ifstream& fin);