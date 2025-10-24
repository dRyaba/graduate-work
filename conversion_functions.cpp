#include <iostream>
#include <string>
#include "DataImporter.h"

// Из edge‑list (строка "u -- v") → KAO/FO+Targets+PArray
void convertEdgeListToKAOFO(const std::string& inPath,
    const std::string& outPath,
    double reliability)
{
    try {
        DataImporter importer;
        importer.convertEdgeListToKAO(inPath, outPath, reliability);
        std::cout << "EdgeList -> KAOFO written to " << outPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error converting EdgeList to KAO: " << e.what() << "\n";
    }
}

// Из KAO/FO(+…) → edge‑list ("u -- v")
void convertKAOFOToEdgeList(const std::string& inPath,
    const std::string& outPath)
{
    try {
        DataImporter importer;
        importer.convertKAOToEdgeList(inPath, outPath);
        std::cout << "KAOFO → EdgeList written to " << outPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error converting KAO to EdgeList: " << e.what() << "\n";
    }
}

