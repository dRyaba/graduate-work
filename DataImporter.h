#pragma once
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "kGraphOperations.h"

/**
 * @brief Централизованная система импорта данных для проекта graduate-work
 * 
 * Этот класс предоставляет единый интерфейс для загрузки графов из различных форматов
 * и управления путями к файлам данных.
 */
class DataImporter {
public:
    /**
     * @brief Конструктор с настройкой базового пути к данным
     * @param baseDataPath Базовый путь к папке с данными (по умолчанию "graphs_data/")
     */
    explicit DataImporter(const std::string& baseDataPath = "graphs_data/");
    
    /**
     * @brief Деструктор
     */
    ~DataImporter() = default;
    
    // Запрещаем копирование и присваивание
    DataImporter(const DataImporter&) = delete;
    DataImporter& operator=(const DataImporter&) = delete;
    
    /**
     * @brief Загружает граф из KAO формата
     * @param filename Имя файла (без пути)
     * @return Умный указатель на загруженный граф
     * @throws std::runtime_error если файл не найден или поврежден
     */
    std::unique_ptr<kGraph> loadKAOGraph(const std::string& filename);
    
    /**
     * @brief Загружает граф из Edge List формата
     * @param filename Имя файла (без пути)
     * @param reliability Значение надежности для всех рёбер (по умолчанию 0.9)
     * @return Умный указатель на загруженный граф
     * @throws std::runtime_error если файл не найден или поврежден
     */
    std::unique_ptr<kGraph> loadEdgeListGraph(const std::string& filename, double reliability = 0.9);
    
    /**
     * @brief Загружает несколько графов из файла GraphsToMerge.txt
     * @return Вектор умных указателей на загруженные графы
     * @throws std::runtime_error если файл не найден или поврежден
     */
    std::vector<std::unique_ptr<kGraph>> loadGraphsToMerge();
    
    /**
     * @brief Получает полный путь к файлу данных
     * @param filename Имя файла
     * @return Полный путь к файлу
     */
    std::string getFullPath(const std::string& filename) const;
    
    /**
     * @brief Проверяет существование файла данных
     * @param filename Имя файла
     * @return true если файл существует, false иначе
     */
    bool fileExists(const std::string& filename) const;
    
    /**
     * @brief Получает список всех доступных файлов данных
     * @return Вектор имен файлов
     */
    std::vector<std::string> getAvailableFiles() const;
    
    /**
     * @brief Устанавливает новый базовый путь к данным
     * @param newBasePath Новый базовый путь
     */
    void setBasePath(const std::string& newBasePath);
    
    /**
     * @brief Получает текущий базовый путь
     * @return Текущий базовый путь
     */
    std::string getBasePath() const { return baseDataPath_; }
    
    /**
     * @brief Конвертирует Edge List в KAO формат и сохраняет
     * @param inputFilename Имя входного файла (Edge List)
     * @param outputFilename Имя выходного файла (KAO)
     * @param reliability Значение надежности
     * @throws std::runtime_error если конвертация не удалась
     */
    void convertEdgeListToKAO(const std::string& inputFilename, 
                             const std::string& outputFilename, 
                             double reliability = 0.9);
    
    /**
     * @brief Конвертирует KAO в Edge List формат и сохраняет
     * @param inputFilename Имя входного файла (KAO)
     * @param outputFilename Имя выходного файла (Edge List)
     * @throws std::runtime_error если конвертация не удалась
     */
    void convertKAOToEdgeList(const std::string& inputFilename, 
                             const std::string& outputFilename);

private:
    std::string baseDataPath_;
    
    /**
     * @brief Внутренняя функция загрузки KAO графа из потока
     * @param fin Поток для чтения
     * @return Умный указатель на загруженный граф
     * @throws std::runtime_error если данные повреждены
     */
    std::unique_ptr<kGraph> loadKAOFromStream(std::ifstream& fin);
    
    /**
     * @brief Валидирует загруженный граф
     * @param graph Указатель на граф для валидации
     * @param filename Имя файла (для сообщений об ошибках)
     * @throws std::runtime_error если граф невалиден
     */
    void validateGraph(const kGraph* graph, const std::string& filename) const;
    
    /**
     * @brief Нормализует путь (заменяет разделители, убирает лишние символы)
     * @param path Путь для нормализации
     * @return Нормализованный путь
     */
    std::string normalizePath(const std::string& path) const;
};
