# Система импорта данных для graduate-work

## Обзор

Новая система импорта данных заменяет старую логику с хардкодом путей на централизованную и гибкую систему управления данными.

## Основные улучшения

### ❌ **Старая система (проблемы):**
- Хардкод путей: `"C:/Users/User/source/repos/graduate work/graduate work/"`
- Отсутствие обработки ошибок
- Дублирование кода загрузки
- Нет валидации данных
- Смешанная логика в разных файлах

### ✅ **Новая система (преимущества):**
- Централизованное управление путями
- Полная обработка ошибок с исключениями
- Валидация загруженных данных
- Единый интерфейс для всех типов файлов
- Кроссплатформенность
- Умные указатели для управления памятью

## Структура файлов

```
graduate-work/
├── DataImporter.h          # Заголовочный файл системы импорта
├── DataImporter.cpp        # Реализация системы импорта
├── graphs_data/            # Папка с данными (новая структура)
│   ├── *.kao.txt          # Графы в KAO формате
│   ├── *.edgelist.txt     # Графы в Edge List формате
│   └── GraphsToMerge.txt  # Файл с графами для объединения
├── graduate work.cpp       # Обновленный основной файл
├── kGraphOperations.cpp    # Обновленные операции с графами
└── example_usage.cpp       # Пример использования
```

## Использование

### Базовое использование

```cpp
#include "DataImporter.h"

// Создаем импортер
DataImporter importer("graphs_data/");

// Загружаем граф из KAO формата
auto graph = importer.loadKAOGraph("K4_kao.txt");

// Загружаем граф из Edge List формата
auto graph2 = importer.loadEdgeListGraph("Geant2009_edgelist.txt", 0.9);

// Загружаем графы для объединения
auto graphs = importer.loadGraphsToMerge();
```

### Обработка ошибок

```cpp
try {
    auto graph = importer.loadKAOGraph("nonexistent.txt");
} catch (const std::exception& e) {
    std::cerr << "Ошибка загрузки: " << e.what() << std::endl;
}
```

### Конвертация форматов

```cpp
// Edge List → KAO
importer.convertEdgeListToKAO("input.edgelist", "output.kao", 0.9);

// KAO → Edge List
importer.convertKAOToEdgeList("input.kao", "output.edgelist");
```

### Получение информации о файлах

```cpp
// Список доступных файлов
auto files = importer.getAvailableFiles();

// Проверка существования файла
bool exists = importer.fileExists("graph.txt");

// Полный путь к файлу
std::string path = importer.getFullPath("graph.txt");
```

## API Reference

### DataImporter

#### Конструктор
```cpp
DataImporter(const std::string& baseDataPath = "graphs_data/")
```

#### Основные методы

| Метод | Описание |
|-------|----------|
| `loadKAOGraph(filename)` | Загружает граф из KAO формата |
| `loadEdgeListGraph(filename, reliability)` | Загружает граф из Edge List формата |
| `loadGraphsToMerge()` | Загружает графы для объединения |
| `fileExists(filename)` | Проверяет существование файла |
| `getAvailableFiles()` | Возвращает список доступных файлов |
| `convertEdgeListToKAO(in, out, rel)` | Конвертирует Edge List в KAO |
| `convertKAOToEdgeList(in, out)` | Конвертирует KAO в Edge List |

## Форматы данных

### KAO формат
```
0,2,5,8,10,13,17,21,24,27,31,35,38,40,43,46,48
2,5,1,3,6,2,4,7,3,8,1,6,9,2,5,7,10,3,6,8,11,4,7,12,5,10,13,6,9,11,14,7,10,12,15,8,11,16,9,14,10,13,15,11,14,16,12,15
0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0
0.9
// Комментарий
```

### Edge List формат
```
1 -- 2
2 -- 3
3 -- 4
// Комментарий
```

## Миграция с старой системы

### Было:
```cpp
std::ifstream fin("C:/Users/User/source/repos/graduate work/graduate work/input.txt");
kGraph G = kGraphFileInput(fin);
```

### Стало:
```cpp
DataImporter importer("graphs_data/");
auto graph = importer.loadKAOGraph("input.txt");
kGraph G = *graph;
```

## Примеры использования

См. файл `example_usage.cpp` для полных примеров использования новой системы импорта.

## Совместимость

- C++11 или выше
- Поддержка std::filesystem (C++17) для автоматического сканирования папок
- Fallback для старых версий C++

## Обработка ошибок

Все методы выбрасывают `std::runtime_error` с описательным сообщением при ошибках:
- Файл не найден
- Неверный формат данных
- Поврежденные данные
- Ошибки валидации

## Производительность

- Использование умных указателей для эффективного управления памятью
- Минимальное копирование данных
- Валидация только при необходимости
- Оптимизированные алгоритмы парсинга
