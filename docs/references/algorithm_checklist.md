# Чеклист соответствия кода алгоритму (ВКР бакалавр и курсовая магистратура)

| Элемент | Описание в работах | Legacy 2f3ee60 | Current HEAD | Статус |
|---------|-------------------|----------------|--------------|--------|
| Отсечение блоков | Только блоки на (s,t)-путях | DecomposeOnBlocksK отсекает | findPathInBlockGraph | Эквивалентно |
| Порядок блоков | Цепь B1…Bn от s к t | 0..BlockNum-1 (порядок декомпозиции) | ordered_block_ids по пути | Эквивалентно |
| Modified Factoring | Algorithm 1: dist>d ∧ dist≤UB → вызов с d+1 | Factoring2VertM | calculateReliabilityMultipleDiameters | Соответствует |
| Свёртка | Algorithm 2: CBR[k]+=R_{i-1}[j]*R_i[k-j], i от BN-1 до 1 | Итеративная convolution | solveRecursiveForBlockChain | Математически эквивалентна |
| Индексация вершин | 1-based (12, 96) | 1-based | 0-based + флаг --1-based | Добавлен --1-based |
| Типы Graph/kGraph | G=(V,E,P), K⊆V терминалы | Graph + kGraph (Targets) | Graph + ReliabilityGraph (target_vertices_) | Соответствует |

## Примечания

- **Индексация:** Флаг `--1-based` в CLI: `--run --1-based file 12 96 8 3` конвертирует 12,96 в 11,95 (0-based).
- **Эталон:** Geant2004, s=12, t=96, d=8 — Decomp Modified Facto ~29 сек (ВКР).
