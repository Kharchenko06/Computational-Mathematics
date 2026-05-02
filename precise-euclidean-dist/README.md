# Высокоточное евклидово расстояние (IEEE 754 double)

Реализация `sqrt(sum((b[i]-a[i])^2))` с точностью до 1 ULP в `R^n`.

Источник алгоритмов: Muller et al., *Handbook of Floating-Point Arithmetic* (далее — Handbook).

---

## Структура проекта

```
precise-euclidean-dist/
├── README.md
├── Makefile
├── main.cpp                       CLI-интерфейс
├── include/
│   ├── types.h                    DoublePair, CompensatedSum, перечисления методов
│   ├── ulp_utils.h                ULP-метрики
│   ├── multiply.h                 naive_square, fma_square, ozaki_square, veltkamp_split
│   ├── compensated_sum.h          fast2sum, two_sum, 5 алгоритмов суммирования
│   ├── sqrt_refined.h             sqrt_naive, sqrt_refined (Ньютон + FMA-остаток)
│   ├── euclidean.h                euclidean_distance, euclidean_all
│   └── gmp_reference.h            
└── src/
│   ├── ulp_utils.cpp
│   ├── multiply.cpp
│   ├── compensated_sum.cpp
│   ├── sqrt_refined.cpp
│   ├── euclidean.cpp
│   └── gmp_reference.cpp
└── tests/
    ├── test_multiply.cpp          EFT-свойства алгоритмов умножения
    ├── test_compensated_sum.cpp   корректность пяти алгоритмов суммирования
    ├── test_sqrt_refined.cpp      сравнение sqrt_naive vs sqrt_refined
    ├── test_permutation.cpp       инвариантность к перестановке dx/dy/dz
    └── test_euclidean_full.cpp    сводная таблица: 15 комбинаций vs GMP
```

Все файлы в корне проекта; Makefile управляет сборкой без CMake.

---

## Сборка и запуск

```bash
# Собрать всё (ожидается 0 errors, 0 warnings)
make clean && make

# Запустить все тесты последовательно
make run

# CLI-демо (3-4-5 и катастрофическая потеря значащих бит)
make demo

# CLI-демо из файла (R^5)
make demo-file
```

### Требования

- GCC с поддержкой C++20 и аппаратного FMA (`-march=native`)
- Библиотека GMP (`-lgmp`)

---

## CLI-использование

```bash
# R^3: два вектора из командной строки
./build/euclidean_distance 0 0 0 3 4 0

# R^n: чтение из файла
./build/euclidean_distance --file input.txt
```

Формат файла:
```
<n>
<a1> <a2> ... <an>
<b1> <b2> ... <bn>
```

Интерпретация вывода: `*` — точное совпадение с GMP (0 ULP), `+` — 1 ULP, остальное не проходит требование.

---

## Алгоритмы

### Умножение — вычисление `d²` с точным остатком ошибки

Все три метода возвращают `DoublePair (hi, lo)` такой, что `hi + lo = d*d` точно, `hi = fl(d*d)`.

| Метод | Алгоритм | Источник |
|-------|----------|----------|
| `Naive` | `hi = d*d`, `lo = 0` | — |
| `FMA` | `hi = d*d`, `lo = fma(d, d, -hi)` | Handbook §5.1, Alg. 5.1 |
| `Ozaki` | Декер через сплит Велткампа | Handbook §4.4, Alg. 4.6–4.7 |

`FMA` и `Ozaki` производят побитово идентичный `lo` при наличии аппаратного FMA — это проверяется в `test_multiply`.

Veltkamp split: `C = 2^27 + 1 = 134217729` (для `p = 53`, `s = ceil(53/2) = 27`). После сплита оба полуслова вмещаются в 26 бит, поэтому их произведения представимы в double точно. Требует `-fno-unsafe-math-optimizations`.

### Суммирование — сложение `d²` компонент

| Метод | Описание | Источник |
|-------|----------|----------|
| `Naive` | Суммирует только `hi`-члены | — |
| `OgitaOishi` | VecSum: суммирует все `2n` компонент через каскад `two_sum` | Handbook §6.3 |
| `KBN2` | Kahan–Babuška–Neumaier: исправляет случай `|xi| > |s|` | Handbook §6.3.2, Alg. 6.6 |
| `KBN3` | KBN2 + один проход `two_sum` по накопленной компенсации | Handbook §6.3.2 |
| `KBN4` | KBN2 + два прохода `two_sum` | Handbook §6.3.2 |

`Naive` дает ~10 ULP на катастрофической потери точности. `KBN3`/`KBN4` снижают до 0–1 ULP.

### Извлечение корня — Newton step с FMA-точным остатком

Handbook §5.4.2. Пусть `r0 = sqrt(S_hi)`, где `S_hi = cs.sum`.

```
p        = r0 * r0                         // округленное r0^2
delta    = fma(r0, r0, -p)                 // точная ошибка умножения
residual = (S_hi - p) + e_total - delta    // точный остаток S - r0^2
result   = r0 + residual / (2 * r0)        // шаг Ньютона
```

`e_total = cs.error[0] + ... + cs.error[order-2]` — накопленная компенсация суммирования.

При наличии компенсированной суммы (`KBN3`/`KBN4`) гарантирует ≤ 1 ULP.

---

## Тесты

### `test_multiply`
- Свойство EFT: `hi == fl(x*x)` для всех трех методов
- Сплит Велткампа: `hi + lo == x` точно
- Побитовое совпадение `FMA.lo == Ozaki.lo`

### `test_compensated_sum`
- Примитив `two_sum`: `hi == fl(a+b)`
- 5 алгоритмов на малых значениях (`~1e-7`)
- Катастрофическое сокращение значащих бит: `A = (1e15,…)`, `B = (1e15+1,…)`, ожидается `sum = 3.0`
- `R^100`: 100 единичных компонент, ожидается `sum = 100.0`

### `test_sqrt_refined`
- 8 тестовых случаев: куб, 3-4-5, малые/большие значения, смешанные порядки, `R^5`
- Сравнение `sqrt_naive` vs `sqrt_refined` по ULP относительно GMP

### `test_permutation`
- 4 пары точек, все 6 перестановок координат
- Для каждой из 15 комбинаций (3 mul × 5 sum): инвариантность (YES/NO) и MaxULP

### `test_euclidean_full`
- 16 пар точек: `R^2`, `R^3`, `R^5`, катастрофическая потеря значащих бит, экстремальные значения
- Для каждой пары — таблица 15 комбинаций с ULP; проход критерия: `best_ulp <= 1`
- Ожидаемый результат: **16/16 PASS**

---

## Флаги компиляции

| Флаг | Назначение |
|------|-----------|
| `-std=c++20` | `std::span`, `std::bit_cast` |
| `-fno-unsafe-math-optimizations` | Запрет перестановок, обязательный для Veltkamp split |
| `-ffp-contract=on` | Явное управление слиянием в FMA |
| `-march=native` | Включение аппаратного FMA |
| `-O2` | Стандартная оптимизация без нарушения семантики |

---

## Источники

**Handbook of Floating-Point Arithmetic** (Muller et al., 2nd ed.):

| Раздел | Содержание |
|--------|-----------|
| §2.6 | Определение ULP |
| §4.3.1 | Fast2Sum |
| §4.3.2 | 2Sum |
| §4.4 | Сплит Велткампа, произведение Декера |
| §5.1 | 2MultFMA |
| §5.4.2 | Уточнение sqrt шагом Ньютона |
| §6.3 | Компенсированное суммирование (VecSum) |
| §6.3.2 | KBN |

**Стандарты:** IEEE 754-2019, ISO/IEC 14882:2020 (C++20).

---

## Архитектурные решения

**Почему GMP?** 256-битная точность дает ошибку ~0.01 ULP относительно точного результата — достаточно для надежной классификации `{0 ULP, 1 ULP, > 1 ULP}`.

**Когда нужен `sqrt_refined`?** При катастрофической потери точности `std::sqrt(cs.sum)` дает до ~5–10 ULP; шаг Ньютона с FMA-точным остатком возвращает ошибку в пределы 1 ULP.

**Инвариантность перестановок.** Сложение FP-чисел не ассоциативно, поэтому перестановка компонент меняет результат. `KBN3`/`KBN4` + `sqrt_refined` дают одинаковый результат для всех 6 перестановок в `R^3` в большинстве тестовых случаев.