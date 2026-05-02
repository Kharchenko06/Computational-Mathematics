#pragma once

// Три алгоритма вычисления x*x с точным членом ошибки.
// Все три возвращают DoublePair (hi, lo): hi + lo = x*x точно, hi = fl(x*x).

#include "types.h"

namespace euclidean {

    // Базовый вариант без компенсации, lo = 0.
    [[nodiscard]] DoublePair naive_square(double x) noexcept;

    // 2MultFMA (Handbook §5.1, Alg. 5.1).
    // Требует аппаратного FMA (включается через -march=native).
    [[nodiscard]] DoublePair fma_square(double x) noexcept;

    // Произведение Деккера через сплит Велткампа (Handbook §4.4, Alg. 4.6 + 4.7).
    // Не требует FMA, ~17 FP-операций.
    [[nodiscard]] DoublePair ozaki_square(double x) noexcept;

    // Вспомогательный сплит, открыт для тестирования.
    // Разбивает x на (hi, lo): hi + lo = x точно, оба по 26 значимых бит.
    // C = 2^27 + 1 (для p=53: s = ceil(53/2) = 27).
    [[nodiscard]] DoublePair veltkamp_split(double x) noexcept;

}