#pragma once

// Коррекция sqrt для CompensatedSum с учетом ошибки

#include "types.h"

namespace euclidean {

    [[nodiscard]] double sqrt_refined(const CompensatedSum& cs) noexcept;

    // Базовый вариант без учета ошибки
    [[nodiscard]] double sqrt_naive(const CompensatedSum& cs) noexcept;

    // Сумма элементов error кроме последнего
    [[nodiscard]] double collect_error(const CompensatedSum& cs) noexcept;

}