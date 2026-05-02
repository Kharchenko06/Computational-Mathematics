#pragma once

// ULP (Unit in the Last Place) по Голдбергу: ulp(x) = 2^(e-p+1),
// e — показатель степени x, p=53 для double (Handbook §2.6).

#include <cstdint>

namespace euclidean {

    // ulp(0) не определен — возвращаем denorm_min.
    // Для NaN/Inf возвращаем сам аргумент.
    [[nodiscard]] double ulp(double x) noexcept;

    // Знаковое расстояние от a до b в единицах ULP.
    // ulp_distance(+0, -0) == 0
    // ulp_distance(NaN, any) == INT64_MAX
    [[nodiscard]] std::int64_t ulp_distance(double a, double b) noexcept;

    [[nodiscard]] std::int64_t ulp_error(double computed, double reference) noexcept;

    // Сырой битовый шаблон IEEE 754 — удобно при отладке округлений.
    [[nodiscard]] std::uint64_t bits(double x) noexcept;

    [[nodiscard]] double next_up(double x) noexcept;
    [[nodiscard]] double next_down(double x) noexcept;

}