#pragma once

// Евклидово расстояние в R^n: sqrt(sum((b[i]-a[i])^2))
// с выбором алгоритмов умножения и суммирования.

#include "types.h"

#include <span>
#include <string_view>

namespace euclidean {

    // Размерности a и b должны совпадать.
    [[nodiscard]] double euclidean_distance(std::span<const double> a,
        std::span<const double> b,
        MultiplyMethod          mul_method,
        SumMethod               sum_method) noexcept;

    // Результаты для всех 15 комбинаций (3 mul x 5 sum), хранятся в data[mul][sum].
    struct AllResults {
        double data[3][5]{};
    };

    [[nodiscard]] AllResults euclidean_all(std::span<const double> a,
        std::span<const double> b) noexcept;

    [[nodiscard]] std::string_view method_name(MultiplyMethod m) noexcept;
    [[nodiscard]] std::string_view method_name(SumMethod s) noexcept;

}