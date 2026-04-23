#pragma once

// Фасад: sqrt(dx^2 + dy^2 + dz^2) через выбранные алгоритмы умножения и суммирования.

#include "types.h"
#include <string_view>

namespace euclidean {

    [[nodiscard]] double euclidean_distance(const Point3D& p1,
        const Point3D& p2,
        MultiplyMethod  mul_method,
        SumMethod       sum_method) noexcept;

    // Запускает все 15 комбинаций (3 mul x 5 sum), результат в data[mul][sum].
    struct AllResults {
        double data[3][5]{};
    };

    [[nodiscard]] AllResults euclidean_all(const Point3D& p1, const Point3D& p2) noexcept;

    [[nodiscard]] std::string_view method_name(MultiplyMethod m) noexcept;
    [[nodiscard]] std::string_view method_name(SumMethod s) noexcept;

} // namespace euclidean