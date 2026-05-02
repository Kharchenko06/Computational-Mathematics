#pragma once

// Коррекция sqrt для CompensatedSum<Order> с учетом ошибки суммирования.
// Handbook of Floating-Point Arithmetic: §5.4.2

#include "types.h"

#include <cmath>

namespace euclidean {

    namespace detail {

        // Сумма всех членов ошибки, от наименьшего (error[Order-2]) к наибольшему (error[0]).
        // Порядок суммирования снижает накопленную ошибку округления.
        template<int Order>
        [[nodiscard]] double collect_error(const CompensatedSum<Order>& cs) noexcept {
            double e = 0.0;
            for (int i = Order - 2; i >= 0; --i) {
                e += cs.error[i];
            }
            return e;
        }

    }

    // Базовый вариант без учета ошибки суммирования.
    template<int Order>
    [[nodiscard]] double sqrt_naive(const CompensatedSum<Order>& cs) noexcept {
        return std::sqrt(cs.sum);
    }

    // Уточненный sqrt: шаг Ньютона с FMA-точным остатком (Handbook §5.4.2).
    // r0 = sqrt(cs.sum); result = r0 + (cs - r0^2) / (2 * r0)
    // Остаток cs - r0^2 вычисляется через fma без потери точности.
    template<int Order>
    [[nodiscard]] double sqrt_refined(const CompensatedSum<Order>& cs) noexcept {
        if (cs.sum == 0.0) {
            return 0.0;
        }

        const double r0 = std::sqrt(cs.sum);

        const double e_total = detail::collect_error(cs);
        if (e_total == 0.0) {
            return r0;
        }

        // точный остаток r0^2 через fma
        const double p = r0 * r0;
        const double delta = std::fma(r0, r0, -p);

        // полный остаток S - r0^2 с учетом накопленной компенсации
        const double residual = (cs.sum - p) + e_total - delta;

        // множитель 0.5 из формулы Ньютона: r1 = r0 + residual / (2 * r0)
        return r0 + (residual * 0.5) / r0;
    }

}