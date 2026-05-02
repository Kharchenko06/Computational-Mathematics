// Handbook of Floating-Point Arithmetic: §5.2 (Теорема 13), §5.4.2

#include "sqrt_refined.h"
#include "compensated_sum.h"

#include <cmath>

namespace euclidean {

double collect_error(const CompensatedSum& cs) noexcept {
    if (cs.order <= 1) {
        return 0.0;
    }
    // error[0] - наибольший член, error[order-2] - наименьший
    // Суммируем от наименьшего к большему, чтобы снизить ошибку округления накопления
    double e = 0.0;
    for (int i = cs.order - 2, i >= 0; --i) {
        e += cs.error[i];
    }
    return e;
}

double sqrt_naive(const CompensatedSum& cs) noexcept {
    return std::sqrt(cs.sum);
}

double sqrt_refined(const CompensatedSum& cs) noexcept {
    if (cs.sum == 0.0) {
        return 0.0;
    }

    const double r0 = std::sqrt(cs.sum);

    const double e_total = collect_error(cs);
    if (e_total == 0.0) {
        return r0;
    }

    // вычисление точного остатка через fma
    const double p = r0 * r0;
    const double delta = std::fma(r0, r0, -p);

    // остаток S - r0^2 с учетом накопленной ошибки
    const double residual = (cs.sum - p) + e_total - delta;

    // множитель 0.5 выбран из формулы Ньютона
    const double correction = (residual * 0.5) / r0;
    
    return r0 + correction;
}

}