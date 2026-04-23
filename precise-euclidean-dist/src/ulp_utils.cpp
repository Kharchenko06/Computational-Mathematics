// Реализация утилит IEEE 754 ULP (Handbook §2.6)

#include "ulp_utils.h"

#include <bit>
#include <cmath>
#include <cstdlib>
#include <limits>

namespace euclidean {

// IEEE 754 binary64 спроектирован так, что целочисленный порядок битовых шаблонов
// совпадает с порядком значений для чисел одного знака. Для отрицательных чисел
// битовое дополнение (~u) инвертирует порядок, выравнивая его с положительными.
// Это позволяет измерять расстояние в ULP через обычное вычитание int64_t.
static std::int64_t to_ordered_int(double x) noexcept {
    const std::uint64_t u = std::bit_cast<std::uint64_t>(x);
    if (u >> 63) {
        return static_cast<std::int64_t>(~u);
    }
    return static_cast<std::int64_t>(u ^ (std::uint64_t{1} << 63));
}

double ulp(double x) noexcept {
    if (std::isnan(x) || std::isinf(x)) {
        return x;
    }
    if (x == 0.0) {
        return std::numeric_limits<double>::denorm_min();
    }
    // frexp: x = frac * 2^exp, frac in [0.5, 1.0)
    // ulp = 2^(exp-1-52) = 2^(exp-53)
    int exp = 0;
    std::frexp(x, &exp);
    const double u = std::scalbn(1.0, exp - 53);
    return (u > 0.0) ? u : std::numeric_limits<double>::denorm_min();
}

std::int64_t ulp_distance(double a, double b) noexcept {
    if (std::isnan(a) || std::isnan(b)) {
        return std::numeric_limits<std::int64_t>::max();
    }
    // Нормализуем -0 в +0, иначе их битовые шаблоны дают ненулевое расстояние
    if (a == 0.0) { a = 0.0; }
    if (b == 0.0) { b = 0.0; }

    return to_ordered_int(b) - to_ordered_int(a);
}

std::int64_t ulp_error(double computed, double reference) noexcept {
    const std::int64_t d = ulp_distance(reference, computed);
    return (d >= 0) ? d : -d;
}

std::uint64_t bits(double x) noexcept {
    return std::bit_cast<std::uint64_t>(x);
}

double next_up(double x) noexcept {
    return std::nextafter(x, std::numeric_limits<double>::infinity());
}

double next_down(double x) noexcept {
    return std::nextafter(x, -std::numeric_limits<double>::infinity());
}

}