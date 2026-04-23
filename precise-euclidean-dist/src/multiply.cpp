// EFT-алгоритмы возведения в квадрат (Handbook §4.4, §5.1)

#include "multiply.h"

#include <cmath>

namespace euclidean {

DoublePair naive_square(double x) noexcept {
    return DoublePair{ x * x, 0.0 };
}

DoublePair fma_square(double x) noexcept {
    const double hi = x * x;
    // fma(x, x, -hi) вычисляет x*x - hi без промежуточного округления,
    // что дает точную ошибку умножения (Handbook §5.1, Alg. 5.1)
    const double lo = std::fma(x, x, -hi);
    return DoublePair{ hi, lo };
}

DoublePair veltkamp_split(double x) noexcept {
    // C = 2^27 + 1 - константа Велткампа для p=53 (s=27, C=2^s+1)
    // После сплита hi и xl вмещаются в 26 бит каждый, поэтому
    // xh*xh, xh*xl, xl*xl - все точно представимы в double (<=52 бита)
    // -fno-unsafe-math-optimizations обязателен: без него компилятор может
    // переставить операции и сломать точность сплита
    constexpr double C = 134217729.0;

    const double gamma = C * x;
    const double delta = x - gamma;
    const double hi    = gamma + delta;
    const double lo    = x - hi;
    return DoublePair{ hi, lo };
}

DoublePair ozaki_square(double x) noexcept {
    const DoublePair s = veltkamp_split(x);
    const double xh = s.hi;
    const double xl = s.lo;

    const double hi = x * x;

    // Последовательность Декера для восстановления ошибки (Handbook §4.4.2):
    // xh*xh точно представимо (26-битные множители),
    // вычитание точно по лемме Штербенца (оба операнда близки по величине)
    const double t1 = xh * xh - hi;
    const double t2 = t1 + (2.0 * xh * xl);
    const double lo = t2 + (xl * xl);

    return DoublePair{ hi, lo };
}

}