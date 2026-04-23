// Handbook of Floating-Point Arithmetic: §4.3.1 (Fast2Sum), §4.3.2 (2Sum), §6.3.2 (Kahan)

#include "compensated_sum.h"

#include <cmath>

namespace euclidean {

DoublePair fast2sum(double a, double b) noexcept {
    const double s = a + b;
    const double z = s - a;
    const double t = b - z;
    return DoublePair{ s, t };
}

DoublePair two_sum(double a, double b) noexcept {
    const double s  = a + b;
    const double a1 = s - b;
    const double b1 = s - a1;
    const double da = a - a1;
    const double db = b - b1;
    const double t  = da + db;
    return DoublePair{ s, t };
}

// Суммирует 6 значений через каскад two_sum, собирая все ошибки в один член.
static CompensatedSum vecsum6(double v0, double v1, double v2,
                               double v3, double v4, double v5) noexcept {
    double s = v0;
    double e = 0.0;

    const double vals[5] = { v1, v2, v3, v4, v5 };
    for (double vi : vals) {
        const DoublePair r = two_sum(s, vi);
        s = r.hi;
        e += r.lo;
    }
    const DoublePair final_r = two_sum(s, e);

    CompensatedSum result;
    result.sum      = final_r.hi;
    result.error[0] = final_r.lo;
    result.order    = 2;
    return result;
}

CompensatedSum naive_sum(DoublePair a, DoublePair b, DoublePair c) noexcept {
    CompensatedSum result;
    result.sum      = (a.hi + b.hi) + c.hi;
    result.error[0] = 0.0;
    result.order    = 1;
    return result;
}

CompensatedSum ogita_oishi_sum(DoublePair a, DoublePair b, DoublePair c) noexcept {
    // hi-члены первыми - они крупнее, это улучшает численную стабильность.
    return vecsum6(a.hi, b.hi, c.hi, a.lo, b.lo, c.lo);
}

CompensatedSum kbn2_sum(DoublePair a, DoublePair b, DoublePair c) noexcept {
    const double inputs[6] = { a.hi, b.hi, c.hi, a.lo, b.lo, c.lo };
    double s = inputs[0];
    double comp = 0.0;

    for (int i = 1; i < 6; ++i) {
        const double xi = inputs[i];
        const double t  = s + xi;
        // Ветка зависит от того, что больше: сумма или новый элемент.
        // Нейманн добавил вторую ветку - оригинальный Каhан её не имел.
        if (std::abs(s) >= std::abs(xi)) {
            comp += (s - t) + xi;
        } else {
            comp += (xi - t) + s;
        }
        s = t;
    }

    CompensatedSum result;
    result.sum      = s;
    result.error[0] = comp;
    result.order    = 2;
    return result;
}

CompensatedSum kbn3_sum(DoublePair a, DoublePair b, DoublePair c) noexcept {
    CompensatedSum pass1 = kbn2_sum(a, b, c);
    const DoublePair r = two_sum(pass1.sum, pass1.error[0]);

    CompensatedSum result;
    result.sum      = r.hi;
    result.error[0] = r.lo;
    result.error[1] = 0.0;
    result.order    = 3;
    return result;
}

CompensatedSum kbn4_sum(DoublePair a, DoublePair b, DoublePair c) noexcept {
    CompensatedSum pass1 = kbn2_sum(a, b, c);
    const DoublePair r2 = two_sum(pass1.sum, pass1.error[0]);
    const DoublePair r3 = two_sum(r2.hi, r2.lo);

    CompensatedSum result;
    result.sum      = r3.hi;
    result.error[0] = r3.lo;
    result.error[1] = r2.lo - (r3.hi - r2.hi);
    result.error[2] = 0.0;
    result.order    = 4;
    return result;
}

}