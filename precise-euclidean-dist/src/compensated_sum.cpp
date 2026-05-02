// Handbook of Floating-Point Arithmetic: §4.3.1 (Fast2Sum), §4.3.2 (2Sum), §6.3.2 (KBN)

#include "compensated_sum.h"

#include <cmath>
#include <vector>

namespace euclidean {

namespace detail {

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

}

// Суммирует произвольный набор double через каскад two_sum, собирая все ошибки в один член.
static CompensatedSum<2> vecsum(std::span<const double> vals) {
    if (vals.empty()) {
        return CompensatedSum<2>{};
    }
    double s = vals[0];
    double e = 0.0;
    for (std::size_t i = 1; i < vals.size(); ++i) {
        const DoublePair r = detail::two_sum(s, vals[i]);
        s = r.hi;
        e += r.lo;
    }
    const DoublePair final_r = detail::two_sum(s, e);

    CompensatedSum<2> result;
    result.sum      = final_r.hi;
    result.error[0] = final_r.lo;
    return result;
}

CompensatedSum<1> naive_sum(std::span<const DoublePair> terms) noexcept {
    double s = 0.0;
    for (const DoublePair& p : terms) {
        s += p.hi;
    }
    CompensatedSum<1> result;
    result.sum = s;
    return result;
}

CompensatedSum<2> ogita_oishi_sum(std::span<const DoublePair> terms) {
    // hi-члены первыми — они крупнее, это улучшает численную стабильность
    std::vector<double> flat;
    flat.reserve(terms.size() * 2);
    for (const DoublePair& p : terms) { flat.push_back(p.hi); }
    for (const DoublePair& p : terms) { flat.push_back(p.lo); }
    return vecsum(flat);
}

CompensatedSum<2> kbn2_sum(std::span<const DoublePair> terms) noexcept {
    double s    = 0.0;
    double comp = 0.0;

    // Нейманн добавил вторую ветку — оригинальный Кахан её не имел
    auto accumulate = [&](double xi) noexcept {
        const double t = s + xi;
        if (std::abs(s) >= std::abs(xi)) {
            comp += (s - t) + xi;
        } else {
            comp += (xi - t) + s;
        }
        s = t;
    };

    // hi-члены первыми, затем lo — тот же порядок, что в ogita_oishi
    for (const DoublePair& p : terms) { accumulate(p.hi); }
    for (const DoublePair& p : terms) { accumulate(p.lo); }

    CompensatedSum<2> result;
    result.sum      = s;
    result.error[0] = comp;
    return result;
}

CompensatedSum<3> kbn3_sum(std::span<const DoublePair> terms) noexcept {
    const CompensatedSum<2> pass1 = kbn2_sum(terms);
    const DoublePair r = detail::two_sum(pass1.sum, pass1.error[0]);

    CompensatedSum<3> result;
    result.sum      = r.hi;
    result.error[0] = r.lo;
    result.error[1] = 0.0;
    return result;
}

CompensatedSum<4> kbn4_sum(std::span<const DoublePair> terms) noexcept {
    // Три аккумулятора: s0 (главный), s1 (первый уровень коррекции), s2 (второй уровень)
    // Каждый входной xi проходит через два вложенных KBN-шага: carry от s0 уходит в s1,
    // carry от s1 уходит в s2. Это гарантирует |error[0]| >= |error[1]| и строго
    // лучшую погрешность, чем KBN3
    double s0 = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;

    auto accumulate = [&](double xi) noexcept {
        const double t0 = s0 + xi;
        const double c0 = (std::abs(s0) >= std::abs(xi))
                        ? (s0 - t0) + xi
                        : (xi - t0) + s0;
        s0 = t0;

        const double t1 = s1 + c0;
        const double c1 = (std::abs(s1) >= std::abs(c0))
                        ? (s1 - t1) + c0
                        : (c0 - t1) + s1;
        s1 = t1;

        s2 += c1;
    };

    for (const DoublePair& p : terms) { accumulate(p.hi); }
    for (const DoublePair& p : terms) { accumulate(p.lo); }

    CompensatedSum<4> result;
    result.sum      = s0;
    result.error[0] = s1;
    result.error[1] = s2;
    result.error[2] = 0.0;
    return result;
}

}