#include "euclidean.h"
#include "multiply.h"
#include "compensated_sum.h"
#include "sqrt_refined.h"

#include <cassert>
#include <vector>

namespace euclidean {

static DoublePair do_square(double v, MultiplyMethod m) noexcept {
    switch (m) {
        case MultiplyMethod::Naive:  return naive_square(v);
        case MultiplyMethod::FMA:    return fma_square(v);
        case MultiplyMethod::Ozaki:  return ozaki_square(v);
    }
    return naive_square(v);
}

double euclidean_distance(std::span<const double> a,
                           std::span<const double> b,
                           MultiplyMethod mul_method,
                           SumMethod sum_method) {
    assert(a.size() == b.size());

    std::vector<DoublePair> squares(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        squares[i] = do_square(b[i] - a[i], mul_method);
    }

    // sqrt_refined инстанцируется под конкретный Order каждой ветки
    switch (sum_method) {
        case SumMethod::Naive:      return sqrt_refined(naive_sum(squares));
        case SumMethod::OgitaOishi: return sqrt_refined(ogita_oishi_sum(squares));
        case SumMethod::KBN2:       return sqrt_refined(kbn2_sum(squares));
        case SumMethod::KBN3:       return sqrt_refined(kbn3_sum(squares));
        case SumMethod::KBN4:       return sqrt_refined(kbn4_sum(squares));
    }
    return sqrt_refined(naive_sum(squares));
}

AllResults euclidean_all(std::span<const double> a,
                          std::span<const double> b) {
    AllResults out;
    constexpr MultiplyMethod mul_methods[3] = {
        MultiplyMethod::Naive, MultiplyMethod::FMA, MultiplyMethod::Ozaki
    };
    constexpr SumMethod sum_methods[5] = {
        SumMethod::Naive, SumMethod::OgitaOishi,
        SumMethod::KBN2,  SumMethod::KBN3, SumMethod::KBN4
    };
    for (int mi = 0; mi < 3; ++mi) {
        for (int si = 0; si < 5; ++si) {
            out.data[mi][si] = euclidean_distance(a, b, mul_methods[mi], sum_methods[si]);
        }
    }
    return out;
}

std::string_view method_name(MultiplyMethod m) noexcept {
    switch (m) {
        case MultiplyMethod::Naive:  return "Naive";
        case MultiplyMethod::FMA:    return "FMA";
        case MultiplyMethod::Ozaki:  return "Ozaki";
    }
    return "Unknown";
}

std::string_view method_name(SumMethod s) noexcept {
    switch (s) {
        case SumMethod::Naive:       return "Naive";
        case SumMethod::OgitaOishi:  return "OgitaOishi";
        case SumMethod::KBN2:        return "KBN2";
        case SumMethod::KBN3:        return "KBN3";
        case SumMethod::KBN4:        return "KBN4";
    }
    return "Unknown";
}

}