#include "euclidean.h"
#include "multiply.h"
#include "compensated_sum.h"
#include "sqrt_refined.h"

namespace euclidean {

static DoublePair do_square(double v, MultiplyMethod m) noexcept {
    switch (m) {
        case MultiplyMethod::Naive:  return naive_square(v);
        case MultiplyMethod::FMA:    return fma_square(v);
        case MultiplyMethod::Ozaki:  return ozaki_square(v);
    }
    return naive_square(v);
}

static CompensatedSum do_sum(DoublePair a, DoublePair b, DoublePair c,
                              SumMethod sm) noexcept {
    switch (sm) {
        case SumMethod::Naive:       return naive_sum(a, b, c);
        case SumMethod::OgitaOishi:  return ogita_oishi_sum(a, b, c);
        case SumMethod::KBN2:        return kbn2_sum(a, b, c);
        case SumMethod::KBN3:        return kbn3_sum(a, b, c);
        case SumMethod::KBN4:        return kbn4_sum(a, b, c);
    }
    return naive_sum(a, b, c);
}

double euclidean_distance(const Point3D&  p1,
                           const Point3D&  p2,
                           MultiplyMethod  mul_method,
                           SumMethod       sum_method) noexcept {
    const double dx = p2.x - p1.x;
    const double dy = p2.y - p1.y;
    const double dz = p2.z - p1.z;

    const DoublePair dx2 = do_square(dx, mul_method);
    const DoublePair dy2 = do_square(dy, mul_method);
    const DoublePair dz2 = do_square(dz, mul_method);

    const CompensatedSum cs = do_sum(dx2, dy2, dz2, sum_method);
    return sqrt_refined(cs);
}

AllResults euclidean_all(const Point3D& p1, const Point3D& p2) noexcept {
    AllResults out;
    const MultiplyMethod mul_methods[3] = {
        MultiplyMethod::Naive, MultiplyMethod::FMA, MultiplyMethod::Ozaki
    };
    const SumMethod sum_methods[5] = {
        SumMethod::Naive, SumMethod::OgitaOishi,
        SumMethod::KBN2,  SumMethod::KBN3, SumMethod::KBN4
    };
    for (int mi = 0; mi < 3; ++mi) {
        for (int si = 0; si < 5; ++si) {
            out.data[mi][si] = euclidean_distance(p1, p2, mul_methods[mi], sum_methods[si]);
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