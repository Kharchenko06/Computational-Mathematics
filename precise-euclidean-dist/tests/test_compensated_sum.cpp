#include "compensated_sum.h"
#include "multiply.h"
#include "ulp_utils.h"
#include "gmp_reference.h"

#include <cmath>
#include <cstdio>
#include <vector>

static double total(const euclidean::CompensatedSum& cs) {
    double s = cs.sum;
    for (int i = 0; i < cs.order - 1; ++i) {
        s += cs.error[i];
    }
    return s;
}

static void test_two_sum_primitive() {
    std::printf("=== two_sum: hi == fl(a+b) ===\n");
    const struct { double a, b; } cases[] = {
        {1.0, 1e-16}, {1e10, 1.0}, {1.0/3.0, 2.0/3.0},
        {1.0, -1.0 + 1e-15}, {1e100, 1e-100}
    };
    int pass = 0;
    for (auto [a, b] : cases) {
        const euclidean::DoublePair r = euclidean::two_sum(a, b);
        const bool ok = (r.hi == a + b);
        std::printf("  two_sum(%+.3e, %+.3e) -> hi=%+.6e lo=%+.6e  %s\n",
            a, b, r.hi, r.lo, ok ? "PASS" : "FAIL");
        if (ok) { ++pass; }
    }
    std::printf("  two_sum: %d/5 passed\n\n", pass);
}

static void test_sum_algorithms() {
    std::printf("=== Алгоритмы суммирования vs GMP ===\n");

    const double dx = 1e-7, dy = 2e-7, dz = 3e-7;
    const euclidean::DoublePair a = euclidean::fma_square(dx);
    const euclidean::DoublePair b = euclidean::fma_square(dy);
    const euclidean::DoublePair c = euclidean::fma_square(dz);

    const euclidean::Point3D p1{0.0, 0.0, 0.0};
    const euclidean::Point3D p2{dx, dy, dz};
    const double ref_dist = euclidean::gmp_euclidean(p1, p2);
    const double ref_sum  = ref_dist * ref_dist;

    const euclidean::CompensatedSum results[5] = {
        euclidean::naive_sum(a, b, c),
        euclidean::ogita_oishi_sum(a, b, c),
        euclidean::kbn2_sum(a, b, c),
        euclidean::kbn3_sum(a, b, c),
        euclidean::kbn4_sum(a, b, c)
    };
    const char* names[5] = {
        "naive      ", "ogita_oishi", "kbn2       ", "kbn3       ", "kbn4       "
    };
    for (int i = 0; i < 5; ++i) {
        const double t = total(results[i]);
        const std::int64_t ulps = euclidean::ulp_error(t, ref_sum);
        std::printf("  %s: sum=%.16e  error[0]=%.6e  ULP=%lld\n",
            names[i], results[i].sum, results[i].error[0],
            static_cast<long long>(ulps));
    }
    std::printf("\n");
}

static void test_catastrophic_cancellation() {
    std::printf("=== Катастрофическая отмена: ожидаем сумму квадратов = 3.0 ===\n");
    const euclidean::Point3D p1{1e15, 1e15, 1e15};
    const euclidean::Point3D p2{1e15 + 1.0, 1e15 + 1.0, 1e15 + 1.0};

    const double dx = p2.x - p1.x;
    const double dy = p2.y - p1.y;
    const double dz = p2.z - p1.z;

    const euclidean::DoublePair a = euclidean::fma_square(dx);
    const euclidean::DoublePair b = euclidean::fma_square(dy);
    const euclidean::DoublePair c = euclidean::fma_square(dz);

    const double ref = euclidean::gmp_euclidean(p1, p2);
    std::printf("  GMP ref: %.16e\n", ref);

    const euclidean::CompensatedSum results[5] = {
        euclidean::naive_sum(a, b, c),      euclidean::ogita_oishi_sum(a, b, c),
        euclidean::kbn2_sum(a, b, c),       euclidean::kbn3_sum(a, b, c),
        euclidean::kbn4_sum(a, b, c)
    };
    const char* names[5] = { "naive", "ogita_oishi", "kbn2", "kbn3", "kbn4" };
    for (int i = 0; i < 5; ++i) {
        const double t = total(results[i]);
        const std::int64_t ulps = euclidean::ulp_error(t, 3.0);
        std::printf("  %s: total=%.16e  ULP vs 3.0=%lld\n",
            names[i], t, static_cast<long long>(ulps));
    }
    std::printf("\n");
}

int main() {
    std::printf("=== test_compensated_sum ===\n\n");
    test_two_sum_primitive();
    test_sum_algorithms();
    test_catastrophic_cancellation();
    return 0;
}