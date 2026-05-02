#include "compensated_sum.h"
#include "multiply.h"
#include "ulp_utils.h"
#include "gmp_reference.h"

#include <cmath>
#include <cstdio>
#include <vector>

// Суммирует все члены CompensatedSum<Order> для сравнения с эталоном
template<int Order>
static double total(const euclidean::CompensatedSum<Order>& cs) {
    double s = cs.sum;
    for (int i = 0; i < Order - 1; ++i) {
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
        const euclidean::DoublePair r = euclidean::detail::two_sum(a, b);
        const bool ok = (r.hi == a + b);
        std::printf("  two_sum(%+.3e, %+.3e) -> hi=%+.6e lo=%+.6e  %s\n",
            a, b, r.hi, r.lo, ok ? "PASS" : "FAIL");
        if (ok) { ++pass; }
    }
    std::printf("  two_sum: %d/5 passed\n\n", pass);
}

static void test_sum_algorithms() {
    std::printf("=== Алгоритмы суммирования vs GMP ===\n");

    const std::vector<double> va = {0.0, 0.0, 0.0};
    const std::vector<double> vb = {1e-7, 2e-7, 3e-7};

    std::vector<euclidean::DoublePair> squares(3);
    for (int i = 0; i < 3; ++i) {
        squares[i] = euclidean::fma_square(vb[i] - va[i]);
    }

    const double ref_dist = euclidean::gmp_euclidean(va, vb);
    const double ref_sum  = ref_dist * ref_dist;

    // Каждый алгоритм возвращает CompensatedSum<N> своего порядка
    auto run = [&](auto cs, const char* name) {
        const double t = total(cs);
        const std::int64_t ulps = euclidean::ulp_error(t, ref_sum);
        std::printf("  %s: sum=%.16e  error[0]=%.6e  ULP=%lld\n",
            name, cs.sum, cs.error.empty() ? 0.0 : cs.error[0],
            static_cast<long long>(ulps));
    };

    run(euclidean::naive_sum(squares),       "naive      ");
    run(euclidean::ogita_oishi_sum(squares),  "ogita_oishi");
    run(euclidean::kbn2_sum(squares),         "kbn2       ");
    run(euclidean::kbn3_sum(squares),         "kbn3       ");
    run(euclidean::kbn4_sum(squares),         "kbn4       ");
    std::printf("\n");
}

static void test_catastrophic_cancellation() {
    std::printf("=== Катастрофическое сокращение: ожидаем сумму квадратов = 3.0 ===\n");
    const std::vector<double> va = {1e15, 1e15, 1e15};
    const std::vector<double> vb = {1e15 + 1.0, 1e15 + 1.0, 1e15 + 1.0};

    std::vector<euclidean::DoublePair> squares(3);
    for (int i = 0; i < 3; ++i) {
        squares[i] = euclidean::fma_square(vb[i] - va[i]);
    }

    const double ref = euclidean::gmp_euclidean(va, vb);
    std::printf("  GMP ref: %.16e\n", ref);

    auto run = [&](auto cs, const char* name) {
        const double t = total(cs);
        const std::int64_t ulps = euclidean::ulp_error(t, 3.0);
        std::printf("  %s: total=%.16e  ULP vs 3.0=%lld\n",
            name, t, static_cast<long long>(ulps));
    };

    run(euclidean::naive_sum(squares),       "naive      ");
    run(euclidean::ogita_oishi_sum(squares),  "ogita_oishi");
    run(euclidean::kbn2_sum(squares),         "kbn2       ");
    run(euclidean::kbn3_sum(squares),         "kbn3       ");
    run(euclidean::kbn4_sum(squares),         "kbn4       ");
    std::printf("\n");
}

static void test_large_dimension() {
    std::printf("=== R^100: sum of 100 ones, expected dist = 10.0 ===\n");
    constexpr int N = 100;
    std::vector<euclidean::DoublePair> squares(N);
    for (int i = 0; i < N; ++i) {
        squares[i] = euclidean::fma_square(1.0);
    }

    const std::vector<double> va(N, 0.0);
    const std::vector<double> vb(N, 1.0);
    const double ref = euclidean::gmp_euclidean(va, vb);
    std::printf("  GMP ref: %.16e\n", ref);

    auto run = [&](auto cs, const char* name) {
        const double t = total(cs);
        const std::int64_t ulps = euclidean::ulp_error(t, ref * ref);
        std::printf("  %s: total=%.16e  ULP=%lld\n",
            name, t, static_cast<long long>(ulps));
    };

    run(euclidean::naive_sum(squares),       "naive      ");
    run(euclidean::ogita_oishi_sum(squares),  "ogita_oishi");
    run(euclidean::kbn2_sum(squares),         "kbn2       ");
    run(euclidean::kbn3_sum(squares),         "kbn3       ");
    run(euclidean::kbn4_sum(squares),         "kbn4       ");
    std::printf("\n");
}

int main() {
    std::printf("=== test_compensated_sum ===\n\n");
    test_two_sum_primitive();
    test_sum_algorithms();
    test_catastrophic_cancellation();
    test_large_dimension();
    return 0;
}