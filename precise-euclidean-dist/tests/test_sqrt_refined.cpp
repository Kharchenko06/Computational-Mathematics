#include "sqrt_refined.h"
#include "compensated_sum.h"
#include "multiply.h"
#include "ulp_utils.h"
#include "gmp_reference.h"

#include <cmath>
#include <cstdio>
#include <vector>

struct TestCase {
    std::vector<double> a;
    std::vector<double> b;
    const char* description;
};

static void run_sqrt_comparison(const TestCase& tc) {
    std::vector<euclidean::DoublePair> squares(tc.a.size());
    for (std::size_t i = 0; i < tc.a.size(); ++i) {
        const double d = tc.b[i] - tc.a[i];
        squares[i] = euclidean::fma_square(d);
    }

    // kbn4_sum возвращает CompensatedSum<4> — тип выводится автоматически
    const auto cs = euclidean::kbn4_sum(squares);

    const double r_naive   = euclidean::sqrt_naive(cs);
    const double r_refined = euclidean::sqrt_refined(cs);
    const double r_ref     = euclidean::gmp_euclidean(tc.a, tc.b);

    const std::int64_t ulp_naive   = euclidean::ulp_error(r_naive,   r_ref);
    const std::int64_t ulp_refined = euclidean::ulp_error(r_refined, r_ref);

    const char* status = (ulp_refined <= 1) ? "PASS" : "FAIL";
    std::printf("  [%s] %s\n", status, tc.description);
    std::printf("    naive:   %.16e  (ULP=%lld)\n", r_naive,   static_cast<long long>(ulp_naive));
    std::printf("    refined: %.16e  (ULP=%lld)\n", r_refined, static_cast<long long>(ulp_refined));
    std::printf("    GMP:     %.16e\n\n", r_ref);
}

int main() {
    std::printf("=== test_sqrt_refined ===\n\n");

    const std::vector<TestCase> cases = {
        { {0,0,0}, {1,1,1},        "unit cube diagonal" },
        { {0,0,0}, {3,4,0},        "3-4-5 triangle" },
        { {1.5, 2.3, -1.7}, {4.1, -0.5, 3.3}, "general floats" },
        { {1e15, 1e15, 1e15}, {1e15+1, 1e15+1, 1e15+1}, "large coords, small diff" },
        { {1e-100, 0, 0}, {2e-100, 0, 0},  "very small values" },
        { {0,0,0}, {1e150, 1e150, 1e150},  "very large values" },
        { {0,0,0}, {1.0, std::sqrt(3.0), 0}, "sum = 4.0 exactly" },
        { {1,2,3,4,5}, {6,7,8,9,10}, "R^5 general" },
    };

    for (const auto& tc : cases) {
        run_sqrt_comparison(tc);
    }
    return 0;
}