#include "sqrt_refined.h"
#include "compensated_sum.h"
#include "multiply.h"
#include "ulp_utils.h"
#include "gmp_reference.h"

#include <cmath>
#include <cstdio>
#include <vector>

struct TestCase {
    euclidean::Point3D p1;
    euclidean::Point3D p2;
    const char* description;
};

static void run_sqrt_comparison(const TestCase& tc) {
    const double dx = tc.p2.x - tc.p1.x;
    const double dy = tc.p2.y - tc.p1.y;
    const double dz = tc.p2.z - tc.p1.z;

    const euclidean::CompensatedSum cs = euclidean::kbn4_sum(
        euclidean::fma_square(dx),
        euclidean::fma_square(dy),
        euclidean::fma_square(dz)
    );

    const double r_naive   = euclidean::sqrt_naive(cs);
    const double r_refined = euclidean::sqrt_refined(cs);
    const double r_ref     = euclidean::gmp_euclidean(tc.p1, tc.p2);

    const std::int64_t ulp_naive   = euclidean::ulp_error(r_naive, r_ref);
    const std::int64_t ulp_refined = euclidean::ulp_error(r_refined, r_ref);

    // Считаем тест пройденным если refined <= 1 ULP
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
    };

    for (const auto& tc : cases) {
        run_sqrt_comparison(tc);
    }
    return 0;
}