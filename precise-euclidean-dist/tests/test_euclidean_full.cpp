#include "euclidean.h"
#include "ulp_utils.h"
#include "gmp_reference.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

struct TestCase {
    std::vector<double> a;
    std::vector<double> b;
    const char* description;
};

static bool run_full_table(const TestCase& tc) {
    const double gmp_ref = euclidean::gmp_euclidean(tc.a, tc.b);
    const euclidean::AllResults all = euclidean::euclidean_all(tc.a, tc.b);

    std::printf("  %s\n", tc.description);
    std::printf("  GMP: %.16e\n", gmp_ref);
    std::printf("  %-12s %-12s  %-22s  %4s\n", "Multiply", "Sum", "Result", "ULP");

    constexpr euclidean::MultiplyMethod mul_methods[3] = {
        euclidean::MultiplyMethod::Naive,
        euclidean::MultiplyMethod::FMA,
        euclidean::MultiplyMethod::Ozaki
    };
    constexpr euclidean::SumMethod sum_methods[5] = {
        euclidean::SumMethod::Naive,    euclidean::SumMethod::OgitaOishi,
        euclidean::SumMethod::KBN2,     euclidean::SumMethod::KBN3,
        euclidean::SumMethod::KBN4
    };

    std::int64_t best_ulp = INT64_MAX;

    for (int mi = 0; mi < 3; ++mi) {
        for (int si = 0; si < 5; ++si) {
            const double result = all.data[mi][si];
            const std::int64_t ulp = euclidean::ulp_error(result, gmp_ref);
            if (ulp < best_ulp) { best_ulp = ulp; }
            const char* mark = (ulp == 0) ? " *" : (ulp <= 1) ? " +" : "";
            std::printf("  %-12s %-12s  %.16e  %4lld%s\n",
                euclidean::method_name(mul_methods[mi]).data(),
                euclidean::method_name(sum_methods[si]).data(),
                result,
                static_cast<long long>(ulp),
                mark);
        }
    }

    const char* verdict = (best_ulp <= 1) ? "PASS" : "FAIL";
    std::printf("  Best ULP: %lld  -> %s\n\n", static_cast<long long>(best_ulp), verdict);
    return (best_ulp <= 1);
}

int main() {
    std::printf("=== test_euclidean_full (* = 0 ULP, + = 1 ULP) ===\n\n");

    const std::vector<TestCase> cases = {
        // R^3 — базовые случаи
        { {0,0,0}, {1,0,0},        "R3: единичное расстояние по X" },
        { {0,0,0}, {1,1,1},        "R3: диагональ куба, sqrt(3)" },
        { {0,0,0}, {3,4,0},        "R3: треугольник 3-4-5" },
        { {1.5, 2.3, -1.7}, {4.1, -0.5, 3.3},  "R3: общий случай" },
        { {0.1, 0.2, 0.3}, {0.4, 0.5, 0.6},    "R3: малые десятичные" },
        { {1e15, 1e15, 1e15}, {1e15+1.0, 1e15+1.0, 1e15+1.0}, "R3: отмена, diff=sqrt(3)" },
        { {1e8, 2e8, 3e8}, {1e8+3.0, 2e8+4.0, 3e8+0.0},       "R3: отмена, diff=5" },
        { {0,0,0}, {1e-100, 1e-100, 1e-100},   "R3: малые значения" },
        { {0,0,0}, {1e100, 1e100, 1e100},       "R3: большие значения" },
        { {0,0,0}, {1e-8, 1.0, 1e8},            "R3: разные порядки" },
        { {0,0,0}, {1,1,0},  "R3: sqrt(2)" },
        // R^2 — плоский случай
        { {0,0}, {3,4},   "R2: треугольник 3-4-5" },
        { {0,0}, {1,1},   "R2: sqrt(2)" },
        // R^5
        { {0,0,0,0,0}, {1,1,1,1,1}, "R5: sqrt(5)" },
        { {1,2,3,4,5}, {6,7,8,9,10}, "R5: diff=(5,5,5,5,5), dist=sqrt(125)" },
        // R^n с катастрофической отменой
        { {1e15,1e15,1e15,1e15,1e15},
          {1e15+1,1e15+1,1e15+1,1e15+1,1e15+1}, "R5: отмена, diff=sqrt(5)" },
    };

    int pass_count = 0;
    for (const auto& tc : cases) {
        if (run_full_table(tc)) { ++pass_count; }
    }

    std::printf("Итог: %d/%zu прошли (<= 1 ULP)\n",
        pass_count, cases.size());

    return (pass_count == static_cast<int>(cases.size())) ? EXIT_SUCCESS : EXIT_FAILURE;
}