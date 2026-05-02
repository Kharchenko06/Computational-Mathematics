#include "euclidean.h"
#include "multiply.h"
#include "compensated_sum.h"
#include "sqrt_refined.h"
#include "ulp_utils.h"
#include "gmp_reference.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

// Все 6 перестановок индексов {0,1,2} для R^3
static const std::array<std::array<int,3>, 6> PERMS3 = {{
    {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}
}};

static double dist_permuted3(const std::vector<double>& deltas,
                              const std::array<int,3>& perm,
                              euclidean::MultiplyMethod mm,
                              euclidean::SumMethod sm) {
    auto sq = [&](int idx) -> euclidean::DoublePair {
        const double d = deltas[perm[idx]];
        switch (mm) {
            case euclidean::MultiplyMethod::Naive:  return euclidean::naive_square(d);
            case euclidean::MultiplyMethod::FMA:    return euclidean::fma_square(d);
            case euclidean::MultiplyMethod::Ozaki:  return euclidean::ozaki_square(d);
        }
        return euclidean::naive_square(d);
    };

    const std::vector<euclidean::DoublePair> terms = { sq(0), sq(1), sq(2) };

    // sqrt_refined инстанцируется под конкретный Order каждой ветки
    switch (sm) {
        case euclidean::SumMethod::Naive:
            return euclidean::sqrt_refined(euclidean::naive_sum(terms));
        case euclidean::SumMethod::OgitaOishi:
            return euclidean::sqrt_refined(euclidean::ogita_oishi_sum(terms));
        case euclidean::SumMethod::KBN2:
            return euclidean::sqrt_refined(euclidean::kbn2_sum(terms));
        case euclidean::SumMethod::KBN3:
            return euclidean::sqrt_refined(euclidean::kbn3_sum(terms));
        case euclidean::SumMethod::KBN4:
            return euclidean::sqrt_refined(euclidean::kbn4_sum(terms));
    }
    return euclidean::sqrt_refined(euclidean::naive_sum(terms));
}

static void test_permutation_invariance(const char* label,
                                         const std::vector<double>& deltas,
                                         double gmp_ref) {
    std::printf("  %s  (GMP ref = %.16e)\n", label, gmp_ref);

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

    std::printf("  %-12s %-12s  %5s  %8s\n", "Multiply", "Sum", "Inv?", "MaxULP");

    for (auto mm : mul_methods) {
        for (auto sm : sum_methods) {
            std::array<double, 6> results{};
            for (int pi = 0; pi < 6; ++pi) {
                results[pi] = dist_permuted3(deltas, PERMS3[pi], mm, sm);
            }
            const bool invariant = std::all_of(results.begin(), results.end(),
                [&](double v){ return v == results[0]; });

            std::int64_t max_ulp = 0;
            for (double r : results) {
                const std::int64_t u = euclidean::ulp_error(r, gmp_ref);
                if (u > max_ulp) { max_ulp = u; }
            }

            std::printf("  %-12s %-12s  %5s  %8lld\n",
                euclidean::method_name(mm).data(),
                euclidean::method_name(sm).data(),
                invariant ? "YES" : "NO",
                static_cast<long long>(max_ulp));
        }
    }
    std::printf("\n");
}

int main() {
    std::printf("=== test_permutation ===\n\n");

    struct { std::vector<double> a, b; const char* label; } cases[] = {
        { {0,0,0}, {1,1,1},                        "R3: Symmetric (1,1,1)" },
        { {0.1, 0.2, 0.3}, {1.4, 0.5, 2.7},       "R3: Asymmetric floats" },
        { {0, 0, 0}, {1e-8, 1.0, 1e8},             "R3: Mixed magnitudes" },
        { {1e15, 2e15, 3e15}, {1e15+1, 2e15+2, 3e15+3}, "R3: Large coords + small diff" },
    };

    for (auto& tc : cases) {
        const double ref = euclidean::gmp_euclidean(tc.a, tc.b);
        std::vector<double> d(tc.a.size());
        for (std::size_t i = 0; i < tc.a.size(); ++i) {
            d[i] = tc.b[i] - tc.a[i];
        }
        test_permutation_invariance(tc.label, d, ref);
    }
    return 0;
}