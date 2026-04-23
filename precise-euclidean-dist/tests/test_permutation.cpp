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

// Все 6 перестановок индексов {0,1,2}
static const std::array<std::array<int,3>, 6> PERMS = {{
    {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}
}};

static double dist_permuted(const std::array<double,3>& deltas,
                             const std::array<int,3>& perm,
                             euclidean::MultiplyMethod mm,
                             euclidean::SumMethod sm) noexcept {
    const double d0 = deltas[perm[0]];
    const double d1 = deltas[perm[1]];
    const double d2 = deltas[perm[2]];

    euclidean::DoublePair sq[3];
    switch (mm) {
        case euclidean::MultiplyMethod::Naive:
            sq[0] = euclidean::naive_square(d0);
            sq[1] = euclidean::naive_square(d1);
            sq[2] = euclidean::naive_square(d2); break;
        case euclidean::MultiplyMethod::FMA:
            sq[0] = euclidean::fma_square(d0);
            sq[1] = euclidean::fma_square(d1);
            sq[2] = euclidean::fma_square(d2); break;
        case euclidean::MultiplyMethod::Ozaki:
            sq[0] = euclidean::ozaki_square(d0);
            sq[1] = euclidean::ozaki_square(d1);
            sq[2] = euclidean::ozaki_square(d2); break;
    }

    euclidean::CompensatedSum cs;
    switch (sm) {
        case euclidean::SumMethod::Naive:      cs = euclidean::naive_sum(sq[0], sq[1], sq[2]); break;
        case euclidean::SumMethod::OgitaOishi: cs = euclidean::ogita_oishi_sum(sq[0], sq[1], sq[2]); break;
        case euclidean::SumMethod::KBN2:       cs = euclidean::kbn2_sum(sq[0], sq[1], sq[2]); break;
        case euclidean::SumMethod::KBN3:       cs = euclidean::kbn3_sum(sq[0], sq[1], sq[2]); break;
        case euclidean::SumMethod::KBN4:       cs = euclidean::kbn4_sum(sq[0], sq[1], sq[2]); break;
    }
    return euclidean::sqrt_refined(cs);
}

static void test_permutation_invariance(const char* label,
                                         const std::array<double,3>& deltas,
                                         double gmp_ref) {
    std::printf("  %s  (GMP ref = %.16e)\n", label, gmp_ref);

    const euclidean::MultiplyMethod mul_methods[3] = {
        euclidean::MultiplyMethod::Naive,
        euclidean::MultiplyMethod::FMA,
        euclidean::MultiplyMethod::Ozaki
    };
    const euclidean::SumMethod sum_methods[5] = {
        euclidean::SumMethod::Naive,    euclidean::SumMethod::OgitaOishi,
        euclidean::SumMethod::KBN2,     euclidean::SumMethod::KBN3,
        euclidean::SumMethod::KBN4
    };

    std::printf("  %-12s %-12s  %5s  %8s\n", "Multiply", "Sum", "Inv?", "MaxULP");

    for (auto mm : mul_methods) {
        for (auto sm : sum_methods) {
            std::array<double, 6> results{};
            for (int pi = 0; pi < 6; ++pi) {
                results[pi] = dist_permuted(deltas, PERMS[pi], mm, sm);
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

    struct { euclidean::Point3D p1, p2; const char* label; } cases[] = {
        { {0,0,0}, {1,1,1},                        "Symmetric (1,1,1)" },
        { {0.1, 0.2, 0.3}, {1.4, 0.5, 2.7},       "Asymmetric floats" },
        { {0, 0, 0}, {1e-8, 1.0, 1e8},             "Mixed magnitudes" },
        { {1e15, 2e15, 3e15}, {1e15+1, 2e15+2, 3e15+3}, "Large coords + small diff" },
    };

    for (auto& tc : cases) {
        const double ref = euclidean::gmp_euclidean(tc.p1, tc.p2);
        const std::array<double,3> d = {
            tc.p2.x - tc.p1.x, tc.p2.y - tc.p1.y, tc.p2.z - tc.p1.z
        };
        test_permutation_invariance(tc.label, d, ref);
    }
    return 0;
}