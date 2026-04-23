#include "multiply.h"
#include "ulp_utils.h"
#include "gmp_reference.h"

#include <cmath>
#include <cstdio>
#include <vector>

static void test_veltkamp_split() {
    std::printf("=== veltkamp_split ===\n");
    const std::vector<double> cases = {
        1.0, -1.0, 1.5, 1.234567890123456789,
        1e100, 1e-100, 3.14159265358979323846
    };
    int pass = 0;
    for (double x : cases) {
        const euclidean::DoublePair s = euclidean::veltkamp_split(x);
        // hi + lo должно быть точно равно x (по построению алгоритма Велткампа)
        const bool ok = (s.hi + s.lo == x);
        std::printf("  split(%+.6e) -> hi=%+.6e lo=%+.6e  hi+lo==x: %s\n",
            x, s.hi, s.lo, ok ? "PASS" : "FAIL");
        if (ok) { ++pass; }
    }
    std::printf("  veltkamp_split: %d/%zu passed\n\n", pass, cases.size());
}

static void test_eft_property(const char* name,
                               euclidean::DoublePair (*fn)(double)) {
    std::printf("=== %s: hi == fl(x*x) ===\n", name);
    const std::vector<double> cases = {
        1.0, 1.5, 1.1, 3.14159265358979,
        1.23456789012345678e10, 1.23456789012345678e-10,
        1.0 / 3.0, std::sqrt(2.0), 1e150, 1e-150
    };
    int pass = 0;
    for (double x : cases) {
        const euclidean::DoublePair r = fn(x);
        const bool hi_ok = (r.hi == x * x);
        const bool lo_small = (r.lo == 0.0 || std::abs(r.lo) < std::abs(r.hi));
        const bool ok = hi_ok && lo_small;
        std::printf("  x=%+.6e  hi=%+.6e  lo=%+.6e  %s\n",
            x, r.hi, r.lo, ok ? "PASS" : "FAIL");
        if (ok) { ++pass; }
    }
    std::printf("  %s: %d/%zu passed\n\n", name, pass, cases.size());
}

static void test_fma_vs_ozaki() {
    std::printf("=== FMA vs Ozaki: lo должны совпадать побитово ===\n");
    const std::vector<double> cases = {
        1.1, 1.0/3.0, 3.14159265358979, 1.23456789012345678
    };
    int pass = 0;
    for (double x : cases) {
        const euclidean::DoublePair fma_r   = euclidean::fma_square(x);
        const euclidean::DoublePair ozaki_r = euclidean::ozaki_square(x);
        const bool ok = (fma_r.hi == ozaki_r.hi) && (fma_r.lo == ozaki_r.lo);
        std::printf("  x=%+.6e  fma.lo=%+.6e  ozaki.lo=%+.6e  %s\n",
            x, fma_r.lo, ozaki_r.lo, ok ? "PASS" : "FAIL");
        if (ok) { ++pass; }
    }
    std::printf("  FMA vs Ozaki: %d/%zu passed\n\n", pass, cases.size());
}

int main() {
    std::printf("=== test_multiply ===\n\n");
    test_veltkamp_split();
    test_eft_property("naive_square", euclidean::naive_square);
    test_eft_property("fma_square",   euclidean::fma_square);
    test_eft_property("ozaki_square", euclidean::ozaki_square);
    test_fma_vs_ozaki();
    return 0;
}