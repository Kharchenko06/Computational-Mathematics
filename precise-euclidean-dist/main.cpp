// Использование: ./build/euclidean_distance x1 y1 z1 x2 y2 z2

#include "euclidean.h"
#include "ulp_utils.h"
#include "gmp_reference.h"

#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 7) {
        std::printf("Использование: %s x1 y1 z1 x2 y2 z2\n", argv[0]);
        return EXIT_FAILURE;
    }

    const euclidean::Point3D p1 = {
        std::strtod(argv[1], nullptr),
        std::strtod(argv[2], nullptr),
        std::strtod(argv[3], nullptr)
    };
    const euclidean::Point3D p2 = {
        std::strtod(argv[4], nullptr),
        std::strtod(argv[5], nullptr),
        std::strtod(argv[6], nullptr)
    };

    std::printf("P1 = (%.17g, %.17g, %.17g)\n", p1.x, p1.y, p1.z);
    std::printf("P2 = (%.17g, %.17g, %.17g)\n\n", p2.x, p2.y, p2.z);

    const double gmp_ref = euclidean::gmp_euclidean(p1, p2);
    std::printf("GMP (256-bit): %.17e\n\n", gmp_ref);

    const euclidean::AllResults all = euclidean::euclidean_all(p1, p2);

    const euclidean::MultiplyMethod mul_methods[3] = {
        euclidean::MultiplyMethod::Naive,
        euclidean::MultiplyMethod::FMA,
        euclidean::MultiplyMethod::Ozaki
    };
    const euclidean::SumMethod sum_methods[5] = {
        euclidean::SumMethod::Naive,
        euclidean::SumMethod::OgitaOishi,
        euclidean::SumMethod::KBN2,
        euclidean::SumMethod::KBN3,
        euclidean::SumMethod::KBN4
    };

    std::printf("  %-10s  %-12s  %-22s  %6s\n", "Multiply", "Sum", "Result", "ULP");
    std::printf("  %s\n", std::string(58, '-').c_str());

    std::int64_t best_ulp = INT64_MAX;
    double       best_val = 0.0;

    for (int mi = 0; mi < 3; ++mi) {
        for (int si = 0; si < 5; ++si) {
            const double result = all.data[mi][si];
            const std::int64_t ulp = euclidean::ulp_error(result, gmp_ref);
            if (ulp < best_ulp) { best_ulp = ulp; best_val = result; }
            // * = точное совпадение с GMP, + = в пределах 1 ULP
            const char* mark = (ulp == 0) ? " *" : (ulp <= 1) ? " +" : "";
            std::printf("  %-10s  %-12s  %.17e  %6lld%s\n",
                euclidean::method_name(mul_methods[mi]).data(),
                euclidean::method_name(sum_methods[si]).data(),
                result,
                static_cast<long long>(ulp),
                mark);
        }
        std::printf("\n");
    }

    std::printf("  %s\n", std::string(58, '-').c_str());
    std::printf("  Лучший результат: %.17e\n", best_val);
    std::printf("  Лучший ULP:       %lld\n", static_cast<long long>(best_ulp));
    std::printf("  GMP эталон:       %.17e\n", gmp_ref);

    return EXIT_SUCCESS;
}