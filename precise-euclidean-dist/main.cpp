// Использование:
//   ./euclidean_distance x1 y1 z1 x2 y2 z2
//   ./euclidean_distance --file input.txt
//
// Формат файла:
//   <n>
//   <a1> <a2> ... <an>
//   <b1> <b2> ... <bn>

#include "euclidean.h"
#include "ulp_utils.h"
#include "gmp_reference.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

static void print_table(std::span<const double> a, std::span<const double> b) {
    const double gmp_ref = euclidean::gmp_euclidean(a, b);
    std::printf("GMP (256-bit): %.17e\n\n", gmp_ref);

    const euclidean::AllResults all = euclidean::euclidean_all(a, b);

    constexpr euclidean::MultiplyMethod mul_methods[3] = {
        euclidean::MultiplyMethod::Naive,
        euclidean::MultiplyMethod::FMA,
        euclidean::MultiplyMethod::Ozaki
    };
    constexpr euclidean::SumMethod sum_methods[5] = {
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
}

static bool read_from_file(const char* path,
                            std::vector<double>& a,
                            std::vector<double>& b) {
    std::ifstream f(path);
    if (!f) {
        std::fprintf(stderr, "Ошибка: не удалось открыть файл '%s'\n", path);
        return false;
    }
    std::size_t n = 0;
    if (!(f >> n) || n == 0) {
        std::fprintf(stderr, "Ошибка: некорректная размерность\n");
        return false;
    }
    a.resize(n);
    b.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (!(f >> a[i])) {
            std::fprintf(stderr, "Ошибка: недостаточно элементов в строке a\n");
            return false;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (!(f >> b[i])) {
            std::fprintf(stderr, "Ошибка: недостаточно элементов в строке b\n");
            return false;
        }
    }
    return true;
}

int main(int argc, char* argv[]) {
    std::vector<double> a, b;

    if (argc == 3 && std::strcmp(argv[1], "--file") == 0) {
        if (!read_from_file(argv[2], a, b)) {
            return EXIT_FAILURE;
        }
    } else if (argc == 7) {
        a = { std::strtod(argv[1], nullptr),
              std::strtod(argv[2], nullptr),
              std::strtod(argv[3], nullptr) };
        b = { std::strtod(argv[4], nullptr),
              std::strtod(argv[5], nullptr),
              std::strtod(argv[6], nullptr) };
    } else {
        std::printf("Использование:\n");
        std::printf("  %s x1 y1 z1 x2 y2 z2\n", argv[0]);
        std::printf("  %s --file input.txt\n", argv[0]);
        std::printf("\nФормат файла:\n  <n>\n  <a1> ... <an>\n  <b1> ... <bn>\n");
        return EXIT_FAILURE;
    }

    std::printf("n = %zu\n", a.size());
    std::printf("A = [");
    for (std::size_t i = 0; i < a.size(); ++i) {
        std::printf("%s%.17g", i ? ", " : "", a[i]);
    }
    std::printf("]\n");
    std::printf("B = [");
    for (std::size_t i = 0; i < b.size(); ++i) {
        std::printf("%s%.17g", i ? ", " : "", b[i]);
    }
    std::printf("]\n\n");

    print_table(a, b);
    return EXIT_SUCCESS;
}