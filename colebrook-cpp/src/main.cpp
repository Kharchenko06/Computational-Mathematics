#include "../include/colebrook.h"
#include "../include/methods.h"
#include "../include/utils.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <string>

// ---------------------------------------------------------------------------
// Проверочные примеры (верифицированы через Wolfram Alpha)
//
// Wolfram Alpha запрос (пример 1):
//   solve 1/sqrt(f) + 2*log10(0.001081 + 0.0000251/sqrt(f)) = 0 for f
//
// Три сценария покрывают типичные инженерные случаи:
//   1. Умеренно шероховатая труба,  Re=100 000
//   2. Практически гладкая труба,   Re=500 000  (высокий Re, малая шероховатость)
//   3. Грубая труба,                Re=10 000   (низкий Re, высокая шероховатость)
// ---------------------------------------------------------------------------
struct Example {
    const char* title;
    double Re, eps, D;
    double ref;     // Wolfram Alpha reference (15 значащих цифр)
};

static const Example EXAMPLES[] = {
    {
        "Re=100 000, eps=0.0002 m, D=0.05 m  (умеренная шероховатость)",
        100000.0, 0.0002, 0.05,
        0.029500688911511   // Wolfram Alpha: f ≈ 0.02950068891
    },
    {
        "Re=500 000, eps=0.00005 m, D=0.2 m  (гладкая труба, высокий Re)",
        500000.0, 0.00005, 0.2,
        0.015870572581033   // Wolfram Alpha: f ≈ 0.01587057258
    },
    {
        "Re=10 000, eps=0.001 m,   D=0.1 m  (шероховатая труба, низкий Re)",
        10000.0, 0.001, 0.1,
        0.043126584706812   // Wolfram Alpha: f ≈ 0.04312658471
    },
};

static void runExamples()
{
    const std::string RESET  = "\033[0m";
    const std::string BOLD   = "\033[1m";
    const std::string CYAN   = "\033[36m";
    const std::string GREEN  = "\033[32m";
    const std::string RED    = "\033[31m";
    const std::string YELLOW = "\033[33m";
    const std::string line72(72, '=');

    std::cout << "\n" << BOLD << CYAN
              << line72 << "\n"
              << "  ВСТРОЕННЫЕ ПРИМЕРЫ  (верификация по Wolfram Alpha)\n"
              << line72 << "\n" << RESET;

    SolverConfig cfg;
    cfg.tol_step = 1e-12;
    cfg.tol_func = 1e-12;
    cfg.max_iter = 200;

    int passed = 0;
    const int N = static_cast<int>(sizeof(EXAMPLES) / sizeof(EXAMPLES[0]));

    for (int i = 0; i < N; ++i) {
        const Example& ex = EXAMPLES[i];
        std::cout << "\n  " << BOLD << "Пример " << (i + 1) << ": "
                  << ex.title << RESET << "\n";

        PipeParams p(ex.Re, ex.eps, ex.D);
        double f0 = initialGuess(p);

        IterResult rN = solveNewton(f0, p, cfg);
        IterResult rH = solveHalley(f0, p, cfg);

        // Берем результат Ньютона (он всегда сходится для Colebrook-White)
        double f_got = rN.converged ? rN.solution : rH.solution;
        double err   = std::fabs(f_got - ex.ref);
        bool   ok    = err < 1e-10;
        if (ok) ++passed;

        std::cout << std::fixed << std::setprecision(15);
        std::cout << "    f* (Newton)    = " << f_got     << "\n";
        std::cout << "    f* (Wolfram)   = " << ex.ref    << "\n";

        std::cout << std::scientific << std::setprecision(3);
        std::cout << "    |∆f|           = " << err << "  ";
        std::cout << (ok ? GREEN + "[OK]" : RED + "[FAIL]") << RESET << "\n";
        std::cout << "    Итераций Newton / Halley : "
                  << rN.iterations << " / " << rH.iterations << "\n";
    }

    std::cout << "\n" << std::string(72, '-') << "\n";
    std::cout << BOLD << "  Результат: " << passed << "/" << N << " примеров OK";
    if (passed == N)
        std::cout << "  " << GREEN << "Все тесты пройдены." << RESET;
    else
        std::cout << "  " << RED << "Есть ошибки!" << RESET;
    std::cout << "\n" << std::string(72, '=') << "\n\n";
}


static double readDouble(const std::string& prompt, double lo, double hi)
{
    double val;
    while (true) {
        std::cout << prompt;
        if (std::cin >> val && val >= lo && val <= hi)
            return val;
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << "  Некорректный ввод. Диапазон: [" << lo << ", " << hi << "]\n";
    }
}

int main()
{
    std::cout << "\n\033[1m\033[36m";
    std::cout << "  +------------------------------------------------------+\n";
    std::cout << "  |   Colebrook-White Equation Solver                   |\n";
    std::cout << "  |   Методы: Newton-Raphson  +  Halley                 |\n";
    std::cout << "  +------------------------------------------------------+\n";
    std::cout << "\033[0m\n";

    // Сначала прогоняем встроенные примеры — удобно для быстрой проверки
    runExamples();

    std::cout << "  Введите параметры трубопровода.\n";
    std::cout << "  (Для демо: Re=100000, eps=0.0002 м, D=0.05 м)\n\n";

    // Re > 4000 - граница турбулентного режима, ниже Colebrook-White неприменим
    double Re  = readDouble("  Число Рейнольдса Re    : ", 2300.0, 1e9);
    double eps = readDouble("  Шероховатость eps, м   : ", 0.0,    1.0);
    double D   = readDouble("  Диаметр трубы D, м     : ", 1e-4,   10.0);

    PipeParams p(Re, eps, D);

    if (Re < 4000)
        std::cout << "\033[33m  Предупреждение: Re < 4000, режим не турбулентный."
                     " Используйте f = 64/Re.\033[0m\n\n";

    // eps/D > 0.05 - за пределами области применимости Colebrook-White
    if (p.relRoughness() > 0.05)
        std::cout << "\033[33m  Предупреждение: eps/D = " << p.relRoughness()
                  << " > 0.05, очень высокая шероховатость.\033[0m\n\n";

    double f0 = initialGuess(p);
    printProblem(p, f0);

    SolverConfig cfg;
    cfg.tol_step = 1e-12;
    cfg.tol_func = 1e-12;
    cfg.max_iter = 200;

    IterResult resNewton = solveNewton(f0, p, cfg);
    printIterations(resNewton);
    printSummary(resNewton);

    IterResult resHalley = solveHalley(f0, p, cfg);
    printIterations(resHalley);
    printSummary(resHalley);

    printComparison(resNewton, resHalley);

    double f_final = resNewton.converged ? resNewton.solution : resHalley.solution;
    std::cout << "\n\033[1m  Результат:\033[0m\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "    f       = " << f_final << "\n";
    std::cout << "    1/sqrt(f) = " << 1.0 / std::sqrt(f_final)
              << "  (левая часть уравнения)\n\n";

    return 0;
}