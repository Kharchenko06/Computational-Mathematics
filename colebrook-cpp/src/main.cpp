#include "../include/colebrook.h"
#include "../include/methods.h"
#include "../include/utils.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <string>

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