#include "../include/utils.h"
#include "../include/colebrook.h"
#include "../include/methods.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

static const std::string RESET  = "\033[0m";
static const std::string BOLD   = "\033[1m";
static const std::string CYAN   = "\033[36m";
static const std::string GREEN  = "\033[32m";
static const std::string YELLOW = "\033[33m";
static const std::string RED    = "\033[31m";
static const std::string MAG    = "\033[35m";

static void line(char c = '-', int n = 72)
{
    std::cout << std::string(n, c) << "\n";
}

double initialGuess(const PipeParams& p)
{
    // Swamee-Jain 1976: явная аппроксимация Colebrook-White
    // 5.74 / Re^0.9 - вязкостный член; показатель 0.9 подобран эмпирически
    double arg = p.eps / (3.7 * p.D) + 5.74 / std::pow(p.Re, 0.9);

    if (arg <= 0.0) return 0.02;

    double logVal = std::log10(arg);

    if (std::fabs(logVal) < 1e-30) return 0.02;

    double f0 = 0.25 / (logVal * logVal);

    // 0.02 - середина типичного диапазона f для турбулентного течения
    if (f0 <= 0.0 || f0 > 1.0 || std::isnan(f0) || std::isinf(f0))
        return 0.02;

    return f0;
}

void printProblem(const PipeParams& p, double f0)
{
    line('=');
    std::cout << BOLD << CYAN
              << "  Уравнение Colebrook-White: численное решение\n"
              << RESET;
    line('=');

    std::cout << BOLD << "\n  Параметры трубопровода:\n" << RESET;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "    Re  = " << std::setw(14) << p.Re      << "\n";
    std::cout << "    eps = " << std::setw(14) << p.eps     << " м\n";
    std::cout << "    D   = " << std::setw(14) << p.D       << " м\n";
    std::cout << "    eps/D = "                << p.relRoughness() << "\n";

    std::cout << "\n  Константы уравнения:\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "    A = eps/(3.7*D) = " << std::setw(16) << p.A << "\n";
    std::cout << "    B = 2.51/Re     = " << std::setw(16) << p.B << "\n";

    std::cout << "\n  " << YELLOW << "Начальное приближение (Swamee-Jain):\n" << RESET;
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "    f0 = " << f0
              << "    (1/sqrt(f0) = " << 1.0 / std::sqrt(f0) << ")\n\n";
    line();
}

void printIterations(const IterResult& res)
{
    std::cout << BOLD << MAG
              << "\n  Метод: " << res.method_name << "\n"
              << RESET;
    line();

    std::cout << BOLD
              << std::setw(5)  << "iter"
              << std::setw(20) << "f"
              << std::setw(18) << "|f_new - f_old|"
              << std::setw(18) << "|F(f)|"
              << RESET << "\n";
    line('-');

    for (const auto& s : res.history) {
        if (s.iter == 0)
            std::cout << YELLOW;
        else if (s.error < 1e-8)
            std::cout << GREEN;

        std::cout << std::setw(5) << s.iter
                  << std::setw(20) << std::fixed << std::setprecision(12) << s.f;

        if (s.error < 0.0)
            std::cout << std::setw(18) << "--";
        else
            std::cout << std::setw(18) << std::scientific << std::setprecision(4) << s.error;

        std::cout << std::setw(18) << std::scientific << std::setprecision(4)
                  << std::fabs(s.F_val);

        if (s.iter == 0)
            std::cout << "  <- начало" << RESET;
        else
            std::cout << RESET;

        std::cout << "\n";
    }
    line('-');
}

void printSummary(const IterResult& res)
{
    std::cout << "\n  " << BOLD << "Итог (" << res.method_name << "):\n" << RESET;

    std::string status = res.converged
        ? (GREEN + "СОШЕЛСЯ" + RESET)
        : (RED   + "НЕ СОШЕЛСЯ" + RESET);

    std::cout << "    Статус         : " << status << "\n";
    std::cout << "    Итераций       : " << res.iterations << "\n";

    std::cout << std::fixed << std::setprecision(15);
    std::cout << "    f*             : " << BOLD << GREEN
              << res.solution << RESET << "\n";
    std::cout << "    1/sqrt(f*)     : " << 1.0 / std::sqrt(res.solution) << "\n";

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "    |F(f*)|        : " << res.final_F    << "\n";
    std::cout << "    Последний шаг  : " << res.final_error << "\n\n";
    line();
}

void printComparison(const IterResult& r1, const IterResult& r2)
{
    line('=');
    std::cout << BOLD << CYAN << "  Сравнение методов\n" << RESET;
    line('=');

    auto yn = [](bool b) -> std::string {
        return b ? (GREEN + std::string("Да") + RESET)
                 : (RED   + std::string("Нет") + RESET);
    };

    std::cout << BOLD
              << std::setw(26) << "Параметр"
              << std::setw(22) << r1.method_name
              << std::setw(22) << r2.method_name
              << RESET << "\n";
    line('-');

    std::cout << std::setw(26) << "Сошелся"
              << std::setw(22) << yn(r1.converged)
              << std::setw(22) << yn(r2.converged) << "\n";

    std::cout << std::setw(26) << "Итераций"
              << std::setw(22) << r1.iterations
              << std::setw(22) << r2.iterations << "\n";

    std::cout << std::fixed << std::setprecision(15);
    std::cout << std::setw(26) << "f*"
              << std::setw(22) << r1.solution
              << std::setw(22) << r2.solution << "\n";

    std::cout << std::scientific << std::setprecision(4);
    std::cout << std::setw(26) << "|F(f*)|"
              << std::setw(22) << r1.final_F
              << std::setw(22) << r2.final_F << "\n";
    line('-');

    double diff = std::fabs(r1.solution - r2.solution);
    std::cout << "\n  |f1* - f2*| = " << std::scientific << std::setprecision(4) << diff << "\n";

    // 1e-10 - практический порог: при разнице меньше этого значения
    // оба метода считаются давшими одинаковый результат
    if (diff < 1e-10)
        std::cout << GREEN << "  Результаты совпадают.\n" << RESET;
    else
        std::cout << YELLOW << "  Есть расхождение.\n" << RESET;

    std::cout << "\n  Быстрее: " << BOLD;
    if      (r1.iterations < r2.iterations) std::cout << r1.method_name;
    else if (r2.iterations < r1.iterations) std::cout << r2.method_name;
    else                                    std::cout << "Ничья";
    std::cout << RESET << "\n";
    line('=');
}