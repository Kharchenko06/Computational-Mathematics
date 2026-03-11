#pragma once
// lorenz.hpp - описание системы Лоренца и общие типы
// Источник: E.N. Lorenz (1963), Wikipedia: Lorenz system

#include <array>
#include <vector>

struct LorenzParams {
    double sigma = 10.0;
    double rho = 28.0;
    double beta = 8.0 / 3.0;
};

using State = std::array<double, 3>;

struct SolverResult {
    std::vector<double> t;  // моменты времени
    std::vector<State> y;   // траектория
    double elapsed_ms;      // время счёта (ms)
    size_t steps;           // число шагов
};

// Правая часть: f(t, s) -> ds / dt
State lorenz_rhs(double t, const State& s, const LorenzParams& p);