#pragma once
// rk4.hpp - метод Рунге-Кутты 4-го порядка (фиксированный шаг)
// Источник: Numerical Recipes in C++, 3rd ed., §17.1

#include "lorenz.hpp"
#include <functional>

State rk4_step(
    std::function<State(double, const State&)> f,
    double t, const State& y, double h
);

SolverResult rk4_integrate(
    std::function<State(double, const State&)> f,
    const State& y0, double t0, double t1, double h
);