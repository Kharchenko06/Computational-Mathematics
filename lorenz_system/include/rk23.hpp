#pragma once
// Bogacki, Shampine — "A 3(2) Pair of Runge-Kutta Formulas", AML 2(4), 1989.
#include "lorenz.hpp"
#include <functional>

struct AdaptiveResult : SolverResult {
    size_t rejected_steps = 0;
};

struct StepOutput {
    State y3;  // решение 3-го порядка
    State y2;  // решение 2-го порядка (для оценки ошибки)
    State k4;  // FSAL: k4 этого шага = k1 следующего
};

StepOutput rk23_step(std::function<State(double, const State&)> f,
    double t, const State& y, const State& k1, double h);

AdaptiveResult rk23_integrate(std::function<State(double, const State&)> f,
    const State& y0, double t0, double t1, double rtol = 1e-6, double atol = 1e-9);