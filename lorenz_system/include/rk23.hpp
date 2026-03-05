#pragma once
// rk23.hpp - адаптивный метод Богацки-Шампайна (BS23)
// Источник: P. Bogacki, L.F. Shampine (1989)
// "A 3(2) Pair of Runge-Kutta Formulas",
// Applied Mathematics Letters, Vol. 2, No. 4, pp. 321-325, 1989.
// Таблица Бутчера: https://en.wikipedia.org/wiki/Bogacki-Shampine_method

#include "lorenz.hpp"
#include <functional>

struct AdaptiveResult : SolverResult {
    size_t rejected_steps = 0;
};

struct StepOutput {
    State y3;   // решение 3-го порядка
    State y2;   // решение 2-го порядка (для оценки ошибки)
    State k4;   // FSAL: переиспользуется как k1 следующего шага
};

StepOutput rk23_step(
    std::function<State(double, const State&)> f,
    double t, const State& y, const State& k1, double h
);

AdaptiveResult rk23_integrate(
    std::function<State(double, const State&)> f,
    const State& y0, double t0, double t1,
    double rtol = 1e-6, double atol = 1e-9
);