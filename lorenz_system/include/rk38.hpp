#pragma once
#include "lorenz.hpp"
#include <functional>

State rk38_step(std::function<State(double, const State&)> f,
    double t, const State& y, double h);

SolverResult rk38_integrate(std::function<State(double, const State&)> f,
    const State& y0, double t0, double t1, double h);