#pragma once
#include <array>
#include <vector>

struct LorenzParams {
    double sigma = 10.0;
    double rho = 28.0;
    double beta = 8.0 / 3.0;
};

using State = std::array<double, 3>;

struct SolverResult {
    std::vector<double> t;
    std::vector<State>  y;
    double elapsed_ms = 0.0;
    size_t steps = 0;
};

State lorenz_rhs(double t, const State& s, const LorenzParams& p);