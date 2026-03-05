#include "lorenz.hpp"

State lorenz_rhs(double /* время здесь не используется */,
                const State& state, const LorenzParams& params) {
    // Классическая система ЛОренца (1963):
    //  dx/dt = sigma * (y - x)
    //  dy/dt = x * (rho - z) - y
    //  dz/dt = x * y - beta * z
    // Для хаоса обычно берут sigma = 10, rho = 28, beta = 8/3

    // Производные
    double dx = params.sigma * (state[1] - state[0]);
    double dy = state[0] * (params.rho - state[2] - state[1]);
    double dz = state[0] * state[1] - params.beta * state[2];

    // Результат
    State dsdt;
    dsdt[0] = dx;
    dsdt[1] = dy;
    dsdt[2] = dz;
    
    return dsdt;
}