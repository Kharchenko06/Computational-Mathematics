#include "lorenz.hpp"

State lorenz_rhs(double, const State& s, const LorenzParams& p) {
    return {
        p.sigma * (s[1] - s[0]),
        s[0] * (p.rho - s[2]) - s[1],
        s[0] * s[1] - p.beta * s[2]
    };
}