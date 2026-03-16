#include "rk4.hpp"
#include <chrono>

State rk4_step(std::function<State(double, const State&)> f,
               double t, const State& y, double h)
{
    State k1 = f(t, y);

    double hh = h * 0.5;
    State mid1;
    for (int i = 0; i < 3; ++i) mid1[i] = y[i] + hh * k1[i];
    State k2 = f(t + hh, mid1);

    State mid2;
    for (int i = 0; i < 3; ++i) mid2[i] = y[i] + hh * k2[i];
    State k3 = f(t + hh, mid2);

    State end_pt;
    for (int i = 0; i < 3; ++i) end_pt[i] = y[i] + h * k3[i];
    State k4 = f(t + h, end_pt);

    State yn;
    double w = h / 6.0;
    for (int i = 0; i < 3; ++i)
        yn[i] = y[i] + w * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    return yn;
}

SolverResult rk4_integrate(std::function<State(double, const State&)> f,
                           const State& y0, double t0, double t1, double h)
{
    SolverResult res;
    if (h <= 0.0 || t1 <= t0) return res;

    int n = int((t1 - t0) / h + 0.5);
    if (n < 1) n = 1;
    double dt = (t1 - t0) / n;

    res.t.push_back(t0);
    res.y.push_back(y0);

    double t = t0;
    State  y = y0;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i) {
        y = rk4_step(f, t, y, dt);
        t += dt;
        res.t.push_back(t);
        res.y.push_back(y);
    }

    auto end = std::chrono::high_resolution_clock::now();
    res.elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
    res.steps = n;
    return res;
}