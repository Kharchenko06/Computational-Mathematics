#include "rk38.hpp"
#include <chrono>

State rk38_step(std::function<State(double, const State&)> f,
                double t, const State& y, double h)
{
    State k1 = f(t, y);

    State tmp;
    for (int i = 0; i < 3; ++i) tmp[i] = y[i] + (h/3.0)*k1[i];
    State k2 = f(t + h/3.0, tmp);

    for (int i = 0; i < 3; ++i) tmp[i] = y[i] + h*(-k1[i]/3.0 + k2[i]);
    State k3 = f(t + 2.0*h/3.0, tmp);

    for (int i = 0; i < 3; ++i) tmp[i] = y[i] + h*(k1[i] - k2[i] + k3[i]);
    State k4 = f(t + h, tmp);

    State yn;
    double w = h / 8.0;
    for (int i = 0; i < 3; ++i)
        yn[i] = y[i] + w*(k1[i] + 3.0*k2[i] + 3.0*k3[i] + k4[i]);
    return yn;
}

SolverResult rk38_integrate(std::function<State(double, const State&)> f,
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
        y = rk38_step(f, t, y, dt);
        t += dt;
        res.t.push_back(t);
        res.y.push_back(y);
    }

    auto end = std::chrono::high_resolution_clock::now();
    res.elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
    res.steps = n;
    return res;
}