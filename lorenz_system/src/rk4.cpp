#include "rk4.hpp"
#include <chrono>

State rk4_step(std::function<State(double, const State&)> f,
               double t, const State& y, double h)
{
    // k1 - производная в начале интервала
    State k1 = f(t, y);

    // пробую дойти до середины с половинным шагом, используя k1
    double hh = h * 0.5;
    State mid1;
    mid1[0] = y[0] + hh * k1[0];
    mid1[1] = y[1] + hh * k1[1];
    mid1[2] = y[2] + hh * k1[2];
    State k2 = f(t + hh, mid1);

    // еще одна попытка в середине, теперь с k2
    State mid2;
    mid2[0] = y[0] + hh * k2[0];
    mid2[1] = y[1] + hh * k2[1];
    mid2[2] = y[2] + hh * k2[2];
    State k3 = f(t + hh, mid2);

    // и наконец, пробую шагнуть на целый шаг с k3
    State end;
    end[0] = y[0] + h * k3[0];
    end[1] = y[1] + h * k3[1];
    end[2] = y[2] + h * k3[2];
    State k4 = f(t + h, end);

    // собираю все вместе с весами (классические коэффициенты)
    State yn;
    double w = h / 6.0;
    yn[0] = y[0] + w * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
    yn[1] = y[1] + w * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
    yn[2] = y[2] + w * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);

    return yn;
}

SolverResult rk4_integrate(std::function<State(double, const State&)> f,
                           const State& y0, double t0, double t1, double h)
{
    SolverResult res;
    if (h <= 0.0 || t1 <= t0) {
        return res;   // невалидные параметры, просто верну пустой результат
    }

    // грубо оценю число шагов, чтобы последний шаг не перескочил t1
    int n = int((t1 - t0) / h + 0.5);
    if (n < 1) n = 1;
    double dt = (t1 - t0) / n;   // подгоняю шаг, чтобы точно попасть в t1

    // сразу запишем начальную точку
    res.t.push_back(t0);
    res.y.push_back(y0);

    double t = t0;
    State y = y0;

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