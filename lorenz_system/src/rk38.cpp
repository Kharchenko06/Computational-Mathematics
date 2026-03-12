#include "rk38.hpp"
#include <chrono>

// один шаг метода Рунге-Кутты 3/8 (четвертого порядка)
State rk38_step(std::function<State(double, const State&)> f,
                double t, const State& y, double h)
{
    // k - коэффициенты по формуле
    State k1 = f(t, y);

    // первая промежуточная точка (на 1/3 шага)
    State tmp1;
    tmp1[0] = y[0] + (h/3.0) * k1[0];
    tmp1[1] = y[1] + (h/3.0) * k1[1];
    tmp1[2] = y[2] + (h/3.0) * k1[2];
    State k2 = f(t + h/3.0, tmp1);

    // вторая промежуточная точка - комбинация k1 и k2
    State tmp2;
    double coeff = -h/3.0;
    tmp2[0] = y[0] + coeff * k1[0] + h * k2[0];
    tmp2[1] = y[1] + coeff * k1[1] + h * k2[1];
    tmp2[2] = y[2] + coeff * k1[2] + h * k2[2];
    State k3 = f(t + 2.0*h/3.0, tmp2);

    // третья промежуточная точка
    State tmp3;
    tmp3[0] = y[0] + h * (k1[0] - k2[0] + k3[0]);
    tmp3[1] = y[1] + h * (k1[1] - k2[1] + k3[1]);
    tmp3[2] = y[2] + h * (k1[2] - k2[2] + k3[2]);
    State k4 = f(t + h, tmp3);

    // финальное приращение - взвешенная сумма
    State y_next;
    double w = h / 8.0;
    y_next[0] = y[0] + w * (k1[0] + 3.0*k2[0] + 3.0*k3[0] + k4[0]);
    y_next[1] = y[1] + w * (k1[1] + 3.0*k2[1] + 3.0*k3[1] + k4[1]);
    y_next[2] = y[2] + w * (k1[2] + 3.0*k2[2] + 3.0*k3[2] + k4[2]);

    return y_next;
}

// интегрирование от t0 до t1 с шагом h (фиксированный, но последний подгоняется)
SolverResult rk38_integrate(
        std::function<State(double, const State&)> f,
        const State& y0, double t0, double t1, double h)
{
    SolverResult res;
    if (h <= 0.0 || t1 <= t0) {
        return res;
    }

    // прикину число шагов, чтобы не промахнуться мимо t1
    int steps = int((t1 - t0) / h + 0.5);
    if (steps < 1) steps = 1;
    double dt = (t1 - t0) / steps;   // точный шаг

    // запишу начальную точку
    res.t.push_back(t0);
    res.y.push_back(y0);

    double t = t0;
    State y = y0;

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < steps; ++i) {
        y = rk38_step(f, t, y, dt);
        t += dt;
        res.t.push_back(t);
        res.y.push_back(y);
        // шаги считаю автоматически
    }

    auto end = std::chrono::high_resolution_clock::now();
    res.elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
    res.steps = steps;

    return res;
}