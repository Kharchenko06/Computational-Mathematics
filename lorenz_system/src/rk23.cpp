#include "rk23.hpp"
#include <cmath>
#include <chrono>
#include <algorithm>

// Таблица Бутчера BS23 (Bogacki-Shampine):
//  c  |  a
// ----+----------------------------
//  0  |  0
// 1/2 | 1/2
// 3/4 |  0   3/4
//  1  | 2/9  1/3  4/9
// ----+----------------------------
//  y3 | 2/9  1/3  4/9   0     (порядок 3, основное решение)
//  y2 | 7/24 1/4  1/3  1/8    (порядок 2, для оценки ошибки)

static constexpr double C2  = 1.0/2.0, C3  = 3.0/4.0;
static constexpr double A21 = 1.0/2.0, A32 = 3.0/4.0;
static constexpr double A41 = 2.0/9.0, A42 = 1.0/3.0, A43 = 4.0/9.0;
static constexpr double B1  = 2.0/9.0, B2  = 1.0/3.0, B3  = 4.0/9.0;
static constexpr double E1  = 7.0/24.0, E2 = 1.0/4.0, E3  = 1.0/3.0, E4 = 1.0/8.0;

// WRMS-норма: e_norm <= 1  <=>  шаг укладывается в допуск
static double error_norm(const State& err, const State& y,
                         const State& yn, double rtol, double atol)
{
    double sum = 0.0;
    for (int i = 0; i < 3; ++i) {
        double sc = atol + rtol * std::max(std::abs(y[i]), std::abs(yn[i]));
        sum += (err[i] / sc) * (err[i] / sc);
    }
    return std::sqrt(sum / 3.0);
}

StepOutput rk23_step(std::function<State(double, const State&)> f,
                     double t, const State& y, const State& k1, double h)
{
    State tmp;

    for (int i = 0; i < 3; ++i) tmp[i] = y[i] + h*A21*k1[i];
    State k2 = f(t + C2*h, tmp);

    for (int i = 0; i < 3; ++i) tmp[i] = y[i] + h*A32*k2[i];
    State k3 = f(t + C3*h, tmp);

    State yn3;
    for (int i = 0; i < 3; ++i)
        yn3[i] = y[i] + h*(B1*k1[i] + B2*k2[i] + B3*k3[i]);

    State k4 = f(t + h, yn3);  // FSAL

    State yn2;
    for (int i = 0; i < 3; ++i)
        yn2[i] = y[i] + h*(E1*k1[i] + E2*k2[i] + E3*k3[i] + E4*k4[i]);

    return {yn3, yn2, k4};
}

AdaptiveResult rk23_integrate(std::function<State(double, const State&)> f,
                              const State& y0, double t0, double t1,
                              double rtol, double atol)
{
    AdaptiveResult res;
    if (t1 <= t0 || rtol <= 0.0 || atol <= 0.0) return res;

    double t     = t0;
    State  y     = y0;
    double h     = (t1 - t0) * 1e-3;
    double h_max = (t1 - t0) * 0.1;

    res.t.push_back(t);
    res.y.push_back(y);
    State k1 = f(t, y);

    auto start = std::chrono::high_resolution_clock::now();

    while (t < t1) {
        if (t + h > t1) h = t1 - t;

        StepOutput out = rk23_step(f, t, y, k1, h);

        State err;
        for (int i = 0; i < 3; ++i) err[i] = out.y3[i] - out.y2[i];
        double e = error_norm(err, y, out.y3, rtol, atol);

        if (e <= 1.0) {
            t  += h;
            y   = out.y3;
            k1  = out.k4;   // FSAL
            res.t.push_back(t);
            res.y.push_back(y);
            res.steps++;
        } else {
            res.rejected_steps++;
        }

        // h_new = h * 0.9 * (1/e)^(1/3)
        // Показатель 1/3 = 1/(p+1), p=2 — порядок схемы оценки ошибки.
        double fac = std::min(std::max(0.9 * std::pow(1.0/(e+1e-15), 1.0/3.0), 0.1), 5.0);
        h = std::min(h * fac, h_max);
    }

    auto end = std::chrono::high_resolution_clock::now();
    res.elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();
    return res;
}