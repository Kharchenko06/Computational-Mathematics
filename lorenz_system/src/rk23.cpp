/* Ключевая идея: на каждом шаге вычисляются два приближения - 3-го и 2-го порядка.
   Их разность даёт оценку локальной ошибки. Если ошибка слишком велика - шаг
   уменьшается и повторяется. Если ошибка мала - шаг можно увеличить */

#include "rk23.hpp"
#include <cmath>
#include <ctime>
#include <algorithm>
#include <functional>

// Норма ошибки: sqrt( mean( (e_i / sc_i)^2 ) )
// где sc_i = atol + rtol * max(|y_i|, |yn_i|)
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

// Один шаг метода Богакки–Шампина (FSAL: k4 первого шага становится k1 следующего)
// Возвращает два приближения (3-го и 2-го порядка) и k4 для следующего шага
StepOutput rk23_step(std::function<State(double, const State&)> f,
                     double t, const State& y,
                     const State& k1, double h)
{
    // Вычисляем k2 (середина шага)
    State y2;
    double half = h * 0.5;
    for (int i = 0; i < 3; ++i)
        y2[i] = y[i] + half * k1[i];
    State k2 = f(t + half, y2);

    // Вычисляем k3 (на 3/4 шага)
    State y3;
    double three_quarter = h * 0.75;
    y3[0] = y[0] + three_quarter * k2[0];
    y3[1] = y[1] + three_quarter * k2[1];
    y3[2] = y[2] + three_quarter * k2[2];
    State k3 = f(t + three_quarter, y3);

    // Приближение 3-го порядка (основное решение) - коэффициенты 2/9, 1/3, 4/9
    double w1 = (2.0/9.0) * h, w2 = (1.0/3.0) * h, w3 = (4.0/9.0) * h;
    State yn3;
    yn3[0] = y[0] + w1 * k1[0] + w2 * k2[0] + w3 * k3[0];
    yn3[1] = y[1] + w1 * k1[1] + w2 * k2[1] + w3 * k3[1];
    yn3[2] = y[2] + w1 * k1[2] + w2 * k2[2] + w3 * k3[2];

    // k4 - производная в конце шага (понадобится для следующего)
    State k4 = f(t + h, yn3);

    // Приближение 2-го порядка (для оценки ошибки) - коэффициенты 7/24, 1/4, 1/3, 1/8
    double v1 = (7.0/24.0) * h, v2 = 0.25 * h, v3 = (1.0/3.0) * h, v4 = 0.125 * h;
    State yn2;
    for (int i = 0; i < 3; ++i) {
        yn2[i] = y[i] + v1*k1[i] + v2*k2[i] + v3*k3[i] + v4*k4[i];
    }

    return {yn3, yn2, k4};
}

// Адаптивное интегрирование от t0 до t1 с заданными точностями
AdaptiveResult rk23_integrate(std::function<State(double, const State&)> f,
                              const State& y0, double t0, double t1,
                              double rtol, double atol)
{
    AdaptiveResult res;
    if (t1 <= t0 || rtol <= 0.0 || atol <= 0.0) {
        return res;     // невалидные параметры, верну пустой результат
    }

    double t = t0;
    State y = y0;
    double h = (t1 - t0) * 1e-3;    // начальный шаг (очень маленький)
    double h_max = (t1 - t0) * 0.1; // ограничу максимальный шаг, чтобы не улететь

    res.t.push_back(t);
    res.y.push_back(y);

    State k1 = f(t, y);             // первая производная (FSAL)

    clock_t start = clock();

    while (t < t1) {
        // последний шаг не должен выходить за t1
        if (t + h > t1) {
            h = t1 - t;
        }

        StepOutput out = rk23_step(f, t, y, k1, h);

        // Ошибка как разность двух приближений
        State err;
        for (int i = 0; i < 3; ++i) {
            err[i] = out.y3[i] - out.y2[i];
        }

        double e = error_norm(err, y, out.y3, rtol, atol);

        if (e <= 1.0) {
            // шаг принят
            t += h;
            y = out.y3;
            k1 = out.k4;            // FSAL: k4 следующего шага
            res.t.push_back(t);
            res.y.push_back(y);
            res.steps++;
        } else {
            res.rejected_steps++;
        }

        // Формула изменения шага (Хайрер, Нерсетт, Ваннер)
        double fac = std::min(std::max(
            0.9 * std::pow(1.0 / (e + 1e-15), 1.0/3.0),
            0.1), 5.0);
        h = std::min(h * fac, h_max);
    }

    clock_t end = clock();
    res.elapsed_ms = double(end - start) * 1000.0 / CLOCKS_PER_SEC;

    return res;
}