#include "../include/methods.h"
#include "../include/colebrook.h"
#include <cmath>

static bool isConverged(double step, double fval, const SolverConfig& cfg)
{
    return std::fabs(step) < cfg.tol_step
        && std::fabs(fval) < cfg.tol_func;
}

// Backtracking: если шаг уводит f <= 0, делим его пополам до 50 раз
// 50 итераций достаточно, чтобы уменьшить шаг в 2^50 ~ 1e15 раз
static double safeStep(double f_curr, double delta)
{
    double f_new = f_curr + delta;
    int tries = 0;
    while (f_new <= 0.0 && tries < 50) {
        delta *= 0.5;
        f_new  = f_curr + delta;
        ++tries;
    }
    return f_new;
}

// Метод Ньютона. Квадратичная сходимость: число верных знаков
// удваивается на каждой итерации вблизи корня
IterResult solveNewton(double f0, const PipeParams& p,
                       const SolverConfig& cfg)
{
    IterResult res;
    res.method_name = "Newton-Raphson";
    res.converged   = false;

    double f = f0;
    res.history.push_back({0, f, F(f, p), -1.0});

    for (int n = 1; n <= cfg.max_iter; ++n) {
        double fval  = F(f, p);
        double dfval = dF(f, p);

        // При |F'| ~ 0 шаг -F/F' становится огромным
        if (std::fabs(dfval) < 1e-30) {
            res.solution    = f;
            res.final_F     = fval;
            res.final_error = 1e30;
            res.iterations  = n - 1;
            return res;
        }

        double delta = -fval / dfval;
        double f_new = safeStep(f, delta);

        if (f_new <= 0.0) {
            res.solution    = f;
            res.final_F     = fval;
            res.final_error = std::fabs(delta);
            res.iterations  = n - 1;
            return res;
        }

        double F_new = F(f_new, p);
        res.history.push_back({n, f_new, F_new, std::fabs(delta)});

        if (isConverged(delta, F_new, cfg)) {
            res.solution    = f_new;
            res.final_F     = F_new;
            res.final_error = std::fabs(delta);
            res.iterations  = n;
            res.converged   = true;
            return res;
        }

        f = f_new;
    }

    res.solution    = f;
    res.final_F     = F(f, p);
    res.final_error = res.history.back().error;
    res.iterations  = cfg.max_iter;
    return res;
}

// Метод Халли. Кубическая сходимость: аппроксимация Паде [1/1] для 1/F,
// что эквивалентно учету кривизны функции через F''
// Формула: delta = -F*F' / (F'^2 - 0.5*F*F'')
IterResult solveHalley(double f0, const PipeParams& p,
                       const SolverConfig& cfg)
{
    IterResult res;
    res.method_name = "Halley";
    res.converged   = false;

    double f = f0;
    res.history.push_back({0, f, F(f, p), -1.0});

    for (int n = 1; n <= cfg.max_iter; ++n) {
        double fval   = F(f,   p);
        double dfval  = dF(f,  p);
        double d2fval = d2F(f, p);

        // Знаменатель Халли: F'^2 - 0.5*F*F''
        double denom = dfval * dfval - 0.5 * fval * d2fval;

        if (std::fabs(denom) < 1e-30) {
            res.solution    = f;
            res.final_F     = fval;
            res.final_error = 1e30;
            res.iterations  = n - 1;
            return res;
        }

        double delta = -(fval * dfval) / denom;
        double f_new = safeStep(f, delta);

        if (f_new <= 0.0) {
            res.solution    = f;
            res.final_F     = fval;
            res.final_error = std::fabs(delta);
            res.iterations  = n - 1;
            return res;
        }

        double F_new = F(f_new, p);
        res.history.push_back({n, f_new, F_new, std::fabs(delta)});

        if (isConverged(delta, F_new, cfg)) {
            res.solution    = f_new;
            res.final_F     = F_new;
            res.final_error = std::fabs(delta);
            res.iterations  = n;
            res.converged   = true;
            return res;
        }

        f = f_new;
    }

    res.solution    = f;
    res.final_F     = F(f, p);
    res.final_error = res.history.back().error;
    res.iterations  = cfg.max_iter;
    return res;
}