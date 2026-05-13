#pragma once

#include "colebrook.h"
#include <vector>
#include <string>

// Одна строка в таблице итераций
struct IterStep {
    int iter;
    double f;
    double F_val;
    // |f_new - f_old|. Для начальной точки (iter=0) равен -1
    double error;
};

struct IterResult {
    std::string method_name;
    double solution;
    double final_F;
    double final_error;
    int iterations;
    bool converged;
    std::vector<IterStep> history;
};

struct SolverConfig {
    double tol_step = 1e-12;
    double tol_func = 1e-12;
    int max_iter = 200;
};

// Метод Ньютона: f_{n+1} = f_n - F(f_n) / F'(f_n)
// Порядок сходимости 2. Требует аналитической F'(f)
IterResult solveNewton(double f0, const PipeParams&   p, const SolverConfig& cfg = SolverConfig{});

// Метод Халли: f_{n+1} = f_n - F*F' / (F'^2 - 0.5*F*F'')
// Порядок сходимости 3. Требует F'(f) и F''(f)
IterResult solveHalley(double f0, const PipeParams&   p, const SolverConfig& cfg = SolverConfig{});