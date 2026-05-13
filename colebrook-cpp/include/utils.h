#pragma once

#include "colebrook.h"
#include "methods.h"

// Формула Swamee-Jain - явное приближение для f
// Погрешность < 3% при 5e3 < Re < 1e8, 1e-6 < eps/D < 0.05
double initialGuess(const PipeParams& p);

void printProblem(const PipeParams& p, double f0);
void printIterations(const IterResult& res);
void printSummary(const IterResult& res);
void printComparison(const IterResult& r1, const IterResult& r2);