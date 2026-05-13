#pragma once

#include <string>

// Физические параметры трубопровода
// A и B вычисляются один раз в конструкторе, чтобы не пересчитывать
// на каждой итерации
struct PipeParams {
    double Re;
    double eps;
    double D;

    // A = eps / (3.7 * D)	- вклад шероховатости в аргумент log10
    // B = 2.51 / Re		- вклад вязкости в аргумент log10
    double A;
    double B;

    PipeParams(double Re, double eps, double D);

    double relRoughness() const { return eps / D; }
};

// F(f) = 1/sqrt(f) + 2*log10(A + B/sqrt(f))
// Возвращает невязку уравнения Colebrook-White. Корень F(f)=0 - это решение.
double F(double f, const PipeParams& p);

// Первая производная F по f. Нужна методу Ньютона и Халли.
double dF(double f, const PipeParams& p);

// Вторая производная F по f. Нужна только методу Халли.
double d2F(double f, const PipeParams& p);