#include "../include/colebrook.h"
#include <cmath>
#include <stdexcept>

PipeParams::PipeParams(double Re_, double eps_, double D_)
    : Re(Re_), eps(eps_), D(D_)
{
    if (Re  <= 0.0) throw std::invalid_argument("Re must be positive");
    if (eps <  0.0) throw std::invalid_argument("eps must be non-negative");
    if (D   <= 0.0) throw std::invalid_argument("D must be positive");

    A = eps / (3.7 * D);
    B = 2.51 / Re;
}

double F(double f, const PipeParams& p)
{
    if (f <= 0.0) return 1e30;

    double sqrtF = std::sqrt(f);
    double arg   = p.A + p.B / sqrtF;

    if (arg <= 0.0) return 1e30;

    return (1.0 / sqrtF) + 2.0 * std::log10(arg);
}

// Вывод F'(f) через правило цепочки:
//   d/df [1/sqrt(f)] = -1 / (2 * f^{3/2})
//   d/df [2*log10(A + B/sqrt(f))] = -B / (ln(10) * f^{3/2} * (A + B/sqrt(f)))
double dF(double f, const PipeParams& p)
{
    if (f <= 0.0) return -1e30;

    double sqrtF = std::sqrt(f);
    double f32   = f * sqrtF;        // f^{3/2}
    double arg   = p.A + p.B / sqrtF;

    if (arg <= 0.0) return -1e30;

    return -1.0 / (2.0 * f32)
           - p.B / (std::log(10.0) * f32 * arg);
}

// F''(f) нужна только для метода Халли. Вывод через дифференцирование F'(f):
//   g  = A + B/sqrt(f)
//   g' = -B / (2 * f^{3/2})
//
//   (f^{3/2} * g)' = (3*A/2)*sqrt(f) + B - используется в числителе term2
double d2F(double f, const PipeParams& p)
{
    if (f <= 0.0) return 1e30;

    double sqrtF = std::sqrt(f);
    double f52   = f * f * sqrtF;    // f^{5/2}
    double f3    = f * f * f;        // f^3
    double g     = p.A + p.B / sqrtF;
    double ln10  = std::log(10.0);

    if (g <= 0.0) return 1e30;

    double term1 = 3.0 / (4.0 * f52);
    double term2 = p.B * (1.5 * p.A * sqrtF + p.B) / (ln10 * f3 * g * g);

    return term1 + term2;
}