#pragma once

// Основные структуры данных, используемые во всех модулях
// Handbook of Floating-Point Arithmetic (Muller et al.):
//   §2.6  - определение ULP
//   §4.3  - Fast2Sum / 2Sum (пары ошибок)
//   §5.1  - 2MultFMA (пары произведений)

#include <array>
#include <cstddef>
#include <vector>   

namespace euclidean {

    // Точка в R^n; владеет координатами
    using Point = std::vector<double>;

    // EFT-результат (Error-Free Transformation): V = hi + lo точно
    // hi = fl(V), lo = ошибка округления
    struct DoublePair {
        double hi{ 0.0 };
        double lo{ 0.0 };
    };

    // Компенсированная сумма: S = sum + err[0] + err[1] + ...
    // Количество учитываемых ошибок зависит от order:
    //   1 - наивное суммирование (err не используется)
    //   2 - KBN-2 / Огита-Оиши (err[0])
    //   3 - KBN-3 (err[0..1])
    //   4 - KBN-4 (err[0..2])
    struct CompensatedSum {
        double sum{ 0.0 };
        std::array<double, 3> error{};
        int order{ 1 };
    };

    // Селекторы методов для верхнеуровневой euclidean()
    enum class MultiplyMethod {
        Naive,  // x*x без компенсации ошибки
        FMA,    // 2MultFMA: (p,e), p+e = x*x точно (Handbook §5.1)
        Ozaki   // произведение Декера через сплит Велткампа (Handbook §4.4)
    };

    enum class SumMethod {
        Naive,
        OgitaOishi,
        KBN2,
        KBN3,
        KBN4
    };

}