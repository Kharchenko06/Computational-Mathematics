#pragma once

// Основные структуры данных, используемые во всех модулях.
// Handbook of Floating-Point Arithmetic (Muller et al.):
//   §2.6  — определение ULP
//   §4.3  — Fast2Sum / 2Sum (пары ошибок)
//   §5.1  — 2MultFMA (пары произведений)

#include <array>
#include <cstddef>
#include <vector>

namespace euclidean {

    // Точка в R^n; владеет координатами
    using Point = std::vector<double>;

    // EFT-результат (Error-Free Transformation): V = hi + lo точно.
    // hi = fl(V), lo — ошибка округления.
    struct DoublePair {
        double hi{ 0.0 };
        double lo{ 0.0 };
    };

    // Компенсированная сумма порядка Order.
    // S = sum + error[0] + ... + error[Order-2] точно.
    // Число членов ошибки = Order - 1:
    //   Order=1 — наивное суммирование (ошибка отсутствует)
    //   Order=2 — один член (KBN2, Огита-Оиши)
    //   Order=3 — два члена (KBN3)
    //   Order=4 — три члена (KBN4)
    // Шаблон допускает любой Order >= 1 без изменения кода алгоритмов.
    template<int Order>
    struct CompensatedSum {
        static_assert(Order >= 1, "Order must be >= 1");
        double sum{ 0.0 };
        std::array<double, Order - 1> error{};
    };

    // Селекторы методов для верхнеуровневой euclidean_distance()
    enum class MultiplyMethod {
        Naive,  // x*x без компенсации ошибки
        FMA,    // 2MultFMA: (p,e), p+e = x*x точно (Handbook §5.1)
        Ozaki   // произведение Деккера через сплит Велткампа (Handbook §4.4)
    };

    enum class SumMethod {
        Naive,
        OgitaOishi,
        KBN2,
        KBN3,
        KBN4
    };

}