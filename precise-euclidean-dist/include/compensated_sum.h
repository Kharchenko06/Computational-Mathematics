#pragma once

// Пять алгоритмов суммирования n пар DoublePair.
// Handbook of Floating-Point Arithmetic: §4.3.1, §4.3.2, §6.3.2

#include "types.h"

#include <span>

namespace euclidean {

	namespace detail {

	    // fast2sum требует |a| >= |b| (Handbook §4.3.1, Alg. 4.3).
	    // Возвращает (hi, lo): hi + lo = a + b точно.
	    [[nodiscard]] DoublePair fast2sum(double a, double b) noexcept;

	    // two_sum без ограничений на порядок, 6 FP-операций (Handbook §4.3.2, Alg. 4.4).
	    // Возвращает (hi, lo): hi + lo = a + b точно.
	    [[nodiscard]] DoublePair two_sum(double a, double b) noexcept;

	}

    // Наивное суммирование только hi-членов, lo игнорируются.
    [[nodiscard]] CompensatedSum<1> naive_sum(std::span<const DoublePair> terms) noexcept;

    // VecSum по Огита-Рамп-Оиши: суммирует все 2n компонент через каскад two_sum
    // (Handbook §6.3). Выделяет временный вектор — не noexcept.
    [[nodiscard]] CompensatedSum<2> ogita_oishi_sum(std::span<const DoublePair> terms);

    // Kahan-Babuka-Neumaier: исправляет случай |xi| > |s|,
    // которого не умеет оригинальный алгоритм Кахана (Handbook §6.3.2, Alg. 6.6).
    [[nodiscard]] CompensatedSum<2> kbn2_sum(std::span<const DoublePair> terms) noexcept;

    // KBN с двумя уровнями компенсации: two_sum применяется к (sum, error[0]).
    [[nodiscard]] CompensatedSum<3> kbn3_sum(std::span<const DoublePair> terms) noexcept;

    // KBN с тремя уровнями компенсации: три аккумулятора s0, s1, s2.
    [[nodiscard]] CompensatedSum<4> kbn4_sum(std::span<const DoublePair> terms) noexcept;

}