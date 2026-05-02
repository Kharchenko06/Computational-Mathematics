#pragma once

// Пять алгоритмов суммирования n пар DoublePair.
// Handbook of Floating-Point Arithmetic: §4.3.1, §4.3.2, §6.3.2

#include "types.h"

#include <span>

namespace euclidean {

	// fast2sum требует |a| >= |b| (Handbook §4.3.1, Alg. 4.3).
	// two_sum без ограничений на порядок, но дороже — 6 FP-операций (Handbook §4.3.2, Alg. 4.4).
	// Оба возвращают (hi, lo): hi + lo = a + b точно.
	[[nodiscard]] DoublePair fast2sum(double a, double b) noexcept;
	[[nodiscard]] DoublePair two_sum(double a, double b) noexcept;

	// Наивное суммирование только hi-членов, lo игнорируются.
	[[nodiscard]] CompensatedSum naive_sum(std::span<const DoublePair> terms) noexcept;

	// VecSum по Огита-Рамп-Оиши: суммирует все 2n компонент (hi и lo)
	// через каскад two_sum, одна компенсация (Handbook §6.3).
	[[nodiscard]] CompensatedSum ogita_oishi_sum(std::span<const DoublePair> terms) noexcept;

	// Kahan-Babuka-Neumaier: исправляет случай когда |xi| > |s|,
	// чего не умеет оригинальный алгоритм Кахана (Handbook §6.3.2, Alg. 6.6).
	[[nodiscard]] CompensatedSum kbn2_sum(std::span<const DoublePair> terms) noexcept;

	// KBN с двумя уровнями компенсации: two_sum применяется к (sum, error[0]).
	[[nodiscard]] CompensatedSum kbn3_sum(std::span<const DoublePair> terms) noexcept;

	// KBN с тремя уровнями компенсации.
	[[nodiscard]] CompensatedSum kbn4_sum(std::span<const DoublePair> terms) noexcept;

}