#pragma once

// Эталонный расчет через GMP mpf_t с точностью 256 бит
// 256 выбрано как запас чтобы ошибка была сильно меньше double

#include <span>

namespace euclidean {

	// Возвращает sqrt(sum((b[i]-a[i])^2)) для векторов произвольной длины.
	// Размерности a и b должны совпадать.
	[[nodiscard]] double gmp_euclidean(std::span<const double> a,
		std::span<const double> b);

}