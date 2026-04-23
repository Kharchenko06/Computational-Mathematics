#pragma once

// Эталонный расчет через GMP mpf_t с точностью 256 бит
// 256 выбрано как запас чтобы ошибка была сильно меньше double

#include "types.h"

namespace euclidean {

	// Возвращает результат sqrt(dx^2+dy^2+dz^2)
	// используется как эталон высокой точности
	[[nodiscard]] double gmp_euclidean(const Point3D& p1, const Point3D& p2);

}