#pragma once

/* 
Вспомогательные средства для анализа по стандарту IEEE 754 ULP(Unit in the Last Place)

Определение ulp Голдберга:
ulp(x) = 2^(e - p + 1), где e — смещенный показатель степени x, а p = 53 для double

Это определение используется во всей книге и соответствует std::numeric_limits<double>::epsilon() = 2^-52

Ссылки:
	Handbook of Floating-Point Arithmetic §2.6
*/

#include <cmath>
#include <cstdint>
#include <limits>

namespace euclidean {

// ulp(x) — размер одного ULP в x (определение Голдберга)
//
// Возвращает 0 для x == 0 (особый случай: ulp(0) некорректно определен)
// Возвращает NaN / Inf для исключительных значений
[[nodiscard]] double ulp(double x) noexcept;

//	ulp_distance(a, b) — знаковое расстояние от a до b в ULP
//
//	Возвращает количество представимых чисел типа double между a и b,
//	со знаком, указывающим направление
//
//	Special cases:
//    ulp_distance(x, x) == 0
//    ulp_distance(+0, -0) == 0
//    ulp_distance(NaN, any) == INT64_MAX (не определено)
[[nodiscard]] std::int64_t ulp_distance(double a, double b) noexcept;

// ulp_error(computed, reference) — абсолютная ошибка ULP
//
// Возвращает |ulp_distance(computed, reference)|
[[nodiscard]] std::int64_t ulp_error(double computed, double reference) noexcept;

// bits(x) — необработанный 754-битный шаблон IEEE в виде uint64_t
// Полезно для отладки и понимания поведения округления
[[nodiscard]] std::uint64_t bits(double x) noexcept;

// next_up / next_down — перемещение на один ULP в заданном направлении
// Тонкие обертки над std::nextafter
[[nodiscard]] double next_up(double x) noexcept;
[[nodiscard]] double next_down(double x) noexcept;

}
