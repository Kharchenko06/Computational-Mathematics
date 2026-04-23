#pragma once

/*
Основные структуры данных, используемые во всех модулях

Ссылки:
	Handbook of Floating-Point Arithmetic (Muller et al.)
	§2.6  – Определение ULP
	§4.3 – Fast2Sum / 2Sum (пары ошибок)
	§5.1 – 2MultFMA (пары товаров)
*/

#include <array>
#include <cstddef>

namespace euclidean {

// Точка в 3D пространстве
struct Point3D {
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
};

// DoublePair — результат преобразования без ошибок (EFT)
// 
// Представляет значение V = hi + lo точно, где:
// hi = fl(V) (округленный результат)
// lo = V - hi (ошибка округления, также тип double)
struct DoublePair {
	double hi{ 0.0 };
	double lo{ 0.0 };
};

// CompensatedSum — результат компенсированного суммирования
//
// Представляет S = sum + err[0] + err[1] + ... точно
// (с точностью до порядка выбранного алгоритма)
//
// порядок == 1 → Наивный (ошибка не используется)
// порядок == 2 → KBN-2 / Огита-Оиши (err[0] используется)
// порядок == 3 → KBN-3 (err[0..1] используется)
// порядок == 4 → KBN-4 (err[0..2] используется)
struct CompensatedSum {
	double sum{ 0.0 };
	std::array<double, 3> error{};
	int order{ 1 };
};

// Селекторы методов - используются euclidean() верхнего уровня
enum class MultiplyMethod {
	Naive,   // p = x*x  (стандарт IEEE 754 умножение)
	FMA,     // 2MultFMA: (p,e) так, что p+e = x*x точно
	Ozaki    // Продукт Dekker через сплит Veltkamp
};

enum class SumMethod {
	Naive,        
	OgitaOishi,
	KBN2,
	KBN3,     
	KBN4       
};

}