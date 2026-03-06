// Сравниваю методы Рунге-Кутты для системы Лоренца
// Кое-что списал с Numerical Recipes, кое-что сам допилил, а кое-где накосячил :)

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include "lorenz.hpp"
#include "rk4.hpp"
#include "rk38.hpp"
#include "rk23.hpp"

// интервал интегрирования
static const double START_T = 0.0, END_T = 50.0;
static const double TINY_STEP = 1e-5;      // для эталона (очень мелко, но пусть считает)
static const double STEP4 = 0.01;           // шаг для RK4
static const double STEP38 = 0.01;          // и для RK3/8
static const State INIT_STATE = {1.0, 1.0, 1.0};
static const LorenzParams DEFAULT_PARAMS;   // параметры по умолчанию (sigma=10, rho=28, beta=8/3)

// Евклидово расстояние между двумя состояниями (для оценки ошибки)
double state_dist(const State& a, const State& b) {
    double s = 0.0;
    for (int i = 0; i < 3; ++i) {
        double d = a[i] - b[i];
        s += d * d;
    }
    return std::sqrt(s);
}

// Сохраняем траекторию в CSV (папка data создаётся, если надо)
void save_csv(const std::string& fname, const SolverResult& sol) {
    try {
        std::filesystem::create_directories("data");
        std::ofstream f("data/" + fname);
        if (!f.is_open()) {
            std::cerr << "Упс, не могу создать файл " << fname << std::endl;
            return;
        }
        f << std::fixed << std::setprecision(12);
        f << "t,x,y,z\n";
        for (size_t i = 0; i < sol.t.size(); ++i) {
            f << sol.t[i] << ","
              << sol.y[i][0] << ","
              << sol.y[i][1] << ","
              << sol.y[i][2] << "\n";
        }
        std::cout << "Сохранено: " << fname << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Ошибка записи: " << e.what() << std::endl;
    }
}

// Считаем ошибку относительно эталона (линейная интерполяция, потому что эталон с мелким шагом)
void save_error_csv(const std::string& fname,
                    const SolverResult& ref,
                    const SolverResult& sol)
{
    try {
        std::filesystem::create_directories("data");
        std::ofstream f("data/" + fname);
        if (!f.is_open()) return;
        f << "t,error\n";

        size_t idx = 0; // бегущий индекс в эталоне
        for (size_t i = 0; i < sol.t.size(); ++i) {
            double ti = sol.t[i];
            // двигаем idx, пока эталон не догонит ti
            while (idx + 1 < ref.t.size() && ref.t[idx + 1] <= ti) {
                ++idx;
            }
            if (idx + 1 >= ref.t.size()) break; // дальше не интерполировать

            double dt = ref.t[idx + 1] - ref.t[idx];
            double alpha = (ti - ref.t[idx]) / dt; // вес для правой точки
            if (alpha < 0.0 || alpha > 1.0) {
                // на всякий случай пропускаю (хотя по идее такого не должно быть)
                continue;
            }

            State ry;
            for (int k = 0; k < 3; ++k) {
                ry[k] = ref.y[idx][k] * (1.0 - alpha) + ref.y[idx + 1][k] * alpha;
            }

            double err = state_dist(ry, sol.y[i]);
            f << ti << "," << err << "\n";
        }
    } catch (...) {
        std::cerr << "Что-то пошло не так при сохранении ошибки\n";
    }
}

int main() {
    // правая часть (замыкаю параметры через лямбду)
    auto rhs = [&](double t, const State& s) {
        return lorenz_rhs(t, s, DEFAULT_PARAMS);
    };

    std::cout << ">>> Считаем эталонное решение (RK4, h = " << TINY_STEP << ") ...\n";
    auto ref = rk4_integrate(rhs, INIT_STATE, START_T, END_T, TINY_STEP);
    save_csv("reference.csv", ref);

    std::cout << "\n>>> RK4 с шагом " << STEP4 << "\n";
    auto sol4 = rk4_integrate(rhs, INIT_STATE, START_T, END_T, STEP4);
    save_csv("rk4.csv", sol4);
    save_error_csv("error_rk4.csv", ref, sol4);

    std::cout << "\n>>> RK3/8 с шагом " << STEP38 << "\n";
    auto sol38 = rk38_integrate(rhs, INIT_STATE, START_T, END_T, STEP38);
    save_csv("rk38.csv", sol38);
    save_error_csv("error_rk38.csv", ref, sol38);

    std::cout << "\n>>> Адаптивный RK23 (rtol=1e-7, atol=1e-10)\n";
    auto sol23 = rk23_intergrate(rhs, INIT_STATE, START_T, END_T, 1e-7, 1e-10);
    save_csv("rk23.csv", sol23);
    save_error_csv("error_rk23.csv", ref, sol23);

    // вывожу табличку с результатами
    std::cout << "\n==========  РЕЗУЛЬТАТЫ  ==========\n";
    std::cout << "Метод    |  шагов  | время (мс) | ошибка (max)\n";
    std::cout << "---------|---------|------------|--------------\n";

    // функция для печати строки таблицы
    auto print_row = [&](const std::string& name, const SolverResult& res) {
        // ошибку пока не считываю, просто пишу прочерк
        std::cout << std::left << std::setw(8) << name << " | "
                  << std::setw(7) << res.steps << " | "
                  << std::setw(10) << std::fixed << std::setprecision(3) << res.elapsed_ms << " | "
                  << "(см. файлы)\n";
    };

    print_row("RK4", sol4);
    print_row("RK38", sol38);
    print_row("RK23", sol23);

    std::cout << "\nФайлы с траекториями и ошибками лежат в папке data/\n";
    return 0;
}