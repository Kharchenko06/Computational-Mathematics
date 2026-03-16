#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <map>
#include <string>
#include "lorenz.hpp"
#include "rk4.hpp"
#include "rk38.hpp"
#include "rk23.hpp"

//  Конфигурация

struct Config {
    LorenzParams params;
    State   y0      = {1.0, 1.0, 1.0};
    double  t_start = 0.0;
    double  t_end   = 50.0;
    double  h_rk4   = 0.01;
    double  h_rk38  = 0.01;
    double  rtol    = 1e-7;
    double  atol    = 1e-10;
};

static void print_config(const Config& c) {
    std::cout << "Параметры:\n"
              << "  σ=" << c.params.sigma << "  ρ=" << c.params.rho
              << "  β=" << c.params.beta << "\n"
              << "  y0 = [" << c.y0[0] << ", " << c.y0[1] << ", " << c.y0[2] << "]\n"
              << "  t  ∈ [" << c.t_start << ", " << c.t_end << "]\n"
              << "  h_rk4=" << c.h_rk4 << "  h_rk38=" << c.h_rk38
              << "  rtol=" << c.rtol << "  atol=" << c.atol << "\n";
}

// Сохраняю фактически использованные параметры — Python читает отсюда
static void save_run_params(const Config& c) {
    std::filesystem::create_directories("data");
    std::ofstream f("data/run_params.cfg");
    if (!f.is_open()) return;
    f << std::setprecision(17);
    f << "sigma = "   << c.params.sigma  << "\n"
      << "rho = "     << c.params.rho    << "\n"
      << "beta = "    << c.params.beta   << "\n"
      << "x0 = "      << c.y0[0]         << "\n"
      << "y0 = "      << c.y0[1]         << "\n"
      << "z0 = "      << c.y0[2]         << "\n"
      << "t_start = " << c.t_start       << "\n"
      << "t_end = "   << c.t_end         << "\n"
      << "h_rk4 = "   << c.h_rk4         << "\n"
      << "h_rk38 = "  << c.h_rk38        << "\n"
      << "rtol = "    << c.rtol           << "\n"
      << "atol = "    << c.atol           << "\n";
}

static double ask(const std::string& prompt, double def) {
    std::cout << "  " << prompt << " [" << def << "]: ";
    std::string s;
    std::getline(std::cin, s);
    return s.empty() ? def : std::stod(s);
}

Config from_console() {
    Config c;
    std::cout << "\nВведите параметры (Enter = оставить значение по умолчанию):\n";
    c.params.sigma = ask("sigma",   c.params.sigma);
    c.params.rho   = ask("rho",     c.params.rho);
    c.params.beta  = ask("beta",    c.params.beta);
    c.y0[0]        = ask("x0",      c.y0[0]);
    c.y0[1]        = ask("y0",      c.y0[1]);
    c.y0[2]        = ask("z0",      c.y0[2]);
    c.t_start      = ask("t_start", c.t_start);
    c.t_end        = ask("t_end",   c.t_end);
    c.h_rk4        = ask("h_rk4",   c.h_rk4);
    c.h_rk38       = ask("h_rk38",  c.h_rk38);
    c.rtol         = ask("rtol",    c.rtol);
    c.atol         = ask("atol",    c.atol);
    return c;
}

Config from_file(const std::string& path) {
    Config c;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "Не могу открыть файл " << path << "\n";
        return c;
    }
    std::map<std::string, double> kv;
    std::string line;
    while (std::getline(f, line)) {
        auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        auto trim = [](std::string s) {
            size_t l = s.find_first_not_of(" \t\r\n");
            size_t r = s.find_last_not_of(" \t\r\n");
            return l == std::string::npos ? std::string{} : s.substr(l, r-l+1);
        };
        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));
        if (!key.empty() && !val.empty()) kv[key] = std::stod(val);
    }
    auto get = [&](const std::string& k, double def) {
        auto it = kv.find(k); return it != kv.end() ? it->second : def;
    };
    c.params.sigma = get("sigma",   c.params.sigma);
    c.params.rho   = get("rho",     c.params.rho);
    c.params.beta  = get("beta",    c.params.beta);
    c.y0[0]        = get("x0",      c.y0[0]);
    c.y0[1]        = get("y0",      c.y0[1]);
    c.y0[2]        = get("z0",      c.y0[2]);
    c.t_start      = get("t_start", c.t_start);
    c.t_end        = get("t_end",   c.t_end);
    c.h_rk4        = get("h_rk4",   c.h_rk4);
    c.h_rk38       = get("h_rk38",  c.h_rk38);
    c.rtol         = get("rtol",    c.rtol);
    c.atol         = get("atol",    c.atol);
    return c;
}

//  Метрики ошибки

struct ErrorMetrics { double abs_max, rel_max; };

ErrorMetrics compute_errors(const SolverResult& ref, const SolverResult& sol) {
    double abs_max = 0.0, rel_max = 0.0;
    size_t idx = 0;
    for (size_t i = 0; i < sol.t.size(); ++i) {
        double ti = sol.t[i];
        while (idx + 1 < ref.t.size() && ref.t[idx + 1] <= ti) ++idx;
        if (idx + 1 >= ref.t.size()) break;

        double alpha = (ti - ref.t[idx]) / (ref.t[idx+1] - ref.t[idx]);
        if (alpha < 0.0 || alpha > 1.0) continue;

        State ry; double ref_norm = 0.0;
        for (int k = 0; k < 3; ++k) {
            ry[k] = ref.y[idx][k]*(1.0-alpha) + ref.y[idx+1][k]*alpha;
            ref_norm += ry[k]*ry[k];
        }
        ref_norm = std::sqrt(ref_norm);

        double diff = 0.0;
        for (int k = 0; k < 3; ++k) { double d = sol.y[i][k]-ry[k]; diff += d*d; }
        diff = std::sqrt(diff);

        abs_max = std::max(abs_max, diff);
        if (ref_norm > 1e-15) rel_max = std::max(rel_max, diff/ref_norm);
    }
    return {abs_max, rel_max};
}

void save_csv(const std::string& fname, const SolverResult& sol) {
    std::filesystem::create_directories("data");
    std::ofstream f("data/" + fname);
    if (!f.is_open()) { std::cerr << "Не могу создать " << fname << "\n"; return; }
    f << std::fixed << std::setprecision(12) << "t,x,y,z\n";
    for (size_t i = 0; i < sol.t.size(); ++i)
        f << sol.t[i] << "," << sol.y[i][0] << ","
          << sol.y[i][1] << "," << sol.y[i][2] << "\n";
    std::cout << "  saved: " << fname << "\n";
}

void save_error_csv(const std::string& fname,
                    const SolverResult& ref, const SolverResult& sol)
{
    std::filesystem::create_directories("data");
    std::ofstream f("data/" + fname);
    if (!f.is_open()) return;
    f << "t,abs_error,rel_error\n";

    size_t idx = 0;
    for (size_t i = 0; i < sol.t.size(); ++i) {
        double ti = sol.t[i];
        while (idx + 1 < ref.t.size() && ref.t[idx+1] <= ti) ++idx;
        if (idx + 1 >= ref.t.size()) break;

        double alpha = (ti - ref.t[idx]) / (ref.t[idx+1] - ref.t[idx]);
        if (alpha < 0.0 || alpha > 1.0) continue;

        State ry; double ref_norm = 0.0;
        for (int k = 0; k < 3; ++k) {
            ry[k] = ref.y[idx][k]*(1.0-alpha) + ref.y[idx+1][k]*alpha;
            ref_norm += ry[k]*ry[k];
        }
        ref_norm = std::sqrt(ref_norm);

        double diff = 0.0;
        for (int k = 0; k < 3; ++k) { double d = sol.y[i][k]-ry[k]; diff += d*d; }
        diff = std::sqrt(diff);

        f << ti << "," << diff << ","
          << (ref_norm > 1e-15 ? diff/ref_norm : 0.0) << "\n";
    }
}

int main(int argc, char* argv[]) {
    Config cfg;

    if (argc > 1) {
        cfg = from_file(argv[1]);
        std::cout << "Параметры из файла " << argv[1] << ":\n";
    } else {
        std::cout << "Источник параметров:\n"
                  << "  1 — ввести вручную\n"
                  << "  2 — загрузить из lorenz.cfg\n"
                  << "Выбор: ";
        std::string choice;
        std::getline(std::cin, choice);
        if (choice == "2") {
            cfg = from_file("lorenz.cfg");
            std::cout << "Параметры из lorenz.cfg:\n";
        } else {
            cfg = from_console();
            std::cout << "Параметры:\n";
        }
    }

    print_config(cfg);

    // Сохраняю фактические параметры — Python скрипты читают отсюда
    save_run_params(cfg);

    std::cout << "\n";

    auto rhs = [&](double t, const State& s) {
        return lorenz_rhs(t, s, cfg.params);
    };

    static const double H_REF = 1e-5;
    std::cout << "Эталон (RK4, h=" << H_REF << ")...\n";
    auto ref = rk4_integrate(rhs, cfg.y0, cfg.t_start, cfg.t_end, H_REF);
    save_csv("reference.csv", ref);

    std::cout << "RK4 (h=" << cfg.h_rk4 << ")...\n";
    auto sol4 = rk4_integrate(rhs, cfg.y0, cfg.t_start, cfg.t_end, cfg.h_rk4);
    save_csv("rk4.csv", sol4);
    save_error_csv("error_rk4.csv", ref, sol4);

    std::cout << "RK3/8 (h=" << cfg.h_rk38 << ")...\n";
    auto sol38 = rk38_integrate(rhs, cfg.y0, cfg.t_start, cfg.t_end, cfg.h_rk38);
    save_csv("rk38.csv", sol38);
    save_error_csv("error_rk38.csv", ref, sol38);

    std::cout << "RK23 (rtol=" << cfg.rtol << ")...\n";
    auto sol23 = rk23_integrate(rhs, cfg.y0, cfg.t_start, cfg.t_end, cfg.rtol, cfg.atol);
    save_csv("rk23.csv", sol23);
    save_error_csv("error_rk23.csv", ref, sol23);

    auto m4  = compute_errors(ref, sol4);
    auto m38 = compute_errors(ref, sol38);
    auto m23 = compute_errors(ref, sol23);

    std::cout << "\n========== РЕЗУЛЬТАТЫ ==========\n";
    std::cout << std::left
              << std::setw(7)  << "Метод"    << " | "
              << std::setw(7)  << "шагов"    << " | "
              << std::setw(10) << "время(мс)"<< " | "
              << std::setw(12) << "abs (max)" << " | "
              << "rel (max)\n"
              << std::string(58, '-') << "\n";

    auto row = [](const std::string& name, const SolverResult& s, ErrorMetrics m) {
        std::cout << std::scientific << std::setprecision(3);
        std::cout << std::left << std::setw(7)  << name        << " | "
                  << std::setw(7)  << s.steps                  << " | "
                  << std::setw(12) << s.elapsed_ms             << " | "
                  << std::setw(12) << m.abs_max                << " | "
                  << m.rel_max << "\n";
    };

    row("RK4",  sol4,  m4);
    row("RK38", sol38, m38);
    row("RK23", sol23, m23);

    std::cout << "\nДанные сохранены в data/\n";
    return 0;
}