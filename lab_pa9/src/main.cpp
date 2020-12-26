#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>

constexpr double In1_T = 1e-4;
constexpr double In1_A = 10.0;

constexpr double L1 = 1e-3;
constexpr double C1 = 1e-9;
constexpr double C2 = 1e-9;

constexpr double It = 1e-12;
constexpr double Cb = 2e-12;
constexpr double MFt = 0.026;
constexpr double Ru = 1e6;
constexpr double RB = 20.0;

constexpr char GNUPLOT_FILENAME[] = "gnuscript.gnuplot";
constexpr char RES_FILENAME[] = "solution_phi4.dat";

constexpr int FIELD_W = 15;
constexpr double CALC_TIME = 1e-3;
constexpr size_t DIMENSION = 13;
constexpr size_t N_MAX = 10;

constexpr double EPS = 1e-4;
constexpr double eps_max = 1e-3;
constexpr double eps_min = 1e-5;

double t_current = 0;

double dt = 1e-8;
double dt_corrected = dt;
double prev_dt = dt;
double dt_min = dt;

enum Base_vars {
    d_phi1 = 0, d_phi2, d_phi3, d_phi4, int_phi1, int_phi2, int_phi3,
    int_phi4, phi1, phi2, phi3, phi4, Ie,
};

double get_E(double t) {
    return In1_A * sin(2 * M_PI / In1_T * t);
}

void fill_M(std::vector<std::vector<double>> &M, const std::vector<double> &cur_X) {
    double I_t_deriv = It / MFt * std::exp((cur_X.at(phi3) - cur_X.at(phi4)) / MFt);

    M.at(d_phi1).at(d_phi1) = 1.0;
    M.at(d_phi1).at(phi1) = -1.0 / dt;

    M.at(d_phi2).at(d_phi2) = 1;
    M.at(d_phi2).at(phi2) = -1.0 / dt;

    M.at(d_phi3).at(d_phi3) = 1.0;
    M.at(d_phi3).at(phi3) = -1.0 / dt;

    M.at(d_phi4).at(d_phi4) = 1.0;
    M.at(d_phi4).at(phi4) = -1.0 / dt;

    M.at(int_phi1).at(int_phi1) = 1.0;
    M.at(int_phi1).at(phi1) = -dt;

    M.at(int_phi2).at(int_phi2) = 1.0;
    M.at(int_phi2).at(phi2) = -dt;

    M.at(int_phi3).at(int_phi3) = 1.0;
    M.at(int_phi3).at(phi3) = -dt;

    M.at(int_phi4).at(int_phi4) = 1.0;
    M.at(int_phi4).at(phi4) = -dt;

    M.at(phi1).at(int_phi1) = -1.0 / L1;
    M.at(phi1).at(Ie) = -1.0;

    M.at(phi2).at(d_phi2) = C1;
    M.at(phi2).at(phi2) = 1.0 / RB;
    M.at(phi2).at(phi3) = -1.0 / RB;
    M.at(phi2).at(Ie) = 1.0;

    M.at(phi3).at(d_phi3) = Cb;
    M.at(phi3).at(d_phi4) = -Cb;
    M.at(phi3).at(phi2) = -1.0 / RB;
    M.at(phi3).at(phi3) = 1.0 / RB + 1 / Ru + I_t_deriv;
    M.at(phi3).at(phi4) = -1.0 / Ru - I_t_deriv;

    M.at(phi4).at(d_phi3) = -Cb;
    M.at(phi4).at(d_phi4) = Cb + C2;
    M.at(phi4).at(phi3) = -1.0 / Ru - I_t_deriv;
    M.at(phi4).at(phi4) = 1.0 / Ru + I_t_deriv;

    M.at(Ie).at(phi1) = -1.0;
    M.at(Ie).at(phi2) = 1.0;
}

void fill_B(std::vector<double> &B, const std::vector<double> &cur_X, const std::vector<double> &prev_X) {
    double I_equation = It * (std::exp((cur_X.at(phi3) - cur_X.at(phi4)) / MFt) - 1);

    B.at(d_phi1) = -(cur_X.at(d_phi1) - (cur_X.at(phi1) - prev_X.at(phi1)) / dt);
    B.at(d_phi2) = -(cur_X.at(d_phi2) - (cur_X.at(phi2) - prev_X.at(phi2)) / dt);
    B.at(d_phi3) = -(cur_X.at(d_phi3) - (cur_X.at(phi3) - prev_X.at(phi3)) / dt);
    B.at(d_phi4) = -(cur_X.at(d_phi4) - (cur_X.at(phi4) - prev_X.at(phi4)) / dt);

    B.at(int_phi1) = -(cur_X.at(int_phi1) - (prev_X.at(int_phi1) - cur_X.at(phi1) * dt));
    B.at(int_phi2) = -(cur_X.at(int_phi2) - (prev_X.at(int_phi2) - cur_X.at(phi2) * dt));
    B.at(int_phi3) = -(cur_X.at(int_phi3) - (prev_X.at(int_phi3) - cur_X.at(phi3) * dt));
    B.at(int_phi4) = -(cur_X.at(int_phi4) - (prev_X.at(int_phi4) - cur_X.at(phi4) * dt));

    B.at(phi1) = -((-cur_X.at(int_phi1)) / L1 - cur_X.at(Ie));
    B.at(phi2) = -(cur_X.at(Ie) + C1 * cur_X.at(d_phi2) + (cur_X.at(phi2) - cur_X.at(phi3)) / RB);
    B.at(phi3) = -(-(cur_X.at(phi2) - cur_X.at(phi3)) / RB + (cur_X.at(phi3) - cur_X.at(phi4)) / Ru +
            Cb * (cur_X.at(d_phi3) - cur_X.at(d_phi4)) + I_equation);
    B.at(phi4) = -(-(cur_X.at(phi3) - cur_X.at(phi4)) / Ru -
            Cb * (cur_X.at(d_phi3) - cur_X.at(d_phi4)) - I_equation + C2 * cur_X.at(d_phi4));

    B.at(Ie) = -(cur_X.at(phi2) - cur_X.at(phi1) - get_E(t_current));
}

std::vector<double> solve_with_gauss(std::vector<std::vector<double> > &A, std::vector<double> &b) {
    size_t max_ind = 0;
    double max_el = 0.0;

    // Прямой ход метода Гаусса
    for (size_t i = 0; i < DIMENSION; ++i) {
        max_el = std::abs(A.at(i).at(i));
        max_ind = i;

        // Поиск максимального элемента в столбце
        for (size_t k = i + 1; k < DIMENSION; ++k) {
            double cur_el = std::abs(A.at(k).at(i));
            if (cur_el > max_el) {
                max_ind = k;
                max_el = cur_el;
            }
        }

        // Если максимальный элемент ноль, то столбец нулевой
        if (max_el < EPS) {
            throw std::domain_error("Нулевой столбец в матрице.");
        }

        // Перестановка строк
        std::swap(A.at(i), A.at(max_ind));
        std::swap(b.at(i), b.at(max_ind));

        // Матричные преобразования
        double diag = A.at(i).at(i);
        for (size_t j = i; j < DIMENSION; ++j) {
            A.at(i).at(j) /= diag;
        }
        b.at(i) /= diag;

        for (size_t k = i + 1; k < DIMENSION; ++k) {
            double coeff = A.at(k).at(i);

            for (size_t j = i; j < DIMENSION; ++j) {
                A.at(k).at(j) -= A.at(i).at(j) * coeff;
            }
            b.at(k) -= b.at(i) * coeff;
        }
    }

    // Обратный ход метода Гаусса
    std::vector<double> X(DIMENSION);
    for (int i = DIMENSION - 1; i >= 0; --i) {
        X.at(i) = b.at(i);
        for (size_t j = i + 1; j < DIMENSION; ++j) {
            X.at(i) -= X.at(j) * A.at(i).at(j);
        }
    }

    return X;
}

bool is_delta_X_le_EPS(const std::vector<double> &delta_X) {
    for (int i = 0; i < DIMENSION; ++i) {
        if (std::abs(delta_X.at(i)) < EPS) {
            continue;
        } else {
            return false;
        }
    }
    return true;
}

void make_gnuplot_script() {
    std::ofstream gnuscript_file(GNUPLOT_FILENAME);

    gnuscript_file << "set term wxt title 'График зависимости phi4(t)'" << std::endl;
    gnuscript_file << "set key right bottom" << std::endl;
    gnuscript_file << "set grid" << std::endl;
    gnuscript_file << "set xrange[" << 0 << ':' << CALC_TIME << ']' << std::endl;
    gnuscript_file << "plot '" << RES_FILENAME << "' using 1:2 with l lc rgb '#4152d1' lw 1.5 title 'phi4(t)'"
                   << std::endl;
    gnuscript_file << "pause -1" << std::endl;

    system((std::string("gnuplot ") + std::string(GNUPLOT_FILENAME)).c_str());
}

int main() {
    std::ofstream phi4_res_file(RES_FILENAME);

    phi4_res_file.setf(std::ios_base::right);

    std::vector<double> delta_X(DIMENSION);

    std::vector<double> cur_X(DIMENSION);
    std::vector<double> prev_X(DIMENSION);
    std::vector<double> prev_prev_X(DIMENSION);

    std::vector<std::vector<double> > M(DIMENSION, std::vector<double>(DIMENSION));
    std::vector<double> B(DIMENSION);

    std::cout << "Идетрасчет..." << std::endl;

    while (t_current < CALC_TIME) {
        dt = dt_corrected;
        bool is_feasible_delta_X = false;
        size_t n = 0;
        while (!is_feasible_delta_X) {
            // Заполняем матрицу узловых проводимостей нулями
            for (auto &row : M) {
                std::fill(row.begin(), row.end(), 0);
            }

            // Заполняем матрицу узловых проводимсотей и вектор невязок значениями
            fill_M(M, cur_X);
            fill_B(B, cur_X, prev_X);

            // Решаем СЛАУ методом Гаусса и находим вектор поправок
            delta_X = solve_with_gauss(M, B);

            // Делаем приращение на величину вектора поправок
            for (size_t i = 0; i < cur_X.size(); ++i) {
                cur_X.at(i) += delta_X.at(i);
            }

            // Проверяем, что полученные приращения меньше EPS
            is_feasible_delta_X = is_delta_X_le_EPS(delta_X);

            // Продолжаем итерации метода Ньютона, если приращение не удовлетворяет по точности (не меньше EPS)
            if (!is_feasible_delta_X) {
                // Следующая итерация метода Ньютона
                n++;

                // Проверяем, не превышено ли маскимальное число итераций для метода Ньютона
                if (n > N_MAX) {
                    // Если превышено, то отбрасываем текущий шаг и меняем шаг по времени
                    n = 0;
                    dt *= 0.5;
                    cur_X = prev_X;

                    // Проверяем сходимость решения
                    if (dt < dt_min) {
                        throw std::domain_error("Решение не сходится.");
                    }
                }
            }
        }
        // Получено допустимое по точности приращение

        // Рассчитываем дельту
        double cur_delta = 0.0;
        for (int i = 0; i < DIMENSION; ++i) {
            double tmp = 0.5 * dt * ((cur_X.at(i) - prev_X.at(i)) / dt -
                                     (prev_X.at(i) - prev_prev_X.at(i)) / prev_dt);
            cur_delta = (tmp > cur_delta) ? tmp : cur_delta;
        }

        // Если интегрирование неудовлетворительно по точности
        if (cur_delta > eps_max && dt_corrected > dt_min) {
            // Уменьшаем шаг по времени и отбрасываем результаты текущего шага
            dt_corrected *= 0.5;
            cur_X = prev_X;
        } else {
            // Удовлетворительная точность, переход к следующему шагу
            // Сохранение значений с предыдущего шага
            prev_prev_X = prev_X;
            prev_X = cur_X;
            prev_dt = dt;

            // Вывод значения phi_4 на текущем временном шаге
            phi4_res_file << std::setw(FIELD_W) << t_current;
            phi4_res_file << std::setw(FIELD_W) << cur_X.at(phi4) << std::endl;

            // Шаг по времени
            t_current += dt;

            if (cur_delta < eps_min || dt_corrected < dt_min) {
                dt_corrected *= 2;
            } else {
                dt_corrected = dt;
            }
        }

        // Проверяем сходимость решения
        if (dt_corrected < dt_min) {
            throw std::domain_error("Решение не сходится.");
        }
    }

    std::cout << "Конец расчета. Построение графика..." << std::endl;
    make_gnuplot_script();
    std::cout << "Заверешениеработы программы..." << std::endl;

    return 0;
}



