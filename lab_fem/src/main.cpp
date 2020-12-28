#include <iostream>
#include <vector>
#include <cmath>

#define CUBEx
constexpr double EPS = 1e-16;

constexpr double X_BEGIN = 0.0;
constexpr double X_END = 90.0;

constexpr double FIRST_BC = 0;
constexpr double SECOND_BC = 2;

constexpr size_t LINEAR_DIM = 2;
constexpr size_t CUBE_DIM = 4;

constexpr size_t ELEMS_NUM = 20;
constexpr double L = (X_END - X_BEGIN) / ELEMS_NUM;

std::vector<double> gauss(std::vector<std::vector<double> > &A, std::vector<double> &b) {
    size_t row_size = A.size();
    size_t col_size = A.back().size();

    // Forward Gauss substitution
    double pivot = 0.0;
    for (size_t i = 0; i < row_size; ++i) {
        for (size_t j = i + 1; j < col_size; ++j) {
            if (std::abs(A.at(j).at(i)) < EPS) {
                continue;
            }

            pivot = A.at(j).at(i) / A.at(i).at(i);
            b.at(j) -= pivot * b.at(i);
            for (size_t k = 0; k < row_size; ++k) {
                A.at(j).at(k) -= pivot * A.at(i).at(k);
            }
        }
    }

    // Backward Gauss substitution
    std::vector<double> x(row_size);
    for (int i = row_size - 1; i >= 0; --i) {
        x.at(i) = b.at(i);
        for (size_t j = i + 1; j < row_size; ++j) {
            x.at(i) -= x.at(j) * A.at(i).at(j);
        }

        x.at(i) /= A.at(i).at(i);
    }

    return x;
}

double analytical_solution(double x) {
    double sqr2_71 = sqrt(2. / 71);
    return ((2 * exp(-3 * sqr2_71 * x) * (exp(3 * sqr2_71 * x) - 1) *
            (5 * exp(3 * sqr2_71 * x) - 2 * exp(3 * sqr2_71 * (x + 90)) - 2 *
            exp(270 * sqr2_71) + 5 * exp(540 * sqr2_71))) / (3 * (exp(540 * sqr2_71) - 1))
    );
}

std::vector<double> build_analytical_solution(std::vector<double> &x_vec) {
    size_t x_vec_size = x_vec.size();

    std::vector<double> y_vec = std::vector<double>(x_vec_size);

    for (size_t i = 0; i < x_vec_size; i++) {
        y_vec.at(i) = analytical_solution(x_vec.at(i));
    }

    return y_vec;
}

std::vector<double> build_linear_solution(size_t elems_num) {
    double L = (X_END - X_BEGIN) / elems_num;

    size_t size = elems_num + 1;
    std::vector<std::vector<double>> A(size, std::vector<double>(size));
    std::vector<double> b(size);

    // Local stiffness matrix for linear FE
    std::vector<std::vector<double> > local_matrix = {
            {-71.0 / L - 6.0 * L, 71.0 / L - 3.0 * L},
            {71.0 / L - 3.0 * L,  -71.0 / L - 6.0 * L}
    };

    // Ansambling and getting global stiffness matrix for linear FE
    for (size_t i = 0; i < elems_num; ++i) {
        for (size_t j = 0; j < LINEAR_DIM; ++j) {
            for (size_t k = 0; k < LINEAR_DIM; ++k) {
                A.at(i + j).at(i + k) += local_matrix.at(j).at(k);
            }
        }
    }

    // Border Conditions
    A.at(0).at(0) = -71.0;
    A.at(1).at(0) = 0.0;
    A.at(size - 1).at(size - 1) = 71.0;
    A.at(size - 2).at(size - 1) = 0.0;

    A.at(size / 2).at(size / 2 - 1) = 0.;
    A.at(size / 2).at(size / 2) = 1.;
    A.at(size / 2).at(size / 2 + 1) = 0.;

    for (size_t i = 2; i < size - 2; ++i) {
        b.at(i) = -60.0 * L;
    }

    b.at(0) = -30.0 * L - FIRST_BC * local_matrix.at(0).at(0);
    b.at(1) = -60.0 * L - FIRST_BC * local_matrix.at(LINEAR_DIM - 1).at(0);
    b.at(size - 2) = -60 * L - SECOND_BC * local_matrix.at(0).at(LINEAR_DIM - 1);
    b.at(size - 1) = -30.0 * L - SECOND_BC * local_matrix.at(LINEAR_DIM - 1).at(LINEAR_DIM - 1);
    b.at(size / 2) = 0.;

    // Solving
    std::vector<double> res = gauss(A, b);
    res.at(0) = FIRST_BC;
    res.at(size - 1) = SECOND_BC;

    return res;
}

std::vector<double> build_cube_solution(size_t elems_num) {
    double L = (X_END - X_BEGIN) / elems_num;

    size_t size = elems_num + 1;
    std::vector<std::vector<double> > A(size, std::vector<double>(size));
    std::vector<double> b(size);

    // Local stiffness matrix for cube FE
    std::vector<std::vector<double> > local_matrix = {
            {
                -2627. / (10. * L) - 48. * L / 35., 13419. / (40. * L) - 297. * L / 280.,
                -1917. / (20. * L) + 27. * L / 70., 923. / (40. * L) - 57. * L / 280.
            },
            {
                13419. / (40. * L) - 297. * L / 280., -3834. / (5. * L) - 243. * L / 35.,
                21087. / (40. * L) + 243. * L / 280., -1917. / (20. * L) + 27. * L / 70.
            },
            {
                -1917. / (20. * L) + 27. * L / 70., 21087. / (40. * L) + 243. * L / 280.,
                -3834. / (5. * L) - 243. * L / 35., 13419. / (40. * L) - 297. * L / 280.
            },
            {
                923. / (40. * L) - 57. * L / 280., -1917. / (20. * L) + 27. * L / 70.,
                13419. / (40. * L) - 297. * L / 280., -2627. / (10. * L) - 48. * L / 35.
            }
    };

    // Local vector of loads
    std::vector<double> local_b = {
        -15.0 * L / 2.,
        -45.0 * L / 2.,
        -45.0 * L / 2.,
        -15.0 * L / 2.
    };

    // Zeroing matrix elements, which refer to internal nodes
    double pivot = 0.0;
    for (size_t i = 1; i < CUBE_DIM - 1; ++i) {
        for (size_t j = 0; j < CUBE_DIM; ++j) {
            if (std::abs(local_matrix.at(j).at(i)) < EPS) {
                continue;
            }

            if (i != j) {
                pivot = local_matrix.at(j).at(i) / local_matrix.at(i).at(i);
                local_b.at(j) -= pivot * local_b.at(i);
                for (size_t k = 0; k < CUBE_DIM; k++) {
                    local_matrix.at(j).at(k) -= pivot * local_matrix.at(i).at(k);
                }
            }
        }
    }

    // Exclude internal nodes
    std::vector<std::vector<double>> local_matrix_mod = {
        {
            local_matrix.at(0).at(0), local_matrix.at(0).at(CUBE_DIM - 1)
        },
        {
            local_matrix.at(CUBE_DIM - 1).at(0), local_matrix.at(CUBE_DIM - 1).at(CUBE_DIM - 1)
        }
    };

    std::vector<double> local_b_mod = {
        local_b.at(0), local_b.at(CUBE_DIM - 1)
    };

    // Ansambling and getting global stiffness matrix for cube FE
    for (size_t i = 0; i < elems_num; ++i) {
        for (size_t j = 0; j < LINEAR_DIM; ++j) {
            for (size_t k = 0; k < LINEAR_DIM; ++k) {
                A.at(i + j).at(i + k) += local_matrix_mod.at(j).at(k);
            }
        }
    }

    // Border Conditions
    A.at(0).at(0) = -71.0;
    A.at(1).at(0) = 0.0;
    A.at(size - 2).at(size - 1) = 0.0;
    A.at(size - 1).at(size - 1) = 71.0;

    b.at(0) = local_b_mod.at(0) - FIRST_BC * local_matrix_mod.at(0).at(0);
    b.at(1) = local_b_mod.at(0) + local_b_mod.at(LINEAR_DIM - 1) - FIRST_BC * local_matrix_mod.at(LINEAR_DIM - 1).at(0);
    b.at(size - 2) = local_b_mod.at(0) + local_b_mod.at(LINEAR_DIM - 1) - SECOND_BC * local_matrix_mod.at(0).at(LINEAR_DIM - 1);
    b.at(size - 1) = local_b_mod.at(LINEAR_DIM - 1) - SECOND_BC * local_matrix_mod.at(LINEAR_DIM - 1).at(LINEAR_DIM - 1);

    for (size_t i = 2; i < size - 2; ++i) {
        b.at(i) = local_b_mod.at(LINEAR_DIM - 1) + local_b_mod.at(0);
    }

    // Solving
    std::vector<double> res = gauss(A, b);
    res.at(0) = FIRST_BC;
    res.at(size - 1) = SECOND_BC;

    return res;
}

double calc_abs_error(const std::vector<double> &y_real, const std::vector<double> &y) {
    double max_err = 0.0;
    for (size_t i = 0; i < y_real.size(); i++) {
        double err = std::abs(y_real.at(i) - y.at(i));

        if (err > max_err) {
            max_err = err;
        }
    }

    return max_err;
}

int main() {
    // Generating evenly spaced nodes
    std::vector<double> x(ELEMS_NUM + 1);
    for (size_t i = 0; i < x.size(); i++) {
        x.at(i) = X_BEGIN + i * L;
    }
    size_t x_size = x.size();

    // Solving
#ifdef CUBE
    std::vector<double> y = build_cube_solution(ELEMS_NUM);
#else
    std::vector<double> y = build_linear_solution(ELEMS_NUM);
#endif
    std::vector<double> y_real = build_analytical_solution(x);

    // Plotting
    FILE *gp = popen("gnuplot -persist", "w");

    fprintf(gp, "$predict << EOD\n");
    for (size_t i = 0; i < x_size; i++) {
        fprintf(gp, "%lf %lf\n", x.at(i), y.at(i));
    }

    fprintf(gp, "EOD\n");
    fprintf(gp, "$real << EOD\n");
    for (size_t i = 0; i < x_size; i++) {
        fprintf(gp, "%lf %lf\n", x.at(i), y_real.at(i));
    }

    fprintf(stdout, "| i  |   X    |   U_real    |    U_FEM    |  Abs error  |\n");
    for (size_t i = 0; i < x_size; i++) {
        fprintf(stdout, "| %2lu | %6.3lf | %6.5e | %6.5e | %6.5e |\n", i, x.at(i), y_real.at(i), y.at(i), std::abs(y.at(i) - y_real.at(i)));
    }

    fprintf(gp, "EOD\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "plot '$predict' using 1:2 with lp lc 'blue' lw 1.5 pt 7 ps 0.5 title 'FEM solution(%zu elements)',"
                "'$real' using 1:2 with lines lc rgb '#23ca5a' lt 1 lw 2 title 'Analytical solution(%zu elements)',\n",
            ELEMS_NUM, ELEMS_NUM
    );
    fprintf(stdout, "Max absolute error: %lf\n", calc_abs_error(y_real, y));

    return 0;
}
