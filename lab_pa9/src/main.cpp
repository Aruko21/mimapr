#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>

// Input signal settings
constexpr double In1_T = 1e-4;
constexpr double In1_A = 10.0;

// Scheme model settings
constexpr double L1 = 1e-3;
constexpr double C1 = 1e-9;
constexpr double C2 = 1e-9;

// Diode settings
constexpr double It = 1e-12;
constexpr double Cb = 2e-12;
constexpr double MFt = 0.026;
constexpr double Ru = 1e6;
constexpr double RB = 20.0;

// Gnuplot settings
constexpr char GNUPLOT_FILENAME[] = "gnuscript.gnuplot";
constexpr char RES_FILENAME[] = "solution_phi4.dat";
constexpr int RESULTS_DIGITS = 13;

// Dimensions settings
constexpr size_t DIMENSION = 13;
constexpr size_t N_MAX = 10;

// Accuracy settings
constexpr double EPS = 1e-4;
constexpr double eps_max = 1e-3;
constexpr double eps_min = 1e-5;

// Time step settings
constexpr double CALC_TIME = 1e-3;
double t_current = 0;
double dt = 1e-8;
double dt_corrected = dt;
double prev_dt = dt;
double dt_min = dt;

// Indices for M matrix and B vector
enum Indices {
    d_phi1 = 0, d_phi2, d_phi3, d_phi4, int_phi1, int_phi2, int_phi3,
    int_phi4, phi1, phi2, phi3, phi4, Ie,
};

double get_E(double t) {
    // Get the voltage value by time
    return In1_A * sin(2 * M_PI / In1_T * t);
}

void M_fill(std::vector<std::vector<double>> &M, const std::vector<double> &cur_X) {
    // M matrix filling according the math model
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

void B_fill(std::vector<double> &B, const std::vector<double> &cur_X, const std::vector<double> &prev_X) {
    // B vector filling according the math model
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

std::vector<double> gauss(std::vector<std::vector<double> > &A, std::vector<double> &b) {
    size_t max_ind = 0;
    double max_el = 0.0;

    // Forward Gauss substitution
    for (size_t i = 0; i < DIMENSION; ++i) {
        max_el = std::abs(A.at(i).at(i));
        max_ind = i;

        // Find max element in row
        for (size_t k = i + 1; k < DIMENSION; ++k) {
            double cur_el = std::abs(A.at(k).at(i));
            if (cur_el > max_el) {
                max_ind = k;
                max_el = cur_el;
            }
        }

        // If max element is zero, then the whole col is zeroed
        if (max_el < EPS) {
            throw std::domain_error("Zero col found");
        }

        // Rows replacing
        std::swap(A.at(i), A.at(max_ind));
        std::swap(b.at(i), b.at(max_ind));

        // Matrix transformations
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

    // Backward Gauss substitution
    std::vector<double> X(DIMENSION);
    for (int i = DIMENSION - 1; i >= 0; --i) {
        X.at(i) = b.at(i);
        for (size_t j = i + 1; j < DIMENSION; ++j) {
            X.at(i) -= X.at(j) * A.at(i).at(j);
        }
    }

    return X;
}

bool is_dX_converge(const std::vector<double> &delta_X) {
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

    gnuscript_file << "set term wxt title 'phi4(t) plot'" << std::endl;
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

    std::cout << "Calculaing started" << std::endl;

    while (t_current < CALC_TIME) {
        dt = dt_corrected;
        bool is_feasible_delta_X = false;
        size_t n = 0;

        // Getting accurate augment
        while (!is_feasible_delta_X) {
            // Matrix initialize
            for (auto &row : M) {
                std::fill(row.begin(), row.end(), 0);
            }

            // Matrix and vector filling
            M_fill(M, cur_X);
            B_fill(B, cur_X, prev_X);

            // Find vector of amendments
            delta_X = gauss(M, B);

            // Augmentation on the amendment value
            for (size_t i = 0; i < cur_X.size(); ++i) {
                cur_X.at(i) += delta_X.at(i);
            }

            // Check the augment accuracy
            is_feasible_delta_X = is_dX_converge(delta_X);

            // If augment not accurate as EPS, then do another Newton iteration
            if (!is_feasible_delta_X) {
                n++;

                // Check the number of iterations for step changing
                if (n > N_MAX) {
                    // Change time step
                    n = 0;
                    dt *= 0.5;
                    cur_X = prev_X;

                    // Check the convergence of solution
                    if (dt < dt_min) {
                        std::cerr << "Solution doesn't converge" << std::endl;
                        return 1;
                    }
                }
            }
        }

        // Delta calculating
        double cur_delta = 0.0;
        for (int i = 0; i < DIMENSION; ++i) {
            double tmp = 0.5 * dt * ((cur_X.at(i) - prev_X.at(i)) / dt -
                                     (prev_X.at(i) - prev_prev_X.at(i)) / prev_dt);
            cur_delta = (tmp > cur_delta) ? tmp : cur_delta;
        }

        // Check accuracy of integration
        if (cur_delta > eps_max && dt_corrected > dt_min) {
            // Reduce time step and rewrite previous step results
            dt_corrected *= 0.5;
            cur_X = prev_X;
        } else {
            // If accuracy is satisfactory, go to the next step
            // Save values from previous step
            prev_prev_X = prev_X;
            prev_X = cur_X;
            prev_dt = dt;

            // Print phi4 value on the current step with fixed number of digits
            phi4_res_file << std::setw(RESULTS_DIGITS) << t_current;
            phi4_res_file << std::setw(RESULTS_DIGITS) << cur_X.at(phi4) << std::endl;

            // Time step
            t_current += dt;

            if (cur_delta < eps_min || dt_corrected < dt_min) {
                dt_corrected *= 2;
            } else {
                dt_corrected = dt;
            }
        }

        // Check solution convergence
        if (dt_corrected < dt_min) {
            std::cerr << "Solution doesn't converge" << std::endl;
            return 1;
        }
    }

    std::cout << "Solution is done. Plotting" << std::endl;
    make_gnuplot_script();

    return 0;
}



