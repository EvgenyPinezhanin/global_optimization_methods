#include <MyMath.h>

// #define DEBUG
#if defined( DEBUG )
    #include <iostream>
#endif

// #define EIGEN
#if defined( EIGEN )
    #include <Eigen/Dense>
#endif

#if defined( MPFR )
    real_number abs_func(real_number x) {
        return mpfr::abs(x, mpfr::mpreal::get_default_rnd());
    }

    real_number sin_func(real_number x) {
        return mpfr::sin(x, mpfr::mpreal::get_default_rnd());
    }

    real_number cos_func(real_number x) {
        return mpfr::cos(x, mpfr::mpreal::get_default_rnd());
    }
#else
    real_number abs_func(real_number x) {
        return std::abs(x);
    }

    real_number sin_func(real_number x) {
        return std::sin(x);
    }

    real_number cos_func(real_number x) {
        return std::cos(x);
    }
#endif

void jordanGaussMethod(const std::vector<std::vector<real_number>> &A, const std::vector<real_number> &B, std::vector<real_number> &X) {
    std::vector<std::vector<real_number>> ATmp(A);
    std::vector<real_number> BTmp(B);
    size_t size = ATmp.size();
    std::vector<size_t> order(size);
    real_number coeff, max_elem;
    size_t max_index_row, max_index_column;

    for (size_t i = 0; i < size; ++i) {
        max_elem = abs_func(ATmp[i][0]);
        for (size_t j = 1; j < size; ++j) {
            if (abs_func(ATmp[i][j]) > max_elem) max_elem = abs_func(ATmp[i][j]);
        }
        if (abs_func(BTmp[i]) > max_elem) max_elem = abs_func(BTmp[i]);
        for (size_t j = 0; j < size; ++j) {
            ATmp[i][j] /= max_elem;
        }
        BTmp[i] /= max_elem;
    }

#if defined( DEBUG )
    std::cout << "SLE after normalization:\n";
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            std::cout << ATmp[i][j] << " ";
        }
        std::cout << BTmp[i] << "\n";
    }
    std::cout << "\n";
#endif

    for (size_t i = 0; i < size; ++i) {
        order[i] = i;
    }

    for (size_t i = 0; i < size; ++i) {
        max_elem = abs_func(ATmp[i][i]);
        max_index_row = i;
        max_index_column = i;
        for (size_t j = i + 1; j < size; ++j) {
            for (size_t k = i + 1; k < size; ++k) {
                if (abs_func(ATmp[j][k]) > max_elem) {
                    max_elem = abs_func(ATmp[j][k]);
                    max_index_row = j;
                    max_index_column = k;
                }
            }
        }
        std::swap(ATmp[i], ATmp[max_index_row]);
        std::swap(BTmp[i], BTmp[max_index_row]);

        for (size_t j = 0; j < size; ++j) {
            std::swap(ATmp[j][i], ATmp[j][max_index_column]);
        }
        std::swap(order[i], order[max_index_column]);

        if (ATmp[i][i] != 0.0) {
            coeff = ATmp[i][i];
            for (size_t k = i; k < size; ++k) {
                ATmp[i][k] /= coeff;
            }
            BTmp[i] /= coeff;
            for (size_t j = i + 1; j < size; ++j) {
                coeff = ATmp[j][i];
                for (size_t k = i; k < size; ++k) {
                    ATmp[j][k] -= coeff * ATmp[i][k];
                }
                BTmp[j] -= coeff * BTmp[i];
            }
        }

    #if defined( DEBUG )
        std::cout << "Iter: " << i << "SLE:\n";
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                std::cout << ATmp[i][j] << " ";
            }
            std::cout << BTmp[i] << "\n";
        }
        std::cout << "\n";
    #endif
    }

    for (size_t i = size - 1; i > 0; --i) {
        for (size_t j = 0; j < i; ++j) {
            BTmp[j] -= ATmp[j][i] * BTmp[i];
            ATmp[j][i] = 0.0;
        }
    }

    for (size_t i = 0; i < size; ++i) {
        X[order[i]] = BTmp[i];
    }
}

#if defined( MPFR )
    void mnk::solve(std::vector<double> &X) const {
        size_t size = X.size();
        std::vector<real_number> X_real(size);
        for (size_t i = 0; i < size; ++i) {
            X_real[i] = X[i];
        }
        solve(X_real);
    }
#endif

void mnk::solve(std::vector<real_number> &X) const {
#if defined( EIGEN )
    size_t size = A.size();
    Eigen::MatrixXd AEigen(size, size);
    Eigen::VectorXd BEigen(size), XEigen(size);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            AEigen(i, j) = A[i][j];
        }
        BEigen(i) = B[i];
    }

    XEigen = AEigen.colPivHouseholderQr().solve(BEigen);

    for (size_t i = 0; i < size; ++i) {
        X[i] = XEigen(i);
    }
#else
    jordanGaussMethod(A, B, X);
#endif
}

double euclideanNormSqr(const std::vector<double> &value) {
    double res = 0.0;
    size_t size = value.size();
    for (int i = 0; i < size; i++) {
        res += value[i] * value[i];
    }
    return res;
}

double euclideanDistance(const std::vector<double> &firstValue, const std::vector<double> &secondValue) {
    double res = 0.0;
    size_t size = firstValue.size();
    for (int i = 0; i < size; i++) {
        res += (firstValue[i] - secondValue[i]) * (firstValue[i] - secondValue[i]);
    }
    return std::sqrt(res);
}

double euclideanDistance(double firstValue, double secondValue) {
    return std::abs(firstValue - secondValue);
}

double chebishevDistance(const std::vector<double> &firstValue, const std::vector<double> &secondValue) {
    double res = std::abs(firstValue[0] - secondValue[0]), resTmp;
    size_t size = firstValue.size();
    for (int i = 1; i < size; i++) {
        resTmp = std::abs(firstValue[i] - secondValue[i]);
        if (resTmp > res) {
            res = resTmp;
        }
    }
    return res;
}

double chebishevDistance(double firstValue, double secondValue) {
    return std::abs(firstValue - secondValue);
}

size_t factorial(size_t num) {
    size_t res = 1;
    for (size_t i = 1; i <= num; ++i) {
        res *= i;
    }
    return res;
}
