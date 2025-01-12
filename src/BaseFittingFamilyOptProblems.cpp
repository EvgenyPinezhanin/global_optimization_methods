#include <opt_problems/BaseFittingFamilyOptProblems.h>

#include <limits>

double OneDimensionalSupportiveOptProblem::computeObjectiveFunction(const double &x) const {
    real_number result = 0.0;

    size_t omegaSize = omega.size();
    for (size_t i = 0; i < omegaSize; ++i) {
        result += coefficients[2 * i] * sin_func(omega[i] * x) +
                  coefficients[2 * i + 1] * cos_func(omega[i] * x);
    }

#if defined( MPFR )
    return typeProblem == TypeProblem::MIN ? result.toDouble() : -result.toDouble();
#else
    return typeProblem == TypeProblem::MIN ? result : -result;
#endif
}

#if defined( TIKHONOV_REGULARIZATION )
double TikhonovFunctional::computeObjectiveFunction(const std::vector<double> &x) const {
    size_t size = B.size();

    double result = 0.0;
    std::vector<double> residual(size, 0.0);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            residual[i] += A[i][j] * x[j];
        }
        residual[i] -= B[i];
    }

    result = euclideanNormSqr(residual) + alpha * euclideanNormSqr(x);

    return result;
}

double TikhonovFunctional::computeConstraintFunction(const std::vector<double> &x, size_t index) const {
    if (index == 0) {
        return computeObjectiveFunction(x);
    } else {
        return std::numeric_limits<double>::quiet_NaN();
    }
}
#endif
