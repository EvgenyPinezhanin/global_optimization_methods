#include <iostream>
#include <cmath>

#include <nlopt.h>

double f_rastrigin(unsigned n, const double* x, double* grad, void* f_data) {
    double sum = 0.0;

    for (size_t i = 0; i < n; ++i) {
        sum += x[i] * x[i] - 10.0 * std::cos(2.0 * M_PI * x[i]);
    }
    
    return 10.0 * n + sum;
}

int main() {
    const unsigned n = 4;
    nlopt_algorithm algorithm = NLOPT_GN_ORIG_DIRECT;
    nlopt_opt nlopt_org_direct = nlopt_create(algorithm, n);
    if (nlopt_org_direct == NULL) std::cerr << "Error of creating DIRECT algorithm";

    // double lb[n]{ -5.12, -5.12, -5.12, -5.12 };
    // double ub[n]{  5.12,  5.12,  5.12,  5.12 };
    // nlopt_set_lower_bound(nlopt_org_direct, lb);
    // nlopt_set_upper_bound(nlopt_org_direct, ub);

    double lb = -5.12, ub = 5.12;
    nlopt_set_lower_bounds1(nlopt_org_direct, lb);
    std::cout << "Lower bound: " << lb << "\n";
    nlopt_set_upper_bounds1(nlopt_org_direct, ub);
    std::cout << "Upper bound: " << ub << "\n";

    nlopt_set_min_objective(nlopt_org_direct, f_rastrigin, nullptr);

    double xtol_rel = 1e-4;
    nlopt_set_xtol_rel(nlopt_org_direct, xtol_rel);

    int maxeval = 1000000;
    nlopt_set_maxeval(nlopt_org_direct, maxeval);

    double x[n]{ 0.0, 0.0, 0.0, 0.0 };
    double opt_f;
    nlopt_result result = nlopt_optimize(nlopt_org_direct, x, &opt_f);

    std::cout << "X* = [" << x[0];
    for (size_t i = 0; i < n; ++i) {
        std::cout << ", " << x[i];
    }
    std::cout << "]\n";

    std::cout << "f(X*) = " << opt_f << "\n";

    std::cout << "Result: " << result << "\n";

    std::cout << "Eval: " << nlopt_get_numevals(nlopt_org_direct) << "\n";

    nlopt_destroy(nlopt_org_direct);
    return 0;
}
