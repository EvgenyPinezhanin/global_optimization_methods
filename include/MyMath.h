#ifndef _MY_MATH_H_
#define _MY_MATH_H_

#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#define MPFR 
#if defined( MPFR )
    #include <mpreal.h>
#endif

#include <vector>

using Matrix = std::vector<std::vector<double>>;

struct point {
    std::vector<double> x;

    point() : x() {};
    point(double _x) : x(_x) {};
    point(double _x, double _y) : x({ _x, _y }) {};
    point(std::vector<double> _x) : x(_x) {};
    point(std::vector<double> _x, double _y) : x(_x) { x.push_back(_y); };

    size_t getDimension() const { return x.size(); };
};

#if defined( MPFR )
    using real_number = mpfr::mpreal;
#else
    using real_number = double;
#endif

real_number abs_func(real_number x);
real_number sin_func(real_number x);
real_number cos_func(real_number x);

void jordanGaussMethod(const std::vector<std::vector<real_number>> &A, const std::vector<real_number> &B, std::vector<double> &X);

class mnk {
private:
    std::vector<std::vector<real_number>> A;
    std::vector<real_number> B;

public:
    mnk() : A(), B() {};
    mnk(const std::vector<std::vector<real_number>>& _A, const std::vector<real_number>& _B)
        : A(_A), B(_B) {};

#if defined( MPFR )
    mnk(const std::vector<std::vector<double>>& _A, const std::vector<double>& _B)
        : A(_A.size(), std::vector<real_number>(_A.size())), B(_B.size())
    {
        size_t size = _A.size();
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                A[i][j] = _A[i][j];
            }
            B[i] = _B[i];
        }
    };
#endif

    void setA(const std::vector<std::vector<real_number>>& _A) { A = _A; };
    void setB(const std::vector<real_number>& _B) { B = _B; };

#if defined( MPFR )
    void solve(std::vector<double> &X) const;
#endif
    void solve(std::vector<real_number> &X) const;
};

double euclideanNormSqr(const std::vector<double> &value);

double euclideanDistance(const std::vector<double> &firstValue, const std::vector<double> &secondValue);
double euclideanDistance(double firstValue, double secondValue);

double chebishevDistance(const std::vector<double> &firstValue, const std::vector<double> &secondValue);
double chebishevDistance(double firstValue, double secondValue);

template<typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

size_t factorial(size_t num);

#endif // _MY_MATH_H_
