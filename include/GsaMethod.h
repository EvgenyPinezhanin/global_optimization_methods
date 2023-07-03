#ifndef GSA_METHOD_H_
#define GSA_METHOD_H_

#include <vector>
#include <functional>

#include <IStronginOptimizationMethod.h>

using std::vector;
using std::function;

class GsaMethod : public IStronginOptimizationMethod<function<double(double)>> {
private:
    double m;

    int insertInSorted(vector<Trial> &trials, Trial trial) override;
    double searchMin() const override;

    void calcCharacteristic() override;

    Trial newTrial(double x) override;
    double newPoint() override;
    double selectNewPoint() override;
    bool stopConditions() override;

    void solveImplementation(vector<double> &X);
    bool solveTestImplementation(const vector<double> &XOpt);
    
public:
    GsaMethod(function<double(double)> _f = nullptr, double _a = 0.0, double _b = 10.0, double _r = 2.0, double _eps = 0.001, 
              int _maxTrials = 1000, int _maxFevals = 1000) : OptimizationMethodNoConstrained(nullptr, 1,
              vector<double>{ _a }, vector<double>{ _b }, _eps, _maxTrials, _maxFevals), f(_f), r(_r), m(0) {};
    
    void setA(double _a) { OptimizationMethod::setA(vector<double>{ _a }); };
    void setB(double _b) { OptimizationMethod::setB(vector<double>{ _b }); };
    void setAB(double _a, double _b) { OptimizationMethod::setAB(vector<double>{ _a }, vector<double>{ _b }); };
    void setR(double _r) { r = _r; };

    double getA() const { return A[0]; };
    double getB() const { return B[0]; };
    double getR() const { return r; };

    double getL() const { return m; };
    
    void solve(int &numberTrials, int &numberFevals, double &x);
    void solve(int &numberTrials, int &numberFevals, vector<double> &X) override;

    bool solveTest(int &numberTrials, int &numberFevals, double xOpt);
    bool solveTest(int &numberTrials, int &numberFevals, vector<double> XOpt) override;
};

#endif // GSA_METHOD_H_
