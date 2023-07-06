#ifndef I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_
#define I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_

#include <vector>

#include <base_classes/opt_methods/IGeneralOptimizationMethod.h>

using std::vector;

template <typename SolutionType, typename TrialType, typename OptimizationTaskType, typename ResultMethodType, typename PointType>
class IGeneralNumericalOptimizationMethod
    : public IGeneralOptimizationMethod<OptimizationTaskType, ResultMethodType, PointType> {
protected:
    vector<TrialType> trialPoints;

    double accuracy;
    int maxTrials, numberTrials;

    virtual TrialType newTrial(const PointType &x) = 0;
    virtual PointType selectNewPoint() = 0;

    virtual SolutionType estimateSolution() const = 0;
    virtual bool stopConditions() = 0;

public:
    IGeneralNumericalOptimizationMethod(const OptimizationTaskType &_task, double _accuracy, int _maxTrials, int _maxFevals)
                                        : IGeneralOptimizationMethod<OptimizationTaskType, ResultMethodType, PointType>(
                                        _task, _maxFevals), trialPoints(0), accuracy(_accuracy), maxTrials(_maxTrials),
                                        numberTrials(0) {};

    void setAccuracy(double _accuracy) { accuracy = _accuracy; };
    double getAccuracy() const { return accuracy; };

    void setMaxTrials(int _maxTrials) { maxTrials = _maxTrials; };
    int getMaxTrials() const { return maxTrials; };

    void getTrialPoints(vector<TrialType> &_trialPoints) const { _trialPoints = trialPoints; };
    int getNumberTrialPoints() const { return trialPoints.size(); };
};

#endif // I_GENERAL_NUMERICAL_OPTIMIZATION_METHOD_H_
