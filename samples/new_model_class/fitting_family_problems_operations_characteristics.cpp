#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(master)
#endif

#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <test_opt_problems/FittingFamilyOptProblems.h>
#include <opt_methods/MggsaMethod.h>
#include <opt_methods/GsaMethod.h>
#include <nlopt.h>
#include <Solver.h>
#include <gnuplot/OutputFile.h>
#include <gnuplot/Script.h>
#include <MyMath.h>
#include <omp.h>

// #define CALC_MGGSA
// #define CALC_DIRECT_ISRES
// #define UPDATE_ADD_INFO
#define DRAW

using OptMethod = GsaMethod<OneDimensionalSupportiveOptProblem>;
using MggsaParameters = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::Parameters;
using GsaParameters = GsaMethod<OneDimensionalSupportiveOptProblem>::Parameters;
using ErrorMetrics = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::ErrorMetrics;
using TypeSolve = MggsaMethod<FittingFamilyOptProblems<OptMethod>>::TypeSolve;
using Task = opt::Task<FittingFamilyOptProblems<OptMethod>>;
using SearchArea = opt::MultiDimensionalSearchArea;

const std::vector<std::string> methodNames{ "mggsa", "direct", "isres" };
const size_t numberMethods = methodNames.size();
const int displayType = 2; // 0 - application, 1 - png, 2 - png(notitle)

void updateAddInfo(const std::string &fileNamePrefix, const std::vector<double> &errors, const size_t familySize) {
    std::vector<std::pair<size_t, double>> operationalCharacteristicsData(familySize);
    std::ifstream operationalCharacteristicsDataFile(fileNamePrefix + "_data");

    for (size_t i = 0; i < familySize; ++i) {
        operationalCharacteristicsDataFile >> operationalCharacteristicsData[i].first
                                           >> operationalCharacteristicsData[i].second;
    }

    operationalCharacteristicsDataFile.close();

    std::ofstream addInfoFile(fileNamePrefix + "_add_info");

    size_t numberSuccessful;

    size_t numberErrors = errors.size();
    for (size_t i = 0; i < numberErrors; ++i) {
        numberSuccessful = std::count_if(operationalCharacteristicsData.begin(), operationalCharacteristicsData.end(),
            [errors, i] (std::pair<size_t, double> elem) {
                return std::abs(elem.second) <= errors[i] && elem.first != 0;
            });

        std::cout << numberSuccessful << "\n";
        
        addInfoFile << "Error: " << errors[i] << ", P = " << (double)numberSuccessful / familySize << "\n";
    }

    addInfoFile.close();
}

void saveOperationalCharacteristics(
    const std::string &fileNamePrefix,
    const std::vector<std::pair<size_t, double>> &operationalCharacteristics,
    const std::vector<std::pair<size_t, double>> &operationalCharacteristicsData,
    const std::vector<double> &errors, const size_t familySize)
{
    OutputFile operationalCharacteristicsFile, operationalCharacteristicsDataFile;

    operationalCharacteristicsFile.open(fileNamePrefix);
    if (!operationalCharacteristicsFile.isOpen()) std::cerr << fileNamePrefix << " opening error\n";

    operationalCharacteristicsDataFile.open(fileNamePrefix + "_data");
    if (!operationalCharacteristicsDataFile.isOpen()) std::cerr << fileNamePrefix << "_data opening error\n";

    size_t operationalCharacteristicsSize = operationalCharacteristics.size();
    for (size_t j = 0; j < operationalCharacteristicsSize; ++j) {
        operationalCharacteristicsFile.addPoint(operationalCharacteristics[j].first,
                                                operationalCharacteristics[j].second, false);
    }
    operationalCharacteristicsFile.close();

    size_t operationalCharacteristicsDataSize = operationalCharacteristicsData.size();
    for (size_t j = 0; j < operationalCharacteristicsDataSize; ++j) {
        operationalCharacteristicsDataFile.addPoint(operationalCharacteristicsData[j].first,
                                                    operationalCharacteristicsData[j].second, false);
    }
    operationalCharacteristicsDataFile.close();

    updateAddInfo(fileNamePrefix, errors, familySize);
}

template<size_t index>
double fittingFamilyFunctions(unsigned n, const double* x, double* grad, void* f_data) {
    FittingFamilyOptProblems<OptMethod> *fitting_data =
        reinterpret_cast<FittingFamilyOptProblems<OptMethod>*>(f_data);

    std::vector<double> x_vector(n);
    for (size_t i = 0; i < n; ++i) {
        x_vector[i] = x[i];
    }

    return fitting_data->computeConstraintFunction(x_vector, index);
}

int main() {
    std::string rootDir = "output_data/new_model_class/fitting_family_problems_operational_characteristics";

    OutputFile varsFile(rootDir + "/vars.txt");
    if (!varsFile.isOpen()) std::cerr << "vars.txt opening error\n";

    double accuracyGsa = 0.01, reliabilityGsa = 2.5;
    size_t maxTrialsGsa = 10000, maxFevalsGsa = 10000;
    GsaParameters gsaParameters(accuracyGsa, 0.0, maxTrialsGsa, maxFevalsGsa, reliabilityGsa);

    GsaMethod<OneDimensionalSupportiveOptProblem> gsa;
    gsa.setParameters(gsaParameters);

    FittingFamilyOptProblems fittingFamilyOptProblems(gsa);

    size_t dimension = fittingFamilyOptProblems.getDimension();
    size_t numberConstraints = fittingFamilyOptProblems.getNumberConstraints();
    size_t familySize = fittingFamilyOptProblems.getFamilySize();
    size_t familyAvailableSize = fittingFamilyOptProblems.getAvailableFamilySize();
    size_t familyNotAvailableSize = familySize - familyAvailableSize;

    std::string familyName;
    fittingFamilyOptProblems.getFamilyName(familyName);

    size_t numberMggsaVariants = 7;
    // size_t numberMggsaVariants = 1;
    std::vector<std::vector<double>> reliability {
        std::vector<double>(numberConstraints + 1, 5.0),
        std::vector<double>(numberConstraints + 1, 4.5),
        std::vector<double>(numberConstraints + 1, 4.0),
        std::vector<double>(numberConstraints + 1, 3.5),
        std::vector<double>(numberConstraints + 1, 3.0),
        std::vector<double>(numberConstraints + 1, 2.5),
        std::vector<double>(numberConstraints + 1, 2.0)
    };
    std::vector<size_t> key{ 3, 3, 3, 3, 3, 3, 3 };
    // std::vector<size_t> key{ 1, 1, 1, 1, 1, 1, 1 };
    // std::vector<size_t> key{ 1 };
    std::vector<double> d(numberMggsaVariants, 0.01);

    std::vector<std::vector<size_t>> K{ { 0, 50000, 25 } };
    
    double error = 1.0;
    size_t maxFevals = 1000000;
    size_t density = 11, increment = 0;
    ErrorMetrics errorMetric = ErrorMetrics::F_ERROR;
    TypeSolve typeSolve = TypeSolve::SOLVE;
    MggsaParameters parameters(0.0, error, 0, maxFevals, errorMetric, std::vector<double>{},
        0.0, density, 0, increment, typeSolve);
        
    MggsaMethod<FittingFamilyOptProblems<OptMethod>> mggsa;
        
    std::vector<double> errors{ error, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0 };

    Solver<FittingFamilyOptProblems<OptMethod>> solver;
    std::vector<std::pair<size_t, double>> operationalCharacteristics;
    std::vector<std::pair<size_t, double>> operationalCharacteristicsData;
    std::stringstream strReport, filesNamePrefix;
    OutputFile operationalCharacteristicsFile, operationalCharacteristicsDataFile;
    double workTime;

    int chunk = 1;

    double totalStartTime = omp_get_wtime();
#if defined( CALC_MGGSA )
#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs() / 2) \
        shared(reliability, key, K) \
        firstprivate(mggsa, numberMggsaVariants, parameters, fittingFamilyOptProblems, familyName, errors, familyAvailableSize)
    for (size_t i = 0; i < numberMggsaVariants; ++i) {
        Task task("FittingFamily", fittingFamilyOptProblems, parameters);
        
        parameters.reliability = reliability[i];
        parameters.key = key[i];
        parameters.maxTrials = K[0][1];
        task.parameters = parameters;
        mggsa.setParameters(parameters);

        solver.calcOperationalCharacteristicsAvailable(mggsa, task, K[0][0], K[0][1], K[0][2],
                                                       operationalCharacteristics, operationalCharacteristicsData, workTime, true);

        strReport << task.name << ", method: " << methodNames[0] << ", error = " << parameters.error
                  << ", max trials = " << parameters.maxTrials << ", max fevals = " << parameters.maxFevals
                  << ", density = " << parameters.density << ", key = " << parameters.key
                  << ", reliability = " << parameters.reliability[0] << ", time: " << workTime << "\n";
        std::cout << strReport.str();
        strReport.str("");

        filesNamePrefix << std::setprecision(2) << rootDir << "/" << methodNames[0] << "/" << task.name << "_"
                        << parameters.density << "_" << parameters.key << "_" << parameters.reliability[0];

        saveOperationalCharacteristics(filesNamePrefix.str(), operationalCharacteristics,
                                       operationalCharacteristicsData, errors, familyAvailableSize);
        filesNamePrefix.str("");
    }
#endif

size_t numberAlgorithms = 2;
std::vector<nlopt_algorithm> algorithms{ NLOPT_GN_ORIG_DIRECT, NLOPT_GN_ISRES };
std::vector<size_t> numberAlgorithmVariants{ 1, 5 };
std::vector<unsigned int> populations{ 15, 20, 25, 30, 35 };
nlopt_result result;

double startTime, endTime;
double resultValue;
int numevals;
std::vector<double> resultPoint(dimension);
std::vector<size_t> numberTrials(familySize);
std::vector<std::vector<double>> optimalPoints;
std::vector<double> coefficients;
std::ostringstream output;
SearchArea searchArea;

for (size_t i = 0; i < familySize; ++i) {
    fittingFamilyOptProblems.setProblemNumber(i);
    if (!fittingFamilyOptProblems.isAvailable()) {
        numberTrials[i] = K[0][1] + 1;
    }
}

size_t numberPoints = (K[0][1] - K[0][0]) / K[0][2] + 1;
operationalCharacteristics.resize(numberPoints);
operationalCharacteristicsData.resize(familySize);

#if defined( CALC_DIRECT_ISRES )
    for (size_t i = 0; i < numberAlgorithms; ++i) {
        nlopt_opt nlopt_opt_algorithm = nlopt_create(algorithms[i], dimension);

        nlopt_set_maxeval(nlopt_opt_algorithm, K[0][1]);

        nlopt_set_min_objective(nlopt_opt_algorithm, fittingFamilyFunctions<8>, &fittingFamilyOptProblems);
        nlopt_add_inequality_constraint(nlopt_opt_algorithm, fittingFamilyFunctions<0>, &fittingFamilyOptProblems, d[0]);
        nlopt_add_inequality_constraint(nlopt_opt_algorithm, fittingFamilyFunctions<1>, &fittingFamilyOptProblems, d[0]);
        nlopt_add_inequality_constraint(nlopt_opt_algorithm, fittingFamilyFunctions<2>, &fittingFamilyOptProblems, d[0]);
        nlopt_add_inequality_constraint(nlopt_opt_algorithm, fittingFamilyFunctions<3>, &fittingFamilyOptProblems, d[0]);
        nlopt_add_inequality_constraint(nlopt_opt_algorithm, fittingFamilyFunctions<4>, &fittingFamilyOptProblems, d[0]);
        nlopt_add_inequality_constraint(nlopt_opt_algorithm, fittingFamilyFunctions<5>, &fittingFamilyOptProblems, d[0]);
        nlopt_add_inequality_constraint(nlopt_opt_algorithm, fittingFamilyFunctions<6>, &fittingFamilyOptProblems, d[0]);
        nlopt_add_inequality_constraint(nlopt_opt_algorithm, fittingFamilyFunctions<7>, &fittingFamilyOptProblems, d[0]);

        for (size_t j = 0; j < numberAlgorithmVariants[i]; ++j) {
            if (algorithms[i] == NLOPT_GN_ISRES) {
                nlopt_set_population(nlopt_opt_algorithm, populations[j]);
            }

            startTime = omp_get_wtime();
            for (size_t k = 0; k < familySize; ++k) {
                fittingFamilyOptProblems.setProblemNumber(k);

                if (!fittingFamilyOptProblems.isAvailable())
                    continue;

                searchArea = fittingFamilyOptProblems.getSearchArea();
                nlopt_set_lower_bounds(nlopt_opt_algorithm, searchArea.lowerBound.data());
                nlopt_set_upper_bounds(nlopt_opt_algorithm, searchArea.upBound.data());

                nlopt_set_stopval(nlopt_opt_algorithm, fittingFamilyOptProblems.getOptimalValue() + error);

                for (size_t i = 0; i < dimension; ++i) {
                    resultPoint[i] = (searchArea.upBound[i] + searchArea.lowerBound[i]) / 2.0;
                }

                result = nlopt_optimize(nlopt_opt_algorithm, resultPoint.data(), &resultValue);

                numevals = nlopt_get_numevals(nlopt_opt_algorithm);
                numberTrials[k] = result == NLOPT_STOPVAL_REACHED ? numevals : K[0][1] + 1;

                operationalCharacteristicsData[k].first = numberTrials[k];
                operationalCharacteristicsData[k].second = fittingFamilyOptProblems.getOptimalValue() - resultValue;

                fittingFamilyOptProblems.getOptimalPoints(optimalPoints);
                auto iter = std::min_element(optimalPoints.begin(), optimalPoints.end(),
                    [&resultPoint] (const std::vector<double> &firstPoint,
                                    const std::vector<double> &secondPoint)
                    {
                        return euclideanDistance(firstPoint, resultPoint) < euclideanDistance(secondPoint, resultPoint);
                    });

                fittingFamilyOptProblems.computeObjectiveFunction(optimalPoints[0]);
                fittingFamilyOptProblems.getCoefficients(coefficients);
                size_t numberCoeff = fittingFamilyOptProblems.getNumberCoefficients();

                output << "Problem number = " << k << ", number trials = " << numevals
                       << ", f(X*) - f(X) = " << fittingFamilyOptProblems.getOptimalValue() - resultValue
                       << ", f(X) = " << resultValue << ", ||X* - X|| = " << euclideanDistance(*iter, resultPoint)
                       << ", X = (" << resultPoint[0];
                for (size_t i = 1; i < dimension; ++i) {
                    output << ", " << resultPoint[i];
                }
                output << "), coeff = (" << coefficients[0];
                for (size_t i = 1; i < numberCoeff; ++i) {
                    output << ", " << coefficients[i];
                }
                output << "), Stopping condition: ";
                if (result == NLOPT_STOPVAL_REACHED) {
                    output << "Error\n";
                } else {
                    output << "Max trials\n";
                }
                std::cout << output.str();
                output.str("");
            }
            endTime = omp_get_wtime();
            workTime = endTime - startTime;

            size_t numberSuccessful, k;
            for (size_t i = 0, k = K[0][0]; k <= K[0][1]; ++i, k += K[0][2]) {
                numberSuccessful = std::count_if(numberTrials.begin(), numberTrials.end(),
                                                 [k] (double elem) { return elem <= k; });
                operationalCharacteristics[i] = std::pair<size_t, double>(
                    k, (double)numberSuccessful / (familyAvailableSize));
            }

            strReport << "FittingFamily" << ", method: " << methodNames[i + 1] << ", error = " << error
                      << ", max trials = " << K[0][1] << ", max fevals = " << K[0][1] << ", time: " << workTime << "\n";
            std::cout << strReport.str();
            strReport.str("");

            if (algorithms[i] == NLOPT_GN_ISRES) {
                filesNamePrefix << rootDir << "/" << methodNames[i + 1] << "/" << "FittingFamily_" << populations[j];
            } else {
                filesNamePrefix << rootDir << "/" << methodNames[i + 1] << "/" << "FittingFamily";
            }
            
            saveOperationalCharacteristics(filesNamePrefix.str(), operationalCharacteristics,
                                           operationalCharacteristicsData, errors, familyAvailableSize);
            filesNamePrefix.str("");
        }

        nlopt_destroy(nlopt_opt_algorithm);
    }
#endif
    double totalEndTime = omp_get_wtime();
    std::cout << "Total time: " << totalEndTime - totalStartTime << std::endl;

#if defined( UPDATE_ADD_INFO )
    std::ostringstream updateFilesNamePrefix;
    // updateFilesNamePrefix << rootDir << "/" << methodNames[2] << "/FittingFamily_" << populations[2];
    updateFilesNamePrefix << rootDir << "/" << methodNames[1] << "/FittingFamily";
    updateAddInfo(updateFilesNamePrefix.str(), errors, familyAvailableSize);
#endif

    varsFile.setVariable("familyName", "FittingFamily");
    varsFile.setVariable("numberMethods", methodNames.size(), false);
    varsFile.initArray("methodNames", methodNames.size());
    for (size_t i = 0; i < methodNames.size(); i++) {
        varsFile.setValueInArray("methodNames", i + 1, methodNames[i]);
    }

    std::vector<size_t> numberVariantsDraw{ numberMggsaVariants - 1,
                                            numberAlgorithmVariants[0],
                                            numberAlgorithmVariants[1] - 2};
    
    varsFile.initArray("numberVariants", methodNames.size());
    for (size_t i = 0; i < methodNames.size(); i++) {
        varsFile.setValueInArray("numberVariants", i + 1, numberVariantsDraw[i], false);
    }
    
    std::vector<std::vector<double>> reliabilityDraw {
        std::vector<double>(numberConstraints + 1, 2.0),
        // std::vector<double>(numberConstraints + 1, 2.5),
        std::vector<double>(numberConstraints + 1, 3.0),
        // std::vector<double>(numberConstraints + 1, 3.5),
        std::vector<double>(numberConstraints + 1, 4.0),
        // std::vector<double>(numberConstraints + 1, 4.5),
        // std::vector<double>(numberConstraints + 1, 5.0),
        std::vector<double>(numberConstraints + 1, 2.0),
        // std::vector<double>(numberConstraints + 1, 2.5),
        std::vector<double>(numberConstraints + 1, 3.0),
        // std::vector<double>(numberConstraints + 1, 3.5),
        std::vector<double>(numberConstraints + 1, 4.0),
        // std::vector<double>(numberConstraints + 1, 4.5),
        // std::vector<double>(numberConstraints + 1, 5.0),
    };

    std::vector<size_t> keyDraw{ 1, 1, 1, 3, 3, 3 };
    // std::vector<size_t> keyDraw{ 3, 3, 3, 3, 3, 3, 3 };
    // std::vector<size_t> keyDraw{ 1, 1, 1, 1, 1, 1, 1 };

    if (numberVariantsDraw[0] > 0) {
        size_t densityDraw = density;
        varsFile.setVariable("density", densityDraw);
        varsFile.initArray("r", numberVariantsDraw[0]);
        varsFile.initArray("key", numberVariantsDraw[0]);
        for (int i = 0; i < numberVariantsDraw[0]; ++i) {
            varsFile.setValueInArray("r", i + 1, reliabilityDraw[i][0]);
            varsFile.setValueInArray("key", i + 1, keyDraw[i]);
        }
    }

    std::vector<unsigned int> populationsDraw{ 15, 25, 35 };

    if (numberVariantsDraw[2] > 0) {
        varsFile.initArray("population", numberVariantsDraw[2]);
        for (size_t i = 0; i < numberVariantsDraw[2]; ++i) {
            varsFile.setValueInArray("population", i + 1, populationsDraw[i]);
        }
    }

    varsFile.close();

#if defined( DRAW )
    Script script("scripts/new_model_class/fitting_family_problems_operational_characteristics.gp");
    script.addArgs(std::vector<int>{ displayType });
    script.start();
    if (script.isError() == 2) std::cerr << "Error gnuplot\n";
    if (script.isError() == 1) std::cerr << "Error chmod\n";
#endif

#if defined( _MSC_VER )
    cin.get();
#endif

	return 0;
}
