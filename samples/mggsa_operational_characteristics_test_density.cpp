#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS
    #define PROC_BIND
#else
    #define PROC_BIND proc_bind(master)
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <mggsa.h>
#include <task.h>
#include <output_results.h>
#include <omp.h>

using namespace std;

// #define CALC

const int familyNumber = 4; // 1 - Grishagin, 2 - GKLS,
                            // 3 - constrained Grishagin, 4 - constrained GKLS
const int displayType = 2; // 0 - application, 1 - png, 2 - png(notitle)

int main() {
    ofstream ofstrOpt("output_data/mggsa_operational_characteristics_test_density_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    const int chunk = 2;

    const int numberFamily = 4;

    TGrishaginProblemFamily grishaginProblems;
    TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
    TGKLSProblemFamily gklsProblems;
    TGKLSConstrainedProblemFamily gklsConstrainedProblems;

    vector<ProblemFamily> problems{ ProblemFamily("GrishaginProblemFamily", &grishaginProblems, TypeConstraints::NoConstraints, "Grishagin"),
                                    ProblemFamily("GKLSProblemFamily", &gklsProblems, TypeConstraints::NoConstraints, "GKLS"),
                                    ProblemFamily("GrishaginProblemConstrainedFamily", &grishaginConstrainedProblems,
                                                  TypeConstraints::Constraints, "GrishaginConstrained"),
                                    ProblemFamily("GKLSProblemConstrainedFamily", &gklsConstrainedProblems,
                                                  TypeConstraints::Constraints, "GKLSConstrained") };

    vector<vector<int>> K{ { 0, 700,  25 },
                           { 0, 1500, 25 },
                           { 0, 2000, 50 },
                           { 0, 3000, 50 } };

    vector<double> r{ 3.0, 4.3, 3.0, 4.5 };
    vector<int> den{ 4, 6, 8, 10, 12, 4, 6, 8, 10, 12 };
    vector<double> d{ 0.0, 0.0, 0.01, 0.01 };

    int sizeDen = den.size();
    vector<vector<vector<double>>> successRate(numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        successRate[i].resize(sizeDen);
        for (int j = 0; j < sizeDen; j++) {
            successRate[i][j].resize((K[i][1] - K[i][0]) / K[i][2] + 1);
        }
    }

    double totalStartTime = omp_get_wtime();
#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics_test_density.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";

    int incr = 0, key = 3;
    double eps = 0.01;
    int maxFevals = 1000000;

    MggsaMethod mggsa(nullptr, -1, -1, vector<double>{}, vector<double>{}, -1.0, -1.0, -1, key, eps, -1, maxFevals, incr);

#pragma omp parallel for schedule(dynamic, chunk) PROC_BIND num_threads(omp_get_num_procs()) collapse(2) \
        shared(numberFamily, r, den, d, successRate) \
        firstprivate(problems, K, mggsa)
    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < sizeDen; j++) {
            vector<double> A, B;
            FunctorFamily functor;
            FunctorFamilyConstrained functorConstrained;

            if (problems[i].type == TypeConstraints::Constraints) {
                functorConstrained.constrainedOptProblemFamily = static_cast<IConstrainedOptProblemFamily*>(problems[i].optProblemFamily);
                (*functorConstrained.constrainedOptProblemFamily)[0]->GetBounds(A, B);
                mggsa.setN((*functorConstrained.constrainedOptProblemFamily)[0]->GetDimension());
                mggsa.setNumberConstraints((*functorConstrained.constrainedOptProblemFamily)[0]->GetConstraintsNumber());
            } else {
                functor.optProblemFamily = static_cast<IOptProblemFamily*>(problems[i].optProblemFamily);
                (*functor.optProblemFamily)[0]->GetBounds(A, B);
                mggsa.setN((*functor.optProblemFamily)[0]->GetDimension());
                mggsa.setNumberConstraints(0);
            }

            mggsa.setAB(A, B);
            mggsa.setR(r[i]);
            mggsa.setD(d[i]);
            mggsa.setDen(den[j]);
            mggsa.setMaxTrials(K[i][1]);

            int numberFunctions = problems[i].optProblemFamily->GetFamilySize();
            vector<int> numberTrialsArray(numberFunctions);
            vector<double> XOpt;
            int numberTrials, numberFevals, numberSuccessful, denTmp;
            TypeSolve type;

            double startTime = omp_get_wtime();
            for (int k = 0; k < numberFunctions; k++) {
                if (problems[i].type == TypeConstraints::Constraints) {
                    functorConstrained.currentFunction = k;
                    XOpt = (*functorConstrained.constrainedOptProblemFamily)[k]->GetOptimumPoint();
                    mggsa.setF(functorConstrained);
                } else {
                    functor.currentFunction = k;
                    XOpt = (*functor.optProblemFamily)[k]->GetOptimumPoint();
                    mggsa.setF(functor);
                }

                if (j < 5) {
                    denTmp = den[j];
                    type = TypeSolve::SOLVE;
                    mggsa.setDen(denTmp);
                    while (true) {
                        if (mggsa.solveTest(XOpt, numberTrials, numberFevals, type)) {
                            numberTrialsArray[k] = numberTrials;
                            break;
                        } else {
                            if (mggsa.getCoincideX()) {
                                denTmp++;
                                mggsa.setDen(denTmp);
                                type = TypeSolve::RESOLVE;
                            } else {
                                numberTrialsArray[k] = K[i][1] + 1;
                                break;
                            }
                        }
                    }
                } else {
                    if (mggsa.solveTest(XOpt, numberTrials, numberFevals)) {
                        numberTrialsArray[k] = numberTrials;
                    } else {
                        numberTrialsArray[k] = K[i][1] + 1;
                    }
                }
            }
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                numberSuccessful = (int)count_if(numberTrialsArray.begin(), numberTrialsArray.end(), [k] (double elem) { return elem <= k; });
                successRate[i][j][k / K[i][2]] = (double)numberSuccessful / numberFunctions;
            }
            double endTime = omp_get_wtime();
            double workTime = endTime - startTime;

            string strOutput = problems[i].name + " den = " + ((j < 5) ? "dynamic(start den = " + to_string(den[j]) + ")" : to_string(den[j])) +
                               " time: " + to_string(workTime) + " thread number: " + to_string(omp_get_thread_num()) + "\n";
            cout << strOutput;
        }
    }

    for (int i = 0; i < numberFamily; i++) {
        for (int j = 0; j < sizeDen; j++) {
            for (int k = K[i][0]; k <= K[i][1]; k += K[i][2]) {
                ofstr << k << " " << successRate[i][j][k / K[i][2]] << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();
#endif
    double totalEndTime = omp_get_wtime();
    double totalWorkTime = totalEndTime - totalStartTime;
    cout << "Total time: " << totalWorkTime << endl;

    setVariableGnuplot(ofstrOpt, "sizeDen", to_string(sizeDen));
    initArrayGnuplot(ofstrOpt, "familyNames", numberFamily);
    for (int i = 0; i < numberFamily; i++) {
        setValueInArrayGnuplot(ofstrOpt, "familyNames", i + 1, problems[i].shortName);
    }
    initArrayGnuplot(ofstrOpt, "den", sizeDen);
    for (int i = 0; i < sizeDen; i++) {
        if (i < 5) {
            setValueInArrayGnuplot(ofstrOpt, "den", i + 1, "dynamic(start m = " + to_string(den[i]) + ")");
        } else {
            setValueInArrayGnuplot(ofstrOpt, "den", i + 1, den[i]);
        }
    }
    ofstrOpt.close();

    vector<int> args{ displayType, familyNumber };
    drawGraphGnuplot("scripts/mggsa_operational_characteristics_test_density.gp", args);

#if defined(_MSC_VER)
    cin.get();
#endif

	return 0;
}
