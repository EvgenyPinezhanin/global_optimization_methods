#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>

#include <Grishagin/GrishaginProblemFamily.hpp>
#include <Grishagin/GrishaginConstrainedProblemFamily.hpp>
#include <GKLS/GKLSProblemFamily.hpp>
#include <GKLS/GKLSConstrainedProblemFamily.hpp>
#include <mggsa.h>

using namespace std;

#define CALC
const int family_number = 3; // 0 - Grishagin, 1 - GKLS,
                             // 2 - Grishagin(constrained), 3 - GKLS(constrained),
                             // 4 - comparison Grishagin and GKLS, 5 - comparison Grishagin and GKLS (constrained)
const int number_family = 4;
int current_func;

TGrishaginProblemFamily grishaginProblems;
double f_grishagin(vector<double> x, int j) {
    switch (j) {
        case 1: return grishaginProblems[current_func]->ComputeFunction(x);
        default: return numeric_limits<double>::quiet_NaN();
    }
};

TGrishaginConstrainedProblemFamily grishaginConstrainedProblems;
double f_constr_grishagin(vector<double> x, int j) {
    int constr = grishaginConstrainedProblems[current_func]->GetConstraintsNumber();
    if (j >= 1 && j <= constr) {
        return grishaginConstrainedProblems[current_func]->ComputeConstraint(j - 1, x);
    } else if (j - 1 == constr) {
        return grishaginConstrainedProblems[current_func]->ComputeFunction(x);
    } else {
        return numeric_limits<double>::quiet_NaN();
    }
};

TGKLSProblemFamily GKLSProblems;
double f_gkls(vector<double> x, int j) {
    switch (j) {
        case 1: return GKLSProblems[current_func]->ComputeFunction({ x });
        default: return numeric_limits<double>::quiet_NaN();
    }
};

TGKLSConstrainedProblemFamily GKLSConstrainedProblems;
double f_constr_gkls(vector<double> x, int j) {
    int constr = GKLSConstrainedProblems[current_func]->GetConstraintsNumber();
    if (j >= 1 && j <= constr) {
        return GKLSConstrainedProblems[current_func]->ComputeConstraint(j - 1, x);
    } else if (j - 1 == constr) {
        return GKLSConstrainedProblems[current_func]->ComputeFunction(x);
    } else {
        return numeric_limits<double>::quiet_NaN();
    }
};

int main() {
#if defined(CALC)
    ofstream ofstr("output_data/mggsa_operational_characteristics.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_operational_characteristics_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    int count_func;

    int start_time, end_time;
    double work_time;

    int K0 = 0, Kmax = 700, Kstep = 25;
    int count_successful = 0;
    int count_trials;

    int den = 10, key = 1, m = 0;
    double eps = 0.01, r = 0.0, d = 0.0;

    vector<double> A, B;
    vector<vector<double>> r_array{ {2.5, 3.0, 3.5},
                                    {3.5, 4.3, 5.0},
                                    {2.5, 3.0, 3.5},
                                    {4.0, 4.5, 4.9} };

    vector<vector<int>> number_trials(number_family);
    number_trials[0].resize(grishaginProblems.GetFamilySize(), 0);
    number_trials[1].resize(GKLSProblems.GetFamilySize(), 0);
    number_trials[2].resize(grishaginConstrainedProblems.GetFamilySize(), 0);
    number_trials[3].resize(GKLSConstrainedProblems.GetFamilySize(), 0);

    mggsa_method mggsa(nullptr, -1, -1, A, B, -1.0, d, den, key, eps);

    ofstr_opt << "Name[1]=\"Grishagin\"" << endl;
    for (int i = 0; i < r_array[0].size(); i++) {
        ofstr_opt << "R[" << (0 * 3) + i + 1 << "]=\"" << r_array[0][i] << "\"" << endl; 
    }

    mggsa.setF(f_grishagin);
    count_func = grishaginProblems.GetFamilySize();
    mggsa.setN(grishaginProblems[0]->GetDimension());
    mggsa.setM(0);

    cout << "GrishaginProblemFamily" << endl;
    ofstr << "# GrishaginProblemFamily" << endl;
    for (int i = 0; i < r_array[0].size(); i++) {
        cout << "r = " << r_array[0][i] << endl;
        mggsa.setR(r_array[0][i]);
        start_time = clock();
        for (int i = 0; i < count_func; i++) {
            current_func = i;
            count_trials = Kmax;
            grishaginProblems[current_func]->GetBounds(A, B);
            mggsa.setAB(A, B);
            if (mggsa.solve_test(grishaginProblems[i]->GetOptimumPoint(), count_trials, Stop::ACCURNUMBER)) {
                number_trials[0][i] = count_trials;
            } else {
                number_trials[0][i] = count_trials + 1;
            }
        }
        for (int i = K0; i <= Kmax; i += Kstep) {
            count_successful = (int)count_if(number_trials[0].begin(), number_trials[0].end(), [i](double elem){ return elem <= i; });
            cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
            ofstr << i << " " << (double)count_successful / count_func << endl;
        }
        ofstr << endl << endl;
        end_time = clock();
        work_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
        cout << "time: " << work_time << endl;
    }

    Kmax = 1500;

    ofstr_opt << "Name[2]=\"GKLS\"" << endl;
    for (int i = 0; i < r_array[1].size(); i++) {
        ofstr_opt << "R[" << (1 * 3) + i + 1 << "]=\"" << r_array[1][i] << "\"" << endl; 
    }

    mggsa.setF(f_gkls);
    count_func = GKLSProblems.GetFamilySize();
    mggsa.setN(GKLSProblems[0]->GetDimension());
    mggsa.setM(0);

    cout << "GKLSProblemFamily" << endl;
    ofstr << "# GKLSProblemFamily" << endl;
    for (int i = 0; i < r_array[1].size(); i++) {
        cout << "r = " << r_array[1][i] << endl;
        mggsa.setR(r_array[1][i]);
        start_time = clock();
        for (int i = 0; i < count_func; i++) {
            current_func = i;
            count_trials = Kmax;
            GKLSProblems[current_func]->GetBounds(A, B);
            mggsa.setAB(A, B);
            if (mggsa.solve_test(GKLSProblems[i]->GetOptimumPoint(), count_trials, Stop::ACCURNUMBER)) {
                number_trials[1][i] = count_trials;
            } else {
                number_trials[1][i] = count_trials + 1;
            }
        }
        for (int i = K0; i <= Kmax; i += Kstep) {
            count_successful = (int)count_if(number_trials[1].begin(), number_trials[1].end(), [i](double elem){ return elem <= i; });
            cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
            ofstr << i << " " << (double)count_successful / count_func << endl;
        }
        ofstr << endl << endl;
        end_time = clock();
        work_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
        cout << "time: " << work_time << endl;
    }

    Kmax = 2500;

    ofstr_opt << "Name[3]=\"GrishaginConstrained\"" << endl;
    for (int i = 0; i < r_array[2].size(); i++) {
        ofstr_opt << "R[" << (2 * 3) + i + 1 << "]=\"" << r_array[2][i] << "\"" << endl; 
    }

    mggsa.setF(f_constr_grishagin);
    count_func = grishaginConstrainedProblems.GetFamilySize();
    mggsa.setN(grishaginConstrainedProblems[0]->GetDimension());
    mggsa.setM(grishaginConstrainedProblems[0]->GetConstraintsNumber());

    cout << "GrishaginProblemConstrainedFamily" << endl;
    ofstr << "# GrishaginProblemConstrainedFamily" << endl;
    for (int i = 0; i < r_array[2].size(); i++) {
        cout << "r = " << r_array[2][i] << endl;
        mggsa.setR(r_array[2][i]);
        start_time = clock();
        for (int i = 0; i < count_func; i++) {
            current_func = i;
            count_trials = Kmax;
            grishaginConstrainedProblems[current_func]->GetBounds(A, B);
            mggsa.setAB(A, B);
            if (mggsa.solve_test(grishaginConstrainedProblems[i]->GetOptimumPoint(), count_trials, Stop::ACCURNUMBER)) {
                number_trials[2][i] = count_trials;
            } else {
                number_trials[2][i] = count_trials + 1;
            }
        }
        for (int i = K0; i <= Kmax; i += Kstep) {
            count_successful = (int)count_if(number_trials[2].begin(), number_trials[2].end(), [i](double elem){ return elem <= i; });
            cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
            ofstr << i << " " << (double)count_successful / count_func << endl;
        }
        ofstr << endl << endl;
        end_time = clock();
        work_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
        cout << "time: " << work_time << endl;
    }

    Kmax = 4000;

    ofstr_opt << "Name[4]=\"GKLSConstrained\"" << endl;
    for (int i = 0; i < r_array[3].size(); i++) {
        ofstr_opt << "R[" << (3 * 3) + i + 1 << "]=\"" << r_array[3][i] << "\"" << endl; 
    }

    mggsa.setF(f_constr_gkls);
    count_func = GKLSProblems.GetFamilySize();
    mggsa.setN(GKLSProblems[0]->GetDimension());
    mggsa.setM(GKLSConstrainedProblems[0]->GetConstraintsNumber());

    cout << "GKLSProblemConstrainedFamily" << endl;
    ofstr << "# GKLSProblemConstrainedFamily" << endl;
    for (int i = 0; i < r_array[3].size(); i++) {
        cout << "r = " << r_array[3][i] << endl;
        mggsa.setR(r_array[3][i]);
        start_time = clock();
        for (int i = 0; i < count_func; i++) {
            current_func = i;
            count_trials = Kmax;
            GKLSProblems[current_func]->GetBounds(A, B);
            mggsa.setAB(A, B);
            if (mggsa.solve_test(GKLSProblems[i]->GetOptimumPoint(), count_trials, Stop::ACCURNUMBER)) {
                number_trials[3][i] = count_trials;
            } else {
                number_trials[3][i] = count_trials + 1;
            }
        }
        for (int i = K0; i <= Kmax; i += Kstep) {
            count_successful = (int)count_if(number_trials[3].begin(), number_trials[3].end(), [i](double elem){ return elem <= i; });
            cout << "K = " << i << " success rate = " << (double)count_successful / count_func << endl;
            ofstr << i << " " << (double)count_successful / count_func << endl;
        }
        ofstr << endl << endl;
        end_time = clock();
        work_time = ((double)end_time - start_time) / CLOCKS_PER_SEC;
        cout << "time: " << work_time << endl;
    }

    ofstr.close();
    ofstr_opt.close();
#endif

    // Plotting operational characteristics(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/mggsa_operational_characteristics.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_operational_characteristics.gp %d", family_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif
	return 0;
}
