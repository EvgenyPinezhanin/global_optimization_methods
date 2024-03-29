#if defined( _MSC_VER )
    #define _CRT_SECURE_NO_WARNINGS    
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <limits>

#include <mggsa.h>
#include <task.h>

using namespace std;

const int sample_func_number = 2; // 0 - f1, 1 - f2, 2 - f3, 3 - f4

double f1(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - x[0] - x[1];
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2(vector<double> x, int j) {
    switch (j) {
        case 1: return (x[0] - 1.0) * (x[0] - 1.0) / 5.0 + (x[1] - 1.0) * (x[1] - 1.0) / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3(vector<double> x, int j) {
    switch (j) {
        case 1: return 1.0 - x[0] - x[1];
        case 2: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f4(vector<double> x, int j) {
    switch (j) {
        case 1: return (x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 2.0) * (x[1] - 2.0) - 2.0;
        case 2: return x[0] * x[0] / 5.0 + x[1] * x[1] / 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/mggsa_sample.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstr_opt("output_data/mggsa_sample_opt.txt");
    if (!ofstr_opt.is_open()) cerr << "File opening error\n";

    double eps = 0.01, r = 2.2, d = 0.0;
    int count_trials, n = 2, den = 10, key = 1, Nmax = 1000;
    vector<double> X(n);
    Stop stop = Stop::ACCURACY;

    vector<task_mggsa> task_array = { task_mggsa(f1, "f1", n, 0, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                                 vector<double>{4.0, 4.0}, vector<double>{}, eps, Nmax, r, d, den, key, stop),
                                      task_mggsa(f2, "f2", n, 0, vector<double>{-4.0, -4.0}, vector<double>{4.0, 4.0},
                                                 vector<double>{1.0, 1.0}, vector<double>{}, eps, Nmax, 1.0, d, den, key, stop),
                                      task_mggsa(f3, "f3", n, 1, vector<double>{-1.0, -1.0}, vector<double>{1.0, 1.0},
                                                 vector<double>{0.5, 0.5}, vector<double>{}, eps, Nmax, r, d, den, key, stop),
                                      task_mggsa(f4, "f4", n, 1, vector<double>{0.0, 0.0}, vector<double>{3.0, 3.0},
                                                 vector<double>{1.0, 1.0}, vector<double>{}, eps, Nmax, r, d, den, key, stop) };

    mggsa_method mggsa(nullptr);

    vector<double> mu;
    vector<vector<double>> points;
    for (int i = 0; i < task_array.size(); i++) {
        if (task_array[i].used) {
            mggsa.setF(task_array[i].f);
            mggsa.setN(task_array[i].n);
            mggsa.setM(task_array[i].m);
            mggsa.setAB(task_array[i].A, task_array[i].B);
            mggsa.setEps(task_array[i].eps);
            mggsa.setNmax(task_array[i].Nmax);
            mggsa.setR(task_array[i].r);
            mggsa.setD(task_array[i].d);
            mggsa.setDen(task_array[i].den);
            mggsa.setKey(task_array[i].key);

            mggsa.solve(count_trials, X, task_array[i].stop);
            mggsa.getMu(mu);

            cout << "Function: " << task_array[i].name << endl;
            cout << "Dimension = " << task_array[i].n << endl;
            cout << "Number of constrained = " << task_array[i].m << endl;
            cout << "[A; B] = [(" << task_array[i].A[0] << ", " << task_array[i].A[1] << "); (" << 
                                     task_array[i].B[0] << ", " << task_array[i].B[1] << ")]"<< endl;
            cout << "X* = (" << task_array[i].X_opt[0] << ", " << task_array[i].X_opt[1] << ")" << endl;
            cout << "f(X*) = " << task_array[i].f(task_array[i].X_opt, task_array[i].m + 1) << endl;
            cout << "Parameters for method:" << endl;
            cout << "eps = " << eps << " r = " << r << " d = " << d << endl;
            cout << "Parameters for constructing the Peano curve:" << endl;
            cout << "m = " << task_array[i].den << " key = " << task_array[i].key << endl;
            cout << "Trials result:" << endl;
            cout << "Number of trials = " << count_trials << endl;
            cout << "Estimation of the Lipschitz constant:" << endl;
            cout << "L(" << task_array[i].name << ") = " << mu[task_array[i].m] << endl;
            for (int j = 0; j < task_array[i].m; j++) {
                cout << "L(g" << j + 1 << ") = " << mu[j] << endl;
            }
            cout << "X = (" << X[0] << ", " << X[1] << ")" << endl;
            cout << "f(X) = " << task_array[i].f(X, task_array[i].m + 1) << endl;
            cout << "||X* - X|| = " << sqrt((task_array[i].X_opt[0] - X[0]) * (task_array[i].X_opt[0] - X[0]) + 
                                            (task_array[i].X_opt[1] - X[1]) * (task_array[i].X_opt[1] - X[1])) << endl;
            cout << "|f(X*) - f(X)| = " << abs(task_array[i].f(task_array[i].X_opt, task_array[i].m + 1) - 
                                               task_array[i].f(X, task_array[i].m + 1)) << endl;
            cout << endl;

            // Saving points for plotting
            ofstr << task_array[i].X_opt[0] << " " << task_array[i].X_opt[1] << " " << 
                     task_array[i].f(task_array[i].X_opt, task_array[i].m + 1) << endl;
            ofstr << endl << endl;
            ofstr << X[0] << " " << X[1] << " " << task_array[i].f(X, task_array[i].m + 1) << endl;
            ofstr << endl << endl;
            mggsa.getPoints(points);
            for (int j = 0; j < points.size(); j++) {
                ofstr << points[j][0] << " " << points[j][1] << " " << 
                         task_array[i].f(points[j], task_array[i].m + 1) << endl;
            }
            ofstr << endl << endl;
        }
    }
    ofstr.close();

    size_t size = task_array.size();
    ofstr_opt << "array AX[" << size << "]" << endl;
    ofstr_opt << "array AY[" << size << "]" << endl;
    ofstr_opt << "array BX[" << size << "]" << endl;
    ofstr_opt << "array BY[" << size << "]" << endl;
    for (int i = 0; i < size; i++) {
        ofstr_opt << "AX[" << i + 1 << "] = " << task_array[i].A[0] << endl;
        ofstr_opt << "BX[" << i + 1 << "] = " << task_array[i].B[0] << endl;
        ofstr_opt << "AY[" << i + 1 << "] = " << task_array[i].A[1] << endl;
        ofstr_opt << "BY[" << i + 1 << "] = " << task_array[i].B[1] << endl;
    }
    ofstr_opt.close();

    // Plotting the function(works with gnuplot)
    int error;
#if defined(__linux__)
    setenv("QT_QPA_PLATFORM", "xcb", false);
    error = system("chmod +x scripts/mggsa_sample.gp");
    if (error != 0) {
        cerr << "Error chmod" << endl;
    }
#endif

    char str[100];
    sprintf(str, "gnuplot -c scripts/mggsa_sample.gp %d", sample_func_number);
    error = system(str);
    if (error != 0) {
        cerr << "Error gnuplot" << endl;
    }

#if defined(_MSC_VER)
    cin.get();
#endif

    return 0;
}
