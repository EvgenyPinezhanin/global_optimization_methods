#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#if defined( _MSC_VER )
    #define _USE_MATH_DEFINES
    #include <math.h>
#else
    #include <cmath>
#endif

#include <imgo.h>
#include <task.h>
#include <output_results.h>

using namespace std;

const int functionNumber = 6; // 1 - f1, 2 - f2, 3 - f3, ...
const int displayType = 1; // 0 - application, 1 - png, 2 - notitle

double f1(double x, int j) {
    switch (j) {
        case 1: return exp(-sin(3.0 * x)) - 1.0 / 10.0 * pow(x - 1.0 / 2.0, 2) - 1.0;
        case 2: return -13.0 / 6.0 * x + sin(13.0 / 4.0 * (2.0 * x + 5.0)) - 53.0 / 12.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f2(double x, int j) {
    switch (j) {
        case 1: return 1.0 / 20.0 - exp(-2.0 / 5.0 * (x + 5.0)) * sin(4.0 / 5.0 * M_PI * (x + 5.0));
        case 2: return (11.0 * x * x - 10.0 * x + 21.0) / (2.0 * (x * x + 1.0));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f3(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1: return 3.0 / 2.0 * (cos(7.0 / 20.0 * (x + 10.0)) - sin(7.0 / 4.0 * (x + 10.0)) + 1.0 / 2.0);
        case 2:
            for (int i = 1; i <= 5; i++) {
                sum += cos(i * x);
            }
            return -sum;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f4(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1:
            for (int i = 1; i <= 5; i++) {
                sum += cos(5.0 / 4.0 * (i + 1.0) * x + i);
            }
            return 6.0 / 25.0 - sum;
        case 2: return 9.0 / 50.0 - 9.0 / 2.0 * exp(-(x - 1.0 / 10.0)) * 
                       sin(2.0 * M_PI * (x - 1.0 / 10.0));
        case 3: return 4.0 * sin(M_PI / 4.0 * x + 1.0 / 20.0) * pow(pow(sin(M_PI / 2.0 * x + 1.0 / 10.0), 3.0) +
                       pow(cos(M_PI / 2.0 * x + 1.0 / 10.0), 3.0), 2.0);
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f5(double x, int j) {
    switch (j) {
        case 1: return 17.0 / 25.0 - 2.0 / 29763.233 * (-1.0 / 6.0 * pow(x, 6) + 52.0 / 25.0 * pow(x, 5) - 39.0 / 80.0 * pow(x, 4) -
                       71.0 / 10.0 * x * x * x + 79.0 / 20.0 * x * x + x - 1.0 / 10.0);
        case 2: return -14.0 / 125.0 * (3.0 * x - 8.0) * sin(252.0 / 125.0 * (x + 3.0 / 2.0)) - 1.0 / 2.0;
        case 3: return sin(0.423531 * x + 3.13531) + sin(10.0 / 3.0 * (0.423531 * x + 3.13531)) +
                       log(0.423531 * x + 3.13531) + 0.36634 - 0.355766 * x;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f6(double x, int j) {
    switch (j) {
        case 1: return 40.0 * (cos(4.0 * x) * (x - sin(x)) * exp( -(x * x) / 2.0));
        case 2: return 2.0 / 25.0 * (x + 4.0) - sin(12.0 / 5.0 * (x + 4.0));
        case 3: return -7.0 / 40.0 * (3.0 * x + 4.0) * sin(63.0 / 20.0 * (x + 4.0));
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f7(double x, int j) {
    switch (j) {
        case 1: return pow(sin(x), 3) * exp(-sin(3.0 * x)) + 1.0 / 2.0;
        case 2: return cos(7.0 / 5.0 * (x + 3.0)) - sin(7.0 * (x + 3.0)) + 3.0 / 10.0;
        case 3: return exp(-cos(4.0 * x - 3.0)) + 1.0 / 250.0 * (4.0 * x - 3.0) * (4.0 * x - 3.0) - 1.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f8(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1:
            return exp(-sin(4.0 * x)) - 1.0 / 10.0 * pow(x - 1.0 / 2.0, 2) - 1.0;
        case 2: 
            for (int i = 1; i <= 5; i++) {
                sum += cos(5.0 * (i + 1.0) * (x + 1.0 / 2.0));
            }
            return 3.0 / 10.0 - sum;
        case 3: return (-21.0 / 20.0 * x - 13.0 / 8.0) * sin(63.0 / 10.0 * x + 63.0 / 4.0) + 1.0 / 5.0;
        case 4: return cos(7.0 / 4.0 * x + 241.0 / 40.0) - sin(35.0 / 4.0 * x + 241.0 / 8.0) - 5.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f9(double x, int j) {
    double sum = 0.0;
    switch (j) {
        case 1: return 1.0 / 40.0 * (x - 4.0) * (x - 32.0 / 5.0) * (x - 9.0) * (x - 11.0) *
                       exp(-1.0 / 10.0 * pow(x - 13.0 / 2.0, 2));
        case 2: return (pow(sin(x + 1.0), 3) + pow(cos(x + 1.0), 3)) * exp(-(x + 1.0) / 10.0);
        case 3: return exp(-cos(3.0 / 5.0 * (x - 5.0 / 2.0))) + 1.0 / 10.0 * pow(3.0 / 25.0 * x - 4.0 / 5.0, 2) - 1.0;
        case 4:
            for (int i = 1; i <= 5; i++) {
                sum += 1.0 / 5.0 * sin((i + 1.0) * x - 1.0) + 2.0;
            }
            return sum;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

double f10(double x, int j) {
    double a, b;
    switch (j) {
        case 1: return 2.0 * exp(-2.0 / M_PI * x) * sin(4.0 * x);
        case 2:
            a = 2.0 / M_PI * x - 1.0 / 2.0; 
            return -a * a * (-a * a + 5.0 * a - 6.0) / (a * a + 1.0) - 1.0 / 2.0;
        case 3: return pow(sin(x), 3) + pow(cos(2.0 * x), 3) - 3.0 / 10.0;
        case 4:
            b = 4.0 / M_PI * (x - 3.0 / 10.0) - 4.0;
            return -1.0 / 500.0 * b * b * b * b * b * b + 3.0 / 100.0 * b * b * b * b - 27.0 / 500.0 * b * b + 3.0 / 2.0;
        default: return numeric_limits<double>::quiet_NaN();
    }
}

int main() {
    ofstream ofstr("output_data/imgo_test.txt");
    if (!ofstr.is_open()) cerr << "File opening error\n";
    ofstream ofstrOpt("output_data/imgo_test_opt.txt");
    if (!ofstrOpt.is_open()) cerr << "File opening error\n";

    double eps = 0.000001, r = 2.0, d = 0.0;
    int maxTrials = 100000, maxFevals = 100000;

    vector<TaskImgo> taskArray = { TaskImgo(f1, "f1(x)", 1, -2.5, 1.5, 1.05738, vector<double>{ 4.640837, 8.666667 }, eps, maxTrials,
                                            maxFevals, r, d),
                                   TaskImgo(f2, "f2(x)", 1, -5.0, 5.0, 1.016, vector<double>{ 2.513269, 6.372595 }, eps, maxTrials,
                                            maxFevals, r, d),
                                   TaskImgo(f3, "f3(x)", 1, -10.0, 10.0, -5.9921, vector<double>{ 3.124504, 13.201241 }, eps, maxTrials,
                                            maxFevals, r, d),
                                   TaskImgo(f4, "f4(x)", 2, 0.0, 4.0, 2.45956, vector<double>{ 29.731102, 35.390605, 12.893183 }, eps,
                                            maxTrials, maxFevals, r, d),
                                   TaskImgo(f5, "f5(x)", 2, -1.5, 11.0, 9.28491, vector<double>{ 5.654617, 0.931981, 2.021595 }, eps, 
                                            maxTrials, maxFevals, r, d),
                                   TaskImgo(f6, "f6(x)", 2, -4.0, 4.0, 2.32396, vector<double>{ 2.48, 25.108154, 8.835339 }, eps, maxTrials,
                                            maxFevals, r, d),
                                   TaskImgo(f7, "f7(x)", 2, -3.0, 2.0, -0.774575, vector<double>{ 8.332010, 5.359309, 6.387862 }, eps,
                                            maxTrials, maxFevals, r, d),
                                   TaskImgo(f8, "f8(x)", 3, -2.5, 1.5, -1.12724, vector<double>{ 20.184982, 90.598898, 6.372137, 10.415012 },
                                            eps, maxTrials, maxFevals, r, d),
                                   TaskImgo(f9, "f9(x)", 3, 0.0, 14.0, 4.0, vector<double>{ 0.873861, 1.682731, 1.254588, 3.843648 }, eps,
                                            maxTrials, maxFevals, r, d),
                                   TaskImgo(f10, "f10(x)", 3, 0.0, 2.0 * M_PI, 4.2250023, vector<double>{ 3.170468, 4.329008, 7.999984, 12.442132 },
                                            eps, maxTrials, maxFevals, r, d) };

    ImgoMethod imgo;

    double x;
    vector<double> L;
    int numberTrials, numberFevals;
    vector<TrialConstrained> trials;

    for (int i = 0; i < taskArray.size(); i++) {
        if (taskArray[i].used) {
            imgo.setF(taskArray[i].f);
            imgo.setNumberConstraints(taskArray[i].numberConstraints);
            imgo.setAB(taskArray[i].A[0], taskArray[i].B[0]);
            imgo.setEps(taskArray[i].eps);
            imgo.setMaxTrials(taskArray[i].maxTrials);
            imgo.setMaxFevals(taskArray[i].maxFevals);
            imgo.setR(taskArray[i].r);
            imgo.setD(taskArray[i].d);

            imgo.solve(numberTrials, numberFevals, x);
            imgo.getL(L);

            printResultImgo(taskArray[i].name, taskArray[i].numberConstraints, taskArray[i].A[0], taskArray[i].B[0], taskArray[i].L,
                            taskArray[i].XOpt[0], taskArray[i].f(taskArray[i].XOpt[0], taskArray[i].numberConstraints + 1),
                            taskArray[i].maxTrials, taskArray[i].maxFevals, taskArray[i].eps, taskArray[i].r, taskArray[i].d,
                            numberTrials, numberFevals, L, x, taskArray[i].f(x, taskArray[i].numberConstraints + 1));

            addPointGnuplot(ofstr, taskArray[i].XOpt[0], taskArray[i].f(taskArray[i].XOpt[0], taskArray[i].numberConstraints + 1));
            addPointGnuplot(ofstr, x, taskArray[i].f(x, taskArray[i].numberConstraints + 1));

            imgo.getTrialPoints(trials);
            addPointsGnuplot(ofstr, trials);
        }
    }
    ofstr.close();

    initArrayGnuplot(ofstrOpt, "numberConstraints", taskArray.size());
    for (int i = 0; i < taskArray.size(); i++) {
        setValueInArrayGnuplot(ofstrOpt, "numberConstraints", i + 1, taskArray[i].numberConstraints, false);
    }
    ofstrOpt.close();

    vector<int> args{ displayType, functionNumber };
    drawGraphGnuplot("scripts/imgo_test.gp", args);

#if defined( _MSC_VER )
    cin.get();
#endif

    return 0;
}
