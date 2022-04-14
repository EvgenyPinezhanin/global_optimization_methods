#include <imgo.h>
#include <iostream>
#include <fstream>
#include <ctime>

// #define DEBUG
// #define TIME_TEST

#if defined(TIME_TEST)
    ofstream ofstr_test;
    vector<int> time_count;
    vector<double> time_mu, time_z_star, time_R, time_trial, time_add_trial, time_add_I;
    int start_time, end_time;
#endif

const double peano_a = 0.0, peano_b = 1.0, peano_random = 0.5;

double dist_vec(vector<double> val1, vector<double> val2) {
    double res = 0.0;
    size_t size = val1.size();
    for (int i = 0; i < size; i++) {
        res += (val1[i] - val2[i]) * (val1[i] - val2[i]);
    }
    res = sqrt(res);
    return res;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void addInSort(vector<trial> &vec, trial tr) {
    vector<trial>::iterator iter = vec.begin();
    vector<trial>::iterator iterEnd = vec.end();
    while(true) {
        if (iter == iterEnd || iter->x > tr.x) break;
        iter++;
    }
    vec.insert(iter, tr);
}

double searchMinXTrial(vector<trial> &trials, int m) {
    double z = 0.0, x = 0.0;
    int k = 0;
    for (int i = 0; i < trials.size(); i++) {
        if (trials[i].nu == m + 1) {
            z = trials[i].z;
            x = trials[i].x;
            k = i + 1;
            break;
        }
    }
    for (int i = k; i < trials.size(); i++) {
        if (trials[i].z < z && trials[i].nu == m + 1) {
            z = trials[i].z;
            x = trials[i].x;
        }
    }
    return x;
}

trial imgo_method::newTrial_single(double x) {
    trial tr(x);
    for (int j = 1; j <= m + 1; j++) {
        if ((f_single(x, j) > 0) || (j == m + 1)) {
            tr.z = f_single(x, j);
            tr.nu = j;
            break;
        }
    }
    return tr;
}

trial imgo_method::newTrial(double x) {
    trial tr(x);
    vector<double> X(n);
    y(x, X);
    for (int j = 1; j <= m + 1; j++) {
        if ((f(X, j) > 0) || (j == m + 1)) {
            tr.z = f(X, j);
            tr.nu = j;
            break;
        }
    }

#if defined(DEBUG)
    std::cout << x << " " << X[0] << " " << X[1] << " " << f(X, m + 1) << " " << tr.nu << std::endl;
#endif

    return tr;
}

double imgo_method::deltaX(double x1, double x2) {
    if (n == 1) {
        return abs(x2 - x1);
    } else {
        return pow(abs(x2 - x1), 1.0 / n);
    }
}

double imgo_method::newPoint(size_t t) {
    if (trial_points[t].nu != trial_points[t - 1].nu) {
        return (trial_points[t].x + trial_points[t - 1].x) / 2.0;
    } else {
        if (n == 1) {
            return (trial_points[t].x + trial_points[t - 1].x) / 2.0 
                   - (trial_points[t].z - trial_points[t - 1].z) / (2.0 * r * mu[trial_points[t].nu - 1]);
        } else {
            return (trial_points[t].x + trial_points[t - 1].x) / 2.0 
                   - sgn(trial_points[t].z - trial_points[t - 1].z) / (2.0 * r) 
                   * pow(abs(trial_points[t].z - trial_points[t - 1].z) / mu[trial_points[t].nu - 1], n);
        }
    }
}

double imgo_method::selectNewPoint(size_t &t, trial last_trial) {
    // Шаг 3
    #if defined(TIME_TEST)
        start_time = clock();
    #endif

    double mu_tmp;
    size_t size_I = I[last_trial.nu - 1].size();
    for (int nu = 0; nu < m + 1; nu++) {
        if (!calc_I[nu]) mu[nu] = 0.0;
    }
    for (int i = 0; i < size_I; i++) {
        if (I[last_trial.nu - 1][i].x != last_trial.x) {
            mu_tmp = abs(I[last_trial.nu - 1][i].z - last_trial.z) / deltaX(last_trial.x, I[last_trial.nu - 1][i].x);
            if (mu_tmp > mu[last_trial.nu - 1]) {
                mu[last_trial.nu - 1] = mu_tmp;
                if (abs(mu[last_trial.nu - 1]) > 6e-13) calc_I[last_trial.nu - 1] = true;
            }
        }
    }
    for (int nu = 0; nu < m + 1; nu++) {
        if (!calc_I[nu] && abs(mu[nu]) < 6e-13) mu[nu] = 1.0;
    }

    // double mu_tmp;
    // size_t size_I = I[last_trial.nu - 1].size();
    // mu[last_trial.nu - 1] = 0.0;
    // for (int i = 1; i < size_I; i++) { // при i = 0 - нет j
    //     for (int j = 0; j < i; j++){
    //         mu_tmp = abs(I[last_trial.nu - 1][i].z - I[last_trial.nu - 1][j].z) / deltaX(I[last_trial.nu - 1][j].x, I[last_trial.nu - 1][i].x);
    //         if (mu_tmp > mu[last_trial.nu - 1]) {
    //             mu[last_trial.nu - 1] = mu_tmp;
    //         }
    //     }
    // }
    // if (abs(mu[last_trial.nu - 1]) < 6e-13) {
    //     mu[last_trial.nu - 1] = 1.0;
    // }

    //for (int nu = 0; nu < m + 1; nu++) {
    //    mu[nu] = 0.0;
    //}
    //double mu_tmp;
    //size_t size_I;
    //for (int nu = 0; nu < m + 1; nu++) {
    //    size_I = I[nu].size();
    //    for (int i = 1; i < size_I; i++) { // при i = 0 - нет j
    //        for (int j = 0; j < i; j++) {
    //            mu_tmp = abs(I[nu][i].z - I[nu][j].z) / deltaX(I[nu][j].x, I[nu][i].x);
    //            if (mu_tmp > mu[nu]) {
    //                mu[nu] = mu_tmp;
    //            }
    //        }
    //    }
    //    if (abs(mu[nu]) < 6e-13) {
    //        mu[nu] = 1.0;
    //    };
    //}

    //for (int i = 1; i < size_I; i++) {
    //    mu_tmp = abs(I[last_I - 1][i].z - I[last_I - 1][i - 1].z) / deltaX(I[last_I - 1][i - 1].x, I[last_I - 1][i].x);
    //    if (mu_tmp > mu[last_I - 1]) {
    //        mu[last_I - 1] = mu_tmp;
    //    }
    //}
    #if defined(TIME_TEST)
        end_time = clock();
        time_mu.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

#if defined(DEBUG)
    cout << "mu: ";
    for (int nu = 0; nu < m + 1; nu++) {
        cout << mu[nu] << " ";
    }
    cout << endl;
#endif

    #if defined(TIME_TEST)
        start_time = clock();
    #endif
    // Шаг 4
    if (n >= 2) {
        for (int nu = 0; nu < m + 1; nu++) {
            if (I[nu].size() != 0) {
                if (nu + 1 != M) {
                    z_star[nu] = -mu[nu] * d;
                } else {
                    z_star[nu] = I[nu][0].z;
                    size_I = I[nu].size();
                    for (int i = 1; i < size_I; i++) {
                        if (I[nu][i].z < z_star[nu]) {
                            z_star[nu] = I[nu][i].z;
                        }
                    }
                }
            }
        }
    } else {
        for (int nu = 0; nu < m + 1; nu++) {
            if (I[nu].size() != 0) {
                z_star[nu] = I[nu][0].z;
                size_I = I[nu].size();
                for (int i = 1; i < size_I; i++) {
                    if (I[nu][i].z < z_star[nu]) {
                        z_star[nu] = I[nu][i].z;
                    }
                }
                if (z_star[nu] <= 0 && nu != m) {
                    z_star[nu] = -mu[nu] * d;
                }
            }
        }
    }
    #if defined(TIME_TEST)
        end_time = clock();
        time_z_star.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

#if defined(DEBUG)
    cout << "z_star: ";
    for (int nu = 0; nu < m + 1; nu++) {
        cout << z_star[nu] << " ";
    }
    cout << endl;
#endif

    #if defined(TIME_TEST)
        start_time = clock();
    #endif
    // Шаг 5, 6
    t = 1;
    double d_x = deltaX(trial_points[0].x, trial_points[1].x); 
    double R, Rtmp = 0.0;
    double mu_v, z_star_v;
    if (trial_points[1].nu == trial_points[0].nu) {
        mu_v = mu[trial_points[1].nu - 1];
        z_star_v = z_star[trial_points[1].nu - 1];
        R = d_x + pow(trial_points[1].z - trial_points[0].z, 2) / (r * r * mu_v * mu_v * d_x) -
            2.0 * (trial_points[1].z + trial_points[0].z - 2.0 * z_star_v) / (r * mu_v);
    } else if (trial_points[0].nu < trial_points[1].nu) {
        mu_v = mu[trial_points[1].nu - 1];
        z_star_v = z_star[trial_points[1].nu - 1];
        R = 2.0 * d_x  - 4.0 * (trial_points[1].z - z_star_v) / (r * mu_v);
    } else {
        mu_v = mu[trial_points[0].nu - 1];
        z_star_v = z_star[trial_points[0].nu - 1];
        R = 2.0 * d_x  - 4.0 * (trial_points[0].z - z_star_v) / (r * mu_v);
    }
#if defined(DEBUG)
    cout << "R[" << trial_points[0].x << ", " << trial_points[1].x << "] = " << R << endl;
#endif
    size_t size_tr_pt = trial_points.size();
    for (size_t i = 2; i < size_tr_pt; i++) {
        d_x = deltaX(trial_points[i - 1].x, trial_points[i].x);
        if (trial_points[i].nu == trial_points[i - 1].nu) {
            mu_v = mu[trial_points[i].nu - 1];
            z_star_v = z_star[trial_points[i].nu - 1];
            Rtmp = d_x + pow(trial_points[i].z - trial_points[i - 1].z, 2) / (r * r * mu_v * mu_v * d_x) -
                   2.0 * (trial_points[i].z + trial_points[i - 1].z - 2.0 * z_star_v) / (r * mu_v);
        } else if (trial_points[i - 1].nu < trial_points[i].nu) {
            mu_v = mu[trial_points[i].nu - 1];
            z_star_v = z_star[trial_points[i].nu - 1];
            Rtmp = 2.0 * d_x  - 4.0 * (trial_points[i].z - z_star_v) / (r * mu_v);
        } else  {
            mu_v = mu[trial_points[i - 1].nu - 1];
            z_star_v = z_star[trial_points[i - 1].nu - 1];
            Rtmp = 2.0 * d_x  - 4.0 * (trial_points[i - 1].z - z_star_v) / (r * mu_v);
        }
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    #if defined(DEBUG)
        cout << "R[" << trial_points[i - 1].x << ", " << trial_points[i].x << "] = " << Rtmp << endl;
    #endif
    }
    #if defined(TIME_TEST)
        end_time = clock();
        time_R.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

    // Шаг 7
    return newPoint(t);
}

void imgo_method::setM(int _m) {
    optimization_method::setM(_m);
    I.resize(m + 1);
    calc_I.resize(m + 1);
    mu.resize(m + 1);
    z_star.resize(m + 1);
}

void imgo_method::getPoints(vector<vector<double>> &points_vec) {
    points_vec.clear();
    vector<double> point(n);
    for (int i = 0; i < trial_points.size(); i++) {
        y(trial_points[i].x, point);
        points_vec.push_back(point);
    }
}

void imgo_method::y(double x, vector<double> &X) {
    mapd(x, den, X.data(), n, key);
    for (int i = 0; i < n; i++) {
        X[i] = X[i] * (B[i] - A[i]) + (A[i] + B[i]) / 2.0;
    }
}

double imgo_method::solve(int &n, Stop stop) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
    }
    trial_points.clear();

    trial tr = newTrial_single(A[0]);
    trial_points.push_back(tr);
    addInSort(I[tr.nu - 1], tr);
    tr = newTrial_single(B[0]);
    trial_points.push_back(tr);
    addInSort(I[tr.nu - 1], tr);
    n = 2;

    double x_k_1;
    size_t t = 1;
    while(true) {
        x_k_1 = selectNewPoint(t, tr);
        tr = newTrial_single(x_k_1);
        // Шаг 1
        addInSort(trial_points, tr);

        // Шаг 2
        addInSort(I[tr.nu - 1], tr);
        n++;
        if (stop == ACCURACY) {
            if (trial_points[t].x - trial_points[t - 1].x <= eps) {
                break;
            }
        } else if (stop == NUMBER){
            if (n >= Nmax) {
                break;
            }
        } else {
            if (trial_points[t].x - trial_points[t - 1].x <= eps || n >= Nmax) {
                break;
            }
        }

    }
    return searchMinXTrial(trial_points, m);
}

void imgo_method::solve(int &count, vector<double> &X, Stop stop) {
    #if defined(TIME_TEST)
        ofstr_test.open("imgo_time_test.txt");
        time_count.clear();
        time_mu.clear();
        time_z_star.clear();
        time_R.clear();
        time_trial.clear();
        time_add_trial.clear();
        time_add_I.clear();
    #endif

    for (int nu = 0; nu < m + 1; nu++) {
        I[nu].clear();
        // mu[nu] = 0.0;
        // mu[nu] = 1.0;
        calc_I[nu] = false;
    }
    trial_points.clear();
    M = 0;

    trial tr{peano_a, -1.0, 0};
    trial_points.push_back(tr);
    tr.x = peano_b;
    trial_points.push_back(tr);

    #if defined(TIME_TEST)
        start_time = clock();
    #endif
    tr = newTrial(peano_random);
    #if defined(TIME_TEST)
        end_time = clock();
    #endif

    #if defined(TIME_TEST)
        time_count.push_back(0);
        time_mu.push_back(0.0);
        time_z_star.push_back(0.0);
        time_R.push_back(0.0);
        time_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

    #if defined(TIME_TEST)
        start_time = clock();
    #endif
    addInSort(trial_points, tr);
    #if defined(TIME_TEST)
        end_time = clock();
        time_add_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

    #if defined(TIME_TEST)
        start_time = clock();
    #endif
    addInSort(I[tr.nu - 1], tr);
    #if defined(TIME_TEST)
        end_time = clock();
        time_add_I.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
    #endif

    if (tr.nu > M) {
        M = tr.nu;
    }
    count = 1;

    double x_k_1;
    size_t t = 1;
    while(true) {
        #if defined(TIME_TEST)
            time_count.push_back(count);
        #endif

        x_k_1 = selectNewPoint(t, tr);

        #if defined(TIME_TEST)
            start_time = clock();
        #endif
        tr = newTrial(x_k_1);
        #if defined(TIME_TEST)
            end_time = clock();
            time_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
        #endif
        
        #if defined(TIME_TEST)
            start_time = clock();
        #endif
        // Шаг 1
        addInSort(trial_points, tr);
        #if defined(TIME_TEST)
            end_time = clock();
            time_add_trial.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
        #endif

        #if defined(TIME_TEST)
            start_time = clock();
        #endif
        // Шаг 2
        addInSort(I[tr.nu - 1], tr);
        #if defined(TIME_TEST)
            end_time = clock();
            time_add_I.push_back((double)(end_time - start_time) / CLOCKS_PER_SEC);
        #endif

        if (tr.nu > M) {
            M = tr.nu;
        }

        count++;
        if (stop == ACCURACY) {
            if (pow(abs(trial_points[t].x - trial_points[t - 1].x), 1.0 / n) <= eps) {
                break;
            }
        } else if (stop == NUMBER){
            if (count >= Nmax) {
                break;
            }
        } else {
            if (pow(abs(trial_points[t].x - trial_points[t - 1].x), 1.0 / n) <= eps || count >= Nmax) {
                break;
            }
        }
    }
    y(searchMinXTrial(trial_points, m), X);
    #if defined(TIME_TEST)
        for (int i = 0; i < count - 1; i++) {
            ofstr_test << time_count[i] << " " << time_mu[i] << " " << time_z_star[i] << " "
                       << time_R[i] << " " << time_trial[i] << " "  << time_add_trial[i] << " " << time_add_I[i] << endl;
        }
        ofstr_test.close();

        #if defined(__linux__)
            int error;
            setenv("QT_QPA_PLATFORM", "xcb", false);
            error = system("chmod +x chart.gp");
            if (error != 0) {
                cerr << "Error chmod" << std::endl;
            }
            error = system("gnuplot -p -c time_test.gp");

            if (error != 0) {
                cerr << "Error gnuplot" << std::endl;
            }
        #endif
    #endif
}

bool imgo_method::solve_test(double x_opt, int k) {
    for (int i = 0; i < I.size(); i++) {
        I[i].clear();
    }
    trial_points.clear();

    trial tr = newTrial_single(A[0]);
    trial_points.push_back(tr);
    addInSort(I[tr.nu - 1], tr);
    tr = newTrial_single(B[0]);
    trial_points.push_back(tr);
    addInSort(I[tr.nu - 1], tr);
    int n = 2;

    double x_k_1;
    size_t t = 1;
    while (true) {
        x_k_1 = selectNewPoint(t, tr);
        tr = newTrial_single(x_k_1);
        // Шаг 1
        addInSort(trial_points, tr);

        // Шаг 2
        addInSort(I[tr.nu - 1], tr);
        n++;
        if (abs(x_k_1 - x_opt) <= eps) {
            return true;
        }
        if (n >= k) {
            return false;
        }
    }
}

void imgo_method::solve_test(vector<double> x_opt, int &count) {
    for (int nu = 0; nu < m + 1; nu++) {
        I[nu].clear();
        calc_I[nu] = false;
    }
    trial_points.clear();
    M = 0;

    trial tr{peano_a, -1.0, 0};
    trial_points.push_back(tr);
    tr.x = peano_b;
    trial_points.push_back(tr);
    tr = newTrial(peano_random);
    addInSort(trial_points, tr);
    addInSort(I[tr.nu - 1], tr);
    if (tr.nu > M) {
        M = tr.nu;
    }
    count = 1;

    double x_k_1;
    size_t t = 1;
    vector<double> X(n);
    while (true) {
        x_k_1 = selectNewPoint(t, tr);
        tr = newTrial(x_k_1);
        // Шаг 1
        addInSort(trial_points, tr);

        // Шаг 2
        addInSort(I[tr.nu - 1], tr);
        if (tr.nu > M) {
            M = tr.nu;
        }
        count++;
        y(x_k_1, X);
        //cout << "x_new = " << x_k_1 << ", dist_vec = " << dist_vec(X, x_opt) << endl;
        if (dist_vec(X, x_opt) <= eps) {
            return ;
        }
    }
}
