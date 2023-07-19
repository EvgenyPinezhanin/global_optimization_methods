#include <opt_methods/PiyavskyMethod.h>

#include <algorithm>
#include <limits>
#include <iterator>

#include <my_math.h>

using std::numeric_limits;
using std::advance;
using std::distance;
using std::min_element;
using std::max;
using std::abs;

const double epsilon = 1e-14;

void PiyavskyMethod::calcCharacteristic() {
    double R = -numeric_limits<double>::infinity(), Rtmp;
    
    int trialPointsSize = (int)trialPoints.size();
    for (int i = 1; i < trialPointsSize; i++) {
        Rtmp = 0.5 * (constantEstimation * (trialPoints[i].x - trialPoints[i - 1].x) -
                      (trialPoints[i].z + trialPoints[i - 1].z));
        if (Rtmp > R) {
            R = Rtmp;
            t = i;
        }
    }
}

double PiyavskyMethod::selectNewPoint() {
    static double M = -1.0;

    // with optimization(const)
    if (trialPoints[t].x == problem.getSearchArea().getUpBound()) {
        M = abs((trialPoints[1].z - trialPoints[0].z) / (trialPoints[1].x - trialPoints[0].x));
    } else {
        M = max({ M, abs((trialPoints[t].z - trialPoints[(size_t)t - 1].z) / 
                         (trialPoints[t].x - trialPoints[(size_t)t - 1].x)), 
                     abs((trialPoints[(size_t)t + 1].z - trialPoints[t].z) / 
                         (trialPoints[(size_t)t + 1].x - trialPoints[t].x)) });
    }

    // without optimization
    // double Mtmp;
    // int trialPointsSize = (int)trialPoints.size();
    // for (int i = 1; i < trialPointsSize; i++) {
    //     Mtmp = abs(trialPoints[i].z - trialPoints[(size_t)i - 1].z) / (trialPoints[i].x - trialPoints[(size_t)i - 1].x);
    //     if (Mtmp > M) M = Mtmp;
    // }

    constantEstimation = (abs(M) <= epsilon) ? 1.0 : reliability * M;

    calcCharacteristic();

    return 0.5 * ((trialPoints[t].x + trialPoints[(size_t)t - 1].x) -
                  (trialPoints[t].z - trialPoints[(size_t)t - 1].z) / constantEstimation);
}