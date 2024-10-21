#ifndef _GENERAL_ZERO_ORDER_TRIAL_H_
#define _GENERAL_ZERO_ORDER_TRIAL_H_

namespace opt {
    struct GeneralZeroOrderTrial {
        double x, z;

        GeneralZeroOrderTrial(double _x = 0.0, double _z = 0.0) : x(_x), z(_z) {};
    };
}

#endif // _GENERAL_ZERO_ORDER_TRIAL_H_