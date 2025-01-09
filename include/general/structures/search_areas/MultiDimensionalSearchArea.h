#ifndef _MULTI_DIMENSIONAL_SEARCH_AREA_H_
#define _MULTI_DIMENSIONAL_SEARCH_AREA_H_

#include <vector>

namespace opt {
    struct MultiDimensionalSearchArea {
        int dimension;
        std::vector<double> lowerBound, upBound;

        MultiDimensionalSearchArea(int _dimension = 2, const std::vector<double> &_lowerBound = std::vector<double>(),
                                   const std::vector<double> &_upBound = std::vector<double>())
            : dimension(_dimension), lowerBound(_lowerBound), upBound(_upBound) {};

        MultiDimensionalSearchArea(int _dimension, double _lowerBound, double _upBound)
            : dimension(_dimension), lowerBound(dimension, _lowerBound), upBound(dimension, _upBound) {};
    };
}

#endif // _MULTI_DIMENSIONAL_SEARCH_AREA_H_
