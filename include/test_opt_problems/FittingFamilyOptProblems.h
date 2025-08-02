#ifndef _FITTING_TEST_FAMILY_OPT_PROBLEMS_H_
#define _FITTING_TEST_FAMILY_OPT_PROBLEMS_H_

#include <opt_methods/ScanningMethod.h>
#include <opt_problems/BaseFittingFamilyOptProblems.h>
#include <general/structures/search_areas/MultiDimensionalSearchArea.h>

const size_t familySizeFitting = 100;
const size_t dimensionFitting = 4;

const std::vector<double> lowerBound{ 0.01, 0.01, 0.01, 0.01 };
const std::vector<double> upperBound{ 2.0, 2.0, 2.0, 2.0 };
const opt::MultiDimensionalSearchArea searchAreaFitting(dimensionFitting, lowerBound, upperBound);

const double alphaFitting = 0.03, deltaFitting = 0.3,
             leftBoundWindowFitting = 1.0, rightBoundWindowFitting = 10.0;

// to generate random values
// std::mt19937_64 gen(30032001);
// std::vector<double> firstValuesFitting(100), secondValuesFitting(100);
// 
// for (size_t i = 0; i < 100; ++i) {
//     firstValuesFitting[i] = ((double)(gen() - gen.min()) / (gen.max() - gen.min()) - 0.5) * 20.0;
//     secondValuesFitting[i] = ((double)(gen() - gen.min()) / (gen.max() - gen.min()) - 0.5) * 20.0;
// }

const double firstPointFitting = 13.0, secondPointFitting = 16.65,
             lastPointFitting = 18.0;
const std::vector<double> firstValuesFitting{
    -7.78107, -4.90685, -3.37658, -2.73023,  -5.94084,
     7.96594, -1.50398, -7.49603, -7.45911,   2.17224,
    -1.14612, -2.40455,  1.03556,  8.00782,   4.49972,
     5.09156,  6.94298, -5.16885,  1.49283,  -2.41119,
     4.15064, -1.94869, -1.83866, -1.58951,   1.09116,
    -9.23451, -9.27416, -7.83039,  5.11324,   6.33287,
     5.89405,  2.29501, -4.59979, -3.95681,  -2.83353,
    -8.73609, -6.50321, -3.74151,  2.46662,  -2.94344,
    -8.67077, -9.0065,   2.89405,  8.98783,   6.6901,
    -8.11935,  4.96405, -6.61755,  0.33684,   0.90863,
    -7.56295, -5.0218,  -7.57752,  4.88265,  -4.20966,
    -4.59785, -7.61615,  2.19408,  5.1424,    9.74696,
    -9.55962, -3.08478,  1.49403, -2.43129,   4.10099,
     2.30374, -7.77526, -6.27269,  9.42625,  -8.02276,
    -1.11275,  3.03413,  5.60052,  8.05824,   6.90708,
    -9.64279, -7.01972, -3.67732,  4.85766,   9.47868,
     3.20472,  1.15996, -8.74756,  0.872564, -4.25306,
    -5.93276,  5.78546, -6.71757, -9.52325,  -4.38575,
    -0.781122, 9.99732, -5.61775,  8.03561,  -7.03563,
     1.69774,  3.79953, -9.02108,  5.0596,   -9.34146
};
const std::vector<double> secondValuesFitting{
    -0.00290043, -4.34861,   9.67357,   2.9791,   -8.32345,
     9.4947,     -5.08374,  -5.55973,  -3.39451,  -0.0823394,
     9.14136,    -2.34577,   4.85999,   7.25819,  -2.04228,
    -9.12765,    -2.78671,   0.892104, -8.20987,  -4.98163,
     2.53174,     7.43376,   6.14034,  -0.800055, -7.53466,
     0.599866,    0.236695,  4.6762,    8.90145,  -3.70739,
    -8.46062,     5.17593,  -1.67315,  -0.729681, -5.54228,
    -7.86698,     6.19509,  -3.52645,  -1.56853,  -0.590776,
     4.08909,    -4.31443,  -0.239009,  7.99223,  -1.62053,
     5.49558,     5.35255,  -6.19067,  -7.54457,   0.99195,
    -9.98292,    -2.90758,  -3.91495,   1.03898,  -5.36399,
     2.12583,     0.630124, -5.813,    -5.08053,  -0.608901,
    -2.47634,    -0.351673,  3.03741,  -5.00285,  -9.86953,
     1.45434,     0.296785, -4.3987,    1.66226,   7.20481,
    -3.55114,    -3.02148,  -6.26793,   0.21012,  -0.0051492,
    -4.42785,     3.84965,  -8.9007,    4.50456,  -4.65675,
     8.16895,     6.88843,   6.25763,  -5.37536,  -1.90303,
     9.42028,    -3.51618,   1.82551,  -0.525289, -9.74838,
     6.00976,     2.57032,  -5.58742,  -5.03271,   7.73093,
     5.69401,    -1.28251,  -9.38067,  -0.453314, -6.41473
};

const std::vector<double> windowPointsFitting{ 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
                                               5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5 };

#if defined( TIKHONOV_REGULARIZATION )
const std::vector<std::vector<std::vector<double>>> optimalPointsFitting {
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0724304, 0.631146, 1.06306, 1.15294 } }, // 0
        std::vector<std::vector<double>>{ std::vector<double>{ 0.197048, 0.684832, 1.33294, 1.44371 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.653981, 0.870179, 1.02759, 1.23116 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0413367, 0.880868, 1.06694, 1.19132 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.342072, 0.515516, 0.673414, 1.48477 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.335756, 0.495111, 0.640377, 1.32833 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.43132, 1.65967, 1.86858, 1.99101 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.313407, 0.459645, 1.28703, 1.41189 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.336241, 0.467904, 1.28169, 1.37059 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.616085, 0.843944, 1.52072, 1.64607 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.416891, 1.11844, 1.23407, 1.37497 } }, // 1
        std::vector<std::vector<double>>{ std::vector<double>{ 1.00913, 1.51489, 1.7923, 1.90842 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.4855, 1.68566, 1.86737, 1.99951 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.296888, 0.447985, 0.649122, 1.35699 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0418225, 0.785643, 1.03877, 1.11601 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.290087, 0.880868, 1.06694, 1.19132 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0413367, 0.69382, 1.06694, 1.17626 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.103038, 0.631632, 1.09804, 1.20687 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.501913, 1.12719, 1.21464, 1.38274 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.25545, 1.6587, 1.81271, 1.98324 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.829612, 1.29505, 1.69732, 1.80809 } }, // 2
        std::vector<std::vector<double>>{ std::vector<double>{ 0.992125, 1.09512, 1.33221, 1.46193 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.03488, 1.14176, 1.27683, 1.99587 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.647179, 0.787587, 1.53384, 1.64218 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.434867, 1.10338, 1.25351, 1.35456 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0680579, 0.63309, 1.06354, 1.1544 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0709729, 0.637462, 1.06597, 1.15634 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0612561, 0.778842, 1.00038, 1.18792 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.971477, 1.50104, 1.61376, 1.74008 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.04328, 0.786615, 0.985323, 1.21415 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.198749, 0.876495, 1.03148, 1.15634 } }, // 3
        std::vector<std::vector<double>>{ std::vector<double>{ 1.25448, 1.66161, 1.81077, 1.9813 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.31365, 0.747991, 1.30865, 1.39221 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.365878, 0.750177, 1.2914, 1.40072 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.17723, 1.63538, 1.80639, 1.9264 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.287172, 0.449928, 0.628717, 1.3225 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0729163, 0.81868, 1.03634, 1.16023 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.856576, 1.39197, 1.68688, 1.81611 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0724304, 0.880868, 1.00524, 1.16071 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.403773, 0.770582, 1.2982, 1.40217 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.0466809, 0.750177, 0.990668, 1.18695 } }, // 4
        std::vector<std::vector<double>>{ std::vector<double>{ 0.341586, 0.438268, 1.27245, 1.36088 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.808478, 1.15634, 1.68785, 1.80056 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.272111, 0.475677, 0.636005, 1.33124 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0413367, 0.662726, 1.06694, 1.17626 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0738879, 0.811879, 0.968319, 1.16314 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.325553, 0.400372, 1.29432, 1.40315 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.29033, 0.450657, 1.28436, 1.40873 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.02759, 1.23359, 1.40072, 1.57853 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.996984, 1.6242, 1.82194, 1.95992 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.336727, 0.513573, 0.648636, 1.48088 } }, // 5
        std::vector<std::vector<double>>{ std::vector<double>{ 0.238102, 0.682159, 1.30598, 1.39537 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.233243, 0.62726, 1.30112, 1.39926 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.291058, 0.683617, 1.28849, 1.41966 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.863621, 1.47675, 1.63028, 1.78089 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0418225, 0.802648, 0.971234, 1.19181 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0641711, 0.636005, 1.07472, 1.168 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.981437, 1.15925, 1.29335, 1.98761 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0729163, 0.81868, 1.06646, 1.16071 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0481384, 0.643778, 1.03828, 1.16023 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.284742, 0.555355, 1.25302, 1.3497 } }, // 6
        std::vector<std::vector<double>>{ std::vector<double>{ 0.784186, 1.19132, 1.67619, 1.78453 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.12427, 1.68688, 1.81854, 1.9813 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.25496, 1.6587, 1.81319, 1.98275 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.591793, 0.830826, 1.04314, 1.18257 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.798761, 1.31569, 1.71505, 1.82874 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0418225, 0.642321, 1.06306, 1.1714 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.290573, 0.493168, 1.30549, 1.391 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.351302, 0.471791, 1.24816, 1.32638 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0855481, 0.846859, 0.997469, 1.14225 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 1.52242, 1.68177, 1.86251, 1.99951 } }, // 7
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0384216, 0.882811, 1.03342, 1.19132 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0724304, 0.880382, 1.00524, 1.19132 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.289358, 0.663455, 1.25812, 1.36501 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0573694, 0.639406, 1.07909, 1.16508 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.344015, 0.441669, 1.23553, 1.34825 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0666003, 0.826454, 0.927994, 1.1714 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.3825, 1.58752, 1.76631, 1.91789 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.311949, 0.460616, 1.31569, 1.42647 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0408508, 0.753578, 0.997955, 1.17966 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 1.36792, 1.60793, 1.79352, 1.92372 } }, // 8
        std::vector<std::vector<double>>{ std::vector<double>{ 1.1267, 1.35165, 1.56444, 1.6791 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0758313, 0.806049, 0.989696, 1.13448 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.44847, 1.12136, 1.24136, 1.38614 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.302961, 0.739246, 1.29699, 1.41068 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.258507, 0.880382, 1.03585, 1.16023 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0491101, 0.814308, 0.976578, 1.17529 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0413367, 0.662726, 1.0514, 1.17626 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0427942, 0.63892, 1.05043, 1.13156 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.35019, 1.55813, 1.78599, 1.90453 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.406202, 1.10581, 1.25059, 1.38128 } }, // 9
        std::vector<std::vector<double>>{ std::vector<double>{ 0.329197, 0.473005, 1.19982, 1.36404 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.343043, 0.407174, 1.34533, 1.43181 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0738879, 0.788558, 1.02128, 1.11116 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.108383, 0.858519, 0.997469, 1.185 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.39294, 1.66258, 1.83651, 1.98081 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.105953, 0.632604, 1.12767, 1.22339 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.316322, 0.474706, 0.627745, 1.32541 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0724304, 0.631632, 1.06743, 1.22241 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.272596, 0.499484, 1.26808, 1.41238 } }
    };
#else
    const std::vector<std::vector<std::vector<double>>> optimalPointsFitting {
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0238464, 0.642321, 1.01739, 1.17626 } }, // 0
        std::vector<std::vector<double>>{ std::vector<double>{ 0.352274, 0.409603, 1.3327, 1.42695 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.750663, 0.816737, 1.02176, 1.12476 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0219031, 0.0379358, 0.106925, 0.263851 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.03745, 0.0816614, 0.153566, 1.10338 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0544543, 0.109354, 0.289601, 1.13642 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.52509, 1.71797, 1.86712, 1.99101 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.266038, 0.498269, 1.29893, 1.40582 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.335756, 0.466447, 1.27877, 1.37205 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.690419, 0.785157, 1.50129, 1.62323 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.0185022, 0.0296765, 0.131703, 0.261422 } }, // 1
        std::vector<std::vector<double>>{ std::vector<double>{ 1.00111, 1.50639, 1.76291, 1.88146 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.48501, 1.64922, 1.77457, 1.98251 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.141906, 0.162797, 0.466447, 1.17577 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0209314, 0.828397, 0.94257, 1.15002 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.230571, 0.911719, 1.04484, 1.13909 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0102429, 0.756493, 1.00524, 1.12962 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0199597, 0.762809, 0.947428, 1.07763 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0243323, 0.0364783, 0.107897, 0.263851 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.2948, 1.67133, 1.80542, 1.9915 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.884755, 1.26711, 1.72429, 1.83166 } }, // 2
        std::vector<std::vector<double>>{ std::vector<double>{ 0.025304, 0.0573694, 0.067572, 0.187089 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0272473, 0.031134, 0.153566, 0.178829 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.686532, 0.77204, 1.53384, 1.63489 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0209314, 0.100123, 0.110326, 0.136075 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0447375, 0.642321, 1.05188, 1.16314 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0471667, 0.641349, 1.05528, 1.16508 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0219031, 0.834713, 0.978522, 1.07229 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.00962, 1.49886, 1.65092, 1.75392 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0223889, 0.861434, 0.966376, 1.06403 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.103524, 0.927509, 1.03585, 1.11359 } }, // 3
        std::vector<std::vector<double>>{ std::vector<double>{ 1.4838, 1.6757, 1.80542, 1.93854 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.287657, 0.719569, 1.27634, 1.40315 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.376566, 0.753578, 1.29578, 1.39343 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.17286, 1.67036, 1.79085, 1.92154 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0301624, 0.0991516, 0.498026, 1.17966 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0413367, 0.92265, 0.987753, 1.07375 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.834227, 1.38177, 1.68882, 1.80639 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0379358, 0.055426, 0.0612561, 0.191461 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.766696, 1.18938, 1.68542, 1.79279 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.016073, 0.733658, 1.04654, 1.15537 } }, // 4
        std::vector<std::vector<double>>{ std::vector<double>{ 0.32361, 0.471791, 1.26468, 1.37108 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0301624, 0.0442517, 0.0879773, 0.137047 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0223889, 0.0641711, 0.515031, 1.19326 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0724304, 0.662726, 1.06694, 1.16071 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0136438, 0.864349, 0.969777, 1.05431 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0500818, 0.0967224, 0.157452, 1.16266 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.331383, 0.418834, 1.31181, 1.41043 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0121863, 0.0505676, 0.0607703, 0.140448 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.991154, 1.67959, 1.78356, 1.90356 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.0126721, 0.0146155, 0.0957507, 1.18743 } }, // 5
        std::vector<std::vector<double>>{ std::vector<double>{ 0.243446, 0.688961, 1.31424, 1.40217 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.369279, 0.41349, 1.288, 1.37108 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.318265, 0.714225, 1.30063, 1.38566 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.805077, 1.49157, 1.61303, 1.73838 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.013158, 0.761351, 1.0616, 1.14565 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0175305, 0.640377, 1.07715, 1.15682 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0121863, 0.0505676, 0.0607703, 0.140448 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0345349, 0.931881, 1.00281, 1.08201 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0607703, 0.636005, 1.05869, 1.1578 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.366364, 0.437296, 1.24476, 1.32784 } }, // 6
        std::vector<std::vector<double>>{ std::vector<double>{ 0.418834, 0.75795, 1.31132, 1.39343 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.09901, 1.71505, 1.82485, 1.94389 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.29189, 1.70242, 1.82291, 1.93951 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.383368, 0.973663, 1.05188, 1.14468 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.776898, 1.3259, 1.68736, 1.80056 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.013158, 0.639891, 1.0752, 1.15537 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.248304, 0.554869, 1.31181, 1.40266 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.288629, 0.567987, 1.25351, 1.34436 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.016073, 0.914391, 0.982894, 1.05286 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 1.55861, 1.72623, 1.84137, 1.96478 } }, // 7
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0583411, 0.0631995, 0.10401, 0.138505 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0277332, 0.927023, 1.00087, 1.10873 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0136438, 0.642806, 1.03925, 1.15197 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0248181, 0.63309, 1.06986, 1.16703 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.356647, 0.421749, 1.19666, 1.41772 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0112146, 0.79536, 1.02176, 1.09075 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.43861, 1.61837, 1.75489, 1.88364 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.372679, 0.400372, 1.32735, 1.41529 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0112146, 0.756493, 1.03585, 1.1131 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 1.46582, 1.6242, 1.76995, 1.88364 } }, // 8
        std::vector<std::vector<double>>{ std::vector<double>{ 1.59165, 1.71943, 1.85352, 1.98032 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0151013, 0.872123, 0.967833, 1.04945 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0257898, 0.0515393, 0.0559119, 0.250248 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.290087, 0.739489, 1.29432, 1.4012 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.181259, 0.94257, 1.03585, 1.11407 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0209314, 0.880382, 0.945485, 1.08881 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0107288, 0.689933, 1.06063, 1.16703 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0398792, 0.628231, 1.05528, 1.13836 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.35116, 1.6174, 1.72526, 1.89919 } },

        std::vector<std::vector<double>>{ std::vector<double>{ 0.0350208, 0.0413367, 0.142391, 0.160367 } }, // 9
        std::vector<std::vector<double>>{ std::vector<double>{ 0.367335, 0.440697, 1.24379, 1.32735 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0141296, 0.0967224, 0.344015, 1.17966 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0107288, 0.751635, 1.07375, 1.15294 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.025304, 0.888156, 1.00524, 1.13739 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 1.59262, 1.71068, 1.8336, 1.95506 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0398792, 0.0583411, 0.0695154, 0.145792 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.0539685, 0.114213, 0.420292, 1.17772 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.025304, 0.632118, 1.09804, 1.19181 } },
        std::vector<std::vector<double>>{ std::vector<double>{ 0.310978, 0.455272, 1.29432, 1.38371 } }
    };
#endif

// mpreal on
// const std::vector<std::vector<std::vector<double>>> optimalPointsFitting {
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.0238464, 0.634061, 1.08152, 1.15002 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.371222, 0.401344, 1.32007, 1.40995 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.741918, 0.87601, 0.982408, 1.08346 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.0583411, 0.0631995, 0.10401, 0.138505 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.0335632, 0.0471667, 0.0515393, 1.28703 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.0194739, 0.0806897, 0.266766, 1.15051 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 1.52509, 1.71797, 1.86712, 1.99101 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.314379, 0.448956, 1.29043, 1.42501 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.347902, 0.456244, 1.27974, 1.36574 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.993583, 1.09658, 1.87004, 1.99393 } },
    // 
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.028219, 0.031134, 0.143363, 0.188546 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 1.00913, 1.51489, 1.7923, 1.90842 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 1.59117, 1.7102, 1.83263, 1.97644 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.353732, 0.389198, 1.25691, 1.49546 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.0136438, 0.830341, 1.00476, 1.0786 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.227899, 0.94257, 1.03634, 1.12962 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.0102429, 0.709366, 1.06694, 1.17626 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.0126721, 0.690905, 1.06111, 1.16703 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 0.0180164, 0.0772888, 0.0894348, 0.144821 } },
    // std::vector<std::vector<double>>{ std::vector<double>{ 1.29189, 1.70242, 1.82291, 1.93951 } },
// }

#if defined( TIKHONOV_REGULARIZATION )
    const std::vector<double> optimalValuesFitting {
        -4.52365, -5.26025, -11.8478, -5.04889, -8.39909, // 0
        -9.40334, -15.0927, -7.48442, -6.25505, -2.65371, 
        
        -14.0928, -6.67006, -15.1666, -7.05359, -4.69209, // 1
        -11.1976, -6.52119, -4.1002, -12.0584, -13.8383, 
        
        -6.96502, -10.3714, -8.5517, -2.98319, -11.5112, // 2
        -5.84953, -5.63254, -8.85234, -19.3343, -7.10207, 
        
        -11.3743, -14.4261, -5.26694, -4.35615, -14.9657, // 3
        -7.56302, -10.0686, -9.04664, -3.30207, -3.65537, 
        
        -8.39889, -7.40059, -2.90611, -7.67385, -5.38046, // 4
        -9.74884, -5.71842, -6.98536, -12.0477, -3.44791, 
        
        -9.95392, -5.05705, -6.61978, -3.92415, -13.3956, // 5
        -4.77029, -5.25417, -8.1197, -8.30461, -5.95772, 
        
        -6.20172, -3.46069, -9.53628, -13.9049, -11.6636, // 6
        -4.28526, -5.1093, -6.20598, -5.3596, -11.9177, 
        
        -10.7845, -5.22585, -9.99107, -4.7345, -4.46234, // 7
        -7.47454, -7.46277, -22.9203, -5.56639, -9.49687, 
        
        -21.5444, -12.6195, -10.8818, -8.27507, -4.92472, // 8
        -11.7985, -6.83026, -5.34834, -4.99816, -24.2445, 
        
        -9.33101, -5.67815, -6.52912, -9.31025, -11.7634, // 9
        -16.5753, -3.84276, -9.17758, -3.78746, -8.47815, 
    };
#else
    const std::vector<double> optimalValuesFitting {
         -4.39812, -5.49008, -12.6529, -7.94221, -9.81698, // 0
        -10.9344, -17.0111,   -7.12513, -6.2141, -2.76506,

        -17.1457, -6.55617, -15.0427, -7.66111, -4.69186, // 1
        -12.0117, -6.58987, -3.41779, -15.8898, -14.057, 

        -5.92744, -15.5127, -13.0735, -3.07574, -14.6724, // 2
        -6.01648, -5.87151, -9.30976, -20.1666, -7.60941, 

        -13.0512, -15.5456, -4.58603, -4.6671, -15.9816, // 3
        -8.1798, -10.9674, -9.0848, -6.2951, -3.68981, 

        -8.95866, -7.40272, -3.6479, -8.22028, -5.11955, // 4
        -10.4913, -6.67487, -7.34365, -14.5174, -3.67739, 

        -12.3376, -5.20478, -6.61751, -4.63993, -13.8956, // 5
        -5.04039, -5.57397, -13.5985, -9.11871, -6.24619, 

        -5.93161, -3.7835, -10.5554, -14.9016, -12.4179, // 6
        -4.35153, -5.43081, -6.14475, -5.63019, -12.8601, 

        -12.0021, -9.44716, -10.7009, -4.61897, -4.70596, // 7
        -6.90997, -8.0072, -25.0219, -5.75251, -9.9702, 

        -23.6709, -24.8198, -11.6794, -11.723, -4.82685, // 8
        -13.545, -7.19689, -5.73958, -5.3355, -25.4979, 

        -11.564, -6.18911, -6.61944, -9.8628, -12.8757, // 9
        -19.4708, -6.63385, -10.2526, -4.39814, -8.93492
    };
#endif

// mpreal on
// const std::vector<double> optimalValuesFitting {
    // -5.29075, -5.64098, -12.8822, -8.98809, -10.2208, 
    // -11.1829, -17.0111, -7.33195, -6.25569, -2.94888,
    // 
    // -17.6818, -6.72019, -17.0377, -7.7694, -4.90138,
    // -12.3503, -6.8715, -4.06631, -16.2776, -14.8417,
// }

#if defined( TIKHONOV_REGULARIZATION )
    std::vector<size_t> badConditioned{};
#else
    std::vector<size_t> badConditioned{ 3, 10, 18, 21, 22, 24, 38, 42, 48, 57, 71, 83, 90, 96 };
#endif

template <typename OptMethod = ScanningMethod<OneDimensionalSupportiveOptProblem>>
class FittingFamilyOptProblems : public BaseFittingFamilyOptProblems<OptMethod> {
public:
    FittingFamilyOptProblems(const OptMethod &_optMethod = OptMethod(), bool _isSortX = false)
    : BaseFittingFamilyOptProblems<OptMethod>(familySizeFitting, dimensionFitting, searchAreaFitting,
      alphaFitting, deltaFitting, leftBoundWindowFitting, rightBoundWindowFitting, windowPointsFitting,
      firstPointFitting, firstValuesFitting, secondPointFitting, secondValuesFitting,
      lastPointFitting, _optMethod, _isSortX, badConditioned, optimalPointsFitting, optimalValuesFitting) {};
};

#endif // _FITTING_TEST_FAMILY_OPT_PROBLEMS_H_
