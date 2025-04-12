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

const std::vector<std::vector<std::vector<double>>> optimalPointsFitting {
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0238464, 0.642321, 1.01739, 1.17626 } }, // 0
    std::vector<std::vector<double>>{ std::vector<double>{ 0.352274, 0.409603, 1.3327, 1.42695 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.750663, 0.816737, 1.02176, 1.12476 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0219031, 0.0379358, 0.106925, 0.263851 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.03745, 0.0816614, 0.153566, 1.10338 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.0544543, 0.109354, 0.289601, 1.13642 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 1.52509, 1.71797, 1.86712, 1.99101 } },
    std::vector<std::vector<double>>{ std::vector<double>{ 0.341343, 0.420049, 1.25375, 1.4306 } },
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

const std::vector<double> optimalValuesFitting {
     -4.39812, -5.49008, -12.6529, -7.94221, -9.81698, // 0
    -10.9344, -17.0111,   -7.12093, -6.2141, -2.76506,

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

// const std::vector<double> optimalValuesFitting {
    // -5.29075, -5.64098, -12.8822, -8.98809, -10.2208, 
    // -11.1829, -17.0111, -7.33195, -6.25569, -2.94888,
    // 
    // -17.6818, -6.72019, -17.0377, -7.7694, -4.90138,
    // -12.3503, -6.8715, -4.06631, -16.2776, -14.8417,
// }

template <typename OptMethod = ScanningMethod<OneDimensionalSupportiveOptProblem>>
class FittingFamilyOptProblems : public BaseFittingFamilyOptProblems<OptMethod> {
public:
    FittingFamilyOptProblems(const OptMethod &_optMethod = OptMethod(), bool _isSortX = false)
    : BaseFittingFamilyOptProblems<OptMethod>(familySizeFitting, dimensionFitting, searchAreaFitting,
      alphaFitting, deltaFitting, leftBoundWindowFitting, rightBoundWindowFitting, windowPointsFitting,
      firstPointFitting, firstValuesFitting, secondPointFitting, secondValuesFitting,
      lastPointFitting, _optMethod, _isSortX, optimalPointsFitting, optimalValuesFitting) {};
};

#endif // _FITTING_TEST_FAMILY_OPT_PROBLEMS_H_
