#ifndef __MEX_FUNCTION_H__
#define __MEX_FUNCTION_H__

#include <vector>
#include <cmath>

#define PI 3.14159265f
#define DEG2RAD(alpha) (alpha * PI / 180.0f)
#define RAD2DEG(alpha) (180.0f * alpha / PI)
#define SIGN(num) (num > 0 ? 1 : -1)
#define MIN(a, b) (a < b ? a : b)
#define MAX(a, b) (a > b ? a : b)
#define GA_MIN_VALUE -15
#define GA_MAX_VALUE 15

#define LAMBDA_A 1
#define LAMBDA_B 1
#define LAMBDA_C 0
#define GA_VERSION "1.4.2"

typedef std::vector<double> vector_d;

typedef struct point {
    double x;
    double y;

    point(double _x, double _y) {
        x = _x;
        y = _y;
    }
} point_t;

typedef struct homo_point {
    double x;
    double y;
    double z;

    homo_point(double _x, double _y, double _z) {
        x = _x;
        y = _y;
        z = _z;
    }
} homo_point_t;

double pdist(const point_t &veca, const point_t &vecb);

void mex_function(const vector_d &match_result_x, const vector_d &match_result_y, const vector_d &compass_reading,
                  vector_d &location_result);

#endif