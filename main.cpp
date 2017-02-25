#include "ibs_matrix.h"
#include "mex_function.h"
#include "mex-it/mex-it.h"
#include <iostream>
#include <python2.7/Python.h>


using namespace mex_binding;
using namespace std;

int main() {
    int nrhs = 3;
    int nlhs = 1;

    int M = 7;

    double A[] = {990.000000, 925.000000, 855.000000, 923.000000, 1343.000000, 1056.000000};
    double B[] = {995.000000, 995.000000, 995.000000, 868.000000, 894.000000, 868.000000};
    double C[] = {253.779806, 228.593566, 213.372202, 170.256359, 358.393997, 26.565051};

    point_t real_loc(1022, 885);

    mxArray *in1 = mxCreateDoubleMatrix(M, 1, mxREAL);
    mxArray *in2 = mxCreateDoubleMatrix(M, 1, mxREAL);
    mxArray *in3 = mxCreateDoubleMatrix(M, 1, mxREAL);

    memcpy(mxGetPr(in1), A, M * sizeof(double));
    memcpy(mxGetPr(in2), B, M * sizeof(double));
    memcpy(mxGetPr(in3), C, M * sizeof(double));

    mxArray *out1 = mxCreateDoubleMatrix(3, 1, mxREAL);

    const mxArray *prhs[nrhs];
    mxArray *plhs[nlhs];

    prhs[0] = in1;
    prhs[1] = in2;
    prhs[2] = in3;
    plhs[0] = out1;

    call_mex_function(mex_function, nlhs, plhs, nrhs, prhs);

//    double result = mxGetScalar(plhs[0]);
//    result = 35;

//    std::vector<double> m(M);
//    memcpy(m.data(), mxGetPr(plhs[0]), M * sizeof(double));
    double result[] = {0, 0, 0, 0, 0, 0};

    memcpy(result, mxGetPr(plhs[0]), 6 * sizeof(double));

    cout << "error :" << pdist(real_loc, point_t(result[0], result[1])) << endl;
//    iBS::Matrix matrix(3);
//    for(int i = 0; i < 3; i++) {
//        for(int j = 0; j < 3; j++) {
//            matrix.data[i][j] = i + j;
//        }
//    }
//    matrix.print();
//
//    iBS::Matrix b(3);
//    for(int i = 0; i < 3; i++) {
//        for(int j = 0; j < 3; j++) {
//            b.data[i][j] = 1;
//        }
//    }
//    b.print();
//
//    iBS::Matrix c = matrix + b;
//
//    matrix.print();
//    b.print();
//    c.print();

//#ifdef DEBUG
//    cout << __FILE__ << " => mex_function: result of (" << x << " + " << y << ") * " << z << " is: " << result << endl;
//#endif

}
