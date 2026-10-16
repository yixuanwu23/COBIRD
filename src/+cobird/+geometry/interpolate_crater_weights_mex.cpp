#include "mex.h"
#include <cmath>

static void require(bool condition, const char* message) {
    if (!condition) {
        mexErrMsgIdAndTxt(
            "cobird:geometry:interpolate_crater_weights_mex:input",
            "%s", message);
    }
}

static const double* dense(const mxArray* array) {
    require(mxIsDouble(array) && !mxIsComplex(array) && !mxIsSparse(array),
            "Inputs must be real, dense double arrays.");
    return mxGetDoubles(array);
}

void mexFunction(int output_count, mxArray* outputs[],
                 int input_count, const mxArray* inputs[]) {
    require(input_count == 3 && output_count == 1,
            "Expected a table, four indices, and four coefficients.");
    const double* table = dense(inputs[0]);
    const double* index = dense(inputs[1]);
    const double* coefficient = dense(inputs[2]);
    const mwSize crater_count = mxGetM(inputs[0]);
    const mwSize bin_count = mxGetN(inputs[0]);
    const mwSize direction_count = mxGetN(inputs[1]);
    require(crater_count > 0 && bin_count > 0 &&
            mxGetM(inputs[1]) == 4 && mxGetM(inputs[2]) == 4 &&
            mxGetN(inputs[2]) == direction_count,
            "Interpolation array sizes are incompatible.");
    for (mwIndex element = 0; element < 4 * direction_count; ++element) {
        require(std::isfinite(index[element]) && index[element] >= 1 &&
                index[element] <= bin_count &&
                std::floor(index[element]) == index[element],
                "A table index is invalid.");
        require(std::isfinite(coefficient[element]) && coefficient[element] >= 0,
                "An interpolation coefficient is invalid.");
    }
    outputs[0] = mxCreateDoubleMatrix(crater_count, direction_count, mxREAL);
    double* output = mxGetDoubles(outputs[0]);
    for (mwIndex direction = 0; direction < direction_count; ++direction) {
        const double* c = coefficient + 4 * direction;
        if (c[0] == 0 && c[1] == 0 && c[2] == 0 && c[3] == 0) {
            continue;
        }
        const double* p0 = table + crater_count *
            (static_cast<mwIndex>(index[4 * direction]) - 1);
        const double* p1 = table + crater_count *
            (static_cast<mwIndex>(index[4 * direction + 1]) - 1);
        const double* p2 = table + crater_count *
            (static_cast<mwIndex>(index[4 * direction + 2]) - 1);
        const double* p3 = table + crater_count *
            (static_cast<mwIndex>(index[4 * direction + 3]) - 1);
        for (mwIndex crater = 0; crater < crater_count; ++crater) {
            output[crater_count * direction + crater] =
                ((p0[crater] * c[0] + p1[crater] * c[1]) +
                 p2[crater] * c[2]) + p3[crater] * c[3];
        }
    }
}
