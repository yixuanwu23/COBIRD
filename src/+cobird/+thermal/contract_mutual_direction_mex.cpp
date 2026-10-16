#include "mex.h"
#include <cmath>

static void require(bool condition, const char* message) {
    if (!condition) {
        mexErrMsgIdAndTxt(
            "cobird:thermal:contract_mutual_direction_mex:input",
            "%s", message);
    }
}

static const mxArray* field(const mxArray* structure, const char* key) {
    require(mxIsStruct(structure) && mxGetNumberOfElements(structure) == 1,
            "An endpoint must be a scalar structure.");
    const mxArray* value = mxGetField(structure, 0, key);
    require(value != nullptr, "An endpoint field is missing.");
    return value;
}

static const double* dense(const mxArray* array, mwSize count) {
    require(mxIsDouble(array) && !mxIsComplex(array) && !mxIsSparse(array) &&
            mxGetNumberOfElements(array) == count,
            "A real array has an incompatible size.");
    return mxGetDoubles(array);
}

static mwSize dimension(const mxArray* structure, const char* key) {
    const double value = *dense(field(structure, key), 1);
    require(std::isfinite(value) && value >= 1 && value == std::floor(value) &&
            value < 1e9, "An endpoint dimension is invalid.");
    return static_cast<mwSize>(value);
}

static mwIndex indexOf(double value, mwSize count) {
    require(std::isfinite(value) && value >= 1 && value <= count &&
            value == std::floor(value), "An endpoint index is invalid.");
    return static_cast<mwIndex>(value) - 1;
}

void mexFunction(int output_count, mxArray* outputs[],
                 int input_count, const mxArray* inputs[]) {
    require(input_count == 7 && output_count == 2,
            "Expected source and receiver endpoints, sources, geometry, and emissivity.");
    const mxArray* source = inputs[0];
    const mxArray* receiver = inputs[1];
    const mwSize source_parent_count = dimension(source, "parent_count");
    const mwSize source_crater_count = dimension(source, "crater_count");
    const mwSize receiver_parent_count = dimension(receiver, "parent_count");
    const mwSize receiver_crater_count = dimension(receiver, "crater_count");
    const mxArray* source_indices = field(source, "indices");
    const mxArray* receiver_indices = field(receiver, "indices");
    const mwSize source_endpoint_count = mxGetNumberOfElements(source_indices);
    const mwSize receiver_endpoint_count = mxGetNumberOfElements(receiver_indices);
    const double* source_index = dense(source_indices, source_endpoint_count);
    const double* receiver_index = dense(receiver_indices, receiver_endpoint_count);
    const double* source_weight = dense(field(source, "weights"),
        source_crater_count * source_endpoint_count);
    const double* receiver_weight = dense(field(receiver, "weights"),
        receiver_crater_count * receiver_endpoint_count);
    const double* rough_fraction = dense(field(source, "fraction"),
        source_endpoint_count);
    const double* receiver_area = dense(field(receiver, "area_ratio"),
        receiver_crater_count);
    const double* smooth_flux = dense(inputs[2], source_parent_count);
    const double* rough_flux = dense(inputs[3],
        source_crater_count * source_parent_count);
    const double* geometry_weight = dense(inputs[4], source_endpoint_count);
    const double* receiver_cosine = dense(inputs[5], receiver_endpoint_count);
    const double emissivity = *dense(inputs[6], 1);
    require(std::isfinite(emissivity) && emissivity >= 0 && emissivity <= 1,
            "Emissivity is invalid.");
    for (mwIndex crater = 0; crater < receiver_crater_count; ++crater) {
        require(std::isfinite(receiver_area[crater]) && receiver_area[crater] > 0,
                "A receiver crater-area ratio is invalid.");
    }

    double total = 0;
    for (mwIndex endpoint = 0; endpoint < source_endpoint_count; ++endpoint) {
        const mwIndex parent = indexOf(source_index[endpoint], source_parent_count);
        const double fraction = rough_fraction[endpoint];
        require(std::isfinite(fraction) && fraction >= 0 && fraction <= 1,
                "A rough fraction is invalid.");
        double directional = 0;
        for (mwIndex crater = 0; crater < source_crater_count; ++crater) {
            directional += source_weight[source_crater_count * endpoint + crater] *
                rough_flux[source_crater_count * parent + crater];
        }
        total += geometry_weight[endpoint] *
            ((1 - fraction) * smooth_flux[parent] + fraction * directional);
    }

    outputs[0] = mxCreateDoubleMatrix(receiver_parent_count, 1, mxREAL);
    outputs[1] = mxCreateDoubleMatrix(
        receiver_crater_count, receiver_parent_count, mxREAL);
    double* smooth_heating = mxGetDoubles(outputs[0]);
    double* rough_heating = mxGetDoubles(outputs[1]);
    for (mwIndex endpoint = 0; endpoint < receiver_endpoint_count; ++endpoint) {
        const mwIndex parent = indexOf(
            receiver_index[endpoint], receiver_parent_count);
        const double incident =
            emissivity * receiver_cosine[endpoint] * total;
        smooth_heating[parent] = incident;
        for (mwIndex crater = 0; crater < receiver_crater_count; ++crater) {
            rough_heating[receiver_crater_count * parent + crater] =
                receiver_weight[receiver_crater_count * endpoint + crater] *
                incident / receiver_area[crater];
        }
    }
}
