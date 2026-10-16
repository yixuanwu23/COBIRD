#include "mex.h"
#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

static void require(bool condition, const char* id, const char* message) {
    if (!condition) {
        mexErrMsgIdAndTxt(id, "%s", message);
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    require(nrhs == 11, "cobird:thermal:calculate_temperature_nonuniform_threads_mex:nrhs",
        "Eleven input arguments are required.");
    require(nlhs <= 2, "cobird:thermal:calculate_temperature_nonuniform_threads_mex:nlhs",
        "At most two output arguments are allowed.");

    const mxArray* temperature_array = prhs[0];
    const mwSize facet_count = mxGetM(temperature_array);
    const mwSize depth_count = mxGetN(temperature_array);
    require(mxIsDouble(temperature_array) && !mxIsComplex(temperature_array) &&
            !mxIsSparse(temperature_array) && mxGetNumberOfDimensions(temperature_array) == 2,
        "cobird:thermal:calculate_temperature_nonuniform_threads_mex:temperature",
        "Temperature must be a dense real double matrix.");
    require(depth_count >= 2,
        "cobird:thermal:calculate_temperature_nonuniform_threads_mex:depth",
        "A nonuniform grid requires at least two depth nodes.");
    require(mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]) && !mxIsSparse(prhs[1]) &&
            mxGetNumberOfElements(prhs[1]) == facet_count,
        "cobird:thermal:calculate_temperature_nonuniform_threads_mex:muSize",
        "Solar-incidence cosines must be dense real doubles matching the facet count.");
    require(mxIsDouble(prhs[2]) && !mxIsComplex(prhs[2]) && !mxIsSparse(prhs[2]) &&
            mxGetNumberOfElements(prhs[2]) == facet_count,
        "cobird:thermal:calculate_temperature_nonuniform_threads_mex:fluxSize",
        "Additional fluxes must be dense real doubles matching the facet count.");
    require(mxIsDouble(prhs[5]) && !mxIsComplex(prhs[5]) && !mxIsSparse(prhs[5]) &&
            mxGetNumberOfElements(prhs[5]) == depth_count,
        "cobird:thermal:calculate_temperature_nonuniform_threads_mex:lowerSize",
        "Lower coefficients must be dense real doubles matching the depth count.");
    require(mxIsDouble(prhs[6]) && !mxIsComplex(prhs[6]) && !mxIsSparse(prhs[6]) &&
            mxGetNumberOfElements(prhs[6]) == depth_count,
        "cobird:thermal:calculate_temperature_nonuniform_threads_mex:upperSize",
        "Upper coefficients must be dense real doubles matching the depth count.");

    const double* temperature = mxGetPr(temperature_array);
    const double* mu = mxGetPr(prhs[1]);
    const double* extra_flux = mxGetPr(prhs[2]);
    const double solar_flux = mxGetScalar(prhs[3]);
    const double epsilon_sigma = mxGetScalar(prhs[4]);
    const double* lower = mxGetPr(prhs[5]);
    const double* upper = mxGetPr(prhs[6]);
    const double surface_beta = mxGetScalar(prhs[7]);
    const double inverse_dt = mxGetScalar(prhs[8]);
    mwSize tile_size = static_cast<mwSize>(mxGetScalar(prhs[9]));
    require(tile_size >= 1,
        "cobird:thermal:calculate_temperature_nonuniform_threads_mex:tileSize",
        "Tile size must be positive.");
    tile_size = std::min(tile_size, facet_count);
    mwSize thread_count = static_cast<mwSize>(mxGetScalar(prhs[10]));
    require(thread_count >= 1,
        "cobird:thermal:calculate_temperature_nonuniform_threads_mex:threadCount",
        "Thread count must be positive.");
    thread_count = std::min(thread_count, facet_count);
    if (facet_count < 32768) {
        thread_count = 1;
    }

    plhs[0] = mxCreateDoubleMatrix(facet_count, depth_count, mxREAL);
    double* result = mxGetPr(plhs[0]);
    double* surface_flux = nullptr;
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(facet_count, 1, mxREAL);
        surface_flux = mxGetPr(plhs[1]);
    }

    auto integrate_range = [&](mwIndex range_begin, mwIndex range_end) {
        std::vector<double> diagonal(static_cast<size_t>(tile_size) * depth_count);

        for (mwIndex begin = range_begin; begin < range_end; begin += tile_size) {
            const mwSize count = std::min<mwSize>(tile_size, range_end - begin);

            for (mwIndex local = 0; local < count; ++local) {
                const mwIndex facet = begin + local;
                const double surface = temperature[facet];
                const double square = surface * surface;
                const double cube = square * surface;
                const double fourth = square * square;
                diagonal[local] = inverse_dt + upper[0]
                    + 4.0 * surface_beta * epsilon_sigma * cube;
                result[facet] = surface * inverse_dt + surface_beta
                    * (solar_flux * mu[facet] + extra_flux[facet]
                       + 3.0 * epsilon_sigma * fourth);
            }

            for (mwSize depth = 1; depth < depth_count; ++depth) {
                const mwIndex global_offset = depth * facet_count + begin;
                const mwIndex tile_offset = depth * tile_size;
                const double center = inverse_dt + lower[depth] + upper[depth];
                for (mwIndex local = 0; local < count; ++local) {
                    diagonal[tile_offset + local] = center;
                    result[global_offset + local] =
                        temperature[global_offset + local] * inverse_dt;
                }
            }

            for (mwSize depth = 1; depth < depth_count; ++depth) {
                const mwIndex current_global = depth * facet_count + begin;
                const mwIndex previous_global = (depth - 1) * facet_count + begin;
                const mwIndex current_tile = depth * tile_size;
                const mwIndex previous_tile = (depth - 1) * tile_size;
                for (mwIndex local = 0; local < count; ++local) {
                    const double multiplier =
                        -lower[depth] / diagonal[previous_tile + local];
                    diagonal[current_tile + local] += upper[depth - 1] * multiplier;
                    result[current_global + local] -=
                        multiplier * result[previous_global + local];
                }
            }

            const mwIndex last_global = (depth_count - 1) * facet_count + begin;
            const mwIndex last_tile = (depth_count - 1) * tile_size;
            for (mwIndex local = 0; local < count; ++local) {
                result[last_global + local] /= diagonal[last_tile + local];
            }

            for (mwSize depth = depth_count - 1; depth-- > 0;) {
                const mwIndex current_global = depth * facet_count + begin;
                const mwIndex next_global = (depth + 1) * facet_count + begin;
                const mwIndex current_tile = depth * tile_size;
                for (mwIndex local = 0; local < count; ++local) {
                    result[current_global + local] =
                        (result[current_global + local]
                         + upper[depth] * result[next_global + local])
                        / diagonal[current_tile + local];
                }
            }

            if (surface_flux != nullptr) {
                for (mwIndex local = 0; local < count; ++local) {
                    const double surface = result[begin + local];
                    const double square = surface * surface;
                    surface_flux[begin + local] =
                        epsilon_sigma * square * square;
                }
            }
        }
    };

    if (thread_count == 1) {
        integrate_range(0, facet_count);
        return;
    }

    const mwSize range_size =
        (facet_count + thread_count - 1) / thread_count;
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(thread_count - 1));
    for (mwSize worker = 0; worker + 1 < thread_count; ++worker) {
        const mwIndex begin = worker * range_size;
        const mwIndex end = std::min<mwIndex>(facet_count, begin + range_size);
        workers.emplace_back(integrate_range, begin, end);
    }
    const mwIndex main_begin = (thread_count - 1) * range_size;
    integrate_range(main_begin, facet_count);
    for (std::thread& worker : workers) {
        worker.join();
    }
}
