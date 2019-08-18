/*
 * This file contains the implementions of functions that are
 * used as smoothing kernels for calculating field quantities
 * with the laws of SPH.
 * */



#include "../simulation-parameters.h"

#include <cmath>

#include "smoothing-kernels.h"



/*
 * cubic_spline_kernel takes in the position of a current particle,
 * a position of a particle to be accumulated, and computes the
 * result of the cubic spline kernel as a function of the difference
 * between these two vectors.
 *
 * Inputs:
 * origin - a 3D position vector for the particle at which we want to calculate
 *          the value of the smoothing kernel
 *
 * r - a 3D position vector for the particle whose fields we are accumulating
 *
 * Outputs:
 * the results of the cubic spline kernel
 * */
__device__ float cubic_spline_kernel(float origin[3], float r[3]) {
    float diff_x_pow2;
    float diff_y_pow2;
    float diff_z_pow2;
    float x;
    float x_norm;

    /* this quantity is 1 / (pi * sl ^ 3) */
    constexpr float inv_pi_slpow3 = 1 / (M_PI * SL * SL * SL);

    /* find the distance between the particle at origin and the
     * particle at r
     * */
    diff_x_pow2 = pow(r[0] - origin[0], 2);
    diff_y_pow2 = pow(r[1] - origin[1], 2);
    diff_z_pow2 = pow(r[2] - origin[2], 2);

    x = sqrt(diff_x_pow2 + diff_y_pow2 + diff_z_pow2);
    x_norm = x / H;

    /* compute the cubic spine kernel */
    if(0 <= x_norm && x_norm < 1)
    {
        return inv_pi_slpow3 * (1 - 1.5 * pow(x_norm, 2) + 0.75 * pow(x_norm, 3));
    }
    else if(1 <= x_norm && x_norm < 2)
    {
        return inv_pi_slpow3 * 0.25 * pow(2 - x_norm, 3);
    }
    else
    {
        return 0;
    }
}
