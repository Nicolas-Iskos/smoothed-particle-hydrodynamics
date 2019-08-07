#include "../simulation-parameters.h"

#include <cmath>

#include "smoothing-kernels.h"


__device__ float cubic_spline_kernel(float origin[3], float r[3]) {
    float diff_x_pow2;
    float diff_y_pow2;
    float diff_z_pow2;
    float x;
    float x_norm;
    constexpr float inv_pi_hpow3 = 1 / (M_PI * H * H * H);

    diff_x_pow2 = pow(r[0] - origin[0], 2);
    diff_y_pow2 = pow(r[1] - origin[1], 2);
    diff_z_pow2 = pow(r[2] - origin[2], 2);

    x = sqrt(diff_x_pow2 + diff_y_pow2 + diff_z_pow2);
    x_norm = x / H;

    if(0 <= x_norm && x_norm < 1)
    {
        return inv_pi_hpow3 * (1 - 1.5 * pow(x_norm, 2) + 0.75 * pow(x_norm, 3));
    }
    else if(1 <= x_norm && x_norm < 2)
    {
        return inv_pi_hpow3 * 0.25 * pow(2 - x_norm, 3);
    }
    else
    {
        return 0;
    }
}
