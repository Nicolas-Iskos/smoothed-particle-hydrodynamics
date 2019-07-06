#include "simulation-parameters.h"
#include "smoothing-kernels.h"
#include <math.h>

__device__ float cubic_spline_kernel(float r) {
    constexpr float inv_pi_hpow3 = 1 / (M_PI * H * H * H);
    float x = r / H;

    if(0 <= x && x < 1)
    {
        return inv_pi_hpow3 * (1 - 1.5 * pow(x, 2) + 0.75 * pow(x, 3));
    }
    else if(1 <= x && x < 2)
    {
        return inv_pi_hpow3 * 0.25 * pow(2 - x, 3);
    }
    else
    {
        return 0;
    }
}
