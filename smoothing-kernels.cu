#include "smoothing-kernels.h"

#include <math.h>

#define H               1
#define INV_PI_HPow3    0.3183

__device__ float cubic_spline_kernel(float r)
{
    float x = r / H;

    if(0 <= x && x < 1)
    {
        return INV_PI_HPow3 * (1 - 1.5 * pow(x, 2) + 0.75 * pow(x, 3));
    }
    else if(1 <= x && x < 2)
    {
        return INV_PI_HPow3 * 0.25 * pow(2 - x, 3);
    }
    else
    {
        return 0;
    }
}
