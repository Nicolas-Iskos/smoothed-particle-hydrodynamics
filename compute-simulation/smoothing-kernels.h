/*
 * This file contains the intereface and definitions for the functions
 * used to compute smoothing kernels. Complete descriptions are in
 * smoothing-kernels.cu.
 * */
__device__ float cubic_spline_kernel(float origin[3], float r[3]);
