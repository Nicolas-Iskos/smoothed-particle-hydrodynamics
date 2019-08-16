/*
 * Units: meters, kilograms, seconds
 * */


#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

/* CUDA-specific parameters */
#define PARTICLES_PER_BLOCK         64
#define GRID_SPACES_PER_BLOCK       64

/* particle parameters */
#define N_PARTICLES                 4096
#define M_PARTICLE                  0.02
#define R_PARTICLE                  0.005

/* SPH parameters */
#define H                           0.0625
#define SL                          0.0625
#define EPSILON                     0.04
#define A_SPH                       2
#define B_SPH                       0.5
#define GAMMA                       1.330

/* experiment evolution parameters */
#define EXP_SPACE_DIM               1
#define DT                          0.004
#define DAMPING_FACTOR              0.5

/* physical constants */
#define G                           9.8
#define C                           1500
#define H_SPEC                      20

/* graphics parameters */
#define WINDOW_SIZE                 600
#define OVERLAP_FACTOR              1.5
#endif

