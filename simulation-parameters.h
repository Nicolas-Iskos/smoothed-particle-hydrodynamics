/*
 * Units: meters, kilograms, seconds
 * */


#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

/* CUDA-specific parameters */
#define PARTICLES_PER_BLOCK         128
#define GRID_SPACES_PER_BLOCK        128
/* particle parameters */
#define N_PARTICLES                 512
#define M_PARTICLE                  0.005
#define R_PARTICLE                  0.01

/* SPH parameters */
#define H                           0.02
#define EPSILON                     0.01
#define A_SPH                       1
#define B_SPH                       2
#define GAMMA                       1.330

/* experiment evolution parameters */
#define EXP_SPACE_DIM               1
#define DT                          0.1

/* physical constants */
#define G                           9.8
#define C                           1500
#define H_SPEC                      4186

#endif

