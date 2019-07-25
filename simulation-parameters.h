/*
 * Units: meters, kilograms, seconds
 * */


#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

/* particle parameters */
#define N_PARTICLES                 512
#define M_PARTICLE                  0.001
#define R_PARTICLE                  0.01

/* SPH parameters */
#define H                           0.02

/* experiment evolution parameters */
#define EXP_SPACE_DIM               1
#define DT                          0.1

/* physical constants */
#define G                           9.8
#define R                           8.134
#define M                           0.018
#define T                           300

#endif

