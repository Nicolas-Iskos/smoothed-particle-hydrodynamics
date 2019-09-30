/*
 * Units: meters, kilograms, seconds
 * */



#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H



/****************************/
/* CUDA-specific parameters */
/****************************/

#define PARTICLES_PER_BLOCK         64
#define GRID_SPACES_PER_BLOCK       64



/***********************/
/* particle parameters */
/***********************/

/* number of particles (must be divisible by two
 * if configuring in with two-block collosion
 * initial conditions)i
 * */
#define N_PARTICLES                 4096

/* particle mass */
#define M_PARTICLE                  0.02

/* radius of particle */
#define R_PARTICLE                  0.005



/******************/
/* SPH parameters */
/******************/

/* radius for which particles are included in SPH calculations */
#define H                           0.0625

/* smoothing radius */
#define SL                          0.0625

/* a generic small number used to prevent division-by-zero errors */
#define EPSILON                     0.01

/* a tuning paramter for the calculation of viscous force */
#define A_SPH                       8

/* a tuning paramter for the calculation of viscous force */
#define B_SPH                       0.3



/***********************************/
/* experiment evolution parameters */
/***********************************/

/* the side length of the cube-shaped experiment space */
#define EXP_SPACE_DIM               1

/* the time step between iterations of the simulation */
#define DT                          0.004

/* the factor by which the particle velocity is damped when
 * they collide with a wall
 * */
#define DAMPING_FACTOR              0.5



/**********************/
/* physical constants */
/**********************/

/* gravitiational constant */
#define G                           12.4

/* speed of sound in the fluid */
#define C                           1500

/* specific heat of the fluid */
#define H_SPEC                      15

/* heat capacity ratio */
#define GAMMA                       1.330



/***********************/
/* graphics parameters */
/***********************/

/* side length of the display window */
#define WINDOW_SIZE                 600

/* a tuning constant that allows the particles to be drawn in an
 * overlapping manner, which can make the mass of particles appear
 * more like a liquid
 *
 * A value of 1 corresponds to no overlap
 * A value of 2 corresponds to 50% overlap
 * */
#define OVERLAP_FACTOR              1.5

/* a tuning constatnt that allows us to slow down the frame rate of the
 * simulation to see it more clearly
 * */
#define SLOWDOWN_FACTOR             1.4



#endif
