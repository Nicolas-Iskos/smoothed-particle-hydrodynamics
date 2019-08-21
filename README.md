## Smoothed particle hydrodynamics simulation

This repository contains the source code and documentation for a 3-D computational
fluid dynamics simulation program using the technique of smoothed
particle hydrodynamics. It was developed using CUDA/C++ and OpenGL. This readme
will contain directions on how to run the simulation, an overview of a few of the unique
features of this simulation and a very brief overview on smoothed particle hydrodynamics.
A detailed explanation of how this program works will be included in the repository wiki.
<br/>
<br/>
<br/>


## How to run this simulation

This repository provides three bash scripts that can be used to run a simulation:

#### play-back.sh
**play-back.sh can be used to view a simulation if your computer does not have nvcc
or an NVIDIA GPU.** Running play-back.sh will graphically play back the results of
whatever simulation-results.csv file is in the play-back-simulation directory.
This repository includes a simulation-results.csv file made from running a
4096-particle simulation for 1.5 seconds starting in the "dam-break" configuration,
which is basically dropping a cube-shaped fluid onto the floor.
After cloning the repository, this simulation can be viewed immediately.
Of course, this script will play back the results of any simulation-results.csv file
placed in the play-back-simulation directory. Once the simulation is running,
- pressing 'e' will exit the simulation playback
- pressing 'r' will restart the simulation playback
<br/>


#### compute.sh
**compute.sh can be used if your computer has nvcc and an NVIDIA GPU.** This script
is to be run with a single command line argument: the number of seconds for which
you would like the simulation to run. It will export the results of the simulation
to a file called simulation-results.csv in the smoothed-particle-hydrodynamics
directory.
<br/>
<br/>
<br/>


#### compute-and-play-back.sh
**compute-and-play-back.sh can be used if your computer has nvcc and and an NVIDIA GPU.**
This script is to be run with a single command line argument: the number of seconds
for which you would like the simulation to run. It computes the results of the simulation,
exports these results to a file called simulation-results.csv in the play-back-simulation
directory and then reads this file as it graphically plays back the simulation. Once the
simulation is running,
- pressing 'e' will exit the simulation playback
- pressing 'r' will restart the simulation playback
<br/>


## A few cool features of this project

Two different forms of parallel execution using CUDA:

- **Particle-wise parallel execution**: All calculations of field values
such as density, pressure and net force for each particle are performed
with one CUDA thread handling each particle.

- **Grid-wise parallel execution**: The cube-shaped experiment space in which
this simulation takes place is divided into many smaller cube-shaped sections
forming a three dimensional grid of grid spaces. During execution,
program data structures store a list of the particles occupying each grid space.
This way, when surrounding particle field values are accumulated to estimate
a field value for a single particle, we can quickly determine in O(1) time which
particles are nearby and should be considered in the calculation. As particles move
during the simulation, the list of particles in each grid space need to be updated.
To accomplish this without the use of critical sections, we employ a set of CUDA kernel
functions that are executed using one CUDA thread per grid space.
<br/>


## About smoothed particle hydrodynamics

Smoothed particle hydrodynamics (SPH) is a computational fluid dynamics
technique used widely in astronomy to model the behavior of fluids. Within
SPH, a fluid is generally modeled as a large number of spherical particles.
The defining feature of SPH is the use of functions known as smoothing kernels
to calculate the values of certain field quantities at the position of each
particle by summing scaled values of certain fields at surrounding particles.
