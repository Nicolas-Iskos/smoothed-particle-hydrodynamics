## Smoothed particle hydrodynamics simulation

This repository contains the source code and documentation for a 3-D smoothed
particle hydrodynamics simulation application written in CUDA/C++. This readme
will contain sections on how to run the simulation, what smoothed particle
hydrodynamics is, and how the code works.



## How to run this simulation

This repository provides three bash scripts that can be used to run a simulation:

#### play-back.sh
This script can be used to view a simulation if your computer does not have nvcc
or an NVIDIA GPU. Running play-back.sh will graphically play back the results of
whatever simulation-results.csv file is in the play-back-simulation directory.
This repository includes a simulation-results.csv file made from running a
4096-particle simulation for 1.5 seconds. After cloning the repository, this
simulation can be viewed immediately. Of course, this script will play back the
results of any simulation-results.csv file placed in the play-back-simulation
directory.


#### compute.sh
This script can be used if your computer has nvcc and an NVIDIA GPU. This script
is to be run with a single command line argument: the number of seconds for which
you would like the simulation to run. It will export the results of the simulation
to a file called simulation-results.csv in the smoothed-particle-hydrodynamics
directory.


#### compute-and-play-back.sh
This script can be used if your computer has nvcc and and an NVIDIA GPU. It
computes the results of the simulation, exports these results to a file called
simulation-results.csv in the play-back-simulation directory and then reads this
file as it graphically plays back the simulation.



## About smoothed particle hydrodynamics

Smoothed particle hydrodynamics (SPH) is a computational fluid dynamics
technique used widely in astronomy to model the behavior of fluids. Within
SPH, a fluid is generally modeled as a large number of spherical particles.
The defining feature of SPH is the use of functions known as smoothing kernels
to calculate the values of certain field quantities at the position of each
particle by summing scaled values of certain fields at surrounding particles.
