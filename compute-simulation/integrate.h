/*
 * This file contains the interface for the set of functions that are used to
 * update the positions and velocities of each particle in the simulation.
 * Each function is described in detail in integrate.cu.
 * */

__global__ void euler_integrate(pi_to_pa_map_t particle_idx_to_addr_map);

__global__ void enforce_boundary_conditions(pi_to_pa_map_t particle_idx_to_addr_map);
