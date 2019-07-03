#include "particle-data-structures.h"

#include <stdint.h>

/*
 * returns a linked list of the grid occupants that are neighbors to
 * the particle target
 * */
__device__ Particle *get_neighbors(particle_id target,
                                          p_to_gr_map_t particle_to_grid_map,
                                          gr_to_p_map_t grid_to_particles_map);

/*
 * uses smooth particle hydrodynamics to calculate density at the position of
 * each particle
 * */
__device__ void calculate_density(particle_id target,
                                  p_to_gr_map_t particle_to_grid_map,
                                  gr_to_p_map_t grid_to_particles_map);

__device__ void calculate_pressure(particle_id target,
                                   p_to_gr_map_t particle_to_grid_map,
                                   gr_to_p_map_t grid_to_particles_map);

__device__ void calculate_internal_energy(particle_id target,
                                          p_to_gr_map_t particle_to_grid_map,
                                          gr_to_p_map_t grid_to_particles_map);

__global__ void calculate_net_force(p_to_gr_map_t particle_to_grid_map,
                                    gr_to_p_map_t grid_to_particles_map);




