#include "particle-data-structures.h"

#include <cstdint>


/*
 * uses smooth particle hydrodynamics to calculate density at the position of
 * each particle
 * */
__device__ void calculate_density(uint32_t particle_idx,
                                  pi_to_gri_map_t particle_to_grid_map,
                                  gri_to_pl_map_t grid_to_particle_list_map);

__device__ void calculate_pressure(uint32_t particle_idx,
                                   pi_to_gri_map_t particle_to_grid_map,
                                   gri_to_pl_map_t grid_to_particle_list_map);

__device__ void calculate_internal_energy(uint32_t particle_idx,
                                          pi_to_gri_map_t particle_to_grid_map,
                                          gri_to_pl_map_t grid_to_particle_list_map);

__global__ void calculate_net_force(pi_to_gri_map_t particle_to_grid_map,
                                    gri_to_pl_map_t grid_to_particle_list_map);




