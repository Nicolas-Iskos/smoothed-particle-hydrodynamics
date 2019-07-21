#include "particle-data-structures.h"

#include <cstdint>


/*
 * uses smooth particle hydrodynamics to calculate density at the position of
 * each particle
 * */
__global__ void calculate_density(gri_to_pl_map_t grid_to_particle_list_map,
                                  pi_to_gri_map_t particle_to_grid_map,
                                  pi_to_pa_map_t particle_idx_to_addr_map,
                                  float (*sm_kernel)(float*, float*));

__global__ void calculate_pressure(gri_to_pl_map_t grid_to_particle_list_map,
                                   pi_to_gri_map_t particle_to_grid_map,
                                   pi_to_pa_map_t particle_idx_to_addr_map);

__global__ void calculate_internal_energy(gri_to_pl_map_t grid_to_particle_list_map,
                                          pi_to_gri_map_t particle_to_grid_map,
                                          pi_to_pa_map_t particle_idx_to_addr_map);

__global__ void calculate_net_force(gri_to_pl_map_t grid_to_particle_list_map,
                                    pi_to_gri_map_t particle_to_grid_map,
                                    pi_to_pa_map_t particle_idx_to_addr_map);




