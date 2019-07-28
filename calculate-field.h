#include "particle-data-structures.h"

#include <cstdint>


/*
 * uses smoothed particle hydrodynamics to calculate density at the position of
 * each particle
 * */
__global__ void calculate_density(gri_to_pl_map_t grid_to_particle_list_map,
                                  pi_to_gri_map_t particle_to_grid_map,
                                  pi_to_pa_map_t particle_idx_to_addr_map);

__global__ void calculate_pressure(pi_to_pa_map_t particle_idx_to_addr_map);

__global__ void calculate_net_force(gri_to_pl_map_t grid_to_particle_list_map,
                                    pi_to_gri_map_t particle_to_grid_map,
                                    pi_to_pa_map_t particle_idx_to_addr_map);

__device__ void add_f_contr_from_pressure(Particle *curr_particle,
                                          Particle *acc_particle,
                                          float *norm_r_curr_acc,
                                          float *total_force);

__device__ void add_f_contr_from_viscosity(Particle *curr_particle,
                                           Particle *acc_particle,
                                           float *norm_r_curr_acc,
                                           float mag_r_curr_acc,
                                           float *total_force);

__device__ void add_f_contr_from_gravity(float *total_force);

__device__ void get_norm_3vector(float *vec, float *norm_vec);

__device__ float get_mag_3vector(float *vec);
