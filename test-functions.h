#include "particle-data-structures.h"
#include <cstdint>

__global__ void device_count_particles_in_grid_slots(
                                     gri_to_pl_map_t grid_to_particle_list_map,
                                     uint32_t *particles_per_grid_slot);


void host_count_particles_in_grid_slots(gri_to_pl_map_t grid_to_particle_list_map,
                                        uint32_t *particles_per_grid_slot_forward,
                                        uint32_t *particles_per_grid_slot_backward);

bool host_grid_consistency_check(gri_to_pl_map_t grid_to_particle_list_map);


void output_particle_idx_to_grid_idx_map(pi_to_gri_map_t curr_particle_to_grid_map);


void insert_particles_test(gri_to_pl_map_t grid_to_particle_list_map,
                                     pi_to_gri_map_t curr_particle_to_grid_map,
                                     pi_to_pa_map_t particle_idx_to_addr_map);

void delete_particles_test(gri_to_pl_map_t grid_to_particle_list_map,
                           pi_to_gri_map_t curr_particle_to_grid_map,
                           pi_to_pa_map_t particle_idx_to_addr_map);
