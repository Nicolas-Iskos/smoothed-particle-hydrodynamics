#include "particle-data-structures.h"
#include <stdint.h>

__global__ void device_count_particles_in_grid_slots(gr_to_p_map_t grid_to_particles_map,
                                              uint32_t *particles_per_grid_slot);


void host_count_particles_in_grid_slots(gr_to_p_map_t grid_to_particles_map,
                                        uint32_t *particles_per_grid_slot);

