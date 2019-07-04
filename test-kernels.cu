#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include <stdint.h>

__global__ void device_count_particles_in_grid_slots(gr_to_p_map_t grid_to_particles_map,
                                              uint32_t *particles_per_grid_slot) {

    uint32_t num_grid_slots = (uint32_t)pow(EXP_SPACE_DIM / H, 3);

    for(size_t i = 0; i < num_grid_slots; i++) {
        uint16_t grid_count = 0;
        for(Particle *p = grid_to_particles_map[i]; p != NULL;
                p = p->next_particle) {
            grid_count++;
        }
        particles_per_grid_slot[i] = grid_count;
    }
}

void host_count_particles_in_grid_slots(gr_to_p_map_t grid_to_particles_map,
                                              uint32_t *particles_per_grid_slot) {

    uint32_t num_grid_slots = (uint32_t)pow(EXP_SPACE_DIM / H, 3);

    for(size_t i = 0; i < num_grid_slots; i++) {
        uint16_t grid_count = 0;
        for(Particle *p = grid_to_particles_map[i]; p != NULL;
                p = p->next_particle) {
            grid_count++;
        }
        particles_per_grid_slot[i] = grid_count;
    }
}
