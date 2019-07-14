#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "test-kernels.h"

#include <stdint.h>
#include <stdio.h>

__global__ void device_count_particles_in_grid_slots(gri_to_pl_map_t grid_to_particle_list_map,
                                              uint32_t *particles_per_grid_slot) {

    uint32_t num_grid_slots = (uint32_t)pow(EXP_SPACE_DIM / H, 3);

    for(size_t i = 0; i < num_grid_slots; i++) {
        uint16_t grid_count = 0;
        for(Particle *p = grid_to_particle_list_map[i]; p != NULL;
                p = p->next_particle) {
            grid_count++;
        }
        particles_per_grid_slot[i] = grid_count;
    }
}

void host_count_particles(gri_to_pl_map_t grid_to_particle_list_map,
                                              uint32_t *particles_per_grid_slot_forward,
                                              uint32_t *particles_per_grid_slot_backward) {

    uint32_t num_grid_slots = (uint32_t)pow(EXP_SPACE_DIM / H, 3);

    for(size_t i = 0; i < num_grid_slots; i++) {
        uint16_t grid_count_forward = 0;
        uint16_t grid_count_backward = 0;
        Particle *prev_particle = NULL;

        for(Particle *p = grid_to_particle_list_map[i]; p != NULL;
                p = p->next_particle) {
            prev_particle = p;
            grid_count_forward++;
        }
        particles_per_grid_slot_forward[i] = grid_count_forward;

        for(Particle *p = prev_particle; p != NULL;
                p = p->prev_particle) {
            grid_count_backward++;
        }
        particles_per_grid_slot_backward[i] = grid_count_backward;
    }
}

bool host_grid_consistency_check(gri_to_pl_map_t grid_to_particle_list_map) {

    uint32_t n_grid_slots = (uint32_t)pow(EXP_SPACE_DIM / H, 3);
    uint32_t *particles_per_grid_slot_forward;
    uint32_t *particles_per_grid_slot_backward;
    cudaMallocManaged(&particles_per_grid_slot_forward, n_grid_slots * sizeof(uint32_t));
    cudaMallocManaged(&particles_per_grid_slot_backward, n_grid_slots * sizeof(uint32_t));

    host_count_particles(grid_to_particle_list_map, particles_per_grid_slot_forward,
                                                    particles_per_grid_slot_backward);

    for(size_t i = 0; i < n_grid_slots; i++) {
        if(particles_per_grid_slot_forward[i] !=
           particles_per_grid_slot_backward[i]) {
            printf("Malformed dll at index %zu\n", i);
            return false;
        }

        if(particles_per_grid_slot_forward[i] != 0) {
            printf("Non-empty cow with %d particles\n",
                    particles_per_grid_slot_forward[i]);
        }
    }

    return true;
}

void output_curr_particle_to_grid_map(pi_to_gri_map_t curr_particle_to_grid_map) {

    for(size_t i = 0; i < N_PARTICLES; i++) {
        printf("particle idx, grid idx: %zu , %d\n", i, curr_particle_to_grid_map[i]);
    }
}
/*
__global__ void insert_particle_test(gri_to_pl_map_t grid_to_particle_list_map,
                                     pi_to_gri_map_t curr_particle_to_grid_map,
                                     pi_to_pa_map_t particle_idx_to_addr_map,
                                     grid_mutex_set_t mutex_set) {

    uint32_t particle_idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t grid_idx = particle_idx % 8;

    Particle *new_particle = particle_idx_to_addr_map[particle_idx];

    device_insert_into_grid(grid_to_particle_list_map,
                            grid_idx,
                            curr_particle_to_grid_map,
                            particle_idx,
                            new_particle,
                            mutex_set);
}


__global__ void delete_particles_test(gri_to_pl_map_t grid_to_particle_list_map,
                                      pi_to_gri_map_t curr_particle_to_grid_map,
                                      pi_to_pa_map_t particle_idx_to_addr_map,
                                      grid_mutex_set_t mutex_set) {

    uint32_t particle_idx = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t grid_idx = curr_particle_to_grid_map[particle_idx];

    Particle *del_particle = particle_idx_to_addr_map[particle_idx];

    device_remove_from_grid(grid_to_particle_list_map,
                            grid_idx,
                            del_particle,
                            mutex_set);
}
*/


