#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "initialize.h"
#include "calculate-field.h"
#include "integrate.h"

#include "test-kernels.h"
#include <stdio.h>

int main() {

    p_to_gr_map_t particle_to_grid_map = gen_particle_to_grid_map();
    gr_to_p_map_t grid_to_particles_map = gen_grid_to_particles_map();

    initialize_dam_break(particle_to_grid_map, grid_to_particles_map);

    uint32_t *grid_counter;
    uint32_t n_grid_spaces = (uint32_t)pow(EXP_SPACE_DIM / H, 3);
    cudaMallocManaged(&grid_counter, n_grid_spaces * sizeof(uint32_t));

    //device_count_particles_in_grid_slots<<<1, 1>>>(grid_to_particles_map, grid_counter);
    //cudaDeviceSynchronize();

    host_count_particles_in_grid_slots(grid_to_particles_map, grid_counter);

    for(int i = 0; i < n_grid_spaces; i++) {
        if(grid_counter[i] != 0) {
            printf("num_particles %d\n", grid_counter[i]);
        }
    }


    for(int i = 0; i < N_PARTICLES; i++) {
        printf("%d\n", particle_to_grid_map[i]);
    }

}

