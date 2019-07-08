#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "calculate-field.h"
#include "integrate.h"

#include "test-kernels.h"
#include <stdio.h>

int main() {

    gri_to_pl_map_t grid_to_particle_list_map = gen_grid_to_particle_list_map();
    pi_to_gri_map_t particle_idx_to_grid_idx_map = gen_particle_idx_to_grid_idx_map();
    pi_to_pa_map_t particle_idx_to_addr_map = gen_particle_idx_to_addr_map();
    grid_mutex_set_t mutex_set = gen_grid_mutex_set();

    initialize_dam_break(grid_to_particle_list_map,
                         particle_idx_to_grid_idx_map,
                         particle_idx_to_addr_map);









    dim3 dimBlock(32, 1);
    dim3 dimGrid(2, 1);

    delete_particles_test<<<dimGrid, dimBlock>>>(grid_to_particle_list_map,
                                                 particle_idx_to_grid_idx_map,
                                                 particle_idx_to_addr_map,
                                                 mutex_set);
    cudaDeviceSynchronize();

    /*
    insert_particles_test<<<dimGrid, dimBlock>>>(grid_to_particle_list_map,
                                                 particle_idx_to_grid_idx_map,
                                                 particle_idx_to_addr_map,
                                                 mutex_set);
    */


    host_grid_consistency_check(grid_to_particle_list_map);
}

