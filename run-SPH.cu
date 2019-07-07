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

    initialize_dam_break(grid_to_particle_list_map,
                         particle_idx_to_grid_idx_map,
                         particle_idx_to_addr_map);







    host_grid_consistency_check(grid_to_particle_list_map);

    output_particle_idx_to_grid_idx_map(particle_idx_to_grid_idx_map);


}

