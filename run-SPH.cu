#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "calculate-field.h"
#include "integrate.h"

#include "test-kernels.h"
#include <stdio.h>
#include <assert.h>

int main() {

    gri_to_pl_map_t grid_to_particle_list_map = gen_grid_to_particle_list_map();

    pi_to_gri_map_t last_particle_to_grid_map = gen_particle_to_grid_map();
    pi_to_gri_map_t curr_particle_to_grid_map = gen_particle_to_grid_map();

    pi_to_pa_map_t particle_idx_to_addr_map = gen_particle_idx_to_addr_map();

    initialize_dam_break(grid_to_particle_list_map,
                         last_particle_to_grid_map,
                         curr_particle_to_grid_map,
                         particle_idx_to_addr_map);





    assert(host_grid_consistency_check(grid_to_particle_list_map));
}

