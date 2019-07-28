#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "calculate-field.h"
#include "integrate.h"

#include "test-functions.h"
#include <cstdio>
#include <cassert>

int main() {
    constexpr uint32_t n_grid_spaces = (EXP_SPACE_DIM * EXP_SPACE_DIM * EXP_SPACE_DIM) /
                                       (H * H * H);

    gri_to_pl_map_t grid_to_particle_list_map = gen_grid_to_particle_list_map();

    pi_to_gri_map_t last_particle_to_grid_map = gen_particle_to_grid_map();
    pi_to_gri_map_t curr_particle_to_grid_map = gen_particle_to_grid_map();

    pi_to_pa_map_t particle_idx_to_addr_map = gen_particle_idx_to_addr_map();

    initialize_dam_break(grid_to_particle_list_map, last_particle_to_grid_map,
                         curr_particle_to_grid_map, particle_idx_to_addr_map);


    for(uint8_t i = 0; i < 1; i++) {

        /* perform SPH calculations to update the forces acting on each
         * particle
         * */
        calculate_density<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                            PARTICLES_PER_BLOCK>>>(grid_to_particle_list_map,
                                                   curr_particle_to_grid_map,
                                                   particle_idx_to_addr_map);

        calculate_pressure<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                             PARTICLES_PER_BLOCK>>>(particle_idx_to_addr_map);

        calculate_net_force<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                              PARTICLES_PER_BLOCK>>>(grid_to_particle_list_map,
                                                     curr_particle_to_grid_map,
                                                     particle_idx_to_addr_map);
        /* adjust the position and velocity of the particles based on the
         * recently-computed net force acting on each particle
         * */
        euler_integrate<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                          PARTICLES_PER_BLOCK>>>(particle_idx_to_addr_map);

        /* ensure that no particle passes into the outer layer of grid spaces
         * in the experimental space, or out of the experimental space
         * entirely
         * */
        enforce_boundary_conditions<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                                      PARTICLES_PER_BLOCK>>>(particle_idx_to_addr_map);

        /* update the particle grid in 3 steps:
         *
         * 1. Keep track of what grid spaces the particles were last in, but also
         *    recompute what their new grid spaces should be and record these
         *    grid spaces as well
         *
         * 2. For each grid space, remove all particles that have left that grid
         *    space.
         *
         * 3. For each grid space, add all particles that have entered that
         *    grid space.
         * */
        update_particle_to_grid_map<<<n_grid_spaces / GRID_SPACES_PER_BLOCK,
                                      GRID_SPACES_PER_BLOCK>>>(last_particle_to_grid_map,
                                                               curr_particle_to_grid_map,
                                                               particle_idx_to_addr_map);

        perform_removals_from_grid<<<n_grid_spaces / GRID_SPACES_PER_BLOCK,
                                     GRID_SPACES_PER_BLOCK>>>(grid_to_particle_list_map,
                                                              last_particle_to_grid_map,
                                                              curr_particle_to_grid_map,
                                                              particle_idx_to_addr_map);

        perform_additions_to_grid<<<n_grid_spaces / GRID_SPACES_PER_BLOCK,
                                     GRID_SPACES_PER_BLOCK>>>(grid_to_particle_list_map,
                                                              last_particle_to_grid_map,
                                                              curr_particle_to_grid_map,
                                                              particle_idx_to_addr_map);


    }
}

