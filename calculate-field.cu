#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "calculate-field.h"

#include "smoothing-kernels.h"


__device__ void calculate_density(uint32_t particle_idx,
                                  pi_to_gri_map_t particle_to_grid_map,
                                  gri_to_pl_map_t grid_to_particle_list_map) {

    uint32_t grid_idx;
    uint32_t grid_pos[3];
    constexpr uint32_t n_grid_spaces_per_dim = (uint32_t)(EXP_SPACE_DIM / H);

    grid_idx = particle_to_grid_map[particle_idx];
    grid_idx_to_grid_pos(grid_idx, grid_pos);

    /* iterate through columns of 3 x 3 x 3 block of grid spaces */
    for(int16_t i = -1; i <= 1; i++) {
        /* iterate through rows of 3 x 3 x 3 block of grid spaces */
        for(int16_t j = -1; i <= 1; j++) {
            /* iterate through layers of 3 x 3 x 3 block of grid spaces */
            for(int16_t k = -1; k <= 1; k++) {
            }
        }
    }
}
