#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "calculate-field.h"

#include "smoothing-kernels.h"


__global__ void calculate_density(gri_to_pl_map_t grid_to_particle_list_map,
                                  pi_to_gri_map_t particle_to_grid_map,
                                  pi_to_pa_map_t particle_idx_to_addr_map) {

    Particle *curr_particle;
    uint32_t curr_particle_idx;
    uint32_t curr_particle_grid_idx;
    uint32_t curr_particle_grid_pos[3];
    float curr_particle_new_density;

    Particle *acc_particle;
    uint32_t acc_particle_grid_idx;
    uint32_t acc_particle_grid_pos[3];
    Particle *acc_particle_grid_list;

    curr_particle_idx = blockDim.x * blockIdx.x + threadIdx.x;
    curr_particle = particle_idx_to_addr_map[curr_particle_idx];
    curr_particle_grid_idx = particle_to_grid_map[curr_particle_idx];
    grid_idx_to_grid_pos(curr_particle_grid_idx, curr_particle_grid_pos);
    curr_particle_new_density = 0;

    /* iterate through columns of 3 x 3 x 3 block of grid spaces */
    for(int16_t layer = -1; layer <= 1; layer++) {
        /* iterate through rows of 3 x 3 x 3 block of grid spaces */
        for(int16_t row = -1; row <= 1; row++) {
            /* iterate through layers of 3 x 3 x 3 block of grid spaces */
            for(int16_t col = -1; col <= 1; col++) {
                /* acquire the particles within the selected grid space */
                acc_particle_grid_pos[0] = curr_particle_grid_pos[0] + col;
                acc_particle_grid_pos[1] = curr_particle_grid_pos[1] + row;
                acc_particle_grid_pos[2] = curr_particle_grid_pos[2] + layer;
                acc_particle_grid_idx = grid_pos_to_grid_idx(acc_particle_grid_pos);
                acc_particle_grid_list = grid_to_particle_list_map[acc_particle_grid_idx];

                for(acc_particle = acc_particle_grid_list;
                    acc_particle != NULL;
                    acc_particle = acc_particle->next_particle) {

                    /*
                     * Sum the density contributions from the particular grid space
                     * selected by the outer 3 loops
                     * */
                    curr_particle_new_density += M_PARTICLE *
                                                 cubic_spline_kernel(
                                                    curr_particle->position,
                                                    acc_particle->position);
                }
            }
        }
    }
    curr_particle->density = curr_particle_new_density;
}
