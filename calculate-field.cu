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
                    curr_particle_new_density +=
                        M_PARTICLE * cubic_spline_kernel(curr_particle->position,
                                                         acc_particle->position);
                }
            }
        }
    }
    curr_particle->density = curr_particle_new_density;
}


__global__ void calculate_pressure(pi_to_pa_map_t particle_idx_to_addr_map) {

    uint32_t particle_idx;
    Particle *particle;

    particle_idx = blockDim.x * blockIdx.x + threadIdx.x;
    particle = particle_idx_to_addr_map[particle_idx];

    /* here we use the ideal gas law as the equation of state
     * relating pressure to density and temperature
     * */
    particle->pressure = particle->density * (R / M) * T;
}


__global__ void calculate_force(gri_to_pl_map_t grid_to_particle_list_map,
                                pi_to_gri_map_t particle_to_grid_map,
                                pi_to_pa_map_t particle_idx_to_addr_map) {


    Particle *curr_particle;
    uint32_t curr_particle_idx;
    uint32_t curr_particle_grid_idx;
    uint32_t curr_particle_grid_pos[3];
    float curr_rho_pow2;
    float curr_p;
    float curr_particle_new_force[3];

    Particle *acc_particle;
    uint32_t acc_particle_grid_idx;
    uint32_t acc_particle_grid_pos[3];
    float acc_rho_pow2;
    float acc_p;
    Particle *acc_particle_grid_list;

    float diff_pos_curr_acc[3];
    float norm_diff_pos_curr_acc[3];

    constexpr float m_particle_pow2 = M_PARTICLE * M_PARTICLE;

    curr_particle_idx = blockDim.x * blockIdx.x + threadIdx.x;

    curr_particle = particle_idx_to_addr_map[curr_particle_idx];
    curr_particle_grid_idx = particle_to_grid_map[curr_particle_idx];
    grid_idx_to_grid_pos(curr_particle_grid_idx, curr_particle_grid_pos);
    curr_rho_pow2 = pow(curr_particle->density, 2);
    curr_p = curr_particle->pressure;
    memset(curr_particle_new_force, 0, 3 * sizeof(float));

    /* iterate through columns of 3 x 3 x 3 block of grid spaces */
    for(int8_t layer = -1; layer <= 1; layer++) {
        /* iterate through rows of 3 x 3 x 3 block of grid spaces */
        for(int8_t row = -1; row <= 1; row++) {
            /* iterate through layers of 3 x 3 x 3 block of grid spaces */
            for(int8_t col = -1; col <= 1; col++) {
                /* acquire the particles within the selected grid space */
                acc_particle_grid_pos[0] = curr_particle_grid_pos[0] + col;
                acc_particle_grid_pos[1] = curr_particle_grid_pos[1] + row;
                acc_particle_grid_pos[2] = curr_particle_grid_pos[2] + layer;
                acc_particle_grid_idx = grid_pos_to_grid_idx(acc_particle_grid_pos);
                acc_particle_grid_list = grid_to_particle_list_map[acc_particle_grid_idx];

                for(acc_particle = acc_particle_grid_list;
                    acc_particle != NULL;
                    acc_particle = acc_particle->next_particle) {

                    acc_rho_pow2 = pow(acc_particle->density, 2);
                    acc_p = acc_particle->pressure;

                    for(uint8_t ax = 0; ax < 3; ax++) {
                        diff_pos_curr_acc[ax] = acc_particle->position[ax] -
                                                curr_particle->position[ax];
                    }

                    get_norm_vector(diff_pos_curr_acc, norm_diff_pos_curr_acc);

                    for(uint8_t ax = 0; ax < 3; ax++) {
                        /*
                         * Sum the density contributions from the particular grid space
                         * selected by the outer 3 loops
                         * */
                        curr_particle_new_force[ax] -=
                            norm_diff_pos_curr_acc[ax] *
                            m_particle_pow2 *
                            (curr_p / curr_rho_pow2 + acc_p / acc_rho_pow2) *
                            cubic_spline_kernel(curr_particle->position,
                                                acc_particle->position);
                    }
                }
            }
        }
    }

    for(uint8_t ax = 0; ax <3; ax++) {
        curr_particle->position[ax] = curr_particle_new_force[ax];
    }
}


__device__ void get_norm_vector(float *vec, float *norm_vec) {

    float x_val;
    float y_val;
    float z_val;
    float mag;

    x_val = vec[0];
    y_val = vec[1];
    z_val = vec[2];
    mag = sqrt(pow(x_val, 2) + pow(y_val, 2) + pow(z_val, 2));

    norm_vec[0] = x_val / mag;
    norm_vec[1] = y_val / mag;
    norm_vec[2] = z_val / mag;
}
