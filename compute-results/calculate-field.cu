#include "../simulation-parameters.h"

#include "particle-data-structures.h"
#include "smoothing-kernels.h"

#include "calculate-field.h"



__global__ void calculate_density(gri_to_pl_map_t grid_to_particle_list_map,
                                  pi_to_gri_map_t particle_to_grid_map,
                                  pi_to_pa_map_t particle_idx_to_addr_map) {

    Particle *curr_particle;
    uint32_t curr_idx;
    uint32_t curr_grid_idx;
    uint32_t curr_grid_pos[3];

    Particle *acc_particle;
    uint32_t acc_grid_idx;
    uint32_t acc_grid_pos[3];
    Particle *acc_grid_list;

    float total_density;

    int8_t col_offset;
    int8_t row_offset;
    int8_t layer_offset;

    curr_idx = blockDim.x * blockIdx.x + threadIdx.x;

    if(curr_idx >= N_PARTICLES) {
        return;
    }

    curr_particle = particle_idx_to_addr_map[curr_idx];
    curr_grid_idx = particle_to_grid_map[curr_idx];
    grid_idx_to_grid_pos(curr_grid_idx, curr_grid_pos);
    total_density = 0;


    for(uint8_t i = 0; i < 27; i++) {

        col_offset = (i % 3) - 1;
        row_offset = ((i % 9) / 3) - 1;
        layer_offset = (i / 9) - 1;

        /* acquire the particles within the selected grid space */
        acc_grid_pos[0] = curr_grid_pos[0] + col_offset;
        acc_grid_pos[1] = curr_grid_pos[1] + row_offset;
        acc_grid_pos[2] = curr_grid_pos[2] + layer_offset;
        acc_grid_idx = grid_pos_to_grid_idx(acc_grid_pos);
        acc_grid_list = grid_to_particle_list_map[acc_grid_idx];


        for(acc_particle = acc_grid_list;
            acc_particle != NULL;
            acc_particle = acc_particle->next_particle) {

            /*
             * Sum the density contributions from the particular grid space
             * selected by the outer 3 loops
             * */
            total_density +=
            M_PARTICLE * cubic_spline_kernel(curr_particle->position,
                                             acc_particle->position);
        }
    }
    curr_particle->density = total_density;
}


__global__ void calculate_pressure(pi_to_pa_map_t particle_idx_to_addr_map) {

    uint32_t particle_idx;
    Particle *particle;

    particle_idx = blockDim.x * blockIdx.x + threadIdx.x;

    if(particle_idx >= N_PARTICLES) {
        return;
    }

    particle = particle_idx_to_addr_map[particle_idx];

    /* here we use the ideal gas law as the equation of state
     * relating pressure to density and temperature
     * */
    particle->pressure = (GAMMA - 1) * H_SPEC * particle->density;
}


__global__ void calculate_net_force(gri_to_pl_map_t grid_to_particle_list_map,
                                    pi_to_gri_map_t particle_to_grid_map,
                                    pi_to_pa_map_t particle_idx_to_addr_map) {

    Particle *curr_particle;
    uint32_t curr_idx;
    uint32_t curr_grid_idx;
    uint32_t curr_grid_pos[3];

    Particle *acc_particle;
    uint32_t acc_grid_idx;
    uint32_t acc_grid_pos[3];
    Particle *acc_grid_list;

    float r_curr_acc[3];
    float r_curr_acc_norm[3];
    float mag_r_curr_acc;
    float total_force[3];

    int8_t col_offset;
    int8_t row_offset;
    int8_t layer_offset;

    curr_idx = blockDim.x * blockIdx.x + threadIdx.x;

    if(curr_idx >= N_PARTICLES) {
        return;
    }

    curr_particle = particle_idx_to_addr_map[curr_idx];
    curr_grid_idx = particle_to_grid_map[curr_idx];
    grid_idx_to_grid_pos(curr_grid_idx, curr_grid_pos);
    memset(total_force, 0, 3 * sizeof(float));

    for(uint8_t i = 0; i < 27; i++) {

        col_offset = (i % 3) - 1;
        row_offset = ((i % 9) / 3) - 1;
        layer_offset = (i / 9) - 1;

        /* acquire the particles within the selected grid space */
        acc_grid_pos[0] = curr_grid_pos[0] + col_offset;
        acc_grid_pos[1] = curr_grid_pos[1] + row_offset;
        acc_grid_pos[2] = curr_grid_pos[2] + layer_offset;
        acc_grid_idx = grid_pos_to_grid_idx(acc_grid_pos);
        acc_grid_list = grid_to_particle_list_map[acc_grid_idx];

        for(acc_particle = acc_grid_list;
            acc_particle != NULL;
            acc_particle = acc_particle->next_particle) {

            for(uint8_t ax = 0; ax < 3; ax++) {
                r_curr_acc[ax] = acc_particle->position[ax] -
                                 curr_particle->position[ax];
            }
            get_norm_3vector(r_curr_acc, r_curr_acc_norm);
            mag_r_curr_acc = get_mag_3vector(r_curr_acc);

            add_f_contr_from_pressure(curr_particle, acc_particle,
                                      r_curr_acc_norm, total_force);

            add_f_contr_from_viscosity(curr_particle, acc_particle,
                                       r_curr_acc_norm, mag_r_curr_acc,
                                       total_force);

        }
    }

    add_f_contr_from_gravity(total_force);

    for(uint8_t ax = 0; ax <3; ax++) {
        curr_particle->force[ax] = total_force[ax];
    }
}


__device__ void add_f_contr_from_pressure(Particle *curr_particle,
                                          Particle *acc_particle,
                                          float *r_curr_acc_norm,
                                          float *total_force) {

    float curr_rho_pow2;
    float curr_p;
    float acc_rho_pow2;
    float acc_p;
    constexpr float m_particle_pow2 = M_PARTICLE * M_PARTICLE;
    constexpr float protection_term = H * H * EPSILON;

    curr_rho_pow2 = pow(curr_particle->density, 2);
    curr_p = curr_particle->pressure;
    acc_rho_pow2 = pow(acc_particle->density, 2);
    acc_p = acc_particle->pressure;

    for(uint8_t ax = 0; ax < 3; ax++) {

        total_force[ax] -=
            r_curr_acc_norm[ax] *
            m_particle_pow2 *
            (curr_p / (curr_rho_pow2 + protection_term) +
            acc_p / (acc_rho_pow2 + protection_term)) *
            cubic_spline_kernel(curr_particle->position,
                                acc_particle->position);
    }
}


__device__ void add_f_contr_from_viscosity(Particle *curr_particle,
                                           Particle *acc_particle,
                                           float *r_curr_acc_norm,
                                           float mag_r_curr_acc,
                                           float *total_force) {

    float v_curr_acc[3];
    float r_curr_acc[3];
    float mag_r_curr_acc_pow2;
    float v_curr_acc_dot_r_curr_acc;
    float grad_v_curr_acc;
    float grad_v_curr_acc_pow2;
    float avg_p_curr_acc;
    float c_curr;
    float c_acc;
    float avg_c_curr_acc;

    constexpr float protection_term = EPSILON * H * H;

    mag_r_curr_acc_pow2 = pow(mag_r_curr_acc, 2);

    for(uint8_t ax = 0; ax < 3; ax++) {
        v_curr_acc[ax] = curr_particle->velocity[ax] - acc_particle->velocity[ax];
        r_curr_acc[ax] = curr_particle->position[ax] - acc_particle->position[ax];
        v_curr_acc_dot_r_curr_acc += (v_curr_acc[ax] * r_curr_acc[ax]);
    }

    /* if the flow is expanding, then the viscous force is zero, so there is no
     * need to compute the rest of the function
     * */
    if(v_curr_acc_dot_r_curr_acc >= 0) {
        return;
    }

    grad_v_curr_acc =
        H * v_curr_acc_dot_r_curr_acc / (mag_r_curr_acc_pow2 + protection_term);
    grad_v_curr_acc_pow2 = pow(grad_v_curr_acc, 2);

    avg_p_curr_acc = (curr_particle->density + acc_particle->density) / 2;

    c_curr = sqrt(GAMMA * curr_particle->pressure /
             (curr_particle->density + protection_term));
    c_acc = sqrt(GAMMA * acc_particle->pressure /
            (acc_particle->density + protection_term));
    avg_c_curr_acc = (c_curr + c_acc) / 2;

    for(uint8_t ax = 0; ax < 3; ax++) {
        total_force[ax] =
        r_curr_acc_norm[ax] *
        (-A_SPH * avg_c_curr_acc * grad_v_curr_acc + B_SPH * grad_v_curr_acc_pow2) /
        avg_p_curr_acc;
    }
}

__device__ void add_f_contr_from_gravity(float *total_force) {

    constexpr float g_force = G * M_PARTICLE;
    total_force[2] -= g_force;
}

__device__ void get_norm_3vector(float *vec, float *norm_vec) {

    float mag;
    constexpr float protection_term = H * H * EPSILON;

    mag = get_mag_3vector(vec);

    norm_vec[0] = vec[0] / (mag + protection_term);
    norm_vec[1] = vec[1] / (mag + protection_term);
    norm_vec[2] = vec[2] / (mag + protection_term);
}


__device__ float get_mag_3vector(float *vec) {

    float x_val_pow2;
    float y_val_pow2;
    float z_val_pow2;
    float mag;

    x_val_pow2 = pow(vec[0], 2);
    y_val_pow2 = pow(vec[1], 2);
    z_val_pow2 = pow(vec[2], 2);
    mag = sqrt(x_val_pow2 + y_val_pow2 + z_val_pow2);

    return mag;
}
