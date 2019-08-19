/*
 * This file contains a set of functions that are designed to
 * calculate the density, pressure and net force acting on each
 * particle within the simulation.
 *
 * These functions work using the computational method of
 * smoothed particle hydrodynamics (SPH). This method involves
 * computing the value of some field at the position of each
 * particle by considering the value of certain fields for surrounding
 * particles within the distance of the smoothing radius H.
 * */



#include "../simulation-parameters.h"

#include "particle-data-structures.h"
#include "smoothing-kernels.h"

#include "calculate-field.h"



/*
 * calculate_density calculates the density at the position of each particle in
 * the simulation and then sets the density parameter of each particle.
 * This function is a CUDA kernel that processes each particle
 * using a single thread.
 *
 * Inputs:
 * grid_to_particle_list_map - an array mapping the index of grid space to a
 *                             doubly-linked-list containing the complete
 *                             set of particles residing in that grid space
 *
 * particle_to_grid_map - an array mapping the index of each particle to
 *                        the index of the grid space it resides in
 *
 * particle_idx_to_addr_map - an array mapping the index of each particle to its
 *                            address in CUDA unified memory
 *
 * Outputs: None
 * */
__global__ void calculate_density(gri_to_pl_map_t grid_to_particle_list_map,
                                  pi_to_gri_map_t particle_to_grid_map,
                                  pi_to_pa_map_t particle_idx_to_addr_map) {
    /* the particle for which we want to calculate the density */
    Particle *curr_particle;
    uint32_t curr_idx;
    uint32_t curr_grid_idx;
    uint32_t curr_grid_pos[3];

    /* the particle whose properties we are accumulating to facilitate
     * the calculation of density
     * */
    Particle *acc_particle;
    uint32_t acc_grid_idx;
    uint32_t acc_grid_pos[3];
    Particle *acc_grid_list;

    /* the final density value calculated by this function */
    float total_density;

    /* the geometrical offsets from the grid space that curr_particle resides
     * in
     * */
    int8_t col_offset;
    int8_t row_offset;
    int8_t layer_offset;



    curr_idx = blockDim.x * blockIdx.x + threadIdx.x;
    if(curr_idx >= N_PARTICLES) {
        return;
    }

    curr_particle = particle_idx_to_addr_map[curr_idx];
    curr_grid_idx = particle_to_grid_map[curr_idx];

    /* convert the row-major grid space index to coordinates in terms
     * of column, row and layer
     * */
    grid_idx_to_grid_pos(curr_grid_idx, curr_grid_pos);

    total_density = 0;



    /* Surrounding the cubic grid space in which the particle resides,
     * there are a total of 26 grid spaces: 9 above, 9 below and 8 on
     * the sides. Because each grid space has side length H, the set
     * of grid spaces that could contain particles that are within
     * curr_particle's smoothing radius H is the 26 surrounding grid
     * spaces in addition to the grid space the particle is currently
     * in, for a total of 27 grid spaces. We accumulate fields from
     * each of the particles in these grid spaces to calculate the
     * density at curr_particle.
     * */
    for(uint8_t i = 0; i < 27; i++) {

        /* choose one of the 27 grid spaces to accumulate particles */
        col_offset = (i % 3) - 1;
        row_offset = ((i % 9) / 3) - 1;
        layer_offset = (i / 9) - 1;

        /* acquire the particles within the selected grid space */
        acc_grid_pos[0] = curr_grid_pos[0] + col_offset;
        acc_grid_pos[1] = curr_grid_pos[1] + row_offset;
        acc_grid_pos[2] = curr_grid_pos[2] + layer_offset;
        acc_grid_idx = grid_pos_to_grid_idx(acc_grid_pos);
        acc_grid_list = grid_to_particle_list_map[acc_grid_idx];

        /* consider the contribution of each particle within the
         * selected grid space
         * */
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



/*
 * calculate_pressure calculates the pressure experienced at the position
 * of each particle and then sets the pressure variable of each particle.
 * This function is a CUDA kernel where each thread performs this task
 * for a single particle.
 *
 * Inputs:
 * particle_idx_to_addr_map - an array mapping the index of each particle to its
 *                            address in CUDA unified memory
 *
 * Outputs: None
 * */
__global__ void calculate_pressure(pi_to_pa_map_t particle_idx_to_addr_map) {

    /* the particle for which we wish to calculate pressure */
    Particle *particle;
    uint32_t particle_idx;



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



/*
 * calculate_net_force computes the net force encountered at the position of
 * each particle and then sets the force variable within each particle.
 *
 * Inputs:
 * grid_to_particle_list_map - an array mapping the index of grid space to a
 *                             doubly-linked-list containing the complete
 *                             set of particles residing in that grid space
 *
 * particle_to_grid_map - an array mapping the index of each particle to
 *                       the index of the grid space it resides in
 *
 * particle_idx_to_addr_map - an array mapping the index of each particle to its
 *                            address in CUDA unified memory
 *
 * Outputs: None
 *
 * */
__global__ void calculate_net_force(gri_to_pl_map_t grid_to_particle_list_map,
                                    pi_to_gri_map_t particle_to_grid_map,
                                    pi_to_pa_map_t particle_idx_to_addr_map) {

    /* the particle for which we wish to calculate the net force */
    Particle *curr_particle;
    uint32_t curr_idx;
    uint32_t curr_grid_idx;
    uint32_t curr_grid_pos[3];

    /* the particle whose fields we are accumulating to calculate
     * pressure
     * */
    Particle *acc_particle;
    uint32_t acc_grid_idx;
    uint32_t acc_grid_pos[3];
    Particle *acc_grid_list;

    /* vectors describing the relative positioning
     * of the current and accumulated particle
     * */
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

    /* Surrounding the cubic grid space in which the particle resides,
     * there are a total of 26 grid spaces: 9 above, 9 below and 8 on
     * the sides. Because each grid space has side length H, the set
     * of grid spaces that could contain particles that are within
     * curr_particle's smoothing radius H is the 26 surrounding grid
     * spaces in addition to the grid space the particle is currently
     * in, for a total of 27 grid spaces. We accumulate fields from
     * each of the particles in these grid spaces to calculate the
     * density at curr_particle.
     * */
    for(uint8_t i = 0; i < 27; i++) {

        /* choose a grid space for which we want to accumulate
         * the fields of constituent particles
         * */
        col_offset = (i % 3) - 1;
        row_offset = ((i % 9) / 3) - 1;
        layer_offset = (i / 9) - 1;

        /* acquire the particles within the selected grid space */
        acc_grid_pos[0] = curr_grid_pos[0] + col_offset;
        acc_grid_pos[1] = curr_grid_pos[1] + row_offset;
        acc_grid_pos[2] = curr_grid_pos[2] + layer_offset;
        acc_grid_idx = grid_pos_to_grid_idx(acc_grid_pos);
        acc_grid_list = grid_to_particle_list_map[acc_grid_idx];

        /* accumulate the contributions from each particle in the
         * grid space to the force due to pressure and the force
         * due to viscosity.
         * */
        for(acc_particle = acc_grid_list;
            acc_particle != NULL;
            acc_particle = acc_particle->next_particle) {

            /* calculate the vectors that define the relative position
             * of the curr_particle and acc_particle
             * */
            for(uint8_t ax = 0; ax < 3; ax++) {
                r_curr_acc[ax] = acc_particle->position[ax] -
                                 curr_particle->position[ax];
            }
            get_norm_3vector(r_curr_acc, r_curr_acc_norm);
            mag_r_curr_acc = get_mag_3vector(r_curr_acc);

            /* add the force contribution from pressure */
            add_f_contr_from_pressure(curr_particle, acc_particle,
                                      r_curr_acc_norm, total_force);

            /* add the force contribution from viscosity */
            add_f_contr_from_viscosity(curr_particle, acc_particle,
                                       r_curr_acc_norm, mag_r_curr_acc,
                                       total_force);
            }
    }

    /* add the contribution due to the gravitational force -
     * this is independent of the accumulated particles
     * */
    add_f_contr_from_gravity(total_force);

    for(uint8_t ax = 0; ax <3; ax++) {
        curr_particle->force[ax] = total_force[ax];
    }
}



/* add_f_contr_from_pressure computes the contribution on the
 * force acting on curr_particle from acc_particle as a result
 * of the fluid pressure exerted and then adds this contribution
 * to total_force.
 *
 * Inputs:
 * curr_particle - the particle for which we are accumulating force
 *                 due to pressure from acc_particle
 *
 * acc_particle - the particle that is contributing to total_force
 *                as a result of the pressure it exerts on curr_particle
 *
 * r_curr_acc_norm - a 3D vector pointing from curr_particle to
 *                   acc_particle
 *
 * total_force - the so far computed force acting on curr_particle
 *
 * Outputs: None
 * */
__device__ void add_f_contr_from_pressure(Particle *curr_particle,
                                          Particle *acc_particle,
                                          float *r_curr_acc_norm,
                                          float *total_force) {

    float curr_rho_pow2;
    float curr_p;
    float acc_rho_pow2;
    float acc_p;
    constexpr float m_particle_pow2 = M_PARTICLE * M_PARTICLE;

    /* used to ensure we never have division by zero */
    constexpr float protection_term = H * H * EPSILON;

    curr_rho_pow2 = pow(curr_particle->density, 2);
    curr_p = curr_particle->pressure;
    acc_rho_pow2 = pow(acc_particle->density, 2);
    acc_p = acc_particle->pressure;

    /* calculate the effect of acc_particle on curr_particle due
     * to pressure in each dimension using the SPH equations
     * */
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



/* add_f_contr_from_viscosity computes the force contribution of
 * acc_particle on curr_particle as a result of the viscous force
 * and then adds this contribution to total_force.
 *
 * Inputs:
 * curr_particle - the particle for which we are accumulating force
 *                 due to pressure from acc_particle
 *
 * acc_particle - the particle that is contributing to total_force
 *                as a result of the pressure it exerts on curr_particle
 *
 * r_curr_acc_norm - a 3D vector pointing from curr_particle to
 *                   acc_particle
 *
 * total_force - the so far computed force acting on curr_particle
 *
 * Outputs: None
 * */
__device__ void add_f_contr_from_viscosity(Particle *curr_particle,
                                           Particle *acc_particle,
                                           float *r_curr_acc_norm,
                                           float mag_r_curr_acc,
                                           float *total_force) {

    float v_curr_acc[3]; /* diff in velocity vector */
    float r_curr_acc[3]; /* diff in position vector */
    float mag_r_curr_acc_pow2;
    float v_curr_acc_dot_r_curr_acc;
    float grad_v_curr_acc;
    float grad_v_curr_acc_pow2;
    float avg_p_curr_acc;
    float c_curr; /* speed of sound for curr_particle */
    float c_acc; /* speed of sound for acc_particle */
    float avg_c_curr_acc;

    /* this term ensures that we never have division by zero */
    constexpr float protection_term = EPSILON * H * H;

    mag_r_curr_acc_pow2 = pow(mag_r_curr_acc, 2);
    v_curr_acc_dot_r_curr_acc = 0;

    /* compute the dot product of the diff in velocity vector and the
     * diff in position vector
     * */
    for(uint8_t ax = 0; ax < 3; ax++) {
        v_curr_acc[ax] = curr_particle->velocity[ax] - acc_particle->velocity[ax];
        r_curr_acc[ax] = curr_particle->position[ax] - acc_particle->position[ax];
        v_curr_acc_dot_r_curr_acc += (v_curr_acc[ax] * r_curr_acc[ax]);
    }

    /* if the flow is expanding (dot product is positive),
     * then the viscous force is zero, so there is no
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

    /* Using the viscosity treatment of INSERT_SCIENTISTS, we use the rules of
     * SPH to calculate the force contribution of acc_particle on curr_particle
     * as a result of the viscous force.
     * */
    for(uint8_t ax = 0; ax < 3; ax++) {
        total_force[ax] -=
        r_curr_acc_norm[ax] *
        (-A_SPH * avg_c_curr_acc * grad_v_curr_acc + B_SPH * grad_v_curr_acc_pow2) /
        (avg_p_curr_acc + protection_term);
    }
}



/* add_f_contr_from_gravity calculates the force on each particle due to
 * gravity and then adds this contribution to total_force.
 *
 * Inputs:
 * total_force - the total force acting on curr_particle calculated so far
 *
 * Outputs: None
 * */
__device__ void add_f_contr_from_gravity(float *total_force) {

    constexpr float g_force = G * M_PARTICLE;

    /* the force of gravity only affects the vertical component
     * of force
     * */
    total_force[2] -= g_force;
}



/*
 * get_norm_3vector computes a normalized version of the
 * 3D vector vec and stores it in the 3D vector norm_vec.
 * */
__device__ void get_norm_3vector(float *vec, float *norm_vec) {

    float mag;
    constexpr float protection_term = H * H * EPSILON;

    mag = get_mag_3vector(vec);

    norm_vec[0] = vec[0] / (mag + protection_term);
    norm_vec[1] = vec[1] / (mag + protection_term);
    norm_vec[2] = vec[2] / (mag + protection_term);
}



/* get_mag_3vector computes and returns the magnitude of vec. */
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
