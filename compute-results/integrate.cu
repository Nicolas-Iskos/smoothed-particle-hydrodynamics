#include "../simulation-parameters.h"

#include "particle-data-structures.h"

#include "integrate.h"

__global__ void euler_integrate(pi_to_pa_map_t particle_idx_to_addr_map) {

    Particle *particle;
    uint32_t particle_idx;

    particle_idx = blockDim.x * blockIdx.x + threadIdx.x;

    if(particle_idx >= N_PARTICLES) {
        return;
    }

    particle = particle_idx_to_addr_map[particle_idx];

    /* calculate velocities and positions after the time step */
    /* rf = ri + vi * dt */
    particle->position[0] += particle->velocity[0] * DT;
    particle->position[1] += particle->velocity[1] * DT;
    particle->position[2] += particle->velocity[2] * DT;

    /* vf = vi + ai * dt */
    particle->velocity[0] += (particle->force[0] / M_PARTICLE) * DT;
    particle->velocity[1] += (particle->force[1] / M_PARTICLE) * DT;
    particle->velocity[2] += (particle->force[2] / M_PARTICLE) * DT;
}


__global__ void leapfrog_integrate(pi_to_pa_map_t particle_idx_to_addr_map) {}


__global__ void enforce_boundary_conditions(pi_to_pa_map_t particle_idx_to_addr_map) {

    uint32_t particle_idx;
    Particle *particle;
    constexpr float max_lim = EXP_SPACE_DIM - H;
    constexpr float protection_term = H * H * EPSILON;

    particle_idx = blockDim.x * blockIdx.x + threadIdx.x;

    if(particle_idx >= N_PARTICLES) {
        return;
    }

    particle = particle_idx_to_addr_map[particle_idx];

    /* for each component of position, ensure that each particle cannot
     * enter the last layer of grid spaces coating the outside of the experiment
     * space. If this condition is violated, produce an inelastic collision between
     * the particle and the wall. Notice that we use protection_term  to
     * ensure that the particle is not updated to be included within the outer
     * layer of grid spaces. This also ensures that the grid list map does not
     * need to be updated after boundary conditions are enforced.
     * */
    for(uint8_t ax = 0; ax < 3; ax++) {
        if(particle->position[ax] >= max_lim) {
            particle->position[ax] = max_lim - protection_term;
            particle->velocity[ax] = - DAMPING_FACTOR * particle->velocity[ax];
        }
        else if(particle->position[ax] < H) {
            particle->position[ax] = H + protection_term;
            particle->velocity[ax] = - DAMPING_FACTOR * particle->velocity[ax];
        }
    }
}

