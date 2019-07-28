#include "simulation-parameters.h"
#include "particle-data-structures.h"

#include "integrate.h"

__global__ void euler_integrate(pi_to_pa_map_t particle_idx_to_addr_map) {

    Particle *particle;
    uint32_t particle_idx;

    particle_idx = blockDim.x * blockIdx.x + threadIdx.x;
    particle = particle_idx_to_addr_map[particle_idx];

    /* calculate velocities and positions after the time step */
    /* rf = ri + vi * dt */
    particle->position[0] = particle->velocity[0] * DT;
    particle->position[1] = particle->velocity[1] * DT;
    particle->position[2] = particle->velocity[2] * DT;

    /* vf = vi + ai * dt */
    particle->velocity[0] = (particle->force[0] / M_PARTICLE) * DT;
    particle->velocity[1] = (particle->force[1] / M_PARTICLE) * DT;
    particle->velocity[2] = (particle->force[2] / M_PARTICLE) * DT;
}


__global__ void leapfrog_integrate(pi_to_pa_map_t particle_idx_to_addr_map) {}


__global__ void enforce_boundary_conditions(pi_to_pa_map_t particle_idx_to_addr_map) {

    uint32_t particle_idx;
    Particle *particle;
    constexpr float max_lim = EXP_SPACE_DIM - H;

    particle_idx = blockDim.x * blockIdx.x + threadIdx.x;
    particle = particle_idx_to_addr_map[particle_idx];

    /* for each component of position, ensure that each particle cannot
     * enter the last layer of grid spaces coating the outside of the experiment
     * space. If this condition is violated, produce an elastic collision between
     * the particle and the wall
     * */
    for(uint8_t ax = 0; ax < 3; ax++) {
        if(particle->position[ax] > max_lim) {
            particle->position[ax] = max_lim;
            particle->velocity[ax] = -particle->velocity[ax];
        }
        else if(particle->position[ax] < H) {
            particle->position[ax] = H;
            particle->velocity[ax] = -particle->velocity[ax];
        }
    }
}

