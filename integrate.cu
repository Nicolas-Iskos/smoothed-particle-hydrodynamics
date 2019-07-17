#include "simulation-parameters.h"
#include "particle-data-structures.h"

#include "integrate.h"

__global__ void euler_integrate(pi_to_pa_map_t particle_idx_to_addr_map){
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



