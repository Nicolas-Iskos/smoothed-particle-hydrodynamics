/*
 * This file contains a set of functions that initialize and manage
 * the particle data structures used in this simulation. The main
 * particle data structures that are managed include:
 *
 * - grid_to_particle_list: an array that maps the index of each grid
 *   space to a doubly-linked list containing each particle residing in
 *   that grid space
 *
 * - particle_to_grid_map: an array mapping the index of each particle to
 *   the index of the grid space that particle resides in
 *
 * - particle_idx_to_addr_map: an array mapping the index of each particle
 *   to its address in CUDA unified memory
 *
 **/



#include "../simulation-parameters.h"

#include "particle-data-structures.h"

#include <cstdint>
#include <stdlib.h>
#include <time.h>



/*
 * returns a pointer to a dynamically allocated
 * particle index to grid index map
 * */
pi_to_gri_map_t gen_particle_to_grid_map() {
    pi_to_gri_map_t particle_to_grid_map;

    cudaMallocManaged(&particle_to_grid_map,
                      N_PARTICLES * sizeof(uint32_t));

    return particle_to_grid_map;
}



/*
 * returns a pointer to a dynamically allocated
 * grid to particle list map
 * */
gri_to_pl_map_t gen_grid_to_particle_list_map() {
    gri_to_pl_map_t grid_to_particle_list_map;
    uint32_t n_grid_spaces;

    n_grid_spaces = (uint32_t)pow(EXP_SPACE_DIM / H, 3);

    cudaMallocManaged(&grid_to_particle_list_map,
            n_grid_spaces * sizeof(Particle*));

    /* ensure that each doubly linked list is NULL-terminated */
    for(size_t i = 0; i < n_grid_spaces; i++) {
        grid_to_particle_list_map[i] = NULL;
    }

    return grid_to_particle_list_map;
}



/*
 * returns a pointer to a dynamically allocated
 * particle index to address map
 * */
pi_to_pa_map_t gen_particle_idx_to_addr_map() {
    pi_to_pa_map_t particle_idx_to_addr_map;

    cudaMallocManaged(&particle_idx_to_addr_map,
            N_PARTICLES * sizeof(Particle*));

    return particle_idx_to_addr_map;
}



/*
 * initialize_dam_break creates a cubic block
 * of particles arranged into a simple cubic lattice centered in
 * the experiment space and sets the input data structures
 * accordingly.
 *
 * Inputs:
 * grid_to_particle_list_map - an array mapping the index of grid space to a
 *                             doubly-linked-list containing the complete
 *                             set of particles residing in that grid space
 *
 * last_particle_to_grid_map - an array mapping the index of each particle to
 *                             the index of the grid space it resided in during
 *                             the previous iteration of the simulation
 *
 * last_particle_to_grid_map - an array mapping the index of each particle to
 *                             the index of the grid space it currently resides
 *                             in
 *
 * particle_idx_to_addr_map - an array mapping the index of each particle to its
 *                            address in CUDA unified memory
 *
 * Outputs: None
 * */
void initialize_dam_break(gri_to_pl_map_t grid_to_particle_list_map,
                          pi_to_gri_map_t last_particle_to_grid_map,
                          pi_to_gri_map_t curr_particle_to_grid_map,
                          pi_to_pa_map_t particle_idx_to_addr_map) {

    /* number of particles forming the side length of the cube lattice */
    uint32_t n_particles_per_dim;
    uint32_t n_particles_per_dim_pow2;
    float cubic_block_rad;
    float particle_spacing;
    float space_center[3];
    float init_particle_pos[3];
    float particle_pos[3];
    uint32_t grid_idx;
    uint32_t particle_idx;
    Particle *new_particle;



    /* initialize random seed to change particle initialization with each
     * run of the program
     * */
    srand(time(NULL));

    /* height_factor describes the fraction by which we move the cubic lattice
     * closer to the ground. A value of 0 means the block is centered, and
     * a value of 1 means the block starts on the ground.
     * */
    constexpr float height_factor = 0.6;
    n_particles_per_dim = (uint32_t)cbrt((float)N_PARTICLES);
    n_particles_per_dim_pow2 = (uint32_t)pow(n_particles_per_dim, 2);

    cubic_block_rad = n_particles_per_dim * R_PARTICLE;
    particle_spacing = 2 * R_PARTICLE;

    space_center[0] = (float)EXP_SPACE_DIM / 2;
    space_center[1] = (float)EXP_SPACE_DIM / 2;
    space_center[2] = (float)EXP_SPACE_DIM / 2;

    init_particle_pos[0] = space_center[0] + cubic_block_rad - R_PARTICLE;
    init_particle_pos[1] = space_center[1] - cubic_block_rad + R_PARTICLE;
    init_particle_pos[2] = space_center[2] + cubic_block_rad - R_PARTICLE -
                           height_factor *
                           ((float)EXP_SPACE_DIM / 2 - cubic_block_rad);

    /*
     * Arrange each particle into its correct grid slot for the
     * simple cubic lattice arrangement
     *
     * Looking down the x-axis towards the origin, we build each
     * slice of the lattice perpendicular to the x-axis
     * from top left to bottom right, proceeding
     * along horizontal rows. The slices of the lattice are built
     * starting at high x-values, going to low x-values.
     * */
    for(particle_idx = 0; particle_idx < N_PARTICLES; particle_idx++) {
        /* compute the position of the particle to be created */
        particle_pos[0] =
            init_particle_pos[0] - particle_spacing *
            (particle_idx / n_particles_per_dim_pow2);
        particle_pos[1] =
            init_particle_pos[1] + particle_spacing *
            (particle_idx % n_particles_per_dim);
        particle_pos[2] =
            init_particle_pos[2] - particle_spacing *
            ((particle_idx % n_particles_per_dim_pow2) / n_particles_per_dim);

        /* initialize the new particle */
        new_particle = new Particle;

        new_particle->position[0] = particle_pos[0];
        new_particle->position[1] = particle_pos[1];
        new_particle->position[2] = particle_pos[2];

        /* we randomly initialize the particle velocity in each direction
         * with the maximum velocity in any direction the velocity due
         * to gravitational acceleration after 3 timesteps
         * */
        new_particle->velocity[0] = (2 * ((float)rand() / RAND_MAX) - 1) *
                                    G * 3 * DT;
        new_particle->velocity[1] = (2 * ((float)rand() / RAND_MAX) - 1) *
                                    G * 3 * DT;
        new_particle->velocity[2] = (2 * ((float)rand() / RAND_MAX) - 1) *
                                    G * 3 * DT;

        new_particle->force[0] = 0;
        new_particle->force[1] = 0;
        new_particle->force[2] = 0;

        new_particle->density = 0;
        new_particle->pressure = 0;

        new_particle->prev_particle = NULL;
        new_particle->next_particle = NULL;

        /* record the address of the new particle */
        particle_idx_to_addr_map[particle_idx] = new_particle;

        /* record the grid index of each particle */
        grid_idx = particle_pos_to_grid_idx(new_particle->position);
        last_particle_to_grid_map[particle_idx] = grid_idx;
        curr_particle_to_grid_map[particle_idx] = grid_idx;

        /*
         * insert the new particle into the correct grid space and
         * record the grid space of the new particle
         * */
        host_insert_into_grid(grid_to_particle_list_map,
                              grid_idx,
                              new_particle);
    }
}



/*
 * host_insert_into_grid takes a given particle and inserts it into
 * the appropriate grid space. It is to be called from the host.
 *
 * Inputs:
 * grid_to_particle_list_map - an array mapping the index of each grid space
 *                             to a doubly linked list containing the
 *                             particles that reside in that grid space
 *
 * grid_idx - the index of the grid space we want to insert into
 *
 * new_particle - a pointer to the particle we wish to insert to the
 *                grid space at grid_idx
 *
 * Outputs: None
 * */
void host_insert_into_grid(gri_to_pl_map_t grid_to_particle_list_map,
                           uint32_t grid_idx,
                           Particle *new_particle) {

    Particle *first_particle_in_grid_slot;

    first_particle_in_grid_slot = grid_to_particle_list_map[grid_idx];

    /* add particle to the correct grid space doubly linked list */
    if(first_particle_in_grid_slot == NULL) {
        new_particle->prev_particle = NULL;
        new_particle->next_particle = NULL;
        grid_to_particle_list_map[grid_idx] = new_particle;
    }
    else {
        first_particle_in_grid_slot->prev_particle = new_particle;
        new_particle->next_particle = first_particle_in_grid_slot;
        new_particle->prev_particle = NULL;
        grid_to_particle_list_map[grid_idx] = new_particle;
    }
}


__global__ void update_particle_to_grid_map(
                                pi_to_gri_map_t last_particle_to_grid_map,
                                pi_to_gri_map_t curr_particle_to_grid_map,
                                pi_to_pa_map_t particle_idx_to_addr_map) {
    uint32_t particle_idx;
    uint32_t pre_update_grid_idx;
    uint32_t updated_grid_idx;
    Particle *particle;

    particle_idx = blockDim.x * blockIdx.x + threadIdx.x;

    if(particle_idx >= N_PARTICLES) {
        return;
    }

    particle = particle_idx_to_addr_map[particle_idx];
    pre_update_grid_idx = curr_particle_to_grid_map[particle_idx];
    updated_grid_idx = particle_pos_to_grid_idx(particle->position);

    /* set the pre-updated grid_idx in the last particle to grid map */
    last_particle_to_grid_map[particle_idx] = pre_update_grid_idx;

    /* set the updated grid idx into the current particle to grid map */
    curr_particle_to_grid_map[particle_idx] = updated_grid_idx;
}



/*
 * Based on the differences between last_particle_to_grid_map and
 * curr_particle_to_grid_map, perform_removals_from_grid will
 * remove from each grid space particles that have
 * exited that grid space since the last iteration. It is a CUDA kernel
 * where each thread performs the work of removing relevant
 * particles from a single grid space.
 *
 * Inputs:
 * grid_to_particle_list_map - an array mapping the index of grid space to a
 *                             doubly-linked-list containing the complete
 *                             set of particles residing in that grid space
 *
 * last_particle_to_grid_map - an array mapping the index of each particle to
 *                             the index of the grid space it resided in during
 *                             the previous iteration of the simulation
 *
 * last_particle_to_grid_map - an array mapping the index of each particle to
 *                             the index of the grid space it currently resides
 *                             in
 *
 * particle_idx_to_addr_map - an array mapping the index of each particle to its
 *                            address in CUDA unified memory
 *
 * Outputs: None
 * */
__global__ void perform_removals_from_grid(
                                gri_to_pl_map_t grid_to_particle_list_map,
                                pi_to_gri_map_t last_particle_to_grid_map,
                                pi_to_gri_map_t curr_particle_to_grid_map,
                                pi_to_pa_map_t particle_idx_to_addr_map) {
    uint32_t grid_idx;
    Particle *del_particle;
    Particle *del_prev_particle;
    Particle *del_next_particle;
    constexpr uint32_t n_grid_spaces = (EXP_SPACE_DIM * EXP_SPACE_DIM * EXP_SPACE_DIM) /
                                       (H * H * H);



    grid_idx = blockDim.x * blockIdx.x + threadIdx.x;
    if(grid_idx >= n_grid_spaces) {
        return;
    }



    /* For the grid space indexed by grid_idx, if we see that a particle
     * indexed by particle_idx was previously in the grid space of grid_idx,
     * but is not anymore, we remove that particle from the grid space
     * indexed by grid_idx.
     * */
    for(uint32_t particle_idx = 0; particle_idx < N_PARTICLES; particle_idx++){
        if((last_particle_to_grid_map[particle_idx] == grid_idx) &&
           (curr_particle_to_grid_map[particle_idx] != grid_idx)) {

            del_particle = particle_idx_to_addr_map[particle_idx];
            del_prev_particle = del_particle->prev_particle;
            del_next_particle = del_particle->next_particle;

            del_particle->next_particle = NULL;
            del_particle->prev_particle = NULL;

            /* remove particle from the grid space at grid_idx */
            if(del_prev_particle == NULL && del_next_particle == NULL) {
                grid_to_particle_list_map[grid_idx] = NULL;
            }
            else if(del_prev_particle == NULL) {
                grid_to_particle_list_map[grid_idx] = del_next_particle;
                del_next_particle->prev_particle = NULL;
            }
            else if(del_next_particle == NULL) {
                del_prev_particle->next_particle = NULL;
            }
            else {
                del_prev_particle->next_particle = del_next_particle;
                del_next_particle->prev_particle = del_prev_particle;
            }
        }
    }
}



/*
 * Based on the differences between last_particle_to_grid_map and
 * curr_particle_to_grid_map, perform_additions_to_grid will
 * add to each grid space particles that have
 * entered that grid space since the last iteration. It is a CUDA kernel
 * where each thread performs the work of adding relevant
 * particles to a single grid space.
 *
 * Inputs:
 * grid_to_particle_list_map - an array mapping the index of grid space to a
 *                             doubly-linked-list containing the complete
 *                             set of particles residing in that grid space
 *
 * last_particle_to_grid_map - an array mapping the index of each particle to
 *                             the index of the grid space it resided in during
 *                             the previous iteration of the simulation
 *
 * last_particle_to_grid_map - an array mapping the index of each particle to
 *                             the index of the grid space it currently resides
 *                             in
 *
 * particle_idx_to_addr_map - an array mapping the index of each particle to its
 *                            address in CUDA unified memory
 *
 * Outputs: None
 * */
__global__ void perform_additions_to_grid(
                                gri_to_pl_map_t grid_to_particle_list_map,
                                pi_to_gri_map_t last_particle_to_grid_map,
                                pi_to_gri_map_t curr_particle_to_grid_map,
                                pi_to_pa_map_t particle_idx_to_addr_map) {
    uint32_t grid_idx;
    Particle *particle;
    Particle *first_particle_in_grid_slot;
    constexpr uint32_t n_grid_spaces = (EXP_SPACE_DIM * EXP_SPACE_DIM * EXP_SPACE_DIM) /
                                       (H * H * H);



    grid_idx = blockDim.x * blockIdx.x + threadIdx.x;
    if(grid_idx >= n_grid_spaces) {
        return;
    }



    /* For the grid space indexed by grid_idx, if we see that a particle
     * indexed by particle_idx was not previously in the grid space of
     * grid_idx but is now, we add that particle to the grid space
     * indexed by grid_idx.
     * */
    for(uint32_t particle_idx = 0; particle_idx < N_PARTICLES; particle_idx++){
        if((last_particle_to_grid_map[particle_idx] != grid_idx) &&
           (curr_particle_to_grid_map[particle_idx] == grid_idx)) {

            particle = particle_idx_to_addr_map[particle_idx];
            first_particle_in_grid_slot = grid_to_particle_list_map[grid_idx];

            /* add particle to the correct grid space doubly linked list */
            if(first_particle_in_grid_slot == NULL) {
                particle->prev_particle = NULL;
                particle->next_particle = NULL;
                grid_to_particle_list_map[grid_idx] = particle;
            }
            else {
                first_particle_in_grid_slot->prev_particle = particle;
                particle->next_particle = first_particle_in_grid_slot;
                particle->prev_particle = NULL;
                grid_to_particle_list_map[grid_idx] = particle;
            }
        }
    }
}



/*
 * particle_pos_to_grid_idx takes a particle index and based in the position
 * of the corresponding particle, returns the grid space that particle
 * resides in.
 * */
__host__ __device__ uint32_t particle_pos_to_grid_idx(float *particle_pos) {
    uint32_t grid_pos[3];

    /* grid space column is related to y coordinate */
    grid_pos[0] = (uint16_t)(particle_pos[1] / H);

    /* grid space row is related to z coordinate */
    grid_pos[1] = (uint16_t)((EXP_SPACE_DIM - particle_pos[2]) / H);

    /* grid space layer is related to x coordinate */
    grid_pos[2] = (uint16_t)((EXP_SPACE_DIM - particle_pos[0]) / H);

    return grid_pos_to_grid_idx(grid_pos);
}



/*
 * Given a column, row and layer of a grid space position,
 * grid_pos_to_grid_idx returns the corresponding grid index
 * */
__host__ __device__ uint32_t grid_pos_to_grid_idx(uint32_t *grid_pos) {

    constexpr uint32_t n_grid_spaces_per_dim = (uint32_t)(EXP_SPACE_DIM / H);
    constexpr uint32_t n_grid_spaces_per_dim_pow2 =
                       (uint32_t)((EXP_SPACE_DIM * EXP_SPACE_DIM) / (H * H));

    /* column + row * num_cols + layer * num_rows * num_cols */
    return grid_pos[0] +
           grid_pos[1] * n_grid_spaces_per_dim +
           grid_pos[2] * n_grid_spaces_per_dim_pow2;
}



/*
 * grid_idx_to_grid_idx takes a grid index and returns the corresponding
 * grid column, row and layer.
 *
 * Inputs:
 * grid_idx - the grid space index we would like to convert to coordinates
 *
 * grid_pos - a 3-element array whose elements correspond to
 *            grid column, grid row and grid layer, respectively
 *
 * Outputs: None
 * */
__host__ __device__ void grid_idx_to_grid_pos(uint32_t grid_idx,
                                              uint32_t *grid_pos) {

    constexpr uint32_t n_grid_spaces_per_dim = (uint32_t)(EXP_SPACE_DIM / H);
    constexpr uint32_t n_grid_spaces_per_dim_pow2 =
                       (uint32_t)((EXP_SPACE_DIM * EXP_SPACE_DIM) / (H * H));

    grid_pos[0] = grid_idx % n_grid_spaces_per_dim;
    grid_pos[1] = (grid_idx % n_grid_spaces_per_dim_pow2) / n_grid_spaces_per_dim;
    grid_pos[2] = grid_idx / n_grid_spaces_per_dim_pow2;
}
