#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include <cstdint>


pi_to_gri_map_t gen_particle_to_grid_map() {
    pi_to_gri_map_t particle_to_grid_map;

    cudaMallocManaged(&particle_to_grid_map,
                      N_PARTICLES * sizeof(uint32_t));

    return particle_to_grid_map;
}


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


pi_to_pa_map_t gen_particle_idx_to_addr_map() {
    pi_to_pa_map_t particle_idx_to_addr_map;

    cudaMallocManaged(&particle_idx_to_addr_map,
            N_PARTICLES * sizeof(Particle*));

    return particle_idx_to_addr_map;
}

/*
 * The dam break initialization function creates a cubic block
 * of particles arranged into a simple cubic lattice centered in
 * the experiment space
 *
 * the experiment space exists in a right handed cartesian
 * coordinate system
 * */
void initialize_dam_break(gri_to_pl_map_t grid_to_particle_list_map,
                          pi_to_gri_map_t last_particle_to_grid_map,
                          pi_to_gri_map_t curr_particle_to_grid_map,
                          pi_to_pa_map_t particle_idx_to_addr_map) {

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



    n_particles_per_dim = (uint32_t)cbrt((float)N_PARTICLES);
    n_particles_per_dim_pow2 = (uint32_t)pow(n_particles_per_dim, 2);

    cubic_block_rad = n_particles_per_dim * R_PARTICLE;
    particle_spacing = 2 * R_PARTICLE;

    space_center[0] = (float)EXP_SPACE_DIM / 2;
    space_center[1] = (float)EXP_SPACE_DIM / 2;
    space_center[2] = (float)EXP_SPACE_DIM / 2;

    init_particle_pos[0] = space_center[0] + cubic_block_rad - R_PARTICLE;
    init_particle_pos[1] = space_center[1] - cubic_block_rad + R_PARTICLE;
    init_particle_pos[2] = space_center[2] + cubic_block_rad - R_PARTICLE;

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

        new_particle->velocity[0] = 0;
        new_particle->velocity[1] = 0;
        new_particle->velocity[2] = 0;

        new_particle->force[0] = 0;
        new_particle->force[1] = 0;
        new_particle->force[2] = 0;

        new_particle->density = 0;
        new_particle->pressure = 0;
        new_particle->internal_energy = 0;

        /* record the address of the new particle */
        particle_idx_to_addr_map[particle_idx] = new_particle;

        /* record the grid index of each particle */
        grid_idx = particle_pos_to_grid_idx(new_particle->position);
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
                                pi_to_pa_map_t particle_idx_to_addr_map,
                                pi_to_gri_map_t last_particle_to_grid_map,
                                pi_to_gri_map_t curr_particle_to_grid_map) {
    uint32_t particle_idx;
    uint32_t pre_update_grid_idx;
    uint32_t updated_grid_idx;
    Particle *particle;

    particle_idx = blockDim.x * blockIdx.x + threadIdx.x;
    particle = particle_idx_to_addr_map[particle_idx];
    pre_update_grid_idx = curr_particle_to_grid_map[particle_idx];
    updated_grid_idx = particle_pos_to_grid_idx(particle->position);

    /* set the pre-updated grid_idx in the last particle to grid map */
    last_particle_to_grid_map[particle_idx] = pre_update_grid_idx;

    /* set the updated grid idx into the current particle to grid map */
    curr_particle_to_grid_map[particle_idx] = updated_grid_idx;
}


__global__ void perform_removals_from_grid(
                                gri_to_pl_map_t grid_to_particle_list_map,
                                pi_to_gri_map_t last_particle_to_grid_map,
                                pi_to_gri_map_t curr_particle_to_grid_map,
                                pi_to_pa_map_t particle_idx_to_addr_map) {
    uint32_t grid_idx;
    Particle *del_particle;
    Particle *del_prev_particle;
    Particle *del_next_particle;

    grid_idx = blockDim.x * blockIdx.x + threadIdx.x;

    for(uint32_t particle_idx = 0; particle_idx < N_PARTICLES; particle_idx++){
        if((last_particle_to_grid_map[particle_idx] == grid_idx) &&
           (curr_particle_to_grid_map[particle_idx] != grid_idx)) {

            del_particle = particle_idx_to_addr_map[particle_idx];
            del_prev_particle = del_particle->prev_particle;
            del_next_particle = del_particle->next_particle;

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


__global__ void perform_additions_to_grid(
                                gri_to_pl_map_t grid_to_particle_list_map,
                                pi_to_gri_map_t last_particle_to_grid_map,
                                pi_to_gri_map_t curr_particle_to_grid_map,
                                pi_to_pa_map_t particle_idx_to_addr_map) {
    uint32_t grid_idx;
    Particle *particle;
    Particle *first_particle_in_grid_slot;

    grid_idx = blockDim.x * blockIdx.x + threadIdx.x;

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


__host__ __device__ uint32_t particle_pos_to_grid_idx(float position[3]) {
    uint32_t grid_pos[3];

    /* grid space column is related to y coordinate */
    grid_pos[0] = (uint16_t)(position[1] / H);

    /* grid space row is related to z coordinate */
    grid_pos[1] = (uint16_t)((EXP_SPACE_DIM - position[2]) / H);

    /* grid space layer is related to x coordinate */
    grid_pos[2] = (uint16_t)((EXP_SPACE_DIM - position[0]) / H);

    return grid_pos_to_grid_idx(grid_pos);
}


__host__ __device__ uint32_t grid_pos_to_grid_idx(uint32_t grid_pos[3]) {

    constexpr uint32_t n_grid_spaces_per_dim = (uint32_t)(EXP_SPACE_DIM / H);
    constexpr uint32_t n_grid_spaces_per_dim_pow2 =
                       (uint32_t)((EXP_SPACE_DIM * EXP_SPACE_DIM) / (H * H));

    /* column + row * num_cols + layer * num_rows * num_cols */
    return grid_pos[0] +
           grid_pos[1] * n_grid_spaces_per_dim +
           grid_pos[2] * n_grid_spaces_per_dim_pow2;
}


__host__ __device__ void grid_idx_to_grid_pos(uint32_t grid_idx,
                                              uint32_t (&grid_coordinates)[3]) {

    constexpr uint32_t n_grid_spaces_per_dim = (uint32_t)(EXP_SPACE_DIM / H);
    constexpr uint32_t n_grid_spaces_per_dim_pow2 =
                       (uint32_t)((EXP_SPACE_DIM * EXP_SPACE_DIM) / (H * H));

    grid_coordinates[0] = grid_idx % n_grid_spaces_per_dim;
    grid_coordinates[1] = (grid_idx % n_grid_spaces_per_dim_pow2) / n_grid_spaces_per_dim;
    grid_coordinates[2] = grid_idx / n_grid_spaces_per_dim_pow2;
}
