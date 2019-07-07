#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include <stdint.h>

pi_to_gri_map_t gen_particle_idx_to_grid_idx_map() {
    pi_to_gri_map_t particle_idx_to_grid_idx_map;

    cudaMallocManaged(&particle_idx_to_grid_idx_map,
                      N_PARTICLES * sizeof(uint32_t));

    return particle_idx_to_grid_idx_map;
}

gri_to_pl_map_t gen_grid_to_particle_list_map() {
    gri_to_pl_map_t grid_to_particle_list_map;

    uint32_t n_grid_spaces = (uint32_t)pow(EXP_SPACE_DIM / H, 3);

    cudaMallocManaged(&grid_to_particle_list_map,
            n_grid_spaces * sizeof(Particle*));

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
                          pi_to_gri_map_t particle_idx_to_grid_idx_map,
                          pi_to_pa_map_t particle_idx_to_addr_map) {

    uint16_t n_grid_spaces_per_dim;
    uint32_t n_grid_spaces_per_dim_pow2;
    uint16_t n_particles_per_dim;
    uint32_t n_particles_per_dim_pow2;

    float cubic_block_rad;
    float particle_spacing;

    float space_center_x;
    float space_center_y;
    float space_center_z;

    float init_particle_pos_x;
    float init_particle_pos_y;
    float init_particle_pos_z;

    float particle_pos_x;
    float particle_pos_y;
    float particle_pos_z;

    uint16_t grid_space_layer;
    uint16_t grid_space_col;
    uint16_t grid_space_row;
    uint32_t grid_idx;

    uint32_t particle_idx;
    Particle *new_particle;



    n_grid_spaces_per_dim = (uint16_t)(EXP_SPACE_DIM / H);
    n_grid_spaces_per_dim_pow2 = (uint32_t)pow(n_grid_spaces_per_dim, 2);
    n_particles_per_dim = (uint16_t)cbrt((float)N_PARTICLES);
    n_particles_per_dim_pow2 = (uint32_t)pow(n_particles_per_dim, 2);

    cubic_block_rad = n_particles_per_dim * PARTICLE_RAD;
    particle_spacing = 2 * PARTICLE_RAD;

    space_center_x = (float)EXP_SPACE_DIM / 2;
    space_center_y = (float)EXP_SPACE_DIM / 2;
    space_center_z = (float)EXP_SPACE_DIM / 2;

    init_particle_pos_x = space_center_x + cubic_block_rad - PARTICLE_RAD;
    init_particle_pos_y = space_center_y - cubic_block_rad + PARTICLE_RAD;
    init_particle_pos_z = space_center_z + cubic_block_rad - PARTICLE_RAD;

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
        particle_pos_x =
            init_particle_pos_x - particle_spacing *
            (particle_idx / n_particles_per_dim_pow2);
        particle_pos_y =
            init_particle_pos_y + particle_spacing *
            (particle_idx % n_particles_per_dim);
        particle_pos_z =
            init_particle_pos_z - particle_spacing *
            ((particle_idx % n_particles_per_dim_pow2) / n_particles_per_dim);


        /* determine grid space into which the new particle goes */
        grid_space_layer = (uint16_t)((EXP_SPACE_DIM - particle_pos_x) / H);
        grid_space_col = (uint16_t)(particle_pos_y / H);
        grid_space_row = (uint16_t)((EXP_SPACE_DIM - particle_pos_z) / H);

        /* determine the index of this grid space in the grid space array */
        grid_idx = grid_space_col +
                   grid_space_row * n_grid_spaces_per_dim +
                   grid_space_layer * n_grid_spaces_per_dim_pow2;


        /* initialize the new particle */
        new_particle = new Particle;

        new_particle->pos_x = particle_pos_x;
        new_particle->pos_y = particle_pos_y;
        new_particle->pos_z = particle_pos_z;

        new_particle->vel_x = 0;
        new_particle->vel_y = 0;
        new_particle->vel_z = 0;

        new_particle->force_x = 0;
        new_particle->force_y = 0;
        new_particle->force_z = 0;

        new_particle->density = 0;
        new_particle->pressure = 0;
        new_particle->internal_energy = 0;

        /* record the address of the new particle */
        particle_idx_to_addr_map[particle_idx] = new_particle;

        /*
         * insert the new particle into the correct grid space and
         * record the grid space of the new particle
         * */
        host_insert_into_grid(grid_to_particle_list_map,
                              grid_idx,
                              particle_idx_to_grid_idx_map,
                              particle_idx,
                              new_particle);
    }
}


void host_insert_into_grid(gri_to_pl_map_t grid_to_particle_list_map,
                           uint32_t grid_idx,
                           pi_to_gri_map_t particle_idx_to_grid_idx_map,
                           uint32_t particle_idx,
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
    /* record grid space in which the new particle is held */
    particle_idx_to_grid_idx_map[particle_idx] = grid_idx;
}


__device__ void device_insert_into_grid(gri_to_pl_map_t grid_to_particle_list_map,
                                        uint32_t grid_idx,
                                        pi_to_gri_map_t particle_idx_to_grid_idx_map,
                                        uint32_t particle_idx,
                                        Particle *new_particle,
                                        grid_mutex_set_t mutex_set) {

    Particle *first_particle_in_grid_slot;

    lock_grid_mutex(mutex_set, grid_idx);

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

    unlock_grid_mutex(mutex_set, grid_idx);

    /* record grid space in which the new particle is held */
    particle_idx_to_grid_idx_map[particle_idx] = grid_idx;
}


__device__ void device_remove_from_grid(gri_to_pl_map_t grid_to_particle_list_map,
                                        uint32_t grid_idx,
                                        Particle *del_particle,
                                        grid_mutex_set_t mutex_set) {
    Particle *del_prev_particle;
    Particle *del_next_particle;

    lock_grid_mutex(mutex_set, grid_idx);



    /* remove the particle from the linked list */
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

    unlock_grid_mutex(mutex_set, grid_idx);
}


__device__ void lock_grid_mutex(grid_mutex_set_t mutex_set, uint32_t mutex_idx) {
    int *mutex_ptr = mutex_set + mutex_idx;

    while(atomicCAS(mutex_ptr, 0, 1) != 0) {
        /* wait until grid space is unlocked */
    }
}


__device__ void unlock_grid_mutex(grid_mutex_set_t mutex_set, uint32_t mutex_idx) {
    int *mutex_ptr = mutex_set + mutex_idx;

    atomicExch(mutex_ptr, 0);
}
