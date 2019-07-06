#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "initialize.h"

#include <stdint.h>


p_to_gr_map_t gen_particle_to_grid_map() {
    p_to_gr_map_t particle_to_grid_map;
    cudaMallocManaged(&particle_to_grid_map, N_PARTICLES * sizeof(uint32_t));

    return particle_to_grid_map;
}

gr_to_p_map_t gen_grid_to_particles_map() {
    gr_to_p_map_t grid_to_particles_map;

    uint32_t n_grid_spaces = (uint32_t)pow(EXP_SPACE_DIM / H, 3);

    cudaMallocManaged(&grid_to_particles_map,
            n_grid_spaces * sizeof(Particle*));

    for(size_t i = 0; i < n_grid_spaces; i++) {
        grid_to_particles_map[i] = NULL;
    }

    return grid_to_particles_map;
}


/*
 * The dam break initialization function creates a cubic block
 * of particles arranged into a simple cubic lattice centered in
 * the experiment space
 *
 * the experiment space exists in a right handed cartesian
 * coordinate system
 * */
void initialize_dam_break(p_to_gr_map_t particle_to_grid_map,
                          gr_to_p_map_t grid_to_particles_map) {

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
    uint32_t grid_index;

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
    for(uint32_t i = 0; i < N_PARTICLES; i++) {
        /* compute the position of the particle to be created */
        particle_pos_x =
            init_particle_pos_x - particle_spacing *
            (i / n_particles_per_dim_pow2);
        particle_pos_y =
            init_particle_pos_y + particle_spacing *
            (i % n_particles_per_dim);
        particle_pos_z =
            init_particle_pos_z - particle_spacing *
            ((i % n_particles_per_dim_pow2) / n_particles_per_dim);


        /* determine grid space into which the new particle goes */
        grid_space_layer = (uint16_t)((EXP_SPACE_DIM - particle_pos_x) / H);
        grid_space_col = (uint16_t)(particle_pos_y / H);
        grid_space_row = (uint16_t)((EXP_SPACE_DIM - particle_pos_z) / H);

        /* determine the index of this grid space in the grid space array */
        grid_index = grid_space_col +
                     grid_space_row * n_grid_spaces_per_dim +
                     grid_space_layer * n_grid_spaces_per_dim_pow2;

        /* record the grid slot of each particle */
        particle_to_grid_map[i] = grid_index;

        /* initialize the new particle */
        new_particle = new Particle;

        new_particle->particle_index = i;

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

        /*
         * insert the newly added particle at the front of the linked list of
         * particles in the grid slot
         * */
        host_insert_into_grid(grid_to_particles_map, grid_index, new_particle);
    }
}

void host_insert_into_grid(gr_to_p_map_t grid_to_particles_map, uint32_t grid_index,
                          Particle *new_particle) {
    Particle *first_particle_in_grid_slot;

    first_particle_in_grid_slot = grid_to_particles_map[grid_index];
    new_particle->next_particle = first_particle_in_grid_slot;
    grid_to_particles_map[grid_index] = new_particle;
}



