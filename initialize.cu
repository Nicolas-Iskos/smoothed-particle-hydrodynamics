#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "initialize.h"

#include <stdint.h>

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

    uint16_t n_grid_spaces_per_dim = EXP_SPACE_DIM / H;

    float cubic_block_rad = CBRT_N_PARTICLES * PARTICLE_RAD;

    float particle_spacing = 2 * PARTICLE_RAD;

    float space_center_x = EXP_SPACE_DIM / 2;
    float space_center_y = EXP_SPACE_DIM / 2;
    float space_center_z = EXP_SPACE_DIM / 2;

    float init_particle_pos_x = space_center_x + cubic_block_rad;
    float init_particle_pos_y = space_center_y - cubic_block_rad;
    float init_particle_pos_z = space_center_z + cubic_block_rad;

    float particle_pos_x;
    float particle_pos_y;
    float particle_pos_z;

    uint16_t grid_slot_x;
    uint16_t grid_slot_y;
    uint16_t grid_slot_z;

    uint32_t grid_index;

    Particle *first_particle_in_grid_slot;
    Particle *new_particle;

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
    for(size_t i = 0; i < N_PARTICLES; i++) {
        particle_pos_x =
            init_particle_pos_x + particle_spacing *
            i / (CBRT_N_PARTICLES * CBRT_N_PARTICLES);
        particle_pos_y =
            init_particle_pos_y + particle_spacing * (i % CBRT_N_PARTICLES);
        particle_pos_z =
            init_particle_pos_z - particle_spacing *
            ((i % CBRT_N_PARTICLES * CBRT_N_PARTICLES) / CBRT_N_PARTICLES);

        grid_slot_x = n_grid_spaces_per_dim - particle_pos_x / H;
        grid_slot_y = particle_pos_y / H;
        grid_slot_z = n_grid_spaces_per_dim - particle_pos_z / H;

        grid_index = grid_slot_y +
                     grid_slot_z * CBRT_N_PARTICLES +
                     grid_slot_x * CBRT_N_PARTICLES * CBRT_N_PARTICLES;


        /* record the grid slot of each particle */
        particle_to_grid_map[i] = grid_index;

        /*
         * add the particle to the list of particles in its corresponding
         * grid slot
         * */
        first_particle_in_grid_slot = grid_to_particles_map[grid_index];

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

        /*
         * insert the newly added particle at the front of the list of
         * particles in the grid slot
         * */
        new_particle->next_particle = first_particle_in_grid_slot;
        grid_to_particles_map[grid_index] = new_particle;
    }

}
