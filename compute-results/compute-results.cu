/*
 * When compiled and run, this file will run the simulation and export its
 * results to a file named "simulation-results.csv."
 * */

#include "../simulation-parameters.h"

#include "particle-data-structures.h"
#include "calculate-field.h"
#include "integrate.h"
#include <fstream>
#include <sstream>
#include <iostream>

#include <cmath>



/*
 * main runs the simulation by first initializing the variables necessary
 * to run the simulation, and then running the simulation by updating
 * the parameters of each particle after successive time steps of length
 * DT.
 *
 * Inputs:
 * argc - the number of command line arguments (should be 3)
 * agrv - an array of the command line arguments: the first argument is
 *        the name of the executable, the second is the number of seconds
 *        for which you would like the simulation to run, and the third
 *        is the file path of the directory to which you would like to
 *        export the simulation results.
 *
 * Outputs: None
 * */
int main(int argc, char **argv) {

    float n_seconds_run_time;
    uint16_t n_iterations;
    std::string save_path;
    std::string full_save_path;
    Particle *particle;




    if(argc != 3) {
        std::cout << "Incorrect number of input arguments" << std::endl;
        return -1;
    }

    /* Determine how many iterations the simulation should complete */
    n_seconds_run_time = std::stof(argv[1]);
    n_iterations = (uint16_t)(n_seconds_run_time / DT);

    save_path = argv[2];
    full_save_path = save_path.append("/simulation-results.csv");
    /* initialize simulation results output file */
    std::ofstream output(full_save_path, std::ios::trunc);



    /* initialize particle data structures with dam break configuration */
    constexpr uint32_t n_grid_spaces = (EXP_SPACE_DIM * EXP_SPACE_DIM * EXP_SPACE_DIM) /
                                       (H * H * H);

    gri_to_pl_map_t grid_to_particle_list_map = gen_grid_to_particle_list_map();

    pi_to_gri_map_t last_particle_to_grid_map = gen_particle_to_grid_map();
    pi_to_gri_map_t curr_particle_to_grid_map = gen_particle_to_grid_map();

    pi_to_pa_map_t particle_idx_to_addr_map = gen_particle_idx_to_addr_map();

    initialize_dam_break(grid_to_particle_list_map, last_particle_to_grid_map,
                         curr_particle_to_grid_map, particle_idx_to_addr_map);




    /* run the simulation */
    for(uint16_t i = 0; i < n_iterations; i++) {
        /* compute the forces acting on each particle using SPH techniques */
        calculate_density<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                            PARTICLES_PER_BLOCK>>>(grid_to_particle_list_map,
                                                   curr_particle_to_grid_map,
                                                   particle_idx_to_addr_map);
        cudaDeviceSynchronize();


        calculate_pressure<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                             PARTICLES_PER_BLOCK>>>(particle_idx_to_addr_map);
        cudaDeviceSynchronize();


        calculate_net_force<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                              PARTICLES_PER_BLOCK>>>(grid_to_particle_list_map,
                                                     curr_particle_to_grid_map,
                                                     particle_idx_to_addr_map);
        cudaDeviceSynchronize();


        /* integrate the position and velocity of the particles based on the
         * recently-computed net force acting on each particle
         * */
        euler_integrate<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                          PARTICLES_PER_BLOCK>>>(particle_idx_to_addr_map);
        cudaDeviceSynchronize();


        /* ensure that no particle passes into the outer layer of grid spaces
         * in the experimental space, or out of the experimental space
         * entirely
         * */
        enforce_boundary_conditions<<<N_PARTICLES / PARTICLES_PER_BLOCK,
                                      PARTICLES_PER_BLOCK>>>(particle_idx_to_addr_map);
        cudaDeviceSynchronize();



        /* update the particle grid in 3 steps:
         *
         * 1. Keep track of what grid spaces the particles were last in, but also
         *    recompute what their new grid spaces should be and record these
         *    grid spaces as well
         *
         * 2. For each grid space, remove all particles that have left that grid
         *    space.
         *
         * 3. For each grid space, add all particles that have entered that
         *    grid space.
         * */
        update_particle_to_grid_map<<<n_grid_spaces / GRID_SPACES_PER_BLOCK,
                                      GRID_SPACES_PER_BLOCK>>>(last_particle_to_grid_map,
                                                               curr_particle_to_grid_map,
                                                               particle_idx_to_addr_map);
        cudaDeviceSynchronize();


        perform_removals_from_grid<<<n_grid_spaces / GRID_SPACES_PER_BLOCK,
                                     GRID_SPACES_PER_BLOCK>>>(grid_to_particle_list_map,
                                                              last_particle_to_grid_map,
                                                              curr_particle_to_grid_map,
                                                              particle_idx_to_addr_map);
        cudaDeviceSynchronize();


        perform_additions_to_grid<<<n_grid_spaces / GRID_SPACES_PER_BLOCK,
                                     GRID_SPACES_PER_BLOCK>>>(grid_to_particle_list_map,
                                                              last_particle_to_grid_map,
                                                              curr_particle_to_grid_map,
                                                              particle_idx_to_addr_map);
        cudaDeviceSynchronize();


        /* create a stream to hold the positions of every particle at the current
         * time step
         * */
        std::stringstream particle_positions_stream;
        /* fill the newly-created stream with the position of each particle */
        for(uint32_t particle_idx = 0; particle_idx < N_PARTICLES; particle_idx++) {
            particle = particle_idx_to_addr_map[particle_idx];
            std::stringstream position_stream;

            position_stream << round(particle->position[0] * 1e4) / 1e4;
            position_stream << " ";
            position_stream << round(particle->position[1] * 1e4) / 1e4;
            position_stream << " ";
            position_stream << round(particle->position[2] * 1e4) / 1e4;
            position_stream << "    ";

            particle_positions_stream << position_stream.str();
        }

        /* add the stream as a line of the output file */
        output << particle_positions_stream.str() << std::endl;
    }
}
