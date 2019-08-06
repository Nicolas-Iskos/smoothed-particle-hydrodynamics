#include "particle-data-structures.h"
#include "simulation-parameters.h"

#include "calculate-field.h"
#include "integrate.h"

#include <GL/gl.h>
#include <GL/glut.h>

#include "test-functions.h"
#include <cstdio>
#include <cassert>

#include <ctime>

/* set up the OpenGL framework */
void initializeGraphics(int *argc_ptr, char **argv) {
    glutInit(argc_ptr, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(300, 300);
    glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
    glutCreateWindow("CUDA-SPH simulation");

    glPointSize((2 * R_PARTICLE / EXP_SPACE_DIM) * WINDOW_SIZE);
}

void renderScene(pi_to_pa_map_t particle_idx_to_addr_map) {

    Particle *particle;
    float gl_pos_x;
    float gl_pos_y;
    float gl_pos_z;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBegin(GL_POINTS);

    for(uint32_t i = 0; i < N_PARTICLES; i++) {
        particle = particle_idx_to_addr_map[i];

        /* the x coordinate in the graphics window is analogous
         * to the y coordinate in a standard orientation rhs system,
         * as is used within the SPH calculations
         * */
        gl_pos_x = 2 * (particle->position[1] / EXP_SPACE_DIM) - 1;

        /* the y coordinate in the graphics window is analogous
         * to the z coordinate in a standard orientation rhs system,
         * as is used within the SPH calculations
         * */
        gl_pos_y = 2 * (particle->position[2] / EXP_SPACE_DIM) - 1;

        /* the z coordinate in the graphics windowis analogous
         * to the x coordinate in a standard orientation rhs system,
         * as is used within the SPH calculations
         * */
        gl_pos_z = 2 * (particle->position[0] / EXP_SPACE_DIM) - 1;

        glVertex3f(gl_pos_x, gl_pos_y, gl_pos_z);
    }

    glEnd();

    glutSwapBuffers();
}

void updateParticles() {

}

int main(int argc, char **argv) {

    constexpr uint32_t n_grid_spaces = (EXP_SPACE_DIM * EXP_SPACE_DIM * EXP_SPACE_DIM) /
                                       (H * H * H);

    gri_to_pl_map_t grid_to_particle_list_map = gen_grid_to_particle_list_map();

    pi_to_gri_map_t last_particle_to_grid_map = gen_particle_to_grid_map();
    pi_to_gri_map_t curr_particle_to_grid_map = gen_particle_to_grid_map();

    pi_to_pa_map_t particle_idx_to_addr_map = gen_particle_idx_to_addr_map();

    initialize_dam_break(grid_to_particle_list_map, last_particle_to_grid_map,
                         curr_particle_to_grid_map, particle_idx_to_addr_map);







    std::clock_t start;
    double duration;

    start = std::clock();


    for(uint8_t i = 0; i < 50; i++) {
        #if 0
        printf("\n");
        for(int i = 0; i < N_PARTICLES; i++)
        {
            printf("position x: %f y: %f z: %f: \n",particle_idx_to_addr_map[i]->position[0],
                                                    particle_idx_to_addr_map[i]->position[1],
                                                    particle_idx_to_addr_map[i]->position[2]);
        }
        for(int i = 0; i < N_PARTICLES; i++)
        {
            printf("velocity x: %f y: %f z: %f: \n",particle_idx_to_addr_map[i]->velocity[0],
                                                    particle_idx_to_addr_map[i]->velocity[1],
                                                    particle_idx_to_addr_map[i]->velocity[2]);
        }
        #endif
        printf("\n");
        for(int i = 0; i < N_PARTICLES; i++)
        {
            uint32_t grid_pos[3];
            grid_idx_to_grid_pos(curr_particle_to_grid_map[i], grid_pos);
            printf("grid pos col: %d row: %d layer: %d \n",grid_pos[0], grid_pos[1], grid_pos[2]);
        }
        /* perform SPH calculations to update the forces acting on each
         * particle
         * */
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

        /* adjust the position and velocity of the particles based on the
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

    }

    duration = ( std::clock() - start) / (double) CLOCKS_PER_SEC;

    printf("\n");
    for(int i = 0; i < N_PARTICLES; i++)
    {
        printf("position x: %f y: %f z: %f: \n",particle_idx_to_addr_map[i]->position[0],
                                                particle_idx_to_addr_map[i]->position[1],
                                                particle_idx_to_addr_map[i]->position[2]);
    }

    for(int i = 0; i < N_PARTICLES; i++)
    {
        printf("velocity x: %f y: %f z: %f: \n",particle_idx_to_addr_map[i]->velocity[0],
                                                particle_idx_to_addr_map[i]->velocity[1],
                                                particle_idx_to_addr_map[i]->velocity[2]);
    }
    printf("\n");



    for(int i = 0; i < N_PARTICLES; i++)
    {
        uint32_t grid_pos[3];
        grid_idx_to_grid_pos(curr_particle_to_grid_map[i], grid_pos);
        printf("grid pos col: %d row: %d layer: %d \n",grid_pos[0], grid_pos[1], grid_pos[2]);
    }













    printf("%f\n", duration);

}

