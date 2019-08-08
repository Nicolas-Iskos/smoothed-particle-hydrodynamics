#include "../simulation-parameters.h"

#include <GL/gl.h>
#include <GL/glut.h>
#include <fstream>
#include <iostream>
#include <sstream>



void InitializeGraphics(int *argc_ptr, char **argv) {
    glutInit(argc_ptr, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(300, 300);
    glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
    glutCreateWindow("CUDA-SPH simulation");

    glPointSize((2 * R_PARTICLE / EXP_SPACE_DIM) * WINDOW_SIZE);
}



/* put every particle on the canvas */
void renderScene(std::ifstream &simulation_results) {

    std::string particle_positions;
    float exp_pos_x;
    float exp_pos_y;
    float exp_pos_z;
    float gl_pos_x;
    float gl_pos_y;
    float gl_pos_z;



    /* read a single line containing the positions of every particle
     *
     * As an example, for Particles p1 and p2 evolving over two
     * timesteps, the simulation_results stream is formatted as:
     *
     * p1_pos_x p1_pos_y p1_pos_z    p2_pos_x p2_pos_y p2_pos_z
     * p1_pos_x p1_pos_y p1_pos_z    p2_pos_x p2_pos_y p2_pos_z
     *
     * where each row corresponds to a different timestep
     *
     * */
    std::getline(simulation_results, particle_positions);

    /* tokenize the line into each particle and its position components */
    std::stringstream position_stream(particle_positions);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBegin(GL_POINTS);

    /* render each particle */
    while(!position_stream.eof()) {

        position_stream >> exp_pos_x;
        position_stream >> exp_pos_y;
        position_stream >> exp_pos_z;

        /* we need to perform these conversions because
         * although OpenGL uses a right handed coordinate system,
         * it is rotated in a silly orientation, with z pointing
         * out of the screen, x pointing to the right and y
         * pointing up.
         * */
        gl_pos_x = 2 * (exp_pos_y / EXP_SPACE_DIM) - 1;
        gl_pos_y = 2 * (exp_pos_z / EXP_SPACE_DIM) - 1;
        gl_pos_z = 2 * (exp_pos_x / EXP_SPACE_DIM) - 1;

        glVertex3f(gl_pos_x, gl_pos_y, gl_pos_z);
    }

    glEnd();
    glutSwapBuffers();
}




void updateParticles() {

}




int main() {}
