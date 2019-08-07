#include <GL/gl.h>
#include <GL/glut.h>

#include "../simulation-parameters.h"

void InitializeGraphics(int *argc_ptr, char **argv) {
    glutInit(argc_ptr, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(300, 300);
    glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
    glutCreateWindow("CUDA-SPH simulation");

    glPointSize((2 * R_PARTICLE / EXP_SPACE_DIM) * WINDOW_SIZE);
}

void renderScene() {

    float gl_pos_x;
    float gl_pos_y;
    float gl_pos_z;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBegin(GL_POINTS);

    for(uint32_t i = 0; i < N_PARTICLES; i++) {
#if 0
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
#endif
    }

    glEnd();

    glutSwapBuffers();
}

void updateParticles() {

}

int main() {}
