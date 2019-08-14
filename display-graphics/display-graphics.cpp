#include "../simulation-parameters.h"

#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <fstream>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <cmath>


/* file containing the positions of each particle and how
 * they evolve with time
 * */
std::ifstream simulation_results("simulation-results.csv");



void delayBetweenFrames()
{
    /* convert from s to ms */
    usleep(DT * 1e6);
    glutPostRedisplay();
}



/* put every particle on the canvas */
void renderFrame() {

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
    
    glMatrixMode(GL_MODELVIEW);
  
    glLoadIdentity();

    glRotatef(-10, 1, 0, 0);
    glRotatef(-20, 0, 1, 0);
    glTranslatef(0,0.2f,0);





    glBegin(GL_QUADS);
    glColor3f(0.9f, 0.9f, 0.9f);
    glVertex3f(-0.8, H - 1, 0.8);
    glVertex3f(0.8, H - 1, 0.8);
    glVertex3f(0.8, H - 1, -0.8);
    glVertex3f(-0.8, H - 1, -0.8);
    glEnd();

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
	/*
	glColor4f(0,0,0,1);
        glPointSize(2 * R_PARTICLE / EXP_SPACE_DIM * WINDOW_SIZE);
	glVertex3f(gl_pos_x, gl_pos_y, gl_pos_z); */

        glColor4f(0.0f, 0.2f, 0.9f, 0.85f);
        glVertex3f(gl_pos_x, gl_pos_y, gl_pos_z);
  
    }
    glEnd();
    glutSwapBuffers();
}



int main(int argc, char **argv) {

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(300, 300);
    glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
    glutCreateWindow("CUDA-SPH simulation");

    glutIdleFunc(delayBetweenFrames);
    glutDisplayFunc(renderFrame);
    
    glClearColor(0.8f,0.8f,0.8f,0.8f);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(2 * R_PARTICLE / EXP_SPACE_DIM * WINDOW_SIZE);
    
    glutMainLoop();
}
