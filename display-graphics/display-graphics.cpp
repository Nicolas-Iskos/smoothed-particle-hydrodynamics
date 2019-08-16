/* we must include this line to prevent warnings about
 * deprecated functions
 * */
#define GL_SILENCE_DEPRECATION

#include "../simulation-parameters.h"

#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
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
    usleep(1.2 * DT * 1e6);
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

    float floor_color[4] = {0.5, 0.5, 0.5, 0.9};
    float particle_color[4] = {0.0, 0.15, 0.82, 1};
    float white[4] = {1.0, 1.0, 1.0, 1.0};

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
    glRotatef(-20, 1, 0, 0);
    glRotatef(-20, 0, 1, 0);
    glTranslatef(0,0.2f,0);



    glBegin(GL_QUADS);

    glMaterialfv(GL_FRONT, GL_AMBIENT, floor_color);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, floor_color);
    glMaterialfv(GL_FRONT, GL_SPECULAR, white);
    glMaterialf(GL_FRONT, GL_SHININESS, 80.0);
    glNormal3f(0.0, 1.0, 0.0);
    glVertex3f(-0.8, H - 1, 0.8);
    glVertex3f(0.8, H - 1, 0.8);
    glVertex3f(0.8, H - 1, -0.8);
    glVertex3f(-0.8, H - 1, -0.8);

    glEnd();



    glBegin(GL_POINTS);

    glMaterialfv(GL_FRONT, GL_AMBIENT, particle_color);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, particle_color);
    glMaterialfv(GL_FRONT, GL_SPECULAR, white);
    glMaterialf(GL_FRONT, GL_SHININESS, 50.0);
    glNormal3f(0.0, 0.0, 1.0);
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



int main(int argc, char **argv) {

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(300, 300);
    glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
    glutCreateWindow("CUDA-SPH simulation");

    glutIdleFunc(delayBetweenFrames);
    glutDisplayFunc(renderFrame);

    glClearColor(0.7 ,0.7 ,0.7 ,0.7);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    float qaAmbientLight[4] = {0.1, 0.1, 0.1, 1.0};
    float qaDiffuseLight[4] = {0.9, 0.9, 0.9, 1.0};
    float qaSpecularLight[4] = {1.0, 1.0, 1.0, 1.0};
    glLightfv(GL_LIGHT0, GL_AMBIENT, qaAmbientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, qaDiffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, qaSpecularLight);

    float qaLightPosition[4] = {0.75, 0.75, 0.5, 1.0};
    glLightfv(GL_LIGHT0, GL_POSITION, qaLightPosition);

    glEnable(GL_POINT_SMOOTH);
    glPointSize(2 * R_PARTICLE * OVERLAP_FACTOR /
                EXP_SPACE_DIM * WINDOW_SIZE);


    sleep(1);

    glutMainLoop();
}
