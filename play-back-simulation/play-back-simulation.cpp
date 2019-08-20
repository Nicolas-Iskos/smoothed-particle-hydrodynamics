/*
 * When compiled and run, this file will graphically display the
 * results of the simulation contained within the file
 * simulation-results.csv using OpenGL. The heart of this program
 * is the function render_frame, which upon each call will read
 * a line containing the position of each particle from
 * simulation-results.csv. Each line i in simulation-results.csv
 * represents the positions of each particle after i time steps
 * of length DT.
 * */



#include "../simulation-parameters.h"

/* we must include this line to prevent warnings about
 * deprecated functions
 * */

#define GL_SILENCE_DEPRECATION
#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif
#include <fstream>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <cmath>



/* global input file stream containing the positions of each
 * particle and how they evolve with time
 * */
std::ifstream simulation_results("simulation-results.csv");



/*
 * delayBetweenFrames places a time delay between each frame so
 * that the simulation runs at a reasonable rate.
 *
 * Inputs: None
 * Outputs: None
 * */
void delay_between_frames()
{
    /* we use the slowdown factor to allow the user to slow down the
     * simulation to make more careful observation
     * */
    usleep(SLOWDOWN_FACTOR * DT * 1e6);
    glutPostRedisplay();
}



/*
 * renderFrame renders the objects visible in the experiment, including
 * the floor of the experiment space and the particles within the simulation.
 *
 * Inputs: None
 * Oututs: None
 * */
void render_frame() {

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



    /* clear the color and depth buffers and then rotate the experiment
     * so that it is visible in an isometric-like view that allows the
     * user to observe all three dimensions
     * */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glRotatef(-20, 1, 0, 0);
    glRotatef(-20, 0, 1, 0);
    glTranslatef(0,0.2f,0);



    /* Construct the floor of the simulation space - this floor exists
     * to give the user an enhanced sense of depth
     * */
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



    /* Render each particle by processing the most recent line read
     * out of the simulation-results file stream
     * */
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



    /* replace the current image buffer with the newly constructed one
     * built above
     * */
    glutSwapBuffers();
}



/*
 * handle_key_press processes user input of the keys 'r' and 'e'.
 * 'r' causes the simulation to restart, while 'e' closes the
 * simulation window by exiting the process.
 *
 * Inputs:
 * key - the name of the key which the user has pressed
 * xmouse - the x coordinate of the location the user has clicked
 * ymouse - the y coordinate of the location the user has clicked
 * */
void handle_key_press(unsigned char key, int xmouse, int ymouse) {

    switch (key) {
        /* pressing the 'e' key will exit the simulation by terminating the
         * process in which it is running
         * */
        case 'e':
            exit(0);
            break;

        /* pressing the 'r' key will reset the simulation, effectively
         * replaying what the user has just seen
         * */
        case 'r':
            simulation_results.clear();
            simulation_results.seekg(0, std::ios::beg);
            break;
    }
}



/*
 * main initializes the window and sets up most of the constant parameters of
 * the simulation, including lighting and particle sizing
 *
 * Inputs: None
 * Outputs: None
 *
 * */
int main(int argc, char **argv) {

    /* the proceeding function calls are made to initialize the simulation */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(300, 300);
    glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
    glutCreateWindow("CUDA-SPH simulation");

    /* register function callbacks used for frame rate pacing, frame rendering
     * and response to user key presses
     * */
    glutIdleFunc(delay_between_frames);
    glutDisplayFunc(render_frame);
    glutKeyboardFunc(handle_key_press);

    /* initialize the graphical properties of the experiment space */
    glClearColor(0.7 ,0.7 ,0.7 ,0.7);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

    /* initialize lighting */
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

    /* initialize the graphical properties of the particles rendered in
     * the simulation
     * */
    glEnable(GL_POINT_SMOOTH);
    glPointSize(2 * R_PARTICLE * OVERLAP_FACTOR /
                EXP_SPACE_DIM * WINDOW_SIZE);

    /* begin the graphical display of the simulation */
    glutMainLoop();
}
