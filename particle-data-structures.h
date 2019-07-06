#ifndef PARTICLE_DATA_STRUCTURES_H
#define PARTICLE_DATA_STRUCTURES_H

#include <stdint.h>

/*
 * classe used to allow new operator
 * to allocate space in unified memory
 * */
class Managed {
    public:
        void *operator new(size_t len) {
            void *ptr;
            cudaMallocManaged(&ptr, len);
            return ptr;
        }

        void operator delete(void *ptr) {
            cudaFree(ptr);
        }
};

/*
 * class containing the attributes of each particle
 * and their allocation and deallocation operations
 * */
class Particle : public Managed {
    public:
        uint32_t particle_index;

        float pos_x;
        float pos_y;
        float pos_z;

        float vel_x;
        float vel_y;
        float vel_z;

        float force_x;
        float force_y;
        float force_z;

        float density;
        float pressure;
        float internal_energy;

        Particle *next_particle;
};





/*****************************************************************************/
/******************************* DATA TYPES **********************************/
/*****************************************************************************/

/*
 * array whose indices correspond to each grid space within the
 * experiment space - an element of the array at index i is a linked list
 * of particles that exist within grid space i.
 * */
typedef Particle** gr_to_p_map_t;

/*
 * array whose indices correspond to each particle - an element
 * of the array corresponds to the grid space in which that
 * particle exists
 * */
typedef uint32_t* p_to_gr_map_t;





/*****************************************************************************/
/******************************** FUNCTIONS **********************************/
/*****************************************************************************/

/*
 * Included below is a set of data structure initialization functions
 * used at the begnning of each simulation. These functions are to
 * be run on the host.
 * */

p_to_gr_map_t gen_particle_to_grid_map();

gr_to_p_map_t gen_grid_to_particles_map();

void initialize_dam_break(p_to_gr_map_t particle_to_grid_map,
                          gr_to_p_map_t grid_to_particles_map);

/*
 * included below is a set of general purpose data structure manipulation
 * functions with versions for both the host and device
 * */

void host_insert_into_grid(gr_to_p_map_t grid_to_particles_map, uint32_t grid_index,
                      Particle *new_particle);



#endif
