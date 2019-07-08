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

        Particle *prev_particle;
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
typedef Particle **gri_to_pl_map_t;

/*
 * array whose indices correspond to each particle - an element
 * of the array corresponds to the grid space in which that
 * particle exists
 * */
typedef uint32_t *pi_to_gri_map_t;



typedef Particle **pi_to_pa_map_t;



typedef int *grid_mutex_set_t;


/*****************************************************************************/
/******************************** FUNCTIONS **********************************/
/*****************************************************************************/

/*
 * Included below is a set of data structure initialization functions
 * used at the begnning of each simulation. These functions are to
 * be run on the host.
 * */


gri_to_pl_map_t gen_grid_to_particle_list_map();

pi_to_gri_map_t gen_particle_idx_to_grid_idx_map();

pi_to_pa_map_t gen_particle_idx_to_addr_map();

grid_mutex_set_t gen_grid_mutex_set();

void initialize_dam_break(gri_to_pl_map_t grid_to_particle_list_map,
                          pi_to_gri_map_t particle_idx_to_grid_idx_map,
                          pi_to_pa_map_t particle_idx_to_addr_map);

/*
 * included below is a set of general purpose data structure manipulation
 * functions with versions for both the host and device
 * */

void host_insert_into_grid(gri_to_pl_map_t grid_to_particle_list_map,
                           uint32_t grid_idx,
                           pi_to_gri_map_t particle_idx_to_grid_idx_map,
                           uint32_t particle_idx,
                           Particle *new_particle);

__device__ void device_insert_into_grid(gri_to_pl_map_t grid_to_particle_list_map,
                                        uint32_t grid_idx,
                                        pi_to_gri_map_t particle_idx_to_grid_idx_map,
                                        uint32_t particle_idx,
                                        Particle *new_particle,
                                        grid_mutex_set_t mutex_set);

__device__ void device_remove_from_grid(gri_to_pl_map_t grid_to_particle_list_map,
                                        uint32_t grid_idx,
                                        Particle *del_particle,
                                        grid_mutex_set_t mutext_set);

__device__ void unlock_grid_mutex(grid_mutex_set_t mutex_set, uint32_t mutex_idx);

__device__ void lock_grid_mutex(grid_mutex_set_t mutex_set, uint32_t mutex_idx);



#endif
