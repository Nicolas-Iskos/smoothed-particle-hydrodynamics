#ifndef PARTICLE_DATA_STRUCTURES_H
#define PARTICLE_DATA_STRUCTURES_H

#include <stdint.h>
#include <array>
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
        float position[3];
        float velocity[3];
        float force[3];

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



/*****************************************************************************/
/******************************** FUNCTIONS **********************************/
/*****************************************************************************/

/*
 * Included below is a set of data structure initialization functions
 * used at the begnning of each simulation. These functions are to
 * be run on the host.
 * */


gri_to_pl_map_t gen_grid_to_particle_list_map();

pi_to_gri_map_t gen_particle_to_grid_map();

pi_to_pa_map_t gen_particle_idx_to_addr_map();

void initialize_dam_break(gri_to_pl_map_t grid_to_particle_list_map,
                          pi_to_gri_map_t last_particle_to_grid_map,
                          pi_to_gri_map_t curr_particle_to_grid_map,
                          pi_to_pa_map_t particle_idx_to_addr_map);






/*
 * included below is a set of general purpose data structure manipulation
 * functions with versions for both the host and device
 * */




/* this function is used for the host executed particle data structure
 * initialization functions
 * */
void host_insert_into_grid(gri_to_pl_map_t grid_to_particle_list_map,
                           uint32_t grid_idx,
                           Particle *new_particle);

/* this function updates the mapping between each particle and its grid space.
 * It also keeps track of the old grid space of each particle, so that the
 * grid-parallel update_grid_to_particle_list_map function can easily remove
 * and add particles to the corresponding grid slots.
 * */
__global__ void update_particle_to_grid_map(
                           pi_to_gri_map_t curr_particle_to_grid_map,
                           pi_to_gri_map_t last_particle_to_grid_map,
                           pi_to_pa_map_t particle_idx_to_addr_map);



__global__ void remove_relevant_particles_from_grid(
                           gri_to_pl_map_t grid_to_particle_list_map,
                           pi_to_gri_map_t last_particle_to_grid_map,
                           pi_to_gri_map_t curr_particle_to_grid_map,
                           pi_to_pa_map_t particle_idx_to_addr_map);

__global__ void add_relevant_particles_to_grid(
                           gri_to_pl_map_t grid_to_particle_list_map,
                           pi_to_gri_map_t last_particle_to_grid_map,
                           pi_to_gri_map_t curr_particle_to_grid_map,
                           pi_to_pa_map_t particle_idx_to_addr_map);

__host__ __device__ uint32_t calculate_grid_idx(float position[]);

/*
__device__ void device_insert_into_grid(gri_to_pl_map_t grid_to_particle_list_map,
                                        uint32_t grid_idx,
                                        pi_to_gri_map_t particle_idx_to_grid_idx_map,
                                        uint32_t particle_idx,
                                        Particle *new_particle);

__device__ void device_remove_from_grid(gri_to_pl_map_t grid_to_particle_list_map,
                                        uint32_t grid_idx,
                                        Particle *del_particle);
*/



#endif
