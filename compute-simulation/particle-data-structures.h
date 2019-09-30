/*
 * This file contains the interface and some of the implementation
 * for the set of functions and data structures used to manage
 * the high-level data structures that allow the simulation to run.
 * Complete descriptions of the functions defined in this file are
 * contained within particle-data-structures.cu.
 * */



#ifndef PARTICLE_DATA_STRUCTURES_H
#define PARTICLE_DATA_STRUCTURES_H



#include <cstdint>



/*
 * classe used to allow new operator to
 * allocate space for each particle in unified memory
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
        /* in each of these vectors,
         * element 0 is the x component,
         * element 1 is the y component,
         * element 2 is the z component
         * */
        float position[3];
        float velocity[3];
        float force[3];

        float density;
        float pressure;

        Particle *prev_particle;
        Particle *next_particle;
};



/*****************************************************************************/
/***************************** DATA STRUCTURES *******************************/
/*****************************************************************************/

/*
 * array whose indices correspond to each grid space within the
 * experiment space - an element of the array at index i is a linked list
 * of particles that exist within grid space i.
 * */
typedef Particle **gri_to_pl_map_t;

/*
 * array whose indices correspond to each particle index.
 * An element of the array corresponds to the index of the
 * grid space in which that particle exists
 * */
typedef uint32_t *pi_to_gri_map_t;

/*
 * array whose indices correspond to each particle index.
 * An element of the array corresponds to the memory address of
 * that particle
 * */
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

void free_grid_to_particle_list_map(gri_to_pl_map_t grid_to_particle_list);

void free_particle_to_grid_map(pi_to_gri_map_t particle_to_grid_map);

void free_particle_idx_to_addr_map(pi_to_pa_map_t particle_idx_to_addr_map);

void initialize_dam_break(gri_to_pl_map_t grid_to_particle_list_map,
                          pi_to_gri_map_t last_particle_to_grid_map,
                          pi_to_gri_map_t curr_particle_to_grid_map,
                          pi_to_pa_map_t particle_idx_to_addr_map);

void initialize_two_block_collision(gri_to_pl_map_t grid_to_particle_list_map,
                                    pi_to_gri_map_t last_particle_to_grid_map,
                                    pi_to_gri_map_t curr_particle_to_grid_map,
                                    pi_to_pa_map_t particle_idx_to_addr_map);

void initialize_block_wall_collision(gri_to_pl_map_t grid_to_particle_list_map,
                                    pi_to_gri_map_t last_particle_to_grid_map,
                                    pi_to_gri_map_t curr_particle_to_grid_map,
                                    pi_to_pa_map_t particle_idx_to_addr_map);
/*
 * included below is a set of general purpose data structure manipulation
 * functions with versions for both the host and device
 * */

void host_insert_into_grid(gri_to_pl_map_t grid_to_particle_list_map,
                           uint32_t grid_idx,
                           Particle *new_particle);

__global__ void update_particle_to_grid_map(
                           pi_to_gri_map_t last_particle_to_grid_map,
                           pi_to_gri_map_t curr_particle_to_grid_map,
                           pi_to_pa_map_t particle_idx_to_addr_map);


__global__ void perform_removals_from_grid(
                           gri_to_pl_map_t grid_to_particle_list_map,
                           pi_to_gri_map_t last_particle_to_grid_map,
                           pi_to_gri_map_t curr_particle_to_grid_map,
                           pi_to_pa_map_t particle_idx_to_addr_map);

__global__ void perform_additions_to_grid(
                           gri_to_pl_map_t grid_to_particle_list_map,
                           pi_to_gri_map_t last_particle_to_grid_map,
                           pi_to_gri_map_t curr_particle_to_grid_map,
                           pi_to_pa_map_t particle_idx_to_addr_map);

/*
 * included below is a set of accessor and utility functions to be used
 * on the data structures described above
 * */

__host__ __device__ uint32_t particle_pos_to_grid_idx(float *position);

__host__ __device__ uint32_t grid_pos_to_grid_idx(uint32_t *grid_pos);

__host__ __device__ void grid_idx_to_grid_pos(uint32_t grid_idx,
                                              uint32_t *grid_pos);




#endif
