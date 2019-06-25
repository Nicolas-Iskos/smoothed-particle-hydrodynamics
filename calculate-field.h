#include <stdint.h>

/*
 * used as the identification number to retreive field quantities for
 * for a given particle
 * */
typedef uint32_t particle_id;

/*
 * used to keep track of which particles exist in each partitioned section
 * (grid slot) of the 3-D experiment space
 * */
typedef struct grid_occupant{
    particle_id id;
    struct grid_occupant *next_grid_occupant;
} grid_occupant_t;

/*
 * a 3-D array of linked lists holding the particle IDs
 * within each grid slot
 * */
typedef grid_occupant_t** grid_t;




/*
 * retreives a scalar quantity for a given particle using the scalar
 * quantity indexing scheme
 * */
__device__ float get_scalar_quantity(particle_id target, float *scalar_field);

/*
 * retreives a 3-D vector quantity for a given particle using the
 * vector quantity indexing scheme
 * */
__device__ float get_vector3_quantity_x(particle_id target,
                                        float *vector3_field);

__device__ float get_vector3_quantity_y(particle_id target,
                                        float *vector3_field);

__device__ float get_vector3_quantity_z(particle_id target,
                                        float *vector3_field);

/*
 * returns a linked list of the grid occupants that are neighbors to
 * the particle target
 * */
__device__ grid_occupant_t *get_neighbors(particle_id target,
                                          float *position_field, grid_t grid);

/*
 * uses smooth particle hydrodynamics to calculate density at the position of
 * each particle
 * */
__global__ void calculate_density(float *position_field, grid_t grid,
                                  float *density_field);

__global__ void calculate_pressure(float *density_field, float *pressure_field);

__global__ void calculate_net_force(float *density_field,
                                    float *pressure_field,
                                    float *position_field, grid_t grid,
                                    float *net_force_field);

__global__ void calculate_internal_energy(float *density_field,
                                          float *pressure_field,
                                          float *velocity_field,
                                          float *position_field, grid_t grid,
                                          float *internal_energy_field);






