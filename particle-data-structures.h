#ifndef PARTICLE_DATA_STRUCTURES_H
#define PARTICLE_DATA_STRUCTURES_H

#include <stdint.h>

typedef uint32_t particle_id;

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

class Particle : public Managed {
    public:
        particle_id id;

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

typedef Particle** gr_to_p_map_t;

typedef int* p_to_gr_map_t;

#endif
