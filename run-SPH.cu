#include "particle-data-structures.h"
#include "calculate-field.h"
#include "integrate.h"
#include "run-SPH.h"

#define NUM_PARTICLES    4096

void initialize_dam_break(p_to_gr_map_t particle_to_grid_map,
                          gr_to_p_map_t grid_to_particles_map) {

}

int main() {

    Particle *a = new Particle;
    a->next_particle = new Particle;


}
