#!/bin/bash

# set thie variable to the file path of your display-graphics folder
savePath="/export/home/phys/nki/private/smoothed-particle-hydrodynamics/display-graphics"

./compute-results/compute "$1" "$savePath"
./display-graphics/display
