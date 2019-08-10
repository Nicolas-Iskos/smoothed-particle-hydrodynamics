#!/bin/bash

# set thie variable to the file path of your display-graphics folder
displayGraphicsPath="/export/home/phys/nki/private/smoothed-particle-hydrodynamics"

./compute-results/compute "$1" "$displayGraphicsPath"
./display-graphics/display
