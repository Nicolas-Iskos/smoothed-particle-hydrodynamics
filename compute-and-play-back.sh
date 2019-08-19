#!/bin/bash

# set thie variable to the file path of your display-graphics folder
savePath="play-back-simulation"

./compute-results/compute "$1" "$savePath"
./display-graphics/display
