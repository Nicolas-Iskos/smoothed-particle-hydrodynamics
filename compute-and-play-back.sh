#!/bin/bash

# set thie variable to the file path of your display-graphics folder
savePath="play-back-simulation"

./compute-simulation/compute "$1" "$savePath"
cd "play-back-simulation"
./play-back
