#!/bin/bash


echo "Starting sea level trough example"

# Copy and rename: 
cp -r examples/sea_level_trough/input_files ./build/input
mv ./build/input/input_files ./build/input/uptrough

# move into build directory and copy psem to build
cd build && cp input/uptrough/uptrough.slurm ./ 

# Run the mesh partitioning for 40 nodes
./bin/partmesh input/uptrough/uptrough.psem 