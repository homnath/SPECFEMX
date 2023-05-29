#!/bin/bash

echo ""
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Finished compiling SPECFEMX: "
echo "  Starting sea level trough example"
echo ""


# Copy and rename: 
cp -r examples/sea_level_trough/input_files ./build/input
mv ./build/input/input_files ./build/input/uptrough

echo "-  Copied input files to build/input directory & renamed"


# move into build directory and copy psem to build
cd build && cp input/uptrough/uptrough.slurm ./ 
echo "-  Copied slurm file to build directory"
echo "-  Running partmesh: "


# Run the mesh partitioning for 40 nodes
./bin/partmesh input/uptrough/uptrough.psem 


