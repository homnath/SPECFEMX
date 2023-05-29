#!/bin/bash

echo ""
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Finished compiling SPECFEMX: "
echo "  Starting sea level trough example"
echo ""


example_path="./examples/sea_level_trough/input_files"
inp_path="input/uptrough"
build_inp_path="./build/${inp_path}"

# Copy and rename: 
cp -r $example_path $build_inp_path
echo "-  Copied input files to ${build_inp_path}" 


# move into build directory and copy psem to build
cd build && cp "${inp_path}/uptrough.slurm" ./ 
echo "-  Copied slurm file to build directory"
echo "-  Running partmesh: "



# Make the output_uptrough dir: 
mkdir output_uptrough 
mkdir tmp 


# Run the mesh partitioning for 40 nodes
./bin/partmesh "${inp_path}/uptrough.psem" > partmesh_output.txt  

# Co
cmp ../examples/sea_level_trough/stable_partmesh_output.txt partmesh_output.txt

