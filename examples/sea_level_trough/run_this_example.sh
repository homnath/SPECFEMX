#!/bin/bash


check_dir_exists() {
    dirname=$1 

    if [ -d "${dirname}" ]
    then
        echo "  -- ${dirname} already exists" 
    else
        mkdir output_uptrough 
        echo "  -- created ${dirname}" 
    fi 

}


function cmp_files(){
    # Compare with stable output from partmesh - should return blank if identical
    if [ -n "$(cmp $1 $2)" ]
    then 
        # Error - not the same
        echo 1
    else 
        echo 0
    fi
    
}




echo ""
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo " - Finished compiling SPECFEMX: "
echo " - Starting sea level trough example"
echo ""


example_path="./examples/sea_level_trough/input_files"
inp_path="input/uptrough"
build_inp_path="./build/${inp_path}"

# Copy and rename: 
cp -r $example_path $build_inp_path
echo "  --  Copied input files to ${build_inp_path}" 


# move into build directory and copy psem to build
cd build && cp "${inp_path}/uptrough.slurm" ./ 
echo "  --  Copied slurm file to build directory"



# Make the output_uptrough dir: 
echo ""
echo " - Checking if directories need to be created"
check_dir_exists "output_uptrough"
check_dir_exists "tmp"
echo ""
 

# Run the mesh partitioning for 40 nodes
echo " - Running partmesh: "
./bin/partmesh "${inp_path}/uptrough.psem" > partmesh_output.txt  

# Before comparison, remove line stating runtime as it will vary 
# slightly for each run
awk '!/total elapsed time/' partmesh_output.txt > tmpfile && mv tmpfile partmesh_output.txt



# Compare with stable output from partmesh - should return 0 if identical
result=$(cmp_files "../examples/sea_level_trough/stable_partmesh_output.txt" "partmesh_output.txt")

if [ $result -eq 1 ]; then 
    echo -e " WARNING: OUTPUT FILES ARE DIFFERENT"
    exit 1
else 
    echo "  -- partmesh output files are the same (yay!)"
fi 