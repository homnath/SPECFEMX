#!/bin/bash


check_dir_exists() {
    dirname=$1 

    if [ -d "${dirname}" ]
    then
        echo "  -- ${dirname} already exists" 
    else
        mkdir ${dirname} 
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



# ____________ COMPARE THE PARTMESH LOG FILE ____________
# Before comparison, remove line stating runtime as it will vary 
# slightly for each run
awk '!/total elapsed time/' partmesh_output.txt > tmpfile && mv tmpfile partmesh_output.txt
# compare
result=$(cmp_files "../examples/sea_level_trough/stable_partmesh_output.txt" "partmesh_output.txt")
# check result
if [ $result -eq 1 ]; then 
    echo -e " WARNING: OUTPUT FILES ARE DIFFERENT"
    exit 1
else 
    echo "  -- partmesh log file is the same!"
fi 




# ___________ COMPARE THE PARTMESH DATA FILES ___________
result=$(cmp_files "../examples/sea_level_trough/stable_partmesh_output.txt" "partmesh_output.txt")

let sum=0

for FILE in ./stable_partition/*  
    do 
        # Extract the file name:
        substr="uptrough"           # Search string
        prefix=${FILE%%$substr*}    
        index=${#prefix}            # Find index in string
        prefix=${FILE:index:100};  
        

        result=$(cmp_files $FILE  "./partition/${prefix}" )

        if [ $result -eq 1 ]; then 
            echo " WARNING: OUTPUT FILES ARE DIFFERENT"
            echo -e "--> FILE:  " $prefix
            exit 1
        fi 

        # Keep a track of partmesh comparisons: 
        sum=$( expr $sum + $result)
    done

    # Double check the sum value: 
    if [ $sum -gt 0 ]; then 
        echo -e " WARNING: PARTMESH OUTPUT FILES ARE DIFFERENT"
        exit 1
    else 
        echo "  -- partmesh output files are the same! "
    fi


 