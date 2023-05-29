#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=40
#SBATCH -t 00:01:40

# load appropriate compiler
source ~/.login

cd $SLURM_SUBMIT_DIR
# problems with mpiexec/mpirun so use srun
mpiexec -n 1 ./bin/pspecfemx ./tests/input_files/test1/test1_trough.psem  -ksp_monitor -ksp_type fgmres -pc_type asm -ksp_gmres_restart 300
