#!/bin/bash
#module load anaconda
#conda activate SPFMX_SL_TESTS



ctr=0;
while [ $ctr = 0 ]; do
   echo "Remove slurms after? (y/n): ";
   read removeslurm
   if [ $removeslurm = 'y' ]
      then
      ctr=1
      echo "Will delete slurms."
   elif [ $removeslurm = 'n' ]
      then
      ctr=2
      echo "Wont delete slurms."
   else
      ctr=0
      echo "must be y/n (case sensitive)"
   fi
done

# Test 1

cp ./tests/slurms/test1_trough.slurm ./ 
sbatch test1_trough.slurm 

if [ $ctr = 1 ] 
  then 
    rm test1_trough.slurm
fi


