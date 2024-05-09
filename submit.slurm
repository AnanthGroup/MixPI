#!/bin/sh

#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -J cp2k_test

module purge
cd $SLURM_SUBMIT_DIR
source /people/john708/cp2k_modules9.sh
echo "Time: `date`"
echo "Using CP2K from location: $CP2K"

##CP2K run command goes here
export IOFILE="argon-in-water"
mpirun -np 1 cp2k_pimd.out
