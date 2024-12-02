#!/bin/sh

#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -J cp2k_test

module purge
source /bin/cp2k-2023.2/tools/toolchain/install/setup

mpirun -np 12 cp2k_pimd.out
