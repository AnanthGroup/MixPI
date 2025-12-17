#!/bin/bash
#SBATCH -N 3
#SBATCH -J water
#SBATCH -t 2:00:00
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=2

source bin/cp2k/tools/toolchain/install/setup

srun --cpu_bind core cp2k_pimd.out

