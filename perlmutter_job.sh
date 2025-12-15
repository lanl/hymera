#!/bin/sh

#SBATCH -t 12:00:00
#SBATCH -C gpu
#SBATCH -N 8
#SBATCH --gpus-per-node 4
#SBATCH -Am3016
#SBATCH -q regular
#SBATCH --mail-type=BEGIN,END,FAIL

# OpenMP settings:
export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun --tasks-per-node=4 --cpu-bind=cores -G 32 --gpu-bind=single:1 ./hybrid -i ../../inputs/avalanche_test.input MHD_Config/EnableRelaxation=1 MHD_Config/EnableReadICFromBinary=0 > hybrid.out
