#!/bin/bash -l
#SBATCH --qos=normal
#SBATCH --nodes=4
#SBATCH --partition=volta-x86
#SBATCH --job-name=hybrid
#SBATCH --time=10:00:00

source ../../volta_sourceme

export OMP_NUM_THREADS=1

mpirun -n 4 ../../bin/hybrid \
  -i hybrid.input \
  parthenon/job/problem_id=hybrid_p \
  > hybrid_p.out

# mpirun -n 4 ../../bin/hybrid \
#   -r hybrid.out8.00004.rhdf \
#   >> hybrid.out
