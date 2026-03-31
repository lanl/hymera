#!/bin/bash -l
#SBATCH --qos=normal
#SBATCH --partition=volta-x86
#SBATCH --job-name=hybrid_array
#SBATCH --time=10:00:00
#SBATCH --array=0-6
#SBATCH --output=slurm-%A_%a.out   # for arrays

VALUES=(32 64 128 256 512 1024 2048)
NN=${VALUES[$SLURM_ARRAY_TASK_ID]}

source ../volta_sourceme

mpirun -n 16 ../bin/hybrid \
  -i hybrid.input \
  ParticleSeed/num_particles_per_block=${NN} \
  parthenon/job/problem_id=hybrid${NN} \
  > hybrid${NN}
