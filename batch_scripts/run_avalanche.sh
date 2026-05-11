#!/bin/bash -l
#SBATCH --qos=normal
#SBATCH --partition=grace-hopper
#SBATCH --job-name=avalanche_no_ps_figure_gamma
#SBATCH --time=10:00:00

# OpenMP settings:
export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export FL="out.txt"
export EXE=../avalanche_test

#run the application:
# applications may perform better with --gpu-bind=none instead of --gpu-bind=single:1

$EXE -header >| $FL

for r in 0.2 0.4 0.6 0.8; do
# for E in 16.5 22.5 28.5; do
for E in 46.5 52.5 58.5 64.5 70.5 ; do
    # compute gamma with two decimal places
    EE=$(awk "BEGIN { printf \"%5.2f\", 4.0 * $r + $E}")
    echo "Processing E = $EE r = $r"
    $EXE -i avalanche_test.input AnalyticField/E_0=$EE ParticleSeed/r_0=$r Simulation/file_path=avalanche_${EE}_${r}_ps.dat
#    srun -n 1 --cpu_bind=cores -G 1 --gpu-bind=single:1 $EXE -N $NP -E_0 $EE -t $NT -r_0 $r -dtLA 1e-5  -method 0 -outfile avalanche_${EE}_${r}_ps.dat >> $FL
  done
done
