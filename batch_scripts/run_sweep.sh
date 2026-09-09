#!/bin/bash -l
#SBATCH --partition=volta-x86
#SBATCH --job-name=cons_sweep
#SBATCH --time=08:00:00
#SBATCH --output=cons_sweep-%j.out
#
# Integrator conservation sweep, run from run/conservation on the cluster:
#     sbatch run_sweep.sh          # or: bash run_sweep.sh
#
# Uses the `avalanche` binary + inputs/fig11a.input (analytic field, no MHD/PETSc).
# Particles are seeded on the magnetic axis (R_a = 3.0, Z = 0) so the poloidal
# drift is minimal and any p_phi / mu drift is the integrator's, not geometry.
# Each run writes its own conservation log to sweep/<name>.dat; big phdf dumps
# are suppressed. Sweeps fixed step hRK (rk4) and tolerance (dopri5).

set -euo pipefail

REPO="/vast/home/obeznosov/git/hymera"
EXE="${EXE:-$REPO/build_volta/bin/conservation}"   # override for CPU: EXE=$REPO/build_cpu/bin/conservation
INPUT="$REPO/run/ava/fig11a.input"
FIELDS="$REPO/run/profile/fields.h5"
NP="${NP:-1}"

# a few seconds of physical time on the magnetic axis
TLIM="${TLIM:-0.3}"          # total pushed time [s]
DT_FORCE="${DT_FORCE:-0.003}" # force/log cadence [s]  -> ~1000 log points

HERE="./"
OUTDIR="$HERE/sweep"
mkdir -p "$OUTDIR"
rm -f "$OUTDIR"/*.dat

source "$REPO"/volta_sourceme

RK4_H=(1e-5 5e-6 1e-6 5e-7 1e-7 5e-8 1e-8)   # fixed step in tau_c
DOPRI_TOL=(1e-3 1e-5 1e-7 1e-9 1e-11)         # rtol = atol

run_case() {
  local name="$1"; shift
  echo "=== $name ==="
  # shellcheck disable=SC2086
  mpirun -n "$NP" "$EXE" -i "$INPUT" \
    parthenon/job/problem_id=cons_${name} \
    parthenon/time/tlim=${TLIM} \
    parthenon/time/dt_force=${DT_FORCE} \
    ParticleSeed/Rseed=3.0 \
    ParticleSeed/Zseed=0.0 \
    Simulation/EnableLargeAngleCollisions=0 \
    Simulation/EnableSmallAngleCollisions=0 \
    Simulation/EnableComputeConservedQuantities=1 \
    Simulation/conservation_log="$OUTDIR/${name}.dat" \
    Simulation/load_fields="$FIELDS" \
    parthenon/output2/dt=-1 \
    parthenon/output8/dt=-1 \
    "$@" \
    > "$OUTDIR/${name}.log" 2>&1
  echo "    -> $OUTDIR/${name}.dat ($(wc -l < "$OUTDIR/${name}.dat" 2>/dev/null || echo 0) rows)"
}

for h in "${RK4_H[@]}"; do
  run_case "rk4_h${h}" Simulation/integrator=rk4 Simulation/hRK=${h}
done

for tol in "${DOPRI_TOL[@]}"; do
  run_case "dopri_tol${tol}" Simulation/integrator=dopri5 Simulation/rtol=${tol} Simulation/atol=${tol}
done

echo "sweep done. logs + .dat in $OUTDIR"
