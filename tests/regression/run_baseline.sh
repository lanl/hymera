#!/usr/bin/env bash
# Run the production MHD configuration and record a regression fingerprint.
#
# This is the tier that exercises the real physics path: the EFIT /
# Grad-Shafranov initial condition (ictype 9) at the production 100x2x200 grid,
# with the full matrix-free SNES + nested fieldsplit + MUMPS stack. It is slow --
# the initial-condition relaxation alone is a complete extra TSSolve -- so it is
# a phase gate, not an inner-loop check.
#
# MHD_Config/monitor=1 is required: it enables the per-step diagnostics (step
# norms, max|div B|, toroidal currents) that the fingerprint is built from.
# Those are off by default because they cost a nested linear solve per step.
#
# Field and restart dumps are suppressed; the fingerprint comes from the log.
#
#   run_baseline.sh                       # 1 step, write baseline if absent
#   run_baseline.sh --steps 3             # longer run
#   run_baseline.sh --update              # overwrite the stored baseline
#   run_baseline.sh --compare             # compare against the stored baseline
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN="$REPO/build/bin/mhd"
# The input deck uses paths relative to build/bin, so run from there.
RUNDIR="$REPO/build/bin"
BASELINE_DIR="$REPO/tests/regression/baselines"

STEPS=1
MODE=check
DT=0.00108   # matches inputs/mhd.input

while [ $# -gt 0 ]; do
  case "$1" in
    --steps)   STEPS="$2"; shift 2 ;;
    --update)  MODE=update; shift ;;
    --compare) MODE=compare; shift ;;
    -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

[ -x "$BIN" ] || { echo "not built: $BIN" >&2; exit 1; }

# shellcheck source=/dev/null
. "$REPO/tests/regression/setup_env.sh" >/dev/null
PY="${HARNESS_PYTHON:-$HOME/.venv-hymera/bin/python}"

TLIM=$($PY -c "print(f'{$DT * $STEPS:.10g}')")
NAME="prod_${STEPS}step_100x02x200"
LOG="$(mktemp -t "${NAME}.XXXXXX.log")"

echo "running $STEPS step(s) at 100x2x200, tlim=$TLIM"
echo "log: $LOG"

cd "$RUNDIR"
set +e
./mhd -i ../../inputs/mhd.input \
      MHD_Config/ic_binary_load=0 \
      MHD_Config/monitor=1 \
      parthenon/time/dt_force="$DT" \
      parthenon/time/tlim="$TLIM" \
      parthenon/output2/dt=-1 \
      parthenon/output8/dt=-1 \
      > "$LOG" 2>&1
RC=$?
set -e
cd - >/dev/null

# PETSc errors do not always reach the exit code, so check the log too.
if [ $RC -ne 0 ]; then
  echo "solver exited $RC; tail of log:" >&2
  tail -20 "$LOG" >&2
  exit $RC
fi
if grep -q 'DIVERGED' "$LOG"; then
  echo "FAIL: solver reported a divergence" >&2
  grep -m5 'DIVERGED' "$LOG" >&2
  exit 1
fi

FP="$(mktemp -t "${NAME}.XXXXXX.json")"
"$PY" "$REPO/tests/regression/fingerprint.py" "$LOG" -o "$FP"

STORED="$BASELINE_DIR/${NAME}.json"
case "$MODE" in
  update)
    mkdir -p "$BASELINE_DIR"
    cp "$FP" "$STORED"
    echo "baseline updated: $STORED"
    ;;
  compare)
    [ -f "$STORED" ] || { echo "no stored baseline at $STORED" >&2; exit 1; }
    "$PY" "$REPO/tests/regression/fingerprint.py" --compare "$STORED" "$FP"
    ;;
  check)
    if [ -f "$STORED" ]; then
      "$PY" "$REPO/tests/regression/fingerprint.py" --compare "$STORED" "$FP"
    else
      mkdir -p "$BASELINE_DIR"
      cp "$FP" "$STORED"
      echo "no baseline existed; recorded one: $STORED"
    fi
    ;;
esac
