#!/usr/bin/env bash
# Tier-0 regression check for the mass-matrix coefficients.
#
# Runs t0_coefficients and compares the SHA-256 of its output with the recorded
# value in baselines/t0_coefficients.sha256. The output itself is ~390 MB of exact
# hexadecimal doubles, too large to commit, so only its hash is stored.
#
#   run_t0.sh <t0_coefficients-binary>            check against the baseline
#   run_t0.sh <t0_coefficients-binary> --update   record a new baseline
#   run_t0.sh <t0_coefficients-binary> --keep F   also keep the full output in F,
#                                                 to diff against another revision
#   run_t0.sh <binary> --baseline NAME            use baselines/NAME.sha256
#                                                 (default t0_coefficients)
#
# A hash mismatch says only that something changed. To see what, run --keep on
# both revisions and diff the two files: each line names the function, the cell
# indices and the stencil location, so a diff localizes the change exactly.
set -euo pipefail

BIN="${1:?usage: run_t0.sh <t0_coefficients> [--update|--keep FILE]}"
shift
MODE=check
KEEP=""
NAME=t0_coefficients
while [ $# -gt 0 ]; do
  case "$1" in
    --update) MODE=update; shift ;;
    --keep)   KEEP="$2"; shift 2 ;;
    --baseline) NAME="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASELINE="$HERE/baselines/${NAME}.sha256"

OUT="${KEEP:-$(mktemp)}"
[ -n "$KEEP" ] || trap 'rm -f "$OUT"' EXIT

"$BIN" > "$OUT"
HASH="$(sha256sum "$OUT" | awk '{print $1}')"
LINES="$(wc -l < "$OUT")"

if [ "$MODE" = update ]; then
  mkdir -p "$(dirname "$BASELINE")"
  printf '%s  %s lines\n' "$HASH" "$LINES" > "$BASELINE"
  echo "recorded baseline: $HASH ($LINES lines)"
  exit 0
fi

[ -f "$BASELINE" ] || { echo "no baseline at $BASELINE; run with --update" >&2; exit 1; }
WANT="$(awk '{print $1}' "$BASELINE")"

if [ "$HASH" = "$WANT" ]; then
  echo "PASS: $NAME bit-identical to baseline ($LINES lines)"
  exit 0
fi

echo "FAIL: $NAME output changed"
echo "  baseline: $WANT"
echo "  now:      $HASH"
echo "  rerun both revisions with --keep and diff the outputs to localize it"
exit 1
