#!/usr/bin/env bash
# Prove that a source change did not alter the generated code of any surviving
# function.
#
# This is the strongest available check for a deletion or a rename, and it needs
# no fixture, no grid, no MPI and no run. It is valid because of three verified
# properties of this build:
#
#   1. no -flto, so there is no cross-translation-unit inlining path
#   2. no statics and no file-scope data anywhere in mhd_core, so removing an
#      extern function cannot change any surviving function's code
#   3. no -ffast-math, so the compiler does not reassociate arithmetic
#
# If any of those change, this proof is void. The script asserts 1 and 3 itself.
#
# LIMITATION: this compares instructions, not data. Changing a string literal --
# for example a PetscLogEventRegister label -- alters .rodata while leaving every
# instruction identical, so it passes. That is correct for behaviour preservation
# (the computed result is unchanged) but means this check alone does not prove a
# diff touched no data. Verify string changes separately with `strings`.
#
# SECOND LIMITATION, and the one that actually bites: references into .rodata
# carry no local symbol, so objdump annotates them with the nearest PRECEDING
# defined function symbol. Deleting the first function in an object therefore
# relabels every such reference in every surviving function, e.g.
#   adrp x0, <FormMaterialPropertiesMatrix>  ->  adrp x0, <FormDiscreteDivergence>
# while the encoded instruction bytes are unchanged. That shows up here as a
# spurious "CHANGED" on many functions at once, all with identical instruction
# counts and all differing only in that annotation.
#
# When that pattern appears, settle it with `rawbytes`, which compares the encoded
# instruction bytes and ignores annotation entirely. Do NOT loosen the normalizer
# to make the symptom go away: `bl <name>` call targets must stay comparable, and
# blurring symbol names would hide a genuinely retargeted call.
#
# Usage:
#   codegen_identity.sh snapshot <name>     # record the current object code
#   codegen_identity.sh compare  <name>     # rebuild and diff against it
#   codegen_identity.sh rawbytes <name> <object> <function>
#                                           # compare encoded bytes of one function
#
# Typical session:
#   ./tests/regression/codegen_identity.sh snapshot before
#   ...make the change...
#   ./tests/regression/codegen_identity.sh compare before
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OBJDIR="$REPO/build/src/CMakeFiles/mhd_core.dir/mhd"
SNAPDIR="${TMPDIR:-/tmp}/hymera-codegen"

MODE="${1:-}"
NAME="${2:-}"
[ -n "$MODE" ] && [ -n "$NAME" ] || { sed -n '2,20p' "$0"; exit 2; }

assert_flags() {
  local db="$REPO/build/compile_commands.json"
  [ -f "$db" ] || { echo "missing $db; configure with CMAKE_EXPORT_COMPILE_COMMANDS=On" >&2; exit 1; }
  local bad
  bad=$(grep -oE '\-flto[^ "]*|\-ffast-math|\-Ofast' "$db" | sort -u || true)
  if [ -n "$bad" ]; then
    echo "REFUSING: this proof requires no LTO and no fast-math, but found:" >&2
    echo "$bad" >&2
    exit 1
  fi
}

# Normalize a disassembly so that only the instruction semantics remain.
# Addresses, file offsets and relocation targets shift whenever anything upstream
# in the object moves, and carry no meaning for this comparison; symbol names do,
# so they are kept.
normalize() {
  sed -E \
    -e 's/^[[:space:]]*[0-9a-f]+:[[:space:]]*//' \
    -e 's/^([0-9a-f]+) </ADDR </' \
    -e 's/\b0x[0-9a-f]+\b/0xADDR/g' \
    -e 's/<([A-Za-z_][A-Za-z0-9_.]*)\+0x[0-9a-f]+>/<\1+OFF>/g' \
    -e 's/[[:space:]][0-9a-f]+[[:space:]]+</ ADDR </g' \
    -e '/^[[:space:]]*$/d' \
    -e '/^Disassembly of section/d' \
    -e '/file format/d'
}

snapshot() {
  assert_flags
  local dest="$SNAPDIR/$NAME"
  rm -rf "$dest"; mkdir -p "$dest"

  shopt -s nullglob
  local objs=("$OBJDIR"/*.o)
  shopt -u nullglob
  [ ${#objs[@]} -gt 0 ] || { echo "no objects in $OBJDIR; build mhd_core first" >&2; exit 1; }

  for o in "${objs[@]}"; do
    local b; b="$(basename "$o" .o)"
    objdump -d --no-show-raw-insn "$o" | normalize > "$dest/$b.asm"
    nm --defined-only "$o" | awk '{print $2, $3}' | sort > "$dest/$b.defined"
    nm -u "$o"            | awk '{print $NF}'     | sort > "$dest/$b.undefined"
  done
  echo "snapshot '$NAME': ${#objs[@]} objects -> $dest"
  wc -l "$dest"/*.asm | tail -1
}

compare() {
  assert_flags
  local base="$SNAPDIR/$NAME"
  [ -d "$base" ] || { echo "no snapshot named '$NAME' (looked in $base)" >&2; exit 1; }

  echo "rebuilding mhd_core"
  cmake --build "$REPO/build" --target mhd_core -j "$(nproc)" >/dev/null

  # Not declared local: the EXIT trap below runs after this function returns, so
  # a local would be out of scope by the time the trap expands it.
  now="$(mktemp -d)"
  trap 'rm -rf "$now"' EXIT
  for o in "$OBJDIR"/*.o; do
    local b; b="$(basename "$o" .o)"
    objdump -d --no-show-raw-insn "$o" | normalize > "$now/$b.asm"
    nm --defined-only "$o" | awk '{print $2, $3}' | sort > "$now/$b.defined"
    nm -u "$o"            | awk '{print $NF}'     | sort > "$now/$b.undefined"
  done

  local fail=0

  echo
  echo "=== symbols removed (expected: exactly the intended deletions) ==="
  for f in "$base"/*.defined; do
    local b; b="$(basename "$f" .defined)"
    if [ -f "$now/$b.defined" ]; then
      comm -23 "$f" "$now/$b.defined" | sed "s/^/  $b: -/" || true
    fi
  done

  echo
  echo "=== symbols added (expected: none, for a pure deletion) ==="
  for f in "$now"/*.defined; do
    local b; b="$(basename "$f" .defined)"
    if [ -f "$base/$b.defined" ]; then
      comm -13 "$base/$b.defined" "$f" | sed "s/^/  $b: +/" || true
    fi
  done

  echo
  echo "=== newly undefined references (expected: none) ==="
  for f in "$now"/*.undefined; do
    local b; b="$(basename "$f" .undefined)"
    if [ -f "$base/$b.undefined" ]; then
      local added
      added=$(comm -13 "$base/$b.undefined" "$f" || true)
      if [ -n "$added" ]; then
        echo "$added" | sed "s/^/  $b: +/"
        fail=1
      fi
    fi
  done

  # The real test. Extract each surviving function's instruction block and
  # compare it in isolation, so that a function merely moving within the object
  # does not register as a change.
  echo
  echo "=== per-function instruction comparison ==="
  for f in "$base"/*.asm; do
    local b; b="$(basename "$f" .asm)"
    [ -f "$now/$b.asm" ] || { echo "  $b: OBJECT DISAPPEARED"; fail=1; continue; }

    local report rc
    set +e
    report=$("$REPO/tests/regression/compare_asm.py" "$f" "$now/$b.asm")
    rc=$?
    set -e
    [ $rc -eq 0 ] || fail=1
    echo "$report" | sed "s/^/  $b: /"
  done

  echo
  if [ $fail -eq 0 ]; then
    echo "PASS: every surviving function has identical generated code"
  else
    echo "FAIL: generated code changed, or a new undefined reference appeared"
  fi
  return $fail
}

# Compare the encoded instruction bytes of one function, ignoring all objdump
# annotation. Use this to settle a suspected annotation artifact: if the byte
# streams match, the machine code is identical regardless of what labels objdump
# printed.
#
# Needs the pre-change object, which `snapshot` does not keep (it stores only
# disassembly), so this reconstructs it by stashing the working tree, rebuilding,
# extracting, then restoring. That means it must be run with a clean index.
rawbytes() {
  local obj="${3:-}" fn="${4:-}"
  [ -n "$obj" ] && [ -n "$fn" ] || {
    echo "usage: codegen_identity.sh rawbytes <name> <object-basename> <function>" >&2
    exit 2
  }

  local target="$OBJDIR/$obj.o"
  [ -f "$target" ] || { echo "no such object: $target" >&2; exit 1; }

  extract_bytes() {
    objdump -d "$1" \
      | awk -v f="<$2>:" '$0 ~ f {go=1; next} go && /^$/ {exit} go' \
      | sed -E 's/^\s*[0-9a-f]+:\s+//' \
      | awk '{print $1}'
  }

  local after; after="$(mktemp)"
  extract_bytes "$target" "$fn" > "$after"

  echo "stashing working tree to rebuild the pre-change object"
  local stashed=0
  if ! git -C "$REPO" diff --quiet || ! git -C "$REPO" diff --cached --quiet; then
    git -C "$REPO" stash push -q --include-untracked -m "codegen_identity rawbytes" \
      && stashed=1
  fi
  cmake --build "$REPO/build" --target mhd_core -j "$(nproc)" >/dev/null

  local before; before="$(mktemp)"
  extract_bytes "$target" "$fn" > "$before"

  if [ $stashed -eq 1 ]; then
    git -C "$REPO" stash pop -q
    cmake --build "$REPO/build" --target mhd_core -j "$(nproc)" >/dev/null
  fi

  local nb na nd
  nb=$(wc -l < "$before"); na=$(wc -l < "$after")
  nd=$(diff "$before" "$after" | grep -c '^[<>]' || true)
  echo "$fn in $obj.o: $nb -> $na instructions, $nd differing encoded bytes"
  rm -f "$before" "$after"
  [ "$nd" -eq 0 ] && echo "IDENTICAL machine code" || echo "MACHINE CODE DIFFERS"
  [ "$nd" -eq 0 ]
}

case "$MODE" in
  snapshot) snapshot ;;
  compare)  compare ;;
  rawbytes) rawbytes "$@" ;;
  *) echo "unknown mode: $MODE" >&2; exit 2 ;;
esac
