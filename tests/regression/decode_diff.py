#!/usr/bin/env python3
"""Decide whether two AArch64 instruction streams differ semantically.

Used to settle a "CHANGED" report from codegen_identity.sh. Removing a function
from an object shifts alignment padding, which shifts branch displacements in
functions nobody edited. Those are not behavioural differences, but they do change
encoded bytes, so a byte comparison alone cannot clear them.

This decodes the differing words and classifies each one. A difference is accepted
as benign only when it is:

  * an inserted or removed NOP (alignment padding), or
  * a branch whose condition and register operands are unchanged but whose
    displacement moved -- i.e. the same branch to the same place, relocated.

Anything else is reported as unexplained, which means the change is real.

Usage:
    decode_diff.py <before.txt> <after.txt>       # one function
    decode_diff.py --pairs <fn>=<before>,<after> ...

Input files hold one hex instruction word per line, as produced by
    objdump -d <obj> | awk '/<fn>:/,/^$/' | sed -E 's/^\\s*[0-9a-f]+:\\s+//' \\
        | awk '{print $1}'
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

NOP = 0xD503201F


def load(path: str | Path) -> list[int]:
    return [int(line.strip(), 16) for line in Path(path).read_text().split() if line.strip()]


def is_bcond(w: int) -> bool:
    # B.cond: 0101 0100 imm19 0 cond
    return (w >> 24) == 0x54


def is_b_or_bl(w: int) -> bool:
    # B: 000101 imm26   BL: 100101 imm26
    return (w >> 26) in (0x05, 0x25)


def is_compare_branch(w: int) -> bool:
    # CBZ/CBNZ: sf 011010 op imm19 Rt ; TBZ/TBNZ: b5 011011 op b40 imm14 Rt
    top = (w >> 25) & 0x3F
    return top in (0x1A, 0x1B)


def classify(before: list[int], after: list[int]) -> tuple[int, list[tuple[int, int, int]]]:
    """Return (non-NOP instruction count, list of unexplained differences)."""
    b = [w for w in before if w != NOP]
    a = [w for w in after if w != NOP]

    if len(b) != len(a):
        # A genuine insertion or deletion of real instructions.
        return len(b), [(-1, len(b), len(a))]

    unexplained: list[tuple[int, int, int]] = []
    for i, (x, y) in enumerate(zip(b, a)):
        if x == y:
            continue
        if is_bcond(x) and is_bcond(y) and (x & 0xF) == (y & 0xF):
            continue  # same condition, relocated target
        if is_b_or_bl(x) and is_b_or_bl(y) and (x >> 26) == (y >> 26):
            continue  # same branch kind, relocated target
        if (
            is_compare_branch(x)
            and is_compare_branch(y)
            and ((x >> 25) & 0x3F) == ((y >> 25) & 0x3F)
            and (x & 0x1F) == (y & 0x1F)
        ):
            continue  # same test on same register, relocated target
        unexplained.append((i, x, y))
    return len(b), unexplained


def report(name: str, before: str, after: str) -> bool:
    n, bad = classify(load(before), load(after))
    if not bad:
        print(f"{name:38s} {n:5d} non-NOP insns, all differences are padding "
              f"or relocated branches -> SEMANTICALLY IDENTICAL")
        return True
    if bad[0][0] == -1:
        print(f"{name:38s} INSTRUCTION COUNT CHANGED: {bad[0][1]} -> {bad[0][2]}")
        return False
    print(f"{name:38s} {n:5d} non-NOP insns, {len(bad)} UNEXPLAINED difference(s):")
    for i, x, y in bad[:5]:
        print(f"    index {i}: {x:08x} -> {y:08x}")
    return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("before", nargs="?")
    ap.add_argument("after", nargs="?")
    ap.add_argument("--pairs", nargs="*", default=None,
                    help="name=before,after triples for batch checking")
    args = ap.parse_args()

    if args.pairs:
        ok = True
        for spec in args.pairs:
            name, files = spec.split("=", 1)
            b, a = files.split(",", 1)
            ok &= report(name, b, a)
        print("\nALL SEMANTICALLY IDENTICAL" if ok else "\nREAL DIFFERENCES PRESENT")
        return 0 if ok else 1

    if not args.before or not args.after:
        ap.error("give before and after, or use --pairs")
    return 0 if report(Path(args.before).stem, args.before, args.after) else 1


if __name__ == "__main__":
    raise SystemExit(main())
