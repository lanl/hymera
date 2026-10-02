#!/usr/bin/env python3
"""Check that a refactor kept every surviving line at its original __LINE__.

Why this matters: PETSc's PetscCall, SETERRQ and PetscCheck macros expand
__LINE__ into the generated code. Deleting or inserting lines above a function
therefore changes the instructions of that function even though nobody edited
it. The refactors in src/mhd compensate with `#line N` directives; this script
proves the compensation is right.

It computes the *effective* line number of every physical line -- what the
preprocessor reports as __LINE__, honouring `#line N` directives -- in both the
old and the new version of a file, then matches every macro line that survives
(identical text, in order) and checks its effective number is unchanged.

Only lines that actually carry __LINE__ into the build are checked: those
containing PetscCall, SETERRQ, PetscCheck or PetscAssert. Comments and plain
statements may move freely.

Usage:
    check_line_numbers.py src/mhd/ts_functions.c            # against HEAD
    check_line_numbers.py src/mhd/ts_functions.c --rev REV
"""

from __future__ import annotations

import argparse
import difflib
import re
import subprocess
import sys
from pathlib import Path

LINE_DIRECTIVE = re.compile(r"^\s*#\s*line\s+(\d+)")
MACRO = re.compile(r"\b(PetscCall|SETERRQ|PetscCheck|PetscAssert)\b")


def effective_lines(text: str) -> list[tuple[int, str]]:
    """Return (effective_line_number, text) for every non-directive line."""
    out: list[tuple[int, str]] = []
    eff = 1
    for line in text.splitlines():
        m = LINE_DIRECTIVE.match(line)
        if m:
            eff = int(m.group(1))
            continue
        out.append((eff, line))
        eff += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--rev", default="HEAD")
    args = ap.parse_args()

    path = Path(args.file)
    root = Path(subprocess.run(["git", "rev-parse", "--show-toplevel"],
                               capture_output=True, text=True, check=True).stdout.strip())
    rel = path.resolve().relative_to(root).as_posix()
    old_text = subprocess.run(["git", "show", f"{args.rev}:{rel}"], cwd=root,
                              capture_output=True, text=True, check=True).stdout
    new_text = path.read_text()

    old = effective_lines(old_text)
    new = effective_lines(new_text)

    old_txt = [t for _, t in old]
    new_txt = [t for _, t in new]
    sm = difflib.SequenceMatcher(None, old_txt, new_txt, autojunk=False)

    checked = 0
    mismatches: list[str] = []
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op != "equal":
            continue
        for k in range(i2 - i1):
            eo, t = old[i1 + k]
            en, _ = new[j1 + k]
            if not MACRO.search(t):
                continue
            checked += 1
            if eo != en:
                mismatches.append(f"  {t.strip()[:80]!r}: line {eo} -> {en}")

    print(f"{path.name}: {checked} surviving macro lines checked against {args.rev}")
    if mismatches:
        print(f"FAIL: {len(mismatches)} changed effective line number")
        for m in mismatches[:15]:
            print(m)
        if len(mismatches) > 15:
            print(f"  ... and {len(mismatches) - 15} more")
        return 1
    print("PASS: every surviving macro line keeps its __LINE__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
