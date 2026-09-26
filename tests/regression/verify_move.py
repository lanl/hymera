#!/usr/bin/env python3
"""Verify that a function-move refactor changed no function body.

Checks the quarantine commits, where functions move out of a compiled source file
into an uncompiled "attic" file and nothing else may change. This proves the move
was byte-exact independently of the compiler, so it is worth running even when the
codegen check passes.

For each source/attic pair:

  * every function now in the attic is byte-identical to what git has for it
  * every function still in the source file is byte-identical to what git has
  * the two files together account for every function the original had
  * no function ended up in both files

Usage:
    verify_move.py <source.c> <attic.c> [--rev REV]
    verify_move.py --all [--rev REV]      # every attic file against its source

REV must be a revision from BEFORE the move.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

# Top-level definitions start at column 0 with one of these return types; a body
# ends at the next line that is exactly "}" at column 0.
DEFINITION = re.compile(
    r"^(?:PetscErrorCode|PetscScalar|PetscInt|PetscBool|PetscReal|void|int)"
    r"[\s*]+(\w+)\s*\("
)


def split_functions(text: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    lines = text.splitlines()
    name: str | None = None
    start = 0
    for i, line in enumerate(lines):
        m = DEFINITION.match(line)
        if m:
            name, start = m.group(1), i
        elif line == "}" and name is not None:
            out[name] = lines[start : i + 1]
            name = None
    return out


def repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "src" / "mhd").is_dir():
            return parent
    raise SystemExit("cannot locate repo root")


def check(source: Path, attic: Path, rev: str, root: Path) -> bool:
    rel = source.relative_to(root).as_posix()
    try:
        original = split_functions(
            subprocess.run(["git", "show", f"{rev}:{rel}"], cwd=root,
                           capture_output=True, text=True, check=True).stdout
        )
    except subprocess.CalledProcessError:
        print(f"{source.name}: cannot read {rev}:{rel}", file=sys.stderr)
        return False

    live = split_functions(source.read_text())
    moved = split_functions(attic.read_text()) if attic.is_file() else {}

    # Only consider functions that existed in the chosen revision; an attic file
    # may also hold functions moved from a different source file in an earlier
    # commit.
    moved_from_here = {k: v for k, v in moved.items() if k in original}

    ok = True
    altered_moved = [f for f in moved_from_here if original[f] != moved_from_here[f]]
    altered_live = [f for f in live if f in original and original[f] != live[f]]
    lost = sorted(set(original) - set(live) - set(moved))
    both = sorted(set(live) & set(moved))

    print(f"{source.name}:")
    print(f"  {rev}: {len(original)} functions -> "
          f"{len(live)} live + {len(moved_from_here)} moved = "
          f"{len(live) + len(moved_from_here)}")

    if altered_moved:
        ok = False
        print(f"  MOVED BODIES ALTERED: {len(altered_moved)} {altered_moved[:5]}")
    else:
        print(f"  moved bodies byte-identical: {len(moved_from_here)}")

    if altered_live:
        ok = False
        print(f"  SURVIVING BODIES ALTERED: {len(altered_live)} {altered_live[:5]}")
    else:
        print(f"  surviving bodies byte-identical: {len(live)}")

    if lost:
        ok = False
        print(f"  LOST (present in {rev}, now nowhere): {lost}")
    if both:
        ok = False
        print(f"  DUPLICATED (in both files): {both}")
    if len(live) + len(moved_from_here) != len(original):
        ok = False
        print("  COUNT MISMATCH: live + moved does not equal the original count")

    # #line directives keep surviving lines at their original numbers, which is
    # what stops the __LINE__ baked into PetscCall/SETERRQ from shifting.
    directives = sum(1 for l in source.read_text().splitlines()
                     if l.startswith("#line"))
    note = "" if directives else "   <-- none present; check codegen carefully"
    print(f"  #line directives: {directives}{note}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", nargs="?")
    ap.add_argument("attic", nargs="?")
    ap.add_argument("--rev", default="HEAD")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()
    root = repo_root()

    if args.all:
        atticdir = root / "src" / "mhd" / "attic"
        if not atticdir.is_dir():
            print("no attic directory")
            return 0
        every_ok = True
        for a in sorted(atticdir.glob("*_attic.c")):
            src = root / "src" / "mhd" / a.name.replace("_attic.c", ".c")
            if not src.is_file():
                print(f"{a.name}: no matching source file")
                every_ok = False
                continue
            every_ok &= check(src, a, args.rev, root)
            print()
        print("ALL MOVES VERIFIED" if every_ok else "VERIFICATION FAILED")
        return 0 if every_ok else 1

    if not args.source or not args.attic:
        ap.error("give both source and attic, or use --all")
    return 0 if check(Path(args.source).resolve(), Path(args.attic).resolve(),
                      args.rev, root) else 1


if __name__ == "__main__":
    raise SystemExit(main())
