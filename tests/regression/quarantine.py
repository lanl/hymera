#!/usr/bin/env python3
"""Move unreachable functions from a compiled source file into its attic file.

The attic, src/mhd/attic/<name>_attic.c, is not compiled. Moving rather than
deleting keeps the code searchable for physics review. This automates the
procedure every quarantine in this repository has followed:

  1. each function body is moved verbatim, from its definition line to the next
     line that is exactly "}" at column 0;
  2. the removed block is replaced by a `#line N` directive, so that every later
     line keeps its __LINE__ -- PETSc's PetscCall/SETERRQ bake __LINE__ into the
     generated code, so without this, functions nobody touched would change;
  3. the function's declaration is removed from the matching header.

It does not decide what is dead. Establish that first (no callers, transitively,
and no undefined reference in any object -- `nm -u`), then verify afterwards with
verify_move.py, check_line_numbers.py and the regression tests.

Usage:
    quarantine.py src/mhd/monitor_functions.c DumpVelocity_Cell DumpEdgeField
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

LINE_DIRECTIVE = re.compile(r"^\s*#\s*line\s+(\d+)")


def effective_number(lines: list[str], idx: int) -> int:
    """__LINE__ the preprocessor assigns to lines[idx], honouring #line."""
    eff = 1
    for k in range(idx):
        m = LINE_DIRECTIVE.match(lines[k])
        eff = int(m.group(1)) if m else eff + 1
    return eff


def find_function(lines: list[str], name: str) -> tuple[int, int]:
    pat = re.compile(
        r"^(?:static\s+)?(?:inline\s+)?(?:PetscErrorCode|PetscScalar|PetscInt|"
        r"PetscBool|PetscReal|void|int)[\s*]+" + re.escape(name) + r"\s*\("
    )
    starts = [i for i, l in enumerate(lines) if pat.match(l) and not l.rstrip().endswith(";")]
    if len(starts) != 1:
        raise SystemExit(f"{name}: expected one definition, found {len(starts)}")
    s = starts[0]
    e = next(i for i in range(s + 1, len(lines)) if lines[i] == "}")
    return s, e


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    src = Path(sys.argv[1])
    names = sys.argv[2:]
    attic = src.parent / "attic" / (src.stem + "_attic.c")
    header = src.with_suffix(".h")

    lines = src.read_text().splitlines()
    spans = sorted((find_function(lines, n) + (n,) for n in names), reverse=True)

    moved: list[tuple[str, list[str]]] = []
    for s, e, n in spans:  # bottom-up, so earlier indices stay valid
        nxt = effective_number(lines, e + 1)
        moved.append((n, lines[s : e + 1]))
        lines[s : e + 1] = [f"#line {nxt}"]
    src.write_text("\n".join(lines) + "\n")

    if not attic.exists():
        raise SystemExit(f"{attic} does not exist; create it with the source's #include block first")
    with attic.open("a") as fh:
        for n, body in reversed(moved):
            fh.write("\n" + "\n".join(body) + "\n")

    if header.exists():
        h = header.read_text().splitlines()
        keep = [l for l in h
                if not any(re.search(r"\b" + re.escape(n) + r"\s*\(", l) and l.rstrip().endswith(";")
                           for n in names)]
        header.write_text("\n".join(keep) + "\n")
        removed_decls = len(h) - len(keep)
    else:
        removed_decls = 0

    total = sum(len(b) for _, b in moved)
    print(f"moved {len(moved)} functions ({total} lines) from {src.name} to {attic.name}; "
          f"removed {removed_decls} declarations from {header.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
