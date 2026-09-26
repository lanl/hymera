#!/usr/bin/env python3
"""Compare two normalized disassemblies function by function.

Called by codegen_identity.sh. Splits each listing into per-symbol instruction
blocks and compares only symbols present in both, so that a function moving
within the object -- which happens whenever anything before it is deleted -- does
not register as a change. A function whose instructions differ does.

Exit status is 0 when every shared function is identical, 1 otherwise.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

# objdump labels a function as "ADDR <name>:" after normalization.
SYMBOL = re.compile(r"^ADDR <([^>]+)>:")


def split_functions(path: Path) -> dict[str, list[str]]:
    """Split a normalized listing into per-symbol instruction blocks.

    Alignment padding is dropped. The linker pads each function out to an
    alignment boundary, so removing a function elsewhere in the object shifts how
    much padding its neighbours receive. That is not a behavioural difference, and
    counting it as one would make every deletion appear to change unrelated code.
    """
    blocks: dict[str, list[str]] = {}
    current: str | None = None
    for line in path.read_text(errors="replace").splitlines():
        m = SYMBOL.match(line)
        if m:
            current = m.group(1)
            blocks[current] = []
            continue
        if current is not None:
            text = line.rstrip()
            if text.strip() in ("nop", "nop\t", "udf\t#0"):
                continue
            blocks[current].append(text)
    return blocks


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: compare_asm.py <before.asm> <after.asm>", file=sys.stderr)
        return 2

    before = split_functions(Path(sys.argv[1]))
    after = split_functions(Path(sys.argv[2]))

    shared = sorted(set(before) & set(after))
    removed = sorted(set(before) - set(after))
    added = sorted(set(after) - set(before))

    changed = []
    for name in shared:
        if before[name] != after[name]:
            # Report the first differing instruction, which is what a reviewer
            # actually needs.
            b, a = before[name], after[name]
            where = next(
                (i for i in range(min(len(b), len(a))) if b[i] != a[i]),
                min(len(b), len(a)),
            )
            detail = f"line {where + 1}: "
            detail += f"{b[where] if where < len(b) else '<end>'!r}"
            detail += f" -> {a[where] if where < len(a) else '<end>'!r}"
            changed.append((name, len(b), len(a), detail))

    parts = [f"{len(shared)} shared functions"]
    if removed:
        parts.append(f"{len(removed)} removed")
    if added:
        parts.append(f"{len(added)} added")

    if not changed:
        print(", ".join(parts) + ", all identical")
        return 0

    print(", ".join(parts) + f", {len(changed)} CHANGED:")
    for name, nb, na, detail in changed[:10]:
        print(f"    {name}: {nb} -> {na} instructions, {detail}")
    if len(changed) > 10:
        print(f"    ... and {len(changed) - 10} more")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
