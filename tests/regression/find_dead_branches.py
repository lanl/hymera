#!/usr/bin/env python3
"""Find if-statements whose condition is a compile-time constant false.

Matches the forms this codebase uses to switch code off:
    if (0)            if(0)
    if (0 && ...)     if (... && 0)
    if (... && 0.0)   (with or without spaces)

For each one it reports the extent of the guarded block (to the matching brace)
and whether there is an `else` -- because deleting `if (0) {A} else {B}` must keep
B, unbraced, and that is a different edit from deleting a plain `if (0) {A}`.

Usage:
    find_dead_branches.py src/mhd/*.c            # list
    find_dead_branches.py src/mhd/*.c --remove   # delete them in place

Removal replaces each deleted block with a `#line` directive so the __LINE__ of
every following line is unchanged; verify with check_line_numbers.py.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

# Condition text, between the parentheses, that is constant false.
FALSE_COND = re.compile(
    r"^\s*(0|0\s*&&.*|.*&&\s*0(\.0*)?)\s*$"
)
IF_START = re.compile(r"^(\s*)(?:\}\s*else\s+)?if\s*\(")
LINE_DIRECTIVE = re.compile(r"^\s*#\s*line\s+(\d+)")


def strip_strings_comments(s: str) -> str:
    s = re.sub(r'"(\\.|[^"\\])*"', '""', s)
    s = re.sub(r"'(\\.|[^'\\])*'", "''", s)
    s = re.sub(r"/\*.*?\*/", "", s)
    return re.sub(r"//.*", "", s)


def find_paren_close(text: str, start: int) -> int:
    """Index of the ')' matching the '(' at text[start]."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "(":
            depth += 1
        elif text[i] == ")":
            depth -= 1
            if depth == 0:
                return i
    return -1


def scan(lines: list[str]):
    """Yield dicts describing each constant-false if-statement."""
    # Work on a comment- and string-stripped copy so braces in comments do not
    # confuse the matcher; indices still line up with the original.
    clean = [strip_strings_comments(l) for l in lines]
    in_block = False
    for i, l in enumerate(lines):
        s = l.strip()
        if in_block:
            if "*/" in s:
                in_block = False
            continue
        if s.startswith("/*") and "*/" not in s:
            in_block = True
            continue
        if s.startswith("//"):
            continue
        m = IF_START.match(clean[i])
        if not m or re.match(r"^\s*\}\s*else", clean[i]):
            # `} else if (0)` chains are left for a human: removing one arm of a
            # chain is a different edit.
            continue

        # Gather the condition, which may span lines.
        joined = clean[i]
        k = i
        p = joined.index("(", m.end() - 1)
        close = find_paren_close(joined, p)
        while close < 0 and k + 1 < len(lines):
            k += 1
            joined += "\n" + clean[k]
            close = find_paren_close(joined, p)
        if close < 0:
            continue
        cond = joined[p + 1 : close]
        if not FALSE_COND.match(cond.replace("\n", " ")):
            continue

        # The guarded statement must be a braced block opening on this line.
        rest = joined[close + 1 :]
        if not rest.lstrip().startswith("{"):
            yield {"line": i, "kind": "unbraced", "cond": cond.strip()}
            continue

        # Find the matching close brace.
        depth = 0
        end = None
        started = False
        for j in range(k, len(lines)):
            seg = clean[j] if j != k else rest
            if j == k:
                seg = rest
            for ch in seg:
                if ch == "{":
                    depth += 1
                    started = True
                elif ch == "}":
                    depth -= 1
                    if started and depth == 0:
                        end = j
                        break
            if end is not None:
                break
        if end is None:
            continue

        tail = clean[end].split("}", 1)[1] if "}" in clean[end] else ""
        has_else = bool(re.match(r"^\s*else\b", tail)) or (
            end + 1 < len(lines) and re.match(r"^\s*else\b", clean[end + 1])
        )
        yield {"line": i, "end": end, "kind": "else" if has_else else "plain",
               "cond": cond.strip()}


def effective_number(lines: list[str], idx: int) -> int:
    eff = 1
    for k in range(idx):
        m = LINE_DIRECTIVE.match(lines[k])
        eff = int(m.group(1)) if m else eff + 1
    return eff


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--remove", action="store_true")
    args = ap.parse_args()

    grand = {"plain": 0, "else": 0, "unbraced": 0}
    for f in args.files:
        path = Path(f)
        lines = path.read_text().splitlines()
        found = list(scan(lines))
        if not found:
            continue
        nested = set()
        # Ignore blocks nested inside another dead block: the outer removal
        # takes them with it.
        spans = [(d["line"], d["end"]) for d in found if d["kind"] == "plain"]
        for d in found:
            for a, b in spans:
                if a < d["line"] <= b:
                    nested.add(d["line"])
        top = [d for d in found if d["line"] not in nested]

        counts = {"plain": 0, "else": 0, "unbraced": 0}
        dead_lines = 0
        for d in top:
            counts[d["kind"]] += 1
            if d["kind"] == "plain":
                dead_lines += d["end"] - d["line"] + 1
        print(f"{path.name}: {counts['plain']} removable blocks ({dead_lines} lines), "
              f"{counts['else']} with else (manual), {counts['unbraced']} unbraced (manual)")
        for d in top:
            if d["kind"] != "plain":
                print(f"    {d['kind']:9s} line {d['line'] + 1}: if ({d['cond'][:60]})")
        for k in grand:
            grand[k] += counts[k]

        if args.remove:
            # Delete bottom-up so earlier indices stay valid.
            for d in sorted((d for d in top if d["kind"] == "plain"),
                            key=lambda d: d["line"], reverse=True):
                a, b = d["line"], d["end"]
                # Effective number of the first line after the block, in the
                # current (already partly edited) text -- that is what the
                # directive must restore.
                nxt = effective_number(lines, b + 1)
                lines[a : b + 1] = [f"#line {nxt}"]
            path.write_text("\n".join(lines) + "\n")

    print(f"\nTOTAL: {grand['plain']} removable, {grand['else']} with else, "
          f"{grand['unbraced']} unbraced")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
