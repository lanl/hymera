#!/usr/bin/env python
"""Extract a regression fingerprint from an MHD solver log.

The solver already prints, unconditionally, everything needed to detect a
changed solution. This script pulls those numbers out of stdout into a stable,
diffable JSON document so that two runs can be compared without re-running
anything.

Sources, all in the existing code:

  AppCtxView          src/mhd/mfd_config.c        ~120 labelled config lines
  mhd_initialize      src/mhd/mhd.c:158-170       derived parameters
  Monitor             src/mhd/ts_functions.c      per-step norms, div B, currents
  mhd_step            src/mhd/mhd.c:660           TSConvergedReason + step count
  PETSc monitors      default_petsc_options.c     SNES/KSP iteration counts

Values are stored at full precision as strings, exactly as printed. The solver
prints with %g (6 significant digits), so the log fingerprint resolves changes
only to that precision -- it is the cheap first-line check. Bit-level
comparison is the state comparator's job (compare_state.py), not this one.

Usage:
    fingerprint.py run.log -o fingerprint.json
    fingerprint.py --compare base.json cand.json [--rtol 1e-12]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

# --- patterns ----------------------------------------------------------------
# Each entry maps a per-step scalar to the regex that captures it. Keep the
# source line cited so a future format change is traceable.

RE_CORRECTOR = re.compile(
    r"Timestep\s+(\d+)\s+\(CORRECTOR\):\s*step size = (\S+), time = (\S+), "
    r"2-norm of X\^\{n\+1\} - X\^\{n\} = (\S+), "
    r"max norm of X\^\{n\+1\} - X\^\{n\} = (\S+)"
)  # ts_functions.c:17593

RE_ERROR = re.compile(
    r"Timestep\s+(\d+):\s*step size = (\S+), time = (\S+), "
    r"2-norm error = (\S+), max norm error = (\S+)"
)  # ts_functions.c:17715  -- only printed when ictype != 9

RE_DIVB = re.compile(
    r"Timestep\s+(\d+):\s*step size = \S+, time = (\S+), "
    r"max norm error of divergence of B = (\S+)"
)  # ts_functions.c:17770

RE_DIVB_RATIO = re.compile(r"max_norm\(div B\) / max_norm\(B\) = (\S+)")  # :17780

RE_CURRENT_IN = re.compile(r"Current intensity inside plasma = (\S+)")  # :17759
RE_CURRENT_OUT = re.compile(r"Current intensity outside plasma = (\S+)")  # :17760
RE_CURRENT_VV = re.compile(r"Current intensity inside vacuum vessel = (\S+)")

RE_CONVERGED = re.compile(r"(\w+) at time (\S+) after (\d+) steps")  # mhd.c:660

RE_SNES = re.compile(r"^\s*(\d+) SNES Function norm (\S+)")
RE_KSP = re.compile(r"^\s*(\d+) KSP Residual norm (\S+)")

# AppCtxView emits "  name = value" lines; capture them as the config block.
RE_CONFIG = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_\[\]]*)\s+= (.*?)\s*$")

# Lines that vary run to run without indicating a changed solution.
NOISE = (
    "PETSC ERROR",
    "Option left:",
    "Configure options:",
    "with 1 MPI process",
    "Reading X vector",
    "Reading from file",
)


def parse_log(path: Path) -> dict:
    fp: dict = {
        "source": path.name,
        "config": {},
        "derived": {},
        "steps": [],
        "div_b": [],
        "div_b_ratio": [],
        "currents": [],
        "mms_error": [],
        "converged": [],
        "snes_counts": [],
        "ksp_counts": [],
    }

    in_config = False
    snes_run: list[str] = []
    ksp_run: list[str] = []

    for line in path.read_text(errors="replace").splitlines():
        if any(n in line for n in NOISE):
            continue

        # AppCtxView block delimits the configuration fingerprint.
        if "AppCtx dump" in line:
            in_config = True
            continue
        if in_config:
            if line.startswith("=====") or "Read initial data" in line:
                in_config = False
            else:
                m = RE_CONFIG.match(line)
                if m and not line.strip().startswith("--"):
                    fp["config"][m.group(1)] = m.group(2)
                continue

        if m := RE_CORRECTOR.search(line):
            fp["steps"].append({
                "step": int(m.group(1)), "dt": m.group(2), "time": m.group(3),
                "dX_l2": m.group(4), "dX_max": m.group(5),
            })
        elif m := RE_ERROR.search(line):
            fp["mms_error"].append({
                "step": int(m.group(1)), "dt": m.group(2), "time": m.group(3),
                "err_l2": m.group(4), "err_max": m.group(5),
            })
        elif m := RE_DIVB.search(line):
            fp["div_b"].append({
                "step": int(m.group(1)), "time": m.group(2), "max_div_b": m.group(3),
            })
        elif m := RE_DIVB_RATIO.search(line):
            fp["div_b_ratio"].append(m.group(1))
        elif m := RE_CURRENT_IN.search(line):
            fp["currents"].append({"where": "plasma", "value": m.group(1)})
        elif m := RE_CURRENT_OUT.search(line):
            fp["currents"].append({"where": "outside", "value": m.group(1)})
        elif m := RE_CURRENT_VV.search(line):
            fp["currents"].append({"where": "vacuum_vessel", "value": m.group(1)})
        elif m := RE_CONVERGED.search(line):
            fp["converged"].append({
                "reason": m.group(1), "time": m.group(2), "steps": int(m.group(3)),
            })
        elif m := RE_SNES.match(line):
            it = int(m.group(1))
            if it == 0 and snes_run:
                fp["snes_counts"].append(len(snes_run) - 1)
                snes_run = []
            snes_run.append(m.group(2))
        elif m := RE_KSP.match(line):
            it = int(m.group(1))
            if it == 0 and ksp_run:
                fp["ksp_counts"].append(len(ksp_run) - 1)
                ksp_run = []
            ksp_run.append(m.group(2))
        else:
            # Derived parameters printed once by mhd_initialize.
            for key, pat in (
                ("lundquist", r"Lundquist number:\s*(\S+)"),
                ("reynolds", r"Reynolds parameter for viscosity:\s*(\S+)"),
                ("resistive_time", r"Characteristic resistive time:\s*(\S+)"),
                ("dt_normalized", r"normalized dt:\s*(\S+)"),
                ("mesh", r"Using a (\d+ x \d+ x \d+) mesh"),
            ):
                if m2 := re.search(pat, line):
                    fp["derived"][key] = m2.group(1)

    if snes_run:
        fp["snes_counts"].append(len(snes_run) - 1)
    if ksp_run:
        fp["ksp_counts"].append(len(ksp_run) - 1)

    return fp


def _num(s: str) -> float | None:
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def compare(base: dict, cand: dict, rtol: float) -> int:
    """Compare two fingerprints. Returns 0 if they agree within rtol."""
    problems: list[str] = []

    # Configuration must match exactly: a differing config means the two runs
    # are not the same run, so any numerical agreement is meaningless.
    for key in sorted(set(base["config"]) | set(cand["config"])):
        b, c = base["config"].get(key), cand["config"].get(key)
        if b != c:
            problems.append(f"config {key}: {b!r} -> {c!r}")

    for key in sorted(set(base["derived"]) | set(cand["derived"])):
        b, c = base["derived"].get(key), cand["derived"].get(key)
        if b != c:
            problems.append(f"derived {key}: {b!r} -> {c!r}")

    # Structural properties must match exactly.
    for key in ("converged",):
        if base[key] != cand[key]:
            problems.append(f"{key}: {base[key]} -> {cand[key]}")

    # Per-step scalars compare within rtol.
    for block, fields in (
        ("steps", ("dX_l2", "dX_max")),
        ("mms_error", ("err_l2", "err_max")),
        ("div_b", ("max_div_b",)),
    ):
        b_rows, c_rows = base[block], cand[block]
        if len(b_rows) != len(c_rows):
            problems.append(
                f"{block}: {len(b_rows)} entries -> {len(c_rows)} entries"
            )
            continue
        for i, (br, cr) in enumerate(zip(b_rows, c_rows)):
            for f in fields:
                bv, cv = _num(br.get(f)), _num(cr.get(f))
                if bv is None or cv is None:
                    if br.get(f) != cr.get(f):
                        problems.append(f"{block}[{i}].{f}: {br.get(f)} -> {cr.get(f)}")
                    continue
                scale = max(abs(bv), abs(cv), 1e-300)
                rel = abs(bv - cv) / scale
                if rel > rtol:
                    problems.append(
                        f"{block}[{i}].{f}: {bv:.17g} -> {cv:.17g}  (rel {rel:.3e})"
                    )

    for i, (bv, cv) in enumerate(zip(base["currents"], cand["currents"])):
        b, c = _num(bv["value"]), _num(cv["value"])
        if b is not None and c is not None:
            scale = max(abs(b), abs(c), 1e-300)
            rel = abs(b - c) / scale
            if rel > rtol:
                problems.append(
                    f"currents[{i}] ({bv['where']}): {b:.17g} -> {c:.17g} "
                    f"(rel {rel:.3e})"
                )

    # Iteration counts are reported but not fatal: a shifted Krylov path is
    # expected for REASSOC-class changes. Surface it, do not fail on it.
    if base["ksp_counts"] != cand["ksp_counts"]:
        print(
            f"note: KSP iteration counts changed "
            f"{base['ksp_counts']} -> {cand['ksp_counts']}"
        )
    if base["snes_counts"] != cand["snes_counts"]:
        print(
            f"note: SNES iteration counts changed "
            f"{base['snes_counts']} -> {cand['snes_counts']}"
        )

    if problems:
        print(f"FAIL: {len(problems)} difference(s) beyond rtol={rtol:g}")
        for p in problems[:40]:
            print(f"  {p}")
        if len(problems) > 40:
            print(f"  ... and {len(problems) - 40} more")
        return 1

    print(f"PASS: fingerprints agree within rtol={rtol:g}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("logs", nargs="+", help="log file, or two .json with --compare")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--rtol", type=float, default=1e-12)
    args = ap.parse_args()

    if args.compare:
        if len(args.logs) != 2:
            ap.error("--compare needs exactly two fingerprint json files")
        base = json.loads(Path(args.logs[0]).read_text())
        cand = json.loads(Path(args.logs[1]).read_text())
        return compare(base, cand, args.rtol)

    fp = parse_log(Path(args.logs[0]))
    text = json.dumps(fp, indent=2, sort_keys=True)
    if args.out:
        Path(args.out).write_text(text + "\n")
        print(f"wrote {args.out}")
    else:
        print(text)

    print(
        f"summary: {len(fp['config'])} config keys, {len(fp['steps'])} corrector "
        f"steps, {len(fp['mms_error'])} MMS error rows, {len(fp['div_b'])} div-B "
        f"rows, {len(fp['currents'])} current readings, "
        f"{len(fp['converged'])} convergence reports",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
