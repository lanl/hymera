#!/usr/bin/env python3
"""Wrap bare PETSc calls in PetscCall(), and add PetscFunctionBeginUser/Return.

The MHD solver discards the return code of roughly 4,400 PETSc calls. PETSc
treats all library errors as fatal, so a discarded code means execution proceeds
on a failed object -- in practice a missing input file becomes a segfault
thousands of lines later instead of a message naming the file.

This applies the transformation mechanically, because doing it by hand across
4,400 sites would introduce more errors than it fixes. It is deliberately
conservative: it only touches statements it can fully recognise, and it reports
everything it declined to touch so the remainder can be reviewed by hand.

What it does, per function:
  1. insert PetscFunctionBeginUser after the opening brace, if absent
  2. rewrite `return (0);` / `return 0;` as PetscFunctionReturn(PETSC_SUCCESS)
  3. wrap whole-statement calls to known error-returning PETSc functions

Ordering matters: PetscCall unwinds through the stack that PetscFunctionBeginUser
pushes, so without step 1 a PETSc traceback still stops at the caller.

What it deliberately does NOT touch:
  * calls whose value is used (assignment, condition, argument, return)
  * calls already inside PetscCall or any other macro
  * multi-line statements, unless the continuation is unambiguous
  * functions that do not return PetscErrorCode -- wrapping those would not
    compile
  * anything inside a comment

Usage:
    add_petsc_checks.py <file.c> [--dry-run] [--report-skipped]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# PETSc calls that return PetscErrorCode. Anything not on this list is left
# alone: an allow-list is the only safe basis here, because several PETSc names
# that look like calls return a value instead (PetscAbsReal returns PetscReal,
# PetscSqr returns a scalar), and PetscAssert is already a checking macro.
ERROR_RETURNING = {
    # DM / DMStag
    "DMCreateGlobalVector", "DMCreateLocalVector", "DMCreateMatrix",
    "DMCreateFieldDecomposition", "DMDestroy", "DMGetCoordinateDM",
    "DMGetCoordinatesLocal", "DMGetLocalVector", "DMGlobalToLocal",
    "DMGlobalToLocalBegin", "DMGlobalToLocalEnd", "DMLocalToGlobal",
    "DMLocalToGlobalBegin", "DMLocalToGlobalEnd", "DMRestoreLocalVector",
    "DMSetFromOptions", "DMSetUp", "DMStagCreate1d", "DMStagCreate2d",
    "DMStagCreate3d", "DMStagCreateCompatibleDMStag", "DMStagGetCorners",
    "DMStagGetDOF", "DMStagGetGlobalSizes", "DMStagGetIsFirstRank",
    "DMStagGetIsLastRank", "DMStagGetLocationSlot", "DMStagGetNumRanks",
    "DMStagMatSetValuesStencil", "DMStagSetUniformCoordinatesExplicit",
    "DMStagVecGetArray", "DMStagVecGetArrayRead", "DMStagVecGetValuesStencil",
    "DMStagVecRestoreArray", "DMStagVecRestoreArrayRead",
    "DMStagVecSetValuesStencil", "DMStagVecSplitToDMDA", "DMDAVecGetArray",
    "DMDAVecGetArrayRead", "DMDAVecRestoreArray", "DMDAVecRestoreArrayRead",
    "DMDAGetCorners", "DMSetMatrixPreallocateOnly",
    # Vec
    "VecAssemblyBegin", "VecAssemblyEnd", "VecAXPBY", "VecAXPBYPCZ", "VecAXPY",
    "VecCopy", "VecDestroy", "VecDuplicate", "VecGetArray", "VecGetArrayRead",
    "VecGetSize", "VecGetSubVector", "VecLoad", "VecNorm", "VecPointwiseMult",
    "VecRestoreArray", "VecRestoreArrayRead", "VecRestoreSubVector", "VecScale",
    "VecSet", "VecShift", "VecView", "VecWAXPY", "VecZeroEntries",
    "VecSetValues", "VecGetOwnershipRange",
    # Mat
    "MatAssemblyBegin", "MatAssemblyEnd", "MatDestroy", "MatDiagonalScale",
    "MatGetDiagonal", "MatMult", "MatMultAdd", "MatSetValues",
    "MatZeroEntries", "MatZeroRowsIS", "MatCreateVecs", "MatScale", "MatAXPY",
    "MatCreateSubMatrix", "MatSetOption", "MatView",
    # TS
    "TSCreate", "TSDestroy", "TSGetConvergedReason", "TSGetDM", "TSGetSNES",
    "TSGetSolution", "TSGetSolveTime", "TSGetStepNumber", "TSGetTime",
    "TSGetTimeStep", "TSMonitorSet", "TSSetDM", "TSSetExactFinalTime",
    "TSSetFromOptions", "TSSetIFunction", "TSSetIJacobian", "TSSetMaxSteps",
    "TSSetMaxTime", "TSSetProblemType", "TSSetRHSFunction", "TSSetSolution",
    "TSSetTime", "TSSetTimeStep", "TSSetType", "TSSetUp", "TSSolve", "TSStep",
    "TSGetMaxTime", "TSSetMaxSNESFailures",
    # KSP / PC / SNES
    "KSPCreate", "KSPDestroy", "KSPGetPC", "KSPSetFromOptions",
    "KSPSetOperators", "KSPSetOptionsPrefix", "KSPSetType", "KSPSetUp",
    "KSPSolve", "PCFieldSplitSetIS", "PCSetFromOptions", "PCSetType",
    "PCSetUp", "PCShellSetApply", "PCShellSetContext", "PCShellSetDestroy",
    "PCShellSetSetUp", "PCFieldSplitGetSubKSP", "PCFieldSplitSetType",
    "SNESSetFromOptions", "SNESSetJacobian", "SNESSetOptionsPrefix",
    "SNESSetType", "SNESSetUseMatrixFree", "SNESGetKSP",
    # IS
    "ISCreateGeneral", "ISCreateStride", "ISDestroy", "ISDifference",
    "ISDuplicate", "ISConcatenate", "ISGetSize", "ISView", "ISSort",
    # Viewer / options / logging / misc
    "PetscViewerASCIIOpen", "PetscViewerBinaryOpen", "PetscViewerDestroy",
    "PetscViewerHDF5Open", "PetscViewerPopFormat", "PetscViewerPushFormat",
    "PetscViewerVTKOpen", "PetscViewerSetUp",
    "PetscOptionsGetBool", "PetscOptionsGetInt", "PetscOptionsGetReal",
    "PetscOptionsInsert", "PetscOptionsInsertFile", "PetscOptionsInsertString",
    "PetscClassIdRegister", "PetscLogEventBegin", "PetscLogEventEnd",
    "PetscLogEventRegister", "PetscObjectSetName", "PetscObjectGetComm",
    "PetscPrintf", "PetscFPrintf", "PetscSNPrintf", "PetscFOpen", "PetscFClose",
    "PetscMalloc1", "PetscFree", "PetscCalloc1", "PetscMemzero",
    "PetscBarrier", "PetscFinalize",
}

# Names that look like calls but are not error-returning. Kept explicit so that a
# future reader sees they were considered rather than missed.
NOT_ERROR_RETURNING = {
    "PetscAbsReal",   # returns PetscReal
    "PetscSqr",       # returns a scalar
    "PetscSqrtReal", "PetscPowReal", "PetscSinReal", "PetscCosReal",
    "PetscExpReal", "PetscLogReal", "PetscMin", "PetscMax",
    "PetscAssert",    # already a checking macro
    "PetscCheck",     # already a checking macro
    "PetscFunctionBeginUser", "PetscFunctionReturn", "PetscCall",
    "PetscError",
}

CALL = re.compile(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*\((.*)\)\s*;\s*$")
DEFINITION = re.compile(
    r"^(?:PetscErrorCode|PetscScalar|PetscInt|PetscBool|PetscReal|void|int)"
    r"[\s*]+(\w+)\s*\("
)
RETURN_ZERO = re.compile(r"^(\s*)return\s*\(?\s*0\s*\)?\s*;\s*$")


def balanced(text: str) -> bool:
    """True if parentheses balance, ignoring those inside string literals."""
    depth = 0
    in_str = False
    esc = False
    for ch in text:
        if esc:
            esc = False
            continue
        if ch == "\\":
            esc = True
            continue
        if ch == '"':
            in_str = not in_str
            continue
        if in_str:
            continue
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                return False
    return depth == 0 and not in_str


def transform(path: Path, dry_run: bool) -> tuple[int, int, int, list[str]]:
    lines = path.read_text().splitlines()
    out: list[str] = []

    wrapped = 0
    begins = 0
    returns = 0
    skipped: list[str] = []

    in_block_comment = False
    # Track whether the current function returns PetscErrorCode; only those may
    # use PetscFunctionBeginUser / PetscFunctionReturn.
    fn_is_errorcode = False
    awaiting_brace = False
    depth = 0

    for i, line in enumerate(lines):
        stripped = line.strip()

        # Track block comments so we never rewrite commented-out code.
        if in_block_comment:
            out.append(line)
            if "*/" in stripped:
                in_block_comment = False
            continue
        if stripped.startswith("/*") and "*/" not in stripped:
            in_block_comment = True
            out.append(line)
            continue
        if stripped.startswith("//") or stripped.startswith("*"):
            out.append(line)
            continue

        m = DEFINITION.match(line)
        if m and line.rstrip().endswith(("{", ")")):
            # PetscErrorCode is a typedef for int, and the public mhd_* entry
            # points in mhd.c are declared `int` for the benefit of the C++ side,
            # so both spellings can carry a PetscCall.
            fn_is_errorcode = line.startswith(("PetscErrorCode", "int "))
            awaiting_brace = True
            depth = 0
            out.append(line)
            if line.rstrip().endswith("{"):
                awaiting_brace = False
                if fn_is_errorcode:
                    # Only add if the function does not already have one.
                    nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
                    if nxt != "PetscFunctionBeginUser;":
                        out.append("  PetscFunctionBeginUser;")
                        begins += 1
                depth = 1
            continue

        if awaiting_brace and stripped == "{":
            out.append(line)
            awaiting_brace = False
            depth = 1
            if fn_is_errorcode:
                nxt = lines[i + 1].strip() if i + 1 < len(lines) else ""
                if nxt != "PetscFunctionBeginUser;":
                    out.append("  PetscFunctionBeginUser;")
                    begins += 1
            continue

        depth += line.count("{") - line.count("}")

        # Rewrite the terminal return of an error-returning function.
        rm = RETURN_ZERO.match(line)
        if rm and fn_is_errorcode:
            out.append(f"{rm.group(1)}PetscFunctionReturn(PETSC_SUCCESS);")
            returns += 1
            continue

        cm = CALL.match(line)
        if cm and "PetscCall" not in line:
            indent, name, args = cm.group(1), cm.group(2), cm.group(3)
            # PetscCall expands to a conditional `return` of a PetscErrorCode, so
            # it only compiles inside a function that returns one. Functions here
            # returning PetscScalar (the whole alpha*/beta* mass-matrix family) or
            # void must be left alone; they need a different treatment and are
            # reported instead.
            if not fn_is_errorcode:
                if name in ERROR_RETURNING:
                    skipped.append(
                        f"{path.name}:{i + 1}: {name} (in non-PetscErrorCode function)"
                    )
                out.append(line)
                continue
            if name in ERROR_RETURNING and balanced(f"({args})"):
                out.append(f"{indent}PetscCall({name}({args}));")
                wrapped += 1
                continue
            if name not in NOT_ERROR_RETURNING and name not in ERROR_RETURNING:
                # An unrecognised call: record it rather than guessing.
                if re.match(r"^(DM|Vec|Mat|TS|KSP|PC|IS|Petsc|SNES)", name):
                    skipped.append(f"{path.name}:{i + 1}: {name}")

        out.append(line)

    if not dry_run:
        path.write_text("\n".join(out) + "\n")

    return wrapped, begins, returns, skipped


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--report-skipped", action="store_true")
    args = ap.parse_args()

    total = [0, 0, 0]
    all_skipped: list[str] = []
    for f in args.files:
        w, b, r, s = transform(Path(f), args.dry_run)
        total[0] += w
        total[1] += b
        total[2] += r
        all_skipped += s
        print(f"{f}: wrapped {w}, added {b} PetscFunctionBeginUser, "
              f"converted {r} returns, skipped {len(s)}")

    print(f"\nTOTAL: wrapped {total[0]}, begins {total[1]}, returns {total[2]}")
    if all_skipped:
        print(f"unrecognised PETSc-looking calls left alone: {len(all_skipped)}")
        if args.report_skipped:
            from collections import Counter
            names = Counter(s.split(": ")[1] for s in all_skipped)
            for n, c in names.most_common():
                print(f"  {c:4d}  {n}")
    if args.dry_run:
        print("\n(dry run: nothing written)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
