# MHD regression harness

Purpose: prove that a change to the C MHD solver does not change the numerical
solution. The solver was inherited with no test coverage; every refactoring step
since has been gated by the checks described here.

Everything here lives in `tests/regression/` and touches no solver source.

## Quick start

```sh
# The inner loop: three exact-hash tests, a few seconds in total.
cd build && ctest

# One-time: provision numpy + matplotlib in a virtualenv (tooling only,
# never linked into the solver).
./tests/regression/setup_env.sh

# Record or check the production baseline (slow: see "Cost" below).
./tests/regression/run_baseline.sh --steps 3 --update    # record
./tests/regression/run_baseline.sh --steps 3 --compare   # check

# Verify the input data decodes and the two copies of the material map agree.
$HOME/.venv-hymera/bin/python tests/regression/viz/plot_materials.py
$HOME/.venv-hymera/bin/python tests/regression/viz/plot_equilibrium.py
```

## What is fingerprinted

`fingerprint.py` scrapes a run log into JSON. Every quantity it reads was
already printed by the solver; no new output code was added.

| Quantity | Source |
|---|---|
| ~88 configuration values | `AppCtxView`, `src/mhd/mfd_config.c` |
| Derived parameters (Lundquist, Reynolds, resistive time, normalized dt, mesh) | `mhd_initialize`, `src/mhd/mhd.c` |
| Per-step `‖Xⁿ⁺¹−Xⁿ‖₂` and max-norm | `Monitor`, `src/mhd/ts_functions.c` |
| Per-step `max‖div B‖` and `max‖div B‖/max‖B‖` | `Monitor` |
| Post-processed-B divergence | `Monitor` |
| Toroidal current inside / outside plasma / in vacuum vessel | `ComputeCurrent` via `Monitor` |
| `TSConvergedReason` and step count | `mhd_step` |
| SNES / KSP iteration counts | PETSc monitors, on by default |

**`MHD_Config/monitor=1` is required.** The per-step diagnostics cost a nested
linear solve, so they default to off. `run_baseline.sh` sets the flag.

Historical note: `Monitor` was registered with `TSMonitorSet` but never ran,
because PETSc only invokes monitors from `TSSolve` while `mhd_step` advances with
`TSStep`. All of the above was therefore missing from every run until it was
called explicitly.

## The sharpest single check: div B

A mimetic discretization satisfies discrete `div·curl = 0` exactly, so
`max‖div B‖/max‖B‖` sits at machine zero and stays there. Measured at the
production configuration:

```
max‖div B‖             = 2.70127e-13   (unchanged across steps)
max‖div B‖/max‖B‖      = 1.34388e-13
post-processed B ratio = 5.94839e-14
```

A discretization error that preserves the solution norms will still move this
number by orders of magnitude. Treat any growth beyond ~2x as a failure.

## Tolerance policy by transformation class

Declare the class before starting a change; do not widen a tolerance to make a
run pass.

| Class | Transformations | Tolerance |
|---|---|---|
| `EXACT` | dead-code deletion; `PetscCall` insertion; renaming; verbatim extraction preserving operation order; magic literal to named constant of identical value | 0 ULP; log values match as strings |
| `REASSOC` | changing sum order or association | `‖ΔX‖∞/‖X‖∞ ≤ 1e-13`, log scalars 1e-12 relative |
| `SOLVER` | anything that can shift the Newton/Krylov path | 1e-8 relative; iteration counts may vary; `TSConvergedReason` and step count must match exactly |

The build uses `-O2 -g -DNDEBUG` with **no LTO and no `-ffast-math`** (verified
in `build/compile_commands.json`). Two consequences:

1. Deleting an `extern` function cannot change the codegen of any surviving
   function, so pure deletion is provable by comparing object code — no run
   required. Assert these flags stay off; the proof is void otherwise.
2. There is no x87 excess precision and no compiler reassociation, so hoisting a
   value recomputed from identical arguments is bit-identical, i.e. `EXACT` and
   not `REASSOC`. In particular, caching a `DMStagGetLocationSlot` result is
   `EXACT`: it is a pure function of the DM.

### Solver-tolerance amplification

`-ts_adapt_type none` fixes the step sequence, so `tⁿ` is reproducible. But the
state at `tⁿ` is only converged to `snes_rtol`/`ksp_rtol`, so a rounding-level
change in the residual can flip a Krylov iteration count and move the state by
`O(ksp_rtol·‖X‖)` — potentially eight orders of amplification in a single step.

Therefore: **validate a `REASSOC` change with a single residual evaluation, not a
multi-step run.** In a one-shot residual evaluation a 1e-16 perturbation stays
1e-16 and is unambiguous.

Since the command-line precedence fix, regression runs can also tighten
`-snes_rtol`/`-ksp_rtol` well below production values to shrink this
amplification. That is a test configuration, not a production change.

## Determinism

Verified: two independent runs of the production configuration at 100x2x200 gave
bit-identical results on every shared step, on all 88 configuration values, and
on every derived parameter. MUMPS is deterministic at fixed rank count, so a
future difference between two identical runs indicates an uninitialized read
rather than expected noise.

## Cost

The production path is expensive. At 100x2x200, one `mhd_step` plus the
initial-condition relaxation is roughly 20 minutes on one core. Almost all of it
is the relaxation — a complete extra `TSSolve` with MUMPS LU — and the implicit
solve itself.

So the full baseline is a phase gate, not an inner-loop check. Newton converges
cleanly when the problem is well posed: `6.91 → 0.13 → 2.9e-3`, with KSP
converging in one iteration per Newton step.

## Input data

`gen_grid_data.py` generates the grid-data files at any resolution, since the
filenames encode the grid (`veceta_grid<NR>x<Nphi>x<NZ>.txt`) and only the
production size exists in-tree.

Verified layouts, all exact functions of the grid:

| File | Field | Shape | At 100x2x200 |
|---|---|---|---|
| `veceta` | `dataC`, cell-centred material tags | `Nr × Nphi × Nz` | 40000 |
| `vecpsi` | `datapsi`, poloidal flux, r-z edge | `(Nr+1) × Nphi × (Nz+1)` | 40602 |
| `vecg` | `datag` = `R·B_φ`, phi-face | `Nr × (Nphi+1) × Nz` | 60000 |

Two cautions:

- `veceta` is bit-identical to `inputs/AxisSymmetricGeometry.dat` transposed
  (zero mismatches over all 20000 cells). The same data is read independently by
  C and by C++, and nothing checks the copies agree. Generate both from one
  source.
- `vecpsi`/`vecg` can be resampled to a coarser grid dimensionally, but the
  result is **not** a discrete Grad-Shafranov solution on the new mesh, and the
  hardcoded cell indices in `betaephi_isolcell` address absolute positions. A
  coarse grid is therefore a different model, not a coarsened one. Use it as a
  crash/NaN canary with its own baseline, never to certify physics.

`datapsi` and `datag` are read only on the `ictype 9/15` path.

## A test that cannot fail is worse than no test

Before trusting a baseline, confirm the comparator detects a real difference.
This has been exercised: comparing a 1-step against a 3-step fingerprint
correctly failed on the differing entry counts. Re-check after any change to the
comparator.

This repository already contains the failure mode being guarded against:
`tests/kinetic/avalanche_test.cpp` is 692 lines that return 0 unconditionally.

## Tiers

| Tier | Test | What it hashes | Time |
|---|---|---|---|
| 0 | `mhd_t0_coefficients` | every live mass-matrix coefficient over all indices, stencil locations and materials, 9.8M values | ~4 s |
| 0 | `mhd_t0_operators` | 13 live mimetic operators applied to a fixed input, 6.3M values | ~2 s |
| 1 | `mhd_t1_residuals` | every live residual and the RHS evaluated once on a fixed state, plus the production residual with `inertia` on, 2.9M values | ~3 s |
| — | `run_baseline.sh` | a full production run, 100×2×200 | ~40 min |

Each tier-0/1 test prints every value as an exact hexadecimal double and compares a
SHA-256 against `tests/regression/baselines/<name>.sha256`. A hash mismatch says only
that something changed; to see what, run `run_t0.sh <binary> --baseline <name>
--keep FILE` on both revisions and diff the files. Every line names the function and
the index, so the diff localizes the change exactly.

Why tier 1 matters: a full time step cannot check a residual bit-exactly, because its
output feeds an inexact Newton-Krylov solve that can amplify a one-ULP change into one
at solver tolerance. A single evaluation has no such amplification. Use the full run
for what only it covers: the composition of everything, the solver stack, restart I/O.

**Re-recording a baseline** is legitimate only when the source tree equals HEAD, so
that the baseline describes committed code. The tests' fixed state is synthetic: a
50×2×186 grid with a tokamak-like nested-shell material layout (all five materials,
four-plasma-neighbour vertices, and the three hardcoded isolated cells), and
asymmetric nowhere-zero fields.

## Tools

| Tool | Use |
|---|---|
| `run_t0.sh` | run a tier-0/1 binary, check or record its hash |
| `run_baseline.sh`, `fingerprint.py` | full production run and its log fingerprint |
| `codegen_identity.sh`, `compare_asm.py` | compare generated code per function before/after; complete proof for pure deletions |
| `decode_diff.py` | decide whether differing AArch64 instructions are only padding and relocated branches |
| `verify_move.py` | prove a function move was byte-exact and that per-file counts reconcile |
| `check_line_numbers.py` | prove every surviving `PetscCall`/`SETERRQ` keeps its `__LINE__` |
| `find_dead_branches.py` | find, and optionally remove, constant-false `if` blocks |
| `add_petsc_checks.py` | wrap bare PETSc calls in `PetscCall`, allow-list based |
| `gen_grid_data.py` | generate grid data files at any resolution |
| `setup_env.sh` | provision numpy/matplotlib in a virtualenv |
| `viz/` | plot the material map and equilibrium inputs |

The Python virtualenv lives in `~/.venv-hymera` and may need recreating with
`setup_env.sh` after the environment is reset.
