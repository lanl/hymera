# MHD regression harness

Purpose: prove that a change to the C MHD solver does not change the numerical
solution. The solver is 49,326 lines of inherited C with no test coverage, so
every refactoring step needs an objective before/after check.

Everything here lives in `tests/regression/` and touches no solver source.

## Quick start

```sh
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

## Known gap

`tests/CMakeLists.txt` still has its test registrations commented out, so `ctest`
reports zero tests. The harness currently runs via the shell scripts above.
Wiring it into `ctest` needs a helper that links `mhd_core` — the existing
`add_parthenon_test` links only Parthenon/HDF5/hflux and so cannot reach the MHD
solver; copy the linking pattern from `src/CMakeLists.txt` instead.
