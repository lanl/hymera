# Refactor changelog

A phase-by-phase record of the refactor of the C MHD solver, `src/mhd`. Each
entry lists its commits, what changed, and how it was shown not to change the
answer. Commit messages carry the full detail; this is the map.

The governing constraint throughout, set by the code owner: **results must stay
consistent with previous runs.** Every change below is either bit-identical to
what came before, or is flagged here as a deliberate, owner-approved exception.

Size of the compiled solver, `src/mhd/*.c`:

| | Lines |
|---|---|
| Before | 49,326 |
| Now | 16,886 (−66%) |
| Quarantined, not compiled (`src/mhd/attic/`) | 28,655 |

## P0 — Regression harness

`833eff2`, `af140f1`

Built before touching any solver code, so that every later step has an objective
check.

- `fingerprint.py` extracts the per-step diagnostics from a run log; `run_baseline.sh`
  runs the production configuration (100×2×200) and compares against a recorded
  baseline. Determinism confirmed: two independent runs agree bit for bit.
- `codegen_identity.sh` compares generated machine code per function before and after
  a change — a complete proof for deletions, needing no run.
- Fixed on the way, because the harness needed them: `Monitor` was registered but
  never ran (PETSc only fires monitors from `TSSolve`, and `mhd_step` uses `TSStep`),
  so no run had ever printed its diagnostics; PETSc options were inserted after the
  command line was consumed, so command-line overrides silently lost; a duplicated
  include guard voided `default_petsc_options.h`; `numC` was discarded; two `%D`
  format specifiers removed from PETSc 3.16.

## P1 — Documentation as discovery

`158e0a3`

The physics existed only in the source. Documented the staggered-grid slot names,
the material model (`dataC` is a five-valued tag field; `|dataC − 1.5| < 0.7` means
"plasma or separatrix-wall"), every field of the `User` struct, and the open physics
questions only the owner can answer.

## P2a — Quarantine of unreachable code

`0823f34`, `d1b2121`, `a408b42`, `c335365`, later `5f8dd63`, `7e27400`

Moved 80+ functions with no callers into `src/mhd/attic/`, which is not compiled, so
the physics they encode stays greppable. Each move proved three ways: byte identity
of every moved and surviving body; identical generated code for every surviving
function; and an unchanged production run.

Two traps found and now documented in the bit-exactness policy:
- PETSc's `PetscCall`/`SETERRQ` expand `__LINE__` into generated code, so deleting
  lines changes functions nobody touched. Every removed block is replaced by a
  `#line` directive.
- Removing a function shifts alignment padding and relabels `.rodata` references in
  `objdump` output. Both look like code changes and are not; `decode_diff.py`
  decodes the instructions to tell them apart.

Reachability is computed transitively: later passes found functions whose only
callers had themselves been quarantined.

Kept by owner decision, not quarantined: `mhd_save_hdf5` / `mhd_load_hdf5`, an
alternative save/load path (`62e1c9b`).

## P3 — Error checking

`76bfc30`

3,816 PETSc calls wrapped in `PetscCall`; 64 `PetscFunctionBeginUser`/`Return`
added; all 16 C++ call sites into `mhd.h` guarded by `MHD_CHECK`, whose handler
emits PETSc's traceback then aborts. A missing input file used to segfault
thousands of lines later; it now names the file. Nothing was failing silently on the
success path: the production run was unchanged.

## P5 — Coefficient and curl deduplication

`f4a442d`, `fee16bd`, `e73791e`, `5f8dd63`, `8234787`

Gated by two tier-0 tests that hash every coefficient and operator over a
9.8-million-point index space in seconds.

| Shared body | Replaces |
|---|---|
| `betae_sum` | 6 `betae*` copies |
| `alphaec_sum` | 5 `alphaec*` copies |
| `edge_average` | `rese`, `condu` |
| `derived_curl` | `FormDerivedCurl` / `_nores` / `_nomp` |

Each public name is a one-line wrapper, so no caller changed.
`mass_matrix_coefficients.c`: 4,231 → 1,271 lines.

## P6 — Residual deduplication and flags

`668e2e2`, `f20302c`, `c2c454f`, `a404205`, `20969e5`, `01778d4`, `9a658e8`,
`6ce3379`, `7973092`

Gated by the tier-1 test, which evaluates every live residual once on a fixed state.
A full time step cannot check the residual bit-exactly: an inexact Newton–Krylov
solve can amplify a one-ULP change.

- **Shared skeleton:** the ~54 slot lookups, 12 edge lengths and cell volume, each
  duplicated across the residuals, now come from `MFD_GetSlots*`/`MFD_UNPACK_SLOTS`,
  `MFD_CellEdgeLengths` and `MFD_CellVolume`.
- **`initialize_ep`** merges `InitializeEP` and `_halo`, which differed only in two
  edge coefficients.
- **`vperp_residual`** merges the production and relaxation residuals behind three
  switches: `f1_inertia`, `f3_resistive`, `f3_jre`.
- **16 never-true `if (0)` blocks deleted**, with the locals they orphaned.
- **The identically-zero `FormRHSFunction_BImplicit`** is no longer registered on
  the implicit path.
- **`MHD_Config/inertia`** (default 0) controls the production advective inertia
  term, which had been multiplied by a literal zero.

Owner decisions recorded here:
- The relaxation is **ideal**, with no resistive term, for consistency with previous
  runs. A probe confirms this was already its behaviour: scaling every resistivity by
  1000 changes none of its outputs.
- `user->jre` and `view3d_t` are kept for the hybrid coupling, and `jre` is its own
  switch.

Floating-point lessons, all because the build contracts to fused multiply-add:
- switch at statement level, never by masking a term inside an expression;
- a shared body taking its switches at run time differed in the last bit, so it is
  `always_inline` and each caller passes constant switches;
- the inertia flag is resolved in the wrapper, which selects between two constant
  specializations;
- deleting the `0.0 *` inertia term changes 30 entries by 1–14 ULP, so with the flag
  off the original statement is kept verbatim.

## Tests now registered with ctest

| Test | Tier | Covers |
|---|---|---|
| `mhd_t0_coefficients` | 0 | every live mass-matrix coefficient |
| `mhd_t0_operators` | 0 | 13 live mimetic operators |
| `mhd_t1_residuals` | 1 | all live residuals, the RHS, and the inertia-on path |

Every harness was negative-tested — shown to fail on a deliberate perturbation —
before being trusted. The tier-0 material layout had to be redesigned for that
reason: the first version had no vertex with four plasma neighbours, so a whole
branch was never evaluated.

## Not yet done

- `EdgeToCellReconstruction_r/_phi/_z` in `geometry.c`: three near-identical copies,
  not yet covered by the operator test.
- With `MHD_Config/inertia=0`, production still builds `niv` and the velocity
  gradient for the multiply-by-zero term. Removing that saves 16 operator calls per
  residual evaluation but moves results by ~1e-15; it waits for a deliberate
  re-baselining.
- Open physics questions remain in [../physics/open-questions.md](../physics/open-questions.md),
  notably Q1, the three hardcoded "isolated" cells.
