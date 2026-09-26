# The C MHD solver

`src/mhd` is the field solver of the Hymera hybrid code: a PETSc `TS` time
integrator over mimetic finite-difference operators on an axisymmetric (r, φ, z)
`DMStag` grid. It is written in C and coupled across a C ABI to the C++/Parthenon
kinetic side in `src/kinetic` and `src/tasks`.

It is 49,326 lines, inherited, and its original author is no longer available.
This documentation exists because most of what the code knows about the physics
is currently recorded nowhere else.

## Reading order

New to this code? Read in this order.

1. **[physics/discretization.md](physics/discretization.md)** — the staggered grid
   and what lives on each stratum. This is the decoder for the ~70 opaque index
   variables (`ivVrmphimzm`, `ivErmzm`, ...) that the residual assembly is written
   in. Without it, `ts_functions.c` is unreadable.
2. **[physics/materials.md](physics/materials.md)** — the material region model,
   and what the `|dataC - 1.5| < 0.7` test scattered across 112 sites actually
   means.
3. **[reference/user-struct.md](reference/user-struct.md)** — the `User` struct is
   the solver's entire mutable state. Which fields are inputs, which are derived,
   which are never written at all.
4. **[reference/function-index.md](reference/function-index.md)** — every function,
   its line range, and whether anything calls it. Generated; regenerate after
   changing the source.
5. **[physics/open-questions.md](physics/open-questions.md)** — questions only the
   code owner can answer. Read before changing anything in the initial-condition
   path.

Before making a change:

6. **[testing/bit-exactness-policy.md](testing/bit-exactness-policy.md)** — the
   normative rules. What must be bit-exact, what may not be, and how to prove it.
7. **[testing/harness.md](testing/harness.md)** — how to run the regression
   harness and what it checks.

## Orientation

**The public API** is eight functions declared in `src/mhd/mhd.h`:
`mhd_PetscInit`, `mhd_initialize`, `mhd_step`, `mhd_getF`, `mhd_resetState`,
`mhd_destroy`, `mhd_savesolution`, `mhd_loadsolution`. Everything the C++ side
does goes through these.

**Structure of the source**, by size:

| File | Lines | Responsibility |
|---|---|---|
| `ts_functions.c` | 22,656 | residual assembly (`FormIFunction*`), initial conditions, `Monitor`, shell preconditioners |
| `geometry.c` | 10,111 | grid metrics, inter-stratum projection/reconstruction operators, boundary index sets, array extraction for the C++ side |
| `mimetic_operators.c` | 6,985 | the discrete operators: curl, divergence, gradient, vector Laplacian, `Delta*` |
| `mass_matrix_coefficients.c` | 4,231 | diagonal mass-matrix coefficients (`alpha*` per cell, `beta*` summed per edge/face/vertex), material-weighted |
| `monitor_functions.c` | 3,352 | diagnostic output: cell averaging, VTK writers, toroidal current |
| `mhd.c` | ~1,000 | the API, object construction, solver setup, save/load |
| `mfd_config.c/.h` | ~340 | the `User` struct and its dump |
| `default_petsc_options.c` | ~150 | the PETSc solver stack configuration |

**A step** is `mhd_step`: one `TSStep` on the implicit system, whose residual is
`FormIFunction_Vperp_viscosity` — the single production residual, registered in
`mhd_initialize`. The initial condition is `FormInitialSolution_psi`, which builds
a tokamak equilibrium from EFIT/Grad-Shafranov data in `inputs/mhd` and relaxes
it with its own nested solve.

## Things that will surprise you

Each of these is verified, and each has cost someone time already.

- **`Monitor` did not run.** It was registered with `TSMonitorSet`, but PETSc
  only invokes monitors from `TSSolve` and `mhd_step` advances with `TSStep`. The
  step norms, divergence-of-B measure and current diagnostics were absent from
  every run. It is now called explicitly, behind `MHD_Config/monitor`
  (default off, because it costs a nested solve per step).
- **The PETSc options in `default_petsc_options.c` used to override your command
  line.** They were inserted after `PetscInitialize` had already consumed `argv`.
  Nothing in the solver stack could be changed without recompiling. The command
  line is now re-applied afterwards and wins.
- **Most of `ts_functions.c` is unreachable.** One of eighteen `FormIFunction_*`
  variants is registered; the rest — roughly 13,200 lines — have no callers. They
  are near-duplicates differing by a handful of tokens. See the function index.
- **Return codes are mostly discarded.** Roughly 91% of PETSc calls in
  `ts_functions.c` ignore their result, and that file contains no
  `PetscFunctionBeginUser`, so error stack traces stop at the caller. This is not
  theoretical: a missing input file produces `SETERRQ`, the caller ignores it, and
  the run continues on garbage or segfaults later.
- **`dataC` is not a resistivity field**, despite living in a file called
  `veceta`. It is a five-valued material tag map, and it is bit-identical to
  `inputs/AxisSymmetricGeometry.dat` transposed — the same data, in two formats,
  read separately by C and by C++, with nothing checking they agree.
- **The grid is effectively pinned to 100x2x200.** The input filenames encode the
  grid, only that resolution exists in-tree, there is no generator for the
  equilibrium data, and `betaephi_isolcell` hardcodes three absolute cell indices.
- **Zero global mutable state.** Verified by `nm`: no `.data`, no `.bss`, no
  file-scope statics, no static functions in 49k lines. Everything is reachable
  from one `User*`. This is the single best thing about the codebase and the
  reason a refactor is tractable at all; `calloc` zero-initialization is
  load-bearing, so do not split or reorder that struct.

## Building and running

See the repository `CLAUDE.md` for the Spack/CMake recipe. In short:

```sh
cmake -B build -DENABLE_CUDA=OFF        # the MHD path is CPU-only regardless
cmake --build build -j $(nproc)
cd build/bin && ./mhd -i ../../inputs/mhd.input
```

Run from `build/bin`: the input decks use relative paths
(`input_folder=../../inputs/mhd`) that only resolve from there. Running from the
repository root fails to open the grid data, and because the error is ignored, the
failure appears later as a segfault rather than a message.
