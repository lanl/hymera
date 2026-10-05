# The C MHD solver

`src/mhd` is the field solver of the Hymera hybrid code: a PETSc `TS` time
integrator over mimetic finite-difference operators on an axisymmetric (r, φ, z)
`DMStag` grid. It is written in C and coupled across a C ABI to the C++/Parthenon
kinetic side in `src/kinetic` and `src/tasks`.

It was inherited and its original author is no longer available. It has since been
refactored from 49,326 lines to 16,886, with every step shown not to change the
answer. [refactor/changelog.md](refactor/changelog.md) records how.

## Reading order

New to this code? Read in this order.

1. **[physics/discretization.md](physics/discretization.md)** — the staggered grid
   and what lives on each stratum. It decodes the opaque slot-index names
   (`ivVrmphimzm`, `ivErmzm`, ...) that the residual assembly is written in.
2. **[physics/materials.md](physics/materials.md)** — the material region model,
   and what the test `|dataC - 1.5| < 0.7` means.
3. **[reference/user-struct.md](reference/user-struct.md)** — the `User` struct is
   the solver's entire mutable state: which fields are inputs, which are derived,
   which are never written.
4. **[reference/function-index.md](reference/function-index.md)** — every function,
   its line range, and what calls it. Generated; regenerate after changing the
   source.
5. **[physics/residual-differences.md](physics/residual-differences.md)** — how the
   production and relaxation residuals differ, and the switches that now express it.
6. **[physics/open-questions.md](physics/open-questions.md)** — physics questions
   only the code owner can answer, and the decisions already taken.

Before making a change:

7. **[testing/bit-exactness-policy.md](testing/bit-exactness-policy.md)** — the
   normative rules: what must stay bit-exact, how to prove it, and the traps that
   make naive edits change results.
8. **[testing/harness.md](testing/harness.md)** — how to run the regression tests.

## Orientation

**The public API** is declared in `src/mhd/mhd.h`: `mhd_PetscInit`,
`mhd_initialize`, `mhd_step`, `mhd_getF`, `mhd_resetState`, `mhd_destroy`, and
two save/load pairs, `mhd_savesolution`/`mhd_loadsolution` (PETSc binary) and
`mhd_save_hdf5`/`mhd_load_hdf5` (HDF5). Every function returns a PETSc error code;
C++ callers wrap each call in `MHD_CHECK`, which reports and aborts on failure.

**Structure of the source:**

| File | Lines | Responsibility |
|---|---|---|
| `ts_functions.c` | 5,797 | residuals, initial conditions, `Monitor`, shell preconditioner, rank-independent vector I/O |
| `geometry.c` | 4,301 | grid metrics, inter-stratum projection and reconstruction, array extraction for the C++ side |
| `mimetic_operators.c` | 2,339 | discrete curl, divergence, gradient, vector Laplacian |
| `monitor_functions.c` | 1,894 | diagnostics: cell averaging, VTK writers, toroidal current |
| `mass_matrix_coefficients.c` | 1,271 | material-weighted mass-matrix coefficients |
| `mhd.c` | 991 | the API, object construction, solver setup |
| `default_petsc_options.c` | 150 | the PETSc solver stack, in commented groups |
| `mfd_config.c/.h` | ~350 | the `User` struct and its dump |

`src/mhd/attic/` holds about 28,000 lines of unreachable code, moved out rather than
deleted so the physics it encodes stays searchable. It is not compiled. Its
[README](../../src/mhd/attic/README.md) says what each function was.

**A step** is `mhd_step`: one `TSStep` on the implicit system. The residual is
`FormIFunction_Vperp_viscosity`. The initial condition is `FormInitialSolution_psi`,
which builds a tokamak equilibrium from EFIT / Grad-Shafranov data in `inputs/mhd`,
relaxes it (`FormIFunction_newequilibrium_Vperp`, ideal — no resistive term), then
solves for EP and tau (`FormIFunction_InitializeEP_halo`).

The production and relaxation residuals share one body, `vperp_residual`, and differ
by three switches; see [physics/residual-differences.md](physics/residual-differences.md).

**Configuration** read from the `<MHD_Config>` deck section that is specific to this
refactor:

| Key | Default | Effect |
|---|---|---|
| `monitor` | 0 | run the per-step diagnostics (step norms, div B, currents); costs a nested solve per step |
| `inertia` | 0 | include the advective inertia term n_i (V·∇)V in the production residual |
| `ictype` | 9 | initial-condition selector; 9 is the production EFIT equilibrium |

PETSc solver options come from `default_petsc_options.c` and can be overridden on
the command line, which takes precedence.

## Things that will surprise you

- **The grid is effectively pinned to 100×2×200.** Input filenames encode the grid,
  only that resolution exists in-tree, there is no generator for the equilibrium
  data, and `alphaecphi_isolcell` hardcodes three absolute cell indices
  (open question Q1).
- **`dataC` is not a resistivity field**, despite living in a file called `veceta`.
  It is a five-valued material tag map, bit-identical to
  `inputs/AxisSymmetricGeometry.dat` — the same data in two formats, read separately
  by C and by C++, with nothing checking they agree.
- **Some results depend on incidental code shape.** The build contracts to fused
  multiply-add, so rewriting an expression — even deleting a term multiplied by
  zero — can move results in the last bit. Several places look improvable and are
  deliberately left as they are; each is commented. Read the bit-exactness policy
  before simplifying arithmetic.
- **Zero global mutable state.** Everything is reachable from one `User*`. The
  struct is `calloc`'d and zero-initialization is load-bearing for fields nobody
  assigns, so do not split or reorder it.
- **Line numbers in compiler and debugger output differ from the file.** The sources
  contain `#line` directives that keep `__LINE__` at its pre-refactor values,
  because `PetscCall` bakes it into the generated code. Use function names, not line
  numbers, to find things.

Fixed during the refactor, recorded so earlier behaviour is explained:

- `Monitor` was registered but never ran, because PETSc only fires monitors from
  `TSSolve` and `mhd_step` uses `TSStep`. It is now called explicitly, behind
  `MHD_Config/monitor`.
- PETSc options in the source overrode the command line. The command line now wins.
- Return codes were almost all discarded, so a missing input file segfaulted
  thousands of lines later. All PETSc calls are now checked and failures name the
  cause.
- Most of the code was unreachable: one of eighteen residual variants was live.
- `FormRHSFunction_BImplicit`, identically zero, was registered on the implicit path
  and built auxiliary fields every Newton iteration only to add zero.

## Building, running, testing

See the repository `CLAUDE.md` for the Spack/CMake recipe. In short:

```sh
cmake -B build -DENABLE_CUDA=OFF        # the MHD path is CPU-only regardless
cmake --build build -j $(nproc)
cd build && ctest                       # three regression tests, a few seconds
cd bin && ./mhd -i ../../inputs/mhd.input
```

Run `mhd` from `build/bin`: the input decks use relative paths
(`input_folder=../../inputs/mhd`) that only resolve from there. From anywhere else
the grid data is not found, and the run now stops with a message naming the file.
