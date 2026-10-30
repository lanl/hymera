# Expanded-domain ITER 15% taper equilibrium evaluator

Frozen model: expanded-domain PINN job 47892285, step 10000. This is inference
for one learned equilibrium, not a new equilibrium solve or a current-conditioned
model. See `model_manifest.json` for checkpoint and executed-source hashes.

## Build and run

Only a C++17 compiler and the standard library are required. Python, JAX, PETSc,
Kokkos and model-loading libraries are not required for inference. CPU host API;
this implementation is not a Kokkos device kernel.

```sh
cmake -S . -B /tmp/iter-pinn-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/iter-pinn-build -j 2
ctest --test-dir /tmp/iter-pinn-build --output-on-failure
/tmp/iter-pinn-build/iter_pinn_evaluate grid /tmp/iter-pinn-grid 100 2 200
```

The output directory must not already exist. The default bounds are normalized
R=[1.525,4.975], Z=[-2.975,2.975], L0=2 m. Custom bounds and length scale can be
appended: `NR NPHI NZ RMIN RMAX ZMIN ZMAX L0_M`. All physical points must lie in
the trained rectangle R=[0.5,14], Z=[-9,9] m. Coordinates outside it are rejected.

For arbitrary point evaluation, write whitespace-separated `R_m Z_m` pairs to
`coordinates.txt`, then run:

```sh
/tmp/iter-pinn-build/iter_pinn_evaluate points coordinates.txt values.txt
```

The output columns are R_m, Z_m, native psi, native g, and Hymera psi. No existing
output file is overwritten.

## Calling from C++

Link `equilibrium.cpp` into the calling executable, or add this directory with
CMake `add_subdirectory` and link the `iter_pinn_equilibrium` target.

```cpp
#include "equilibrium.hpp"

auto native = iter_pinn::evaluate(6.2, 0.6); // physical R,Z in metres
// native.psi: poloidal flux / (2*pi), Wb/rad (numerically V s)
// native.g: R*B_phi, T m

auto hymera = iter_pinn::evaluate_hymera(3.1, 0.3, 2.0);
// normalized input coordinates, dimensional output fields for Hymera

iter_pinn::Grid grid{100, 2, 200, 1.525, 4.975, -2.975, 2.975, 2.0};
auto fields = iter_pinn::evaluate_grid(grid);
// fields.psi and fields.g can supply the existing datapsi/datag arrays.
```

This directory provides a callable library and grid exporter. Hymera's production
initialization and input files are unchanged. For a subsequent integration,
select the PINN before the psi/g file-loading branch of `mhd_initialize`.
Retain the existing magnetic reconstruction and material tags.

## Conventions and grid layout

Native PINN fields obey B_R=-psi_Z/R, B_Z=psi_R/R, B_phi=g/R. Hymera uses the
opposite poloidal-flux convention in `FormInitialSolution_psi`. The adapter returns
`psi_hymera=-psi_native+psi_offset`, preserving physical poloidal B. It returns
`g_native` unchanged by default, preserving physical toroidal B. An explicit
`g_sign=-1` reverses the physical toroidal field; this is not a gauge change.

Both interfaces return dimensional fields. Do not divide psi or g by B0 or L0
again: Hymera already forms B in physical units and then divides B by B0.
The model's psi is flux/(2*pi), not the toroidally integrated poloidal flux.

The ASCII grid files have a count header followed by one double per line:

- psi: (NR+1)*NPHI*(NZ+1), R/Z vertices and phi cell centers.
  Offset `i+j*(NR+1)+k*NPHI*(NR+1)`.
- g: NR*(NPHI+1)*NZ, R/Z cell centers and phi faces.
  Offset `i+j*NR+k*(NPHI+1)*NR`.

Every toroidal plane repeats the same axisymmetric equilibrium, including the
last g face. The supplied historical g file had a background-only last plane;
this exporter deliberately evaluates that face consistently with all others.
`veceta` is a separate material map and is not regenerated or changed.

**Initialization caveat:** Hymera can overwrite these fields with a binary
restart (`ic_binary_load=1`), or relax the imported equilibrium when constructing
a new state. Evaluator parity does not establish that a subsequently relaxed
Hymera state is identical. No MHD evolution is run by this tool.

## What is evaluated

- Four hidden tanh layers of width 32 and a linear output, with the archived
  Fortran-order weight packing, original coordinate map, Gaussian seed and
  20 V s flux scale.
- All eleven finite rectangular coils at their learned currents: the same
  10x10 far quadrature and four-triangle 32x32 Duffy quadrature near/inside coils.
  The Green function uses the executed JAX evaluator's eight AGM iterations.
  This is the saved numerical model's coil quadrature, not a new exact coil integral.
- Analytic 15% tapered Taylor G profile with the saved axis/X-point fluxes,
  alpha, connected-support mask (80x128) and wall polygon. G=32.86 T m outside
  that support. The saved selector is reproduced rather than inferred from psi
  alone, which would incorrectly activate disconnected/external flux regions.

The neural self-flux is smooth. The finite-coil source has sharp boundaries,
and G uses the archived discrete support selector; the full evaluator does
not promise global C-infinity smoothness or differentiability of its branch choices.
It exposes field values, not automatic derivatives.

## Validation and provenance

`test_evaluator.cpp` checks unit/sign conversion, invalid coordinates, domain
corners, and every entry of a small staggered grid. `validate.py` compares the
C++ result to the original frozen Python/JAX source, including coil interiors,
coil boundaries, quadrature transitions, source-mask transitions, and the cached
fields used in the figures. It also exports and checks a full 100x2x200 grid.

`validate.py` is a developer check and requires the qsmhd-pinn workspace and its
Python dependencies; it is not an inference dependency. Run it on HiPerGator:

```sh
python validate.py --workspace /blue/cmcdevitt/j.arnaud/qsmhd-pinn \
  --executable /absolute/path/to/iter_pinn_evaluate --output /new/validation/path
```

`export_model.py` is maintained in qsmhd-pinn under
`free_boundary_equilibrium/exports/iter_taper15_cpp_20261027/`; running it there
regenerates `model_data.hpp` and `model_manifest.json` from the archived checkpoint.
Do not hand-edit those generated files.
