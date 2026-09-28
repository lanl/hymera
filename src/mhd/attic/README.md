# attic/

This directory holds source that is **not built** and **not part of any CMake
target**. `ts_functions_attic.c` contains 22 functions moved verbatim (byte-for-byte
identical bodies) out of `../ts_functions.c`. Each had zero callers anywhere in
`src/`, `tests/`, or the built objects (verified with `nm`) as of commit
`0823f34`. They were moved rather than deleted so the physics they encode stays
in-tree and greppable for review before a later, separate deletion commit.

These functions are retained for physics review only. They contain known
defects and must **not** be re-enabled (re-linked, re-wired to a caller, or
added to a CMake target) without first regenerating the regression baselines
under `tests/regression/`. See `docs/mhd/physics/open-questions.md` for the
open physics questions that motivated leaving this code unreachable rather
than fixing it in place.

## What each function was

All of these are alternate `IFunction`/`FormInitialSolution`/`SampleShellPC*`
variants that were superseded by the live code path. Where useful, the
description below is relative to the residual the live driver actually calls,
`FormIFunction_Vperp_viscosity` (constraint: `[(curl B) x B = -lambda*Lap(V)].e_{R/Z}`,
`V.B = 0`, with the `eta_phi`/`eta_perp` split resistive term in the tau/EP
constraints).

- **FormInitialSolution_LargeData** -- initial-solution setup (`X_0`) reading
  B-field/levelset data via `ReadALine` line-by-line instead of the bulk-array
  `ReadInitialData` path used by the live `FormInitialSolution`.
- **FormIFunction_Vperp_viscosity_halo** -- same eta_phi/eta_perp resistive
  split as the live residual, but only applied inside the halo/wall region.
- **FormIFunction_Vperp_viscosity_halo_isolcell** -- as above, but with
  eta_phi additionally lowered on three isolated cells.
- **FormIFunction_Inertia** -- adds ion inertia `n_i m_i (V.grad)V` to the
  momentum constraint (in place of the live residual's viscosity-only RHS),
  still with `V.B = 0`.
- **FormIFunction_Inertia_viscosity** -- same inertia term as above, plus the
  `-lambda*Lap(V)` viscosity term together (no `V.B=0` projection split).
- **FormIFunction_Inertia_V_ni** -- inertia term `n_i m_i (V.grad)V` with no
  viscosity and no `V.B=0` constraint at all.
- **FormIFunction_Inertia_V** -- inertia term using a frozen `n_i(t=0)`
  instead of the evolving density, with `[V.grad V].B = 0`.
- **FormIFunction_DampingV** -- replaces the momentum constraint with an
  artificial damping term `dampV*V - (curl B) x B = 0`.
- **FormInitialpsi** -- sets up an initial poloidal flux `psi` forced to be
  nonzero only inside the vacuum vessel.
- **FormInitialSolution_psi_fromNphi2** -- initial-solution setup from
  `psi`/`G(psi)` data generated on a coarse 2-cell-in-phi mesh, interpolated
  onto the full `Nr x Nphi x Nz` mesh.
- **FormIFunction_BImplicit2** -- vertex-centered variant of the momentum
  equation, `n_i dV/dt - curl(B) x B`, as an implicit left-hand side.
- **FormIFunction_newequilibrium** -- momentum residual
  `(curl B) x B = (1/Re)*Lap(V) + n_i*dV/dt` used to relax toward a new
  equilibrium, no `V.B=0` split.
- **FormIFunction_InitializeEPV** -- initializer for `V`, `EP`, and `tau`
  together, solving `[(curl B0) x B0 = -lambda*Lap(V0)].e_{R/Z}`, `V0.B0=0`.
- **FormIFunction2** -- treats `tau` and `Phi` as auxiliary variables with
  `tau = -V x B + (eta/mu_0)(curl B)` (resistive Ohm's law), `Phi = 0`.
- **FormIFunction** -- treats `tau` and `Phi` as auxiliary variables with both
  set to zero (ideal, no-resistivity variant of FormIFunction2).
- **SampleShellPCSetUp_ApproximateDiag** -- shell-preconditioner setup using
  an approximate diagonal in place of an exact one.
- **SampleShellPCSetUp_SuperLU** -- shell-preconditioner setup using a
  SuperLU direct factorization instead of the live field-split approach.
- **SampleShellPCSetUp_Diag** -- shell-preconditioner setup using an exact
  diagonal approximation.
- **SampleShellPCApply_ApproximateDiag** -- apply step paired with
  `SampleShellPCSetUp_ApproximateDiag`.
- **SampleShellPCApply_Diag** -- apply step paired with
  `SampleShellPCSetUp_Diag`.
- **SampleShellPCDestroy_ApproximateDiag** -- destroy/cleanup paired with
  `SampleShellPCSetUp_ApproximateDiag`.
- **SampleShellPCDestroy_Diag** -- destroy/cleanup paired with
  `SampleShellPCSetUp_Diag`.

## geometry.c

`geometry_attic.c` contains 19 functions moved verbatim (byte-for-byte
identical bodies) out of `../geometry.c`. The first 17 had zero callers
anywhere in `src/`, `tests/`, or the built objects (verified transitively --
several were only reachable from other already-dead functions), confirmed
with the same codegen-identity check used for `ts_functions.c`.

The remaining two, `ComputeIsEBoundary` and `VertexToEdgeReconstruction_scalar`,
were held back at that time because each still had a caller compiled into
`mhd_core` (`ComputeIsEBoundary` was called from `FormDerivedGradDivergence`
and `FormDerivedGradEtaDivergence` in `mimetic_operators.c`;
`VertexToEdgeReconstruction_scalar` was called from geometry.c itself, inside
`ComputeIsEBoundary`). The later quarantine pass on `mimetic_operators.c`
(see below) removed both callers as unreachable, so this commit re-verified
(via `nm -u`) that neither function has any remaining caller and appended
them to this file.

- **ComputeIsEBoundary** -- builds an `IS` of edge-located boundary degrees of
  freedom, via a dummy Jacobian/PC/KSP probe (the edge-DOF analogue of the
  already-quarantined `ComputeIsBBoundary`/`ComputeIsCBoundary`).
- **VertexToEdgeReconstruction_scalar** -- scalar-field variant of the live
  `VertexToEdgeReconstruction` (vertex-to-edge reconstruction operator).

## What each function was

Most of these are inter-stratum projection/reconstruction operators for the
mimetic staggered grid (moving field data between vertex/edge/face/cell
locations). Every `*Mat` variant listed below is the matrix-assembling twin
of a live Vec-applying operator of the same base name (e.g.
`VertexToEdgeReconstructionMat` assembles the matrix form of the live
`VertexToEdgeReconstruction`); `EdgeToVertexProjection_Original` is a
superseded earlier implementation of the live `EdgeToVertexProjection`.

- **EBoundaryAdjusters** -- assembles Jacobian adjustment entries (`LM`,
  `Offset`) for edge-located degrees of freedom on domain boundary faces.
- **ComputeIsBBoundary** -- builds an `IS` of the B-field (face-located)
  boundary degrees of freedom, via a dummy Jacobian/PC/KSP probe.
- **ComputeIsCBoundary** -- builds an `IS` of the cell-centered boundary
  degrees of freedom, via the same dummy Jacobian/PC/KSP probe pattern.
- **ReadDataInVec** -- reads face/cell field data (`F_r2`, `F_phi2`, `F_z2`,
  `C2`) from a secondary DM layout into vectors, for cross-mesh data import.
- **CellToVertexProjectionVector** -- projects a cell-centered vector field
  onto vertex locations (vector-field counterpart of the live
  `CellToVertexProjectionScalar`).
- **VertexToCellReconstruction** -- reconstructs a cell-centered field from
  vertex-located values.
- **VertexToEdgeReconstructionMat** -- matrix form of the live
  `VertexToEdgeReconstruction` (vertex-to-edge reconstruction operator).
- **VertexToFaceReconstructionMat** -- matrix form of the live
  `VertexToFaceReconstruction` (vertex-to-face reconstruction operator).
- **EdgeToCellReconstructionMat** -- matrix form of the live
  `EdgeToCellReconstruction_{r,phi,z}` family (edge-to-cell reconstruction).
- **FaceToCellReconstructionMat** -- matrix-assembling face-to-cell
  reconstruction operator (no live Vec-applying analogue of this exact name
  remains; superseded by other reconstruction paths).
- **FaceToVertexProjectionMat** -- matrix form of the live
  `FaceToVertexProjection` (face-to-vertex projection operator).
- **EdgeToVertexProjection_Original** -- superseded earlier implementation of
  the live `EdgeToVertexProjection` (edge-to-vertex projection).
- **EdgeToVertexProjectionMat** -- matrix form of the live
  `EdgeToVertexProjection` (edge-to-vertex projection operator).
- **CellToFaceProjectionMat** -- matrix form of the live
  `CellToFaceProjection` (cell-to-face projection operator).
- **FromPetscVecToArray** -- unpacks a PETSc `Vec` (B-field and current
  components) into flat C arrays (`gf_BR/BP/BZ`, `g_R/P/Z`) for interop
  outside the DMStag machinery.
- **CellCoordArrays** -- fills flat R/Z coordinate arrays (`vecCR`, `vecCZ`)
  for cell centers.
- **ScatterTest** -- diagnostic exercising a PETSc `VecScatter` between two
  layouts; a standalone correctness probe, not part of any solve path.

## mimetic_operators.c

`mimetic_operators_attic.c` contains 12 functions moved verbatim (byte-for-byte
identical bodies) out of `../mimetic_operators.c`. Each had zero callers
anywhere in `src/`, `tests/`, or the built objects (verified with `nm`,
several only reachable from other functions in this same group), confirmed
with the same codegen-identity check used for `ts_functions.c` and
`geometry.c`.

One candidate from the original task list, `ApplyDeltastar2`, was **not**
moved: it is still called from `FormIFunction_Initializepsi` in
`../ts_functions.c` (a live, currently-compiled function that this task was
not scoped to touch), so moving it breaks the `mhd_core` build. It remains
declared in `mimetic_operators.h` and defined in `mimetic_operators.c`.

The divergence/gradient family below (`FormGradDerivedDivergence`,
`FormDerivedGradDivergence`, `FormDerivedGradEtaDivergence`,
`FormDerivedDivergence`, `FormDerivedGradient`, `FormDiscreteGradient`) are
alternative discrete mimetic divergence/gradient operator formulations,
distinct from the live `FormDiscreteDivergence`/`ApplyDerivedDivergence`/
`ApplyVectorLaplacian` path. `FormDerivedCurlExt` is a variant of the live
`FormDerivedCurl`/`FormDerivedCurlnores`/`FormDerivedCurlnomp` family (note
its different signature: it takes `Vec,Vec,void*` with no `TS` argument).

- **FormMaterialPropertiesMatrix** -- assembles a diagonal matrix of material
  property coefficients (`MPdiag`) on cell centers.
- **FormFaceMassMatrix** -- assembles the diagonal face-located mass matrix
  (`Mfdiag`) for the mimetic B-field discretization.
- **FormEdgeMassMatrix** -- assembles the diagonal edge-located mass matrix
  (`Mediag`) for the mimetic E/tau-field discretization.
- **FormGradDerivedDivergence** -- assembles a "gradient of derived
  divergence" operator matrix (`GD`); calls the now-attic
  `FormDerivedDivergence` internally.
- **FormDerivedGradDivergence** -- assembles a "derived gradient of
  divergence" operator matrix (`GD`), a different discrete composition than
  `FormGradDerivedDivergence`.
- **FormDerivedGradEtaDivergence** -- as `FormDerivedGradDivergence`, but with
  an `eta`-weighted (resistivity-weighted) divergence term.
- **FormDerivedDivergence** -- assembles the derived discrete divergence
  operator matrix (`D`); an alternate matrix-assembling formulation next to
  the live `ApplyDerivedDivergence`.
- **FormDerivedGradient** -- assembles both the derived divergence (`D`) and
  derived gradient (`G`) operator matrices together in one pass.
- **FormDiscreteGradient** -- assembles the discrete primary gradient
  operator matrix (`G`).
- **FormDerivedCurlExt** -- vector-field variant of the live
  `FormDerivedCurl`/`FormDerivedCurlnores`/`FormDerivedCurlnomp` family,
  taking `X`/`F` directly with no `TS`/`user` context argument.
- **FormDiscreteGradientEP_tilde** -- computes `(1/r) \tilde{\nabla}(EP)`, a
  variant of the live `FormDiscreteGradientEP`/`FormDiscreteGradientEP_noMat`
  used only by the now-attic `ApplyDeltastar`.
- **ApplyDeltastar** -- applies the discrete `Delta*` operator to `EP` via
  `FormDiscreteGradientEP_tilde` followed by `ApplyDerivedDivergence`; a
  composed-operator alternative to the live `ApplyDeltastar2`.

## monitor_functions.c

`monitor_functions_attic.c` contains 6 functions moved verbatim (byte-for-byte
identical bodies) out of `../monitor_functions.c`. `DumpSolution` and
`Dump1stVertexField` are referenced only from commented-out call sites in
`../ts_functions.c`; `getHermiteDataFD`'s only caller is `createHermiteFD`,
which is itself dead; `createHermiteFD` and `multiplybyR` have zero callers
anywhere. Confirmed with the same codegen-identity check used elsewhere in
this directory.

- **DumpSolution** -- writes a full solution dump (B, V, EP, tau, n_i, and
  derived fields) to disk; superseded by the live `DumpSolution_Cell`.
- **Dump1stVertexField** -- writes a single vertex-located vector field to
  disk, e.g. for the `CurlBxB` diagnostic referenced (commented out) in
  `FormDummyIJacobian4`.
- **DumpPsi_Cell** -- writes the poloidal flux function `psi` (cell-centered)
  to disk.
- **getHermiteDataFD** -- fills a per-cell Hermite finite-difference stencil
  data block from a flat field array; called only by `createHermiteFD`.
- **createHermiteFD** -- builds a Hermite finite-difference interpolant table
  over the full mesh, cell by cell, by calling `getHermiteDataFD`.
- **multiplybyR** -- scales a flat field array by `R` element-wise.

