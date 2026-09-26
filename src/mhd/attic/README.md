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
