# Open questions

Decisions that only the code owner can make, because the answers are not in the
code and the original author is unavailable. Each entry gives what the code does
now, the evidence, what changing it would cost, and the risk of leaving it.

Nothing here is a bug to be fixed unilaterally.

Status as of the refactor: Q2 and Q8 are resolved; Q3 is now a single argument
away from either answer; the rest remain open. Line numbers in the bodies refer
to the pre-refactor source and are kept for history -- use function names to find
the current code. Several of these change
production physics results.

Severity: **HIGH** = affects production physics output; **MEDIUM** = affects
diagnostics or reproducibility; **LOW** = hygiene.

---

## Q1 — What are the three hardcoded "isolated cells"? (HIGH)

> **Status.** Still open. After the coefficient merge (`fee16bd`) the three
> indices appear once, in `alphaecphi_isolcell`, which now chooses between
> `etawallphi_isol_cell` and `etawallphi` and calls the shared `alphaec_sum`. The
> commented-out `z > 0` alternative was not carried over; it is quoted below.
> Line numbers below are pre-refactor.

`alphaecphi_isolcell` (`src/mhd/mass_matrix_coefficients.c:3361`) selects a
special resistivity, `user->etawallphi_isol_cell`, for three cells named by
absolute index. The test appears three times, at lines 3426, 3451 and 3469:

```c
if((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176)){
//if(arrCoord[ez][ephi][er][icp[2]] > 0.0){
  alpha = cellvolume * (user -> eta0 / (user -> etawallphi_isol_cell)) / 4.0;
```

The commented-out line directly beneath each is the interesting part: it is a
test on `z > 0`, i.e. a *positional* criterion, suggesting the index list
replaced something geometric.

**This is on the production path.** `betaephi_isolcell`
(`mass_matrix_coefficients.c:1326`) calls `alphaecphi_isolcell`, and has two call
sites: `ts_functions.c:5582` in a dead residual variant, and
`ts_functions.c:10383` inside `FormIFunction_InitializeEP_halo` — which is live
as a sub-solve of `FormInitialSolution_psi`, the production initial condition.

**The obvious hypothesis is false.** These cells are not isolated in the material
field. Verified: all three have tag 0 (blanket wall), all four r/z neighbours of
each also have tag 0, and the entire tag-0 region is a *single* connected
component of 3980 cells. Geometrically they sit in the blanket wall just above
the top of the plasma, at `ez` 176 and 183 out of 200.

So "isolated cell" means something else — perhaps isolated in the toroidal
current path, or in a conductor topology not expressed by the tag field.

**Question:** what property do these three cells have? If it is computable from
the tag field or the geometry, the treatment should be computed, and would then
survive a mesh change.

**Risk of leaving it:** any change to `NR`/`NZ` silently relocates this physics to
different cells, or drops it entirely when the indices fall out of range. No
error is raised. This is the single largest obstacle to ever running this code at
another resolution.

---

## Q2 — Was the `jre` runaway-current source meant to be dropped? (RESOLVED)

> **Resolved.** `jre` and `view3d_t` are kept for the hybrid coupling (owner
> decision). The halo variant that lacked `jre` was unreachable and is in the
> attic. In the merged residual `jre` is its own switch, `f3_jre`: on for
> production, off for relaxation, as before. The relaxation is ideal by decision
> (`6ce3379`).

The production residual `FormIFunction_Vperp_viscosity`
(`ts_functions.c:2910-3885`) reads `user->jre` — the runaway-electron current
handed over from the C++ kinetic side — at 12 sites, contributing a source term to
four residual rows.

Its near-twin `FormIFunction_Vperp_viscosity_halo`
(`ts_functions.c:3886-4866`) is 98% identical but omits `jre` entirely. The
remaining differences are three mass-matrix coefficient substitutions
(`betae2` becomes `betaephi2` on the r-z edge and `betaeperp2` on the phi-z and
r-phi edges).

**Question:** was the halo variant intended to supersede the production residual?
If so, dropping `jre` is a physics regression, because the whole point of the
hybrid coupling is that the kinetic current feeds back into the field solve. If
not, why is `jre` absent there?

**Consequence for refactoring:** when these variants are unified, `jre` must
become an explicit parameter rather than being silently harmonised in either
direction.

**Status (P6-MERGE-RES):** `FormIFunction_Vperp_viscosity` and
`FormIFunction_newequilibrium_Vperp` now share one body, `vperp_residual`, which
takes an `MFD_ResidualTerms` switch struct. `jre` is the explicit switch
`f3_jre`, separate from the resistive switch `f3_resistive`. It is ON for
production and OFF for the relaxation, which matches today's behaviour. See
`residual-differences.md`.

**Note: the relaxation is ideal.** The owner decided that the initial-condition
relaxation (`FormIFunction_newequilibrium_Vperp`) uses ideal Ohm's law, with no
resistive `curl2(B)` term (`f3_resistive = FALSE`), for consistency with previous
runs. This is deliberate, not drift. The log line after the relaxation solve in
`FormInitialSolution_psi` used to read "Saving solution after relaxation with
resitivity %le". That implied the relaxation used `user->etaplasma`, which it does
not. It now reads "Saving relaxed solution (run etaplasma %le, not used by the
ideal relaxation), MPI size = %d". The format arguments are unchanged and no
computed value is affected.

---

## Q3 — Is the "plain halo" initial condition meant to be unreachable? (OPEN, now a one-word change)

> **Status.** `InitializeEP` and `_halo` now share one body, `initialize_ep`
> (`a404205`), taking the r-z and other-edge coefficients as arguments. Selecting
> plain halo -- `betaephi2` instead of `betaephi_isolcell` on r-z edges -- is now a
> single argument. The question of which is intended remains; answering it
> changes the production initial condition. Line references below are
> pre-refactor.

For the `Vperp` residual family, the halo and isolated-cell treatments are two
separate functions (`_halo` and `_halo_isolcell`), differing in a single token.

For the initialization family that distinction was collapsed:
`FormIFunction_InitializeEP_halo` hardcodes `betaephi_isolcell` at
`ts_functions.c:10383`. There is no way to select halo treatment *without* the
isolated-cell treatment in the initial condition.

**Question:** is the isolated-cell treatment intended unconditionally during
initialization?

**Risk:** fixing this changes the production initial condition, and therefore
every number downstream of it. It cannot be done as a behaviour-preserving
refactor.

---

## Q4 — Which definition of the reported toroidal current is correct? (OPEN)

> **Status.** The disabled resistive block was deleted with the other never-true
> branches (`20969e5`); the reported current is unchanged. The question stands: if
> the resistive definition is wanted, it must be reinstated deliberately. Line
> references below are pre-refactor.

`ComputeCurrent` (`src/mhd/monitor_functions.c:2475`) is live and its output is
printed each step as "Current intensity inside plasma". It contains an `if(0)`
block beginning near line 2511 that disables a resistive evaluation,
`J = (1/eta)(tau - grad EP)`, leaving only the `J = curl(B)/mu0` path from
`FormDerivedCurlnores`.

Roughly 100 lines are dead inside a function that runs.

**Question:** which definition is intended? In a resistive MHD calculation the
two agree only where Ohm's law is satisfied exactly, so the disabled branch may
have been a consistency check — or the active branch may be the fallback.

**Risk:** this is reported physics, not an internal detail. Any published current
value depends on the answer. Re-enabling the block would change the reported
numbers.

---

## Q5 — Dead features, or unfinished ones? (MEDIUM)

> **Status.** Partly overtaken. `DiagMe1`, the latent crash listed below, is gone:
> its only reader was in quarantined code. Several other fields named here are
> now read nowhere at all; see [../reference/user-struct.md](../reference/user-struct.md)
> for the current list.

Several `User` fields are read but never assigned. Since the struct is `calloc`'d
at `src/mhd/mhd.c:64`, they are permanently zero, so the code they guard never
runs:

- `savesol` — guards a `SaveSolution` call in `mhd.c`. Permanently 0, so the call
  is dead. Note `savesol` *does* appear in input decks, where it is silently
  ignored (see Q6).
- `prestep`, `prev_current`, `present_current` — declared, never written, never
  meaningfully read.
- `dataz`, `dataphi`, `datar` with `numz`, `numphi`, `numr` — never populated,
  yet `AppCtxView` prints them as set/NULL.
- `isE_boundary`, `isB_boundary`, `isni_boundary` — set to `NULL` in
  `src/kinetic/kinetic.cpp:225-227` and never assigned again. `ComputeIsEBoundary`
  and `ComputeIsBBoundary` take an `IS*` out-parameter and write to the caller's
  local, not to these fields. `ComputeIsBBoundary` has no callers at all.
- `testSpGD`, `testSpGDsamerhs` — hardcoded to 0 in `kinetic.cpp`, read nowhere.
- `mi`, `density` — written from the C++ side, read only by the debug print.
- `DiagMe1` — never written, but *read* at `src/mhd/mimetic_operators.c:1264` via
  `VecGetSubVector`, under the condition `etawall == 1 && etaplasma == 1`. If that
  branch is ever taken it operates on an unwritten `Vec` handle.

**Question:** for each, remove or wire up? `DiagMe1` deserves separate attention:
it is a latent crash, not merely dead weight.

---

## Q6 — Should unknown deck parameters be an error? (MEDIUM)

These appear in input files but are read nowhere in `src/`:

- `MHD_Config/EnableRelaxation` and `MHD_Config/EnableReadICFromBinary` — both
  passed by `perlmutter_job.sh`, and both silently doing nothing. The live
  parameter is `MHD_Config/ic_binary_load`.
- `pred_loop`, `delay_kinetic`, `savesol`, `nPR` — present in various
  `inputs/*.input`, read nowhere.

So the production GPU job script has been passing two overrides that have no
effect, and nothing reported it.

**Question:** should an unrecognised `MHD_Config` key be a hard error? That is
the cheapest permanent fix for this class of bug. The cost is that any stale deck
in someone's working directory stops running until cleaned up.

---

## Q7 — Restore ARKIMEX, or make it an error? (LOW)

In `mhd_initialize`'s `switch (user->tstype)`, `case 4` has its entire body
commented out and falls through to `break` with no `TSSetType` call, leaving
whatever `TSSetFromOptions` later selects. Setting `tstype=4` therefore does
something undefined rather than failing.

Related: `FormIJacobian_BImplicit` (`ts_functions.c:21-25`) has an **empty body**
and returns success, and is registered when `jtype == 0`. That path is currently
unreachable because `-snes_mf_operator` is checked first, but a registered
Jacobian that silently computes nothing is a trap.

**Question:** restore the ARKIMEX case, or make both an explicit
`SETERRQ(..., "not implemented")`? The latter is behaviour-preserving today, since
both paths are unreachable.

---

## Q8 — What are the magic constants in the diagnostic path? (RESOLVED)

> **Resolved by removal.** Both constants sat in never-true diagnostic branches and
> were deleted with them (`20969e5`). Neither appears in the current code.

Two unexplained literals near `Monitor` in `ts_functions.c`:

- `VecAXPBY(LapV, -10.0, 10.0, user->X0)` — the `10.0` appears to be an inlined
  `1/dt`, implying a hardcoded timestep of 0.1 that no longer matches any deck.
- `(int)(step + 800)` — an unexplained step-number offset in a dump filename.

Both are in diagnostic code. **Question:** what were they for, and is either
still meaningful?

---

## Q9 — Should the duplicated material map be deduplicated? (MEDIUM)

`inputs/mhd/veceta_grid100x02x200.txt` and `inputs/AxisSymmetricGeometry.dat`
contain the same 20000 material tags — verified bit-identical, zero mismatches.
The first is read by C (`ReadInitialData`), the second independently by C++
(`src/kinetic/kinetic.cpp`, via a bare `ifstream >>` loop with no open check and
no element count check). Nothing verifies the copies agree.

**Question:** should one be derived from the other at build or load time, or
should a consistency check be added at startup? As it stands, editing one and not
the other produces a run in which the C and C++ halves disagree about where the
plasma is, with no diagnostic.

---

## Q10 — Is `stag_vec_io`'s format acceptable for long-term restarts? (MEDIUM)

`mhd_savesolution`/`mhd_loadsolution` achieve genuine rank-independence via the
DMDA permutation trick documented above `stag_vec_io`, now in
`src/mhd/ts_functions.c` (it moved from `mhd.c` in `7999840`, and `188f772`
reworked the load path). That property is
valuable and was clearly deliberate.

But the file itself has no header, no magic number, no version, and no dimension
record. `PetscObjectSetName` is called on each record and then ignored by the
binary viewer, so the 11 records are matched **purely by loop order**. Nothing
records `Nr/Nphi/Nz`, the DOF layout, the TS time, or the step number — a restart
restores the fields with the wrong clock, relying on Parthenon's own restart for
the time.

**Question:** is a version header and a dimension check wanted? The cost is that
existing restart files become unreadable, or need a migration path.
