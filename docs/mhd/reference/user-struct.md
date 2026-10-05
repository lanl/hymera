# The `User` struct (`src/mhd/mfd_config.h`)

> **Line numbers.** File:line references in this document were taken before the
> refactor shrank `ts_functions.c` from 22,656 to ~5,800 lines; many no longer
> point at the right place. Function names are still accurate -- grep for them.
> The field list itself was re-checked against `mfd_config.h` and is current: 99
> fields, including `numC`, `monitor` and `inertia`, added during the refactor.

`User` (defined at `src/mhd/mfd_config.h:97-200`) is the single mutable-state
struct of the entire `mhd_core` C library. It is passed by pointer (`User*`,
aliased as `void*` in most PETSc callback signatures) to essentially every
function in `src/mhd/*.c`.

**Load-bearing fact:** `nm` over all 8 `mhd_core` object files shows no
globals, no file-scope statics, and no static functions anywhere in the
49k-line solver. `User` is calloc'd once, at `src/mhd/mhd.c:64`:

```c
*user = (User*) calloc(1, sizeof(User));
```

Zero-initialization from that `calloc` is therefore the *only* initializer
for any field nobody explicitly assigns. Any field below marked
`NEVER WRITTEN` sits at `0` / `NULL` / `'\0'` for the process lifetime unless
a PETSc/MPI call writes into it as an out-parameter. **Because of this, the
struct must never be split, reordered, or converted to designated
initializers** — that would either break the `calloc` zero-init guarantee for
untouched trailing fields, or (if reordered) do nothing bad by itself, but
any refactor that replaces the blanket `calloc` with a partial/designated
initializer silently reintroduces uninitialized reads.

Total field count in the struct: **99** (verified by parsing
`mfd_config.h:97-200`, correctly handling the unterminated comment described
in Hazards below).

Producer categories used in the tables:
- **DECK** — read via `pin->GetOrAdd*("MHD_Config", ...)` (or another deck
  section) in `src/kinetic/kinetic.cpp`; parameter name and default given.
- **HARDCODED** — assigned a literal in `kinetic.cpp`.
- **DERIVED** — computed in `kinetic.cpp` from other inputs.
- **C-INTERNAL** — assigned inside `src/mhd/*.c`; `file:line` given.
- **NEVER WRITTEN** — no assignment anywhere in `src/mhd/*.c` or
  `src/kinetic/kinetic.cpp`; flagged prominently.

All `mhd_context->FIELD` assignments referenced as DECK/HARDCODED/DERIVED are
from `src/kinetic/kinetic.cpp:165-244` (the `if (mhd_context != nullptr)`
block inside `Kinetic::Initialize`), unless noted otherwise.

---

## 1. Physical / normalization parameters

One-line meaning: reference scales and resistivities used to non-dimensionalize
the MHD equations and to select material resistivity by region.

| Field | Type | Producer | Value / source | Read sites | Meaning |
|---|---|---|---|---|---|
| `density` | `PetscReal` | DERIVED | `nD0` (Deuterium density, `Reference/nD0` default `1e20`) | none in `src/mhd/*.c` | Reference density; set into the struct but not read by any `mhd_core` `.c` file (`grep -c` on all `*.c` = 0 hits besides the `AppCtxView` dump). |
| `L0` | `PetscReal` | DERIVED | `a` (`Reference/a`, default `2.0`, minor radius/reference length) | `mhd.c` (dt/resistivity prints, coordinate setup), `mass_matrix_coefficients.c`, `ts_functions.c` (B0/L0 normalizations), `monitor_functions.c` (`I *= L0*L0`) | Reference length used to convert normalized geometry/coefficients back to physical units. |
| `B0` | `PetscReal` | DECK | `Reference/B0`, default `5.3` (on-axis field, Tesla) | `ts_functions.c` (normalizing `datar/dataphi/dataz`, and `datag`), `monitor_functions.c` (`J *= B0/L0`) | On-axis magnetic field magnitude; normalizes the raw B-field IC data. |
| `V_A` | `PetscReal` | DERIVED | `VA = B0/sqrt(mi*nD0*mu0)` (`Derived/VA`) | `mhd.c:201` (dt print), `ts_functions.c:15185` (`VecScale(Xcopy,1.0/V_A)`) | Alfvén speed; used for one legacy diagnostic print and one rescale. |
| `mu0` | `PetscReal` | DECK/DERIVED | `mu0` constant (`PhysicalConstants::mu0`) | `mass_matrix_coefficients.c` (resistivity ratios `mu0/eta*`), `mimetic_operators.c:1888,1893,1903` (saved/overwritten to `1.0` during `FormGradDerivedDivergence`) | Vacuum permeability; also temporarily forced to `1.0` (see Hazards-adjacent note below) while computing a normalized derived gradient. |
| `mi` | `PetscReal` | HARDCODED | `mi = PhysicalConstants::amu` (atomic mass unit) | **none** — `grep -rn "user *-> *mi\b" *.c` finds zero hits outside `AppCtxView` | Set but never consumed inside `mhd_core`; only reported by `AppCtxView`. |
| `eta0` | `PetscReal` | DERIVED | `eta0 = a*VA*mu0` (`Derived/eta0`) | `mhd.c` (prints), `mass_matrix_coefficients.c` (numerator of every `eta0/eta*` resistivity ratio) | Reference resistivity; used as normalization numerator for the region-tag resistivity lookup. |
| `eta` | `PetscReal` | DECK | `Derived/eta`, default `1.0` | `mass_matrix_coefficients.c:2481,2587,...` (only inside `/* ... */` **commented-out** debug lines, plus one live use at `2481`: `alpha = cellvolume*(eta0/eta)/4.0`) | Resistivity scale factor `[-]`; mostly referenced from dead/commented code — verify at call site before relying on it. |
| `etawall` | `PetscReal` | DECK | `Geometry/etawall`, default `4.4e-2` Ohm·m | `mass_matrix_coefficients.c` (`dataC==0` branch), `mimetic_operators.c:1890,1895,1905` (saved/forced to `1.0`/restored) | Resistivity of the blanket wall region (`dataC==0`, per the material tag map). |
| `etawallperp` | `PetscReal` | DECK | `Geometry/etawallperp`, default = `etawall` | `mass_matrix_coefficients.c` (poloidal-face resistivity branches) | Poloidal-direction wall resistivity, distinct from `etawallphi`. |
| `etawallphi` | `PetscReal` | DECK | `Geometry/etawallphi`, default = `etawall` | `mass_matrix_coefficients.c` (toroidal-face resistivity branches) | Toroidal-direction wall resistivity. |
| `etawallphi_isol_cell` | `PetscReal` | DECK | `Geometry/etawallphi_isol_cell`, default = `etawall` | `mass_matrix_coefficients.c` (isolated-cell toroidal branches, e.g. line 3428) | Lower toroidal resistivity used only in isolated wall cells. |
| `etaplasma` | `PetscReal` | DERIVED | `1.0/sigmapar` (Spitzer parallel conductivity; `Plasma/etaplasma` override) | `mhd.c` (prints, Lundquist number), `mass_matrix_coefficients.c` (`dataC==1` branch), `mimetic_operators.c:1889,1894,1904`, `ts_functions.c:20678` | Resistivity inside the separatrix (`dataC==1`, "plasma"). |
| `etasepwal` | `PetscReal` | DECK | `Geometry/etasepwal`, default = `etaplasma` | `mass_matrix_coefficients.c` (`dataC==2` branch, ~33 sites) | Resistivity inside the plasma chamber but outside the separatrix (`dataC==2`, "separatrix-wall"). |
| `etaVV` | `PetscReal` | DECK | `Geometry/etaVV`, default `1.30288e-6` Ohm·m | `mass_matrix_coefficients.c` (`dataC==-1` branch, ~34 sites), `mhd.c` (print) | Resistivity inside the vacuum-vessel wall (`dataC==-1`). |
| `etaout` | `PetscReal` | DECK | `Geometry/etaout`, default `1.30288e-3` Ohm·m | `mass_matrix_coefficients.c` (default/else branch, `dataC==-2`), `mimetic_operators.c:1891,1896,1906` | Resistivity outside the outer vacuum vessel / exterior (`dataC==-2`), also the fallback for any unmatched tag value. |
| `dampV` | `PetscScalar` | DECK | `Numerical/dampV`, default `0.01` | `ts_functions.c` (~29 sites; stabilization term in the velocity equation) | Stabilization coefficient in the constraint `((∇×B)×B − dampV·V)·e_{R/Z} = 0`. |
| `Re` | `PetscScalar` | DECK | `Derived/Re`, default `200.0` | `ts_functions.c` (viscous Laplacian term `(1/Re)∇²V`, ~5 sites), `mhd.c:203` (print) | Reynolds parameter in the velocity constraint equation. Several commented-out lines in `ts_functions.c:17572-17582` show it was once meant to be time-varying (`user->Re = ...`); those assignments are dead code (all commented out), so `Re` is effectively read-only after deck load. |

---

## 2. Geometry and mesh

One-line meaning: domain extents, grid resolution, cell spacing, and the two
PETSc `DMStag` objects (and their coordinate-array view) that discretize the
axisymmetric (R,φ,Z) mesh.

| Field | Type | Producer | Value / source | Read sites | Meaning |
|---|---|---|---|---|---|
| `rmin` | `PetscReal` | DERIVED | `Rmin * a` (`Geometry/rmin` default `1.525`, times `L0`) | `mimetic_operators.c:5884-5896` (manufactured-solution BCs), `monitor_functions.c` (`Dump1stVertexField` coordinate setup) | Minimum radius (physical units, since multiplied by `a`). |
| `rmax` | `PetscReal` | DERIVED | `Rmax * a` (`Geometry/rmax` default `4.975`, times `L0`) | same as `rmin`, plus several manufactured-solution formulas in `ts_functions.c`/`mimetic_operators.c` | Maximum radius. |
| `phimin` | `PetscReal` | HARDCODED | `0.0` | `mhd.c` (coordinate setup), `monitor_functions.c` | Minimum azimuth. |
| `phimax` | `PetscReal` | HARDCODED | `2.0*M_PI` | `mhd.c`, `monitor_functions.c`; also used to derive `dphi` in `kinetic.cpp:208` | Maximum azimuth. |
| `zmin` | `PetscReal` | DERIVED | `Zmin * a` (`Geometry/zmin` default `-2.975`, times `L0`) | `mhd.c`, `monitor_functions.c` | Minimum height. |
| `zmax` | `PetscReal` | DERIVED | `Zmax * a` (`Geometry/zmax` default `2.975`, times `L0`) | `mhd.c`, `ts_functions.c` (manufactured-solution exponential decay `exp(-Z/zmax)`, e.g. lines 14075-14077, 14934) | Maximum height. |
| `Nr` | `PetscInt` | DERIVED | `NR` (`Numerical/NR`, default `100`) | `mhd.c` (DMStag sizes, `NR=user->Nr`), `geometry.c`, `monitor_functions.c`, plus the `numC`/`numpsi`/`numg` size checks in `mhd.c:110,124,133` | Global grid points in the r direction. |
| `Nphi` | `PetscInt` | DERIVED | `Nphi` (`Numerical/Nphi`, default `2`) | same call sites as `Nr` | Global grid points in the φ direction. |
| `Nz` | `PetscInt` | DERIVED | `NZ` (`Numerical/NZ`, default `200`) | same call sites as `Nr` | Global grid points in the z direction. |
| `coorda` | `DM` | C-INTERNAL | `mhd.c:169,171` (`DMStagCreate3d(... &user->coorda)`) | ~150 sites across `geometry.c`, `mimetic_operators.c`, `mass_matrix_coefficients.c`, `ts_functions.c` (every stencil/coordinate lookup uses `user->coorda`) | Auxiliary `DMStag` used only to fetch cell/face/edge/vertex-center coordinates (`arrCoord`); destroyed implicitly with the coordinate-array restore at `mhd.c:832` (the `DM` handle itself is not explicitly `DMDestroy`'d in `mhd_destroy`, unlike `da`). |
| `da` | `DM` | C-INTERNAL | `mhd.c:159` or `mhd.c:161` (`DMStagCreate3d(... &user->da)`, branch chosen by `phibtype`) | `mhd.c` (TS setup, `DMCreateGlobalVector`), destroyed at `mhd.c:839` | Primary solution `DMStag` (dof layout: 4 vertex dofs [3×V + EP], 1 edge dof, 1 face dof, 1 cell dof). |
| `arrCoord` | `PetscScalar****` | C-INTERNAL | `mhd.c:184` (`DMStagVecGetArrayRead(dmCoorda, coordaLocal, &user->arrCoord)`) | pervasive — every mimetic-operator/mass-matrix/ts function indexes `arrCoord[ez][ephi][er][component]` for coordinate values | Read-only 4D view into `coorda`'s local coordinate vector; restored (invalidated) at `mhd.c:832`. |
| `ictype` | `PetscInt` | DECK | `MHD_Config/ictype`, default `9` | dozens of sites in `ts_functions.c`/`mhd.c`; `9` and `15` select the EFIT/Grad-Shafranov equilibrium path (production), `1-8,10-13` select manufactured analytic solutions (regression-harness convergence tests) | Initial-condition selector; gates whether `datapsi`/`datag` files are read (`mhd.c:121`) and which `FormInitialSolution*`/`FormExactSolution*` branch runs. |
| `phibtype` | `PetscInt` | DECK | `MHD_Config/phibtype`, default `1` | `mhd.c:158,168` (periodic vs. non-periodic φ boundary for `DMStagCreate3d`), `mimetic_operators.c:2781` | Boundary-condition type for φ: nonzero selects `DM_BOUNDARY_PERIODIC` in φ. |
| `dr` | `PetscReal` | DERIVED | `dR = (Rmax-Rmin)/NR` | `mhd.c` (print), `geometry.c`/`monitor_functions.c` (cell-volume/edge-length formulas) | Radial cell size. |
| `dphi` | `PetscReal` | DERIVED | `(phimax-phimin)/Nphi` (computed in `kinetic.cpp:208` *after* `phimax`/`phimin`/`Nphi` are set on the same struct) | `mhd.c` (print), `geometry.c:9182` (cell-volume formula), commented-out debug prints in `mimetic_operators.c` | Azimuthal cell size. |
| `dz` | `PetscReal` | DERIVED | `dZ = (Zmax-Zmin)/NZ` | `mhd.c` (print), `geometry.c:9182` | Vertical cell size. |

---

## 3. Time control

One-line meaning: normalized timestep, simulation time window, integrator
selection, and step counters used across restarts.

| Field | Type | Producer | Value / source | Read sites | Meaning |
|---|---|---|---|---|---|
| `dt` | `PetscReal` | DERIVED | `dt_mhd/tauA` (`parthenon/time.dt_force` divided by Alfvén time) | `mhd.c` (dt prints, `t += dt`/`t -= dt` in `mhd_step`/`mhd_resetState`), `ts_functions.c` (e.g. `VecScale(ni, dt)`), also used to derive `dumpfreq` | Normalized MHD timestep length. |
| `itime` | `PetscReal` | DERIVED | `itime/tauA` (`Numerical/itime`, default `0.0`, divided by Alfvén time) | `ts_functions.c` — gates ~11 `ictype==9/15 && itime==0.0` branches that pick the one-time equilibrium-loading path vs. restart path | Initial time for the MHD counters; used purely as a "is this a fresh start" flag (`itime==0.0`) rather than as a running clock. |
| `ftime` | `PetscReal` | DERIVED | `final_time/tauA` (`parthenon/time.tlim` divided by Alfvén time) | `mhd.c:307` (`ftime = user->ftime` local copy for `TSSetMaxTime` region), `ts_functions.c:17605` (`time >= user->ftime` dump-trigger check); also used to derive `dumpfreq` in `kinetic.cpp:237` | Final simulation time (normalized). |
| `tstype` | `PetscInt` | DECK | `MHD_Config/tstype`, default `2` | `mhd.c:229-263` (`switch(tstype)` selecting `TSEULER`/`TSBEULER`/`TSCN`/.../`TSTHETA`), `mhd.c:266` (`tstype > 1` gates matrix-free/Jacobian setup) | Timestepping method: 1=Forward Euler, 2=Backward Euler, 3=Crank-Nicholson, 4=(ARKIMEX, currently commented out — falls through to no-op), 5=ROSW, 6=BDF2, 7=Theta. |
| `oldstep` | `PetscInt` | HARDCODED | `kinetic.cpp:240` sets it to `0` | `ts_functions.c` (`step + user->oldstep` used pervasively as the "absolute" step counter for logging/dump filenames; `ictype==9 && oldstep==0` gates one-time-only dump logic at line 17623) | Last step number from a previous run, added to the current step counter so restart logs/filenames continue numbering. Always `0` in the current pipeline since nothing ever loads a nonzero value (no restart path sets it) — see Hazards for why editing nearby lines is risky. |
| `n_record` | `PetscInt` | DECK-adjacent (HARDCODED) | `kinetic.cpp:241` sets it to `0` | **none** — `grep -n "n_record\b" *.c` (all of `src/mhd/*.c`) hits only the `AppCtxView` print in `mfd_config.c:56` | Documented as "Counter for successful TS steps" but never incremented anywhere in `mhd_core`; effectively dead. |
| `n_record_Steady_jRE` | `PetscInt` | HARDCODED | `kinetic.cpp:242` sets it to `0` | **none** — same as above, only in `AppCtxView` (`mfd_config.c:57`) | Documented as "counter for TS step where j_RE change is less than 10%" but never incremented; dead. |

---

## 4. Live PETSc objects (Vec / IS / Mat / KSP / TS)

One-line meaning: solution vectors, index sets partitioning the DOF layout,
and the block-preconditioner matrices/solvers built once per run inside
`ts_functions.c`'s Jacobian/preconditioner setup routines.

| Field | Type | Producer | Value / source | Read sites | Meaning |
|---|---|---|---|---|---|
| `X` | `Vec` | C-INTERNAL | `mhd.c:188` (`DMCreateGlobalVector(user->da, &user->X)`) | `mhd.c` (`TSSetSolution`, `FormInitialSolution*`, `SaveSolution`), destroyed `mhd.c:835` | Current global solution vector (all fields packed per the `da` dof layout). |
| `X0` | `Vec` | C-INTERNAL | `mhd.c:187` | `mhd.c:684,743` (`mhd_step`/`mhd_resetState` copy `X`↔`X0` for the reset/rollback path), `ts_functions.c` (`Monitor`'s corrector-step diff, `FormIFunction_BImplicit`'s `LapV` term) | Solution vector saved at the start of a step; used to compute the corrector difference and to support `mhd_resetState`'s rollback. |
| `isV` | `IS` | C-INTERNAL | `mhd.c:387` (`ISDifference(islist[0], user->isEP, &user->isV)`) | `mhd.c`, `ts_functions.c` (Schur-complement/preconditioner assembly, ~10+ sites) | Index set for the velocity DOFs. |
| `isni` | `IS` | C-INTERNAL | `mhd.c:391` (`ISDuplicate`) | `mhd.c`, `ts_functions.c` (field-decomposition lists) | Index set for ion number density DOFs. |
| `isB` | `IS` | C-INTERNAL | `mhd.c:390` (`ISDuplicate`) | `mhd.c`, `ts_functions.c` (Schur complement for the B-field block, divergence checks in `Monitor`) | Index set for magnetic-field DOFs. |
| `isEP` | `IS` | C-INTERNAL | `mhd.c:362` (`ISDuplicate(isEPdup, ...)`) | `mhd.c`, `ts_functions.c` (`MatCreateSubMatrix(J,isEP,isEP,...,&DiagBlock_EP)`) | Index set for electrostatic-potential DOFs. |
| `istau` | `IS` | C-INTERNAL | `mhd.c:389` (`ISDuplicate`) | `mhd.c`, `ts_functions.c` (Schur complement for B block) | Index set for the divergence-free `tau` field. |
| `isE_boundary` | `IS` | C-INTERNAL | `geometry.c:1019` inside `ComputeIsEBoundary` (`ISDuplicate(isdup, isE_boundary)`) — but see read-nowhere note | HARDCODED to `NULL` in `kinetic.cpp:226`; `mfd_config.c:97` only reads it for the debug dump | `ComputeIsEBoundary` (`geometry.c:556`) writes through an output *pointer parameter* named `isE_boundary`, not through `user->isE_boundary` — **`grep -rn "user *-> *isE_boundary\s*="` finds zero hits**. The function is called from `mimetic_operators.c:1359,1922` with a *local* variable `isC_boundary` as the actual argument, so `user->isE_boundary` itself is never written. Confirmed **read-nowhere / effectively dead** (only ever `NULL`). |
| `isB_boundary` | `IS` | C-INTERNAL (same caveat) | Same pattern via `ComputeIsBBoundary` (`geometry.c:1026,1208`), called with local `isC_boundary`, not `user->isB_boundary` | Same as above | Same situation as `isE_boundary`: the struct field is never actually assigned; only `NULL`-checked in `AppCtxView`. |
| `isni_boundary` | `IS` | NEVER WRITTEN | — | `mfd_config.c:99` (`AppCtxView` NULL-check only) | No producer anywhere; stays `NULL`. |
| `DiagMe` | `Vec` | C-INTERNAL | `mimetic_operators.c:2920` (`VecCopy(Mediag, user->DiagMe)` inside `FormDerivedDivergence`) | `mimetic_operators.c:1272` (`VecGetSubVector(user->DiagMe, isBE, &user->MeVec)`) | Diagonal of the mass matrix `M_e` **with** material properties, on the extended DM. Note: the `Vec` handle itself is never explicitly created (no `VecDuplicate`/`DMCreateGlobalVector` targets `user->DiagMe` directly) — `VecCopy` into an uncreated `Vec` would fail, so in practice this relies on `DiagMe` having been created as a side effect elsewhere (likely aliasing/reuse not visible from a single-file grep); flagged as worth double-checking against the full call graph if this path is ever touched. |
| `MeVec` | `Vec` | C-INTERNAL | `mimetic_operators.c:1272` (`VecGetSubVector` output) | `mimetic_operators.c:1275-1276` (`VecAssemblyBegin/End`) | Sub-vector of `DiagMe` restricted to edge/face DOFs (`isBE`), i.e. `M_e` diagonal in the regular (non-extended) DM. |
| `DiagMe1` | `Vec` | NEVER WRITTEN | — | `mimetic_operators.c:1264` (read as input to `VecGetSubVector`, only when `etawall==1 && etaplasma==1`), `mfd_config.c:103` | Documented as "`M_e` diagonal **without** material properties, extended DM" — no code path assigns it, so the `VecGetSubVector(user->DiagMe1, ...)` call at `mimetic_operators.c:1264` reads an uninitialized/NULL `Vec` handle whenever `etawall==1 && etaplasma==1` (which only occurs transiently inside `FormGradDerivedDivergence`'s "set resistivities to 1" trick, `mimetic_operators.c:1893-1896`). **NEVER WRITTEN — flagged.** |
| `MeVec1` | `Vec` | C-INTERNAL (conditionally) | `mimetic_operators.c:1264` (`VecGetSubVector` output, same guarded branch) | `mimetic_operators.c:1267-1268` | Sub-vector of `DiagMe1`; only populated in the same guarded branch, and since `DiagMe1` is never written, this is downstream of the same problem. |
| `ts` | `TS` | C-INTERNAL | `mhd.c:216` (`TSCreate`) | `mhd.c` (`TSStep`, `TSSetDM`, `TSMonitorSet`, ...), `ts_functions.c` (`TSGetSNES`, ×4) | The PETSc time integrator handle; destroyed `mhd.c:837`. |
| `DiagBlock_V` | `Mat` | C-INTERNAL | `ts_functions.c:18004` etc. (`MatCreate`, then `MatGetSchurComplement` output, 3 separate preconditioner-setup functions) | `ts_functions.c` (KSP operator setup, ~12 sites), destroyed `ts_functions.c:18205/18310/18451` | Schur-complement diagonal block for the velocity block of the block preconditioner. |
| `DiagBlock_EP` | `Mat` | C-INTERNAL | `ts_functions.c:18352` (`MatCreateSubMatrix(J,isEP,isEP,...)`) | `ts_functions.c:18372` (KSP operator), destroyed `18452` | Diagonal block for the electrostatic-potential block. |
| `DiagBlock_B` | `Mat` | C-INTERNAL | `ts_functions.c:18331` (`MatCreate`), `18355` (`MatGetSchurComplement` output) | `ts_functions.c:18382`, destroyed `18453` | Diagonal block for the magnetic-field block. |
| `KSP_V` | `KSP` | C-INTERNAL | `ts_functions.c:18069` etc. (`KSPCreate`, 3 call sites) | `ts_functions.c` (tolerances, `KSPGetPC`, `KSPSetUp`), destroyed `18209/18312/18456` | Solver for the velocity block of the user-provided preconditioner. |
| `KSP_B` | `KSP` | C-INTERNAL | `ts_functions.c:18381` | `ts_functions.c`, destroyed `18457` | Solver for the magnetic-field block. |
| `KSP_EP` | `KSP` | C-INTERNAL | `ts_functions.c:18371` | `ts_functions.c`, destroyed `18455` | Solver for the electrostatic-potential block. |
| `isALL_V` | `IS` | C-INTERNAL | `ts_functions.c:18017` (`ISDifference(isALL, user->isV, &user->isALL_V)`) | `ts_functions.c` (Schur complement, `MatCreateSubMatrix` for off-diagonal blocks), destroyed `ts_functions.c:18211/18314` and also `mhd.c:824` (belt-and-suspenders destroy in `mhd_destroy`) | Index set for `{EP, tau, B, ni}` — i.e. everything except velocity. |
| `OffDiagBlock_U` | `Mat` | C-INTERNAL | `ts_functions.c:18006` (`MatCreate`), `18031` (`MatCreateSubMatrix` output) | `ts_functions.c:18182` (`MatMult`), destroyed `18207` | Upper off-diagonal Jacobian block in the `{ETBN,V}` partitioning. |
| `OffDiagBlock_L` | `Mat` | C-INTERNAL | `ts_functions.c:18005`, `18032` | `ts_functions.c:18177` (`MatMultAdd`), destroyed `18206` | Lower off-diagonal Jacobian block. |

---

## 5. Input data arrays (raw IC data read from `inputs/mhd`)

One-line meaning: flat arrays read once at startup from text files under
`user->input_folder`, indexed by hand-rolled offset arithmetic into the
(R,φ,Z) grid; and their lengths.

| Field | Type | Producer | Value / source | Read sites | Meaning |
|---|---|---|---|---|---|
| `dataC` | `PetscReal*` | C-INTERNAL | `mhd.c:109` (`ReadInitialData(&user->dataC, &user->numC, "veceta_grid...")`), unconditional (every `ictype`) | `geometry.c:10093`, `mass_matrix_coefficients.c` (region-tag branches, ~30+ sites), `monitor_functions.c` (current-integral region masks) | **Verified:** a material-region *tag* field, not a resistivity value. Layout is `Nr*Nphi*Nz`, indexed `[er + ephi*Nr + ez*Nphi*Nr]` (cell-centered). The tag→material map, read from `mass_matrix_coefficients.c` near line 2432 (`alphaec2`-style branches): `1`=plasma (`etaplasma`), `2`=separatrix-wall (`etasepwal`), `0`=blanket wall (`etawall`), `-1`=vacuum vessel (`etaVV`), `-2`/else=exterior (`etaout`). **Verified bit-identical to `inputs/AxisSymmetricGeometry.dat` transposed**, which the C++ side reads independently near `src/kinetic/kinetic.cpp:303-325` into its own `indicator` view — **two separate copies of the same geometry data, with no cross-consistency check between them.** |
| `dataz` | `PetscReal*` | NEVER WRITTEN | — | `ts_functions.c:14737,14752` (only inside the `ictype==9 && itime==0.0` branch of `FormInitialSolution`, under the initial-B-field-loading `if(1){...}` block) | Documented as "initial B_z component" data, layout `Nr*Nphi*Nz` (element-ish); **no `ReadInitialData` (or any other) call ever assigns `user->dataz`** — confirmed via `grep -rnE "&\s*user\s*->\s*dataz|user\s*->\s*dataz\s*=" *.c` returning zero hits. The read sites at `ts_functions.c:14737/14752` therefore dereference a `NULL` pointer whenever that code path executes with `ictype==9`. **This is a live bug risk, not dead code** — it is only masked by the fact that `ictype==9`'s IC data actually comes through `datar`/`dataphi`/`dataz`-style filenames that were apparently superseded by the `datapsi`/`datag` (Grad-Shafranov) path added later (see `mhd.c:99-145`), which is what current input decks (`ictype 9/15`) actually rely on. Flagged as **NEVER WRITTEN**, confirming the task's expectation. |
| `dataphi` | `PetscReal*` | NEVER WRITTEN | — | `ts_functions.c:14733,14747` (same branch as `dataz`) | Same situation as `dataz`: initial B_φ component, never assigned. **NEVER WRITTEN**, confirmed. |
| `datar` | `PetscReal*` | NEVER WRITTEN | — | `ts_functions.c:14729,14742` (same branch) | Same situation: initial B_r component, never assigned. **NEVER WRITTEN**, confirmed. |
| `numz` | `PetscInt` | NEVER WRITTEN | — | `ts_functions.c:14759` (bounds-check print against a `countz` accumulator) | Length of `dataz`; since `dataz` is never populated via any `ReadInitialData` call, `numz` is never set either. **NEVER WRITTEN**, confirmed. |
| `numphi` | `PetscInt` | NEVER WRITTEN | — | `ts_functions.c:14758` | Length of `dataphi`. **NEVER WRITTEN**, confirmed. |
| `numr` | `PetscInt` | NEVER WRITTEN | — | `ts_functions.c:14757` | Length of `datar`. **NEVER WRITTEN**, confirmed. |
| `datag` | `PetscReal*` | C-INTERNAL | `mhd.c:132` (`ReadInitialData(&user->datag, &user->numg, "vecg_grid...")`), only when `ictype==9 \|\| ictype==15` | `ts_functions.c:20326,20330,21118,21122` (`FormInitialSolution_psi`'s B_φ = G(ψ)/R reconstruction) | R·B_φ, phi-face layout `Nr*(Nphi+1)*Nz`. Freed at `ts_functions.c:20362` (and a commented-out duplicate free at `21157`). Consumed only on the `ictype` 9/15 equilibrium path, as expected. |
| `datapsi` | `PetscReal*` | C-INTERNAL | `mhd.c:123` (`ReadInitialData(&user->datapsi, &user->numpsi, "vecpsi_grid...")`), only when `ictype==9 \|\| ictype==15` | `ts_functions.c:20337-20348,21129-21140` (poloidal-flux edge values feeding the E-field IC) | Poloidal flux ψ, r-z edge layout `(Nr+1)*Nphi*(Nz+1)`. Freed at `ts_functions.c:20363`. |
| `numg` | `PetscInt` | C-INTERNAL | `mhd.c:132` (`ReadInitialData` out-param), checked at `mhd.c:133` against `Nr*(Nphi+1)*Nz` | `ts_functions.c:20334` (bounds check) | Length of `datag`; validated to match grid size at load time. |
| `numpsi` | `PetscInt` | C-INTERNAL | `mhd.c:123`, checked at `mhd.c:124` against `(Nr+1)*Nphi*(Nz+1)` | `ts_functions.c:20351` (bounds check) | Length of `datapsi`; validated at load time. |
| `numC` | `PetscInt` | C-INTERNAL | `mhd.c:109` (`ReadInitialData` out-param), checked at `mhd.c:110` against `Nr*Nphi*Nz` | none besides the check itself | Length of `dataC`. **Recently added** field per the task context — it now holds a real value where previously (per the task description) the count was written to a discarded local; confirmed by the explicit `PetscCheck(user->numC == user->Nr*user->Nphi*user->Nz, ...)` guard immediately after the read, which is exactly the kind of consistency check that a discarded local could never support. |

---

## 6. Flags (deck-controlled run-mode switches)

One-line meaning: integer toggles read from the deck (or hardcoded) that gate
optional/legacy code paths — Jacobian strategy, debug verbosity, dump/monitor
diagnostics, and legacy save/test paths.

| Field | Type | Producer | Value / source | Read sites | Meaning |
|---|---|---|---|---|---|
| `jtype` | `PetscInt` | DECK | `MHD_Config/jtype`, default `2` | `mhd.c:279-286` (`switch`-like `if` chain selecting analytical vs. FD-with/without-coloring Jacobian) | Jacobian strategy: 0=user-provided analytical, 1=slow FD, 2=FD with coloring. |
| `debug` | `PetscInt` | DECK | `MHD_Config/debug`, default `0` | `mass_matrix_coefficients.c` (7 sites gating verbose `PetscPrintf`/`MatView`/`VecView` dumps), `mimetic_operators.c` (multiple `if(user->debug)` blocks) | Enables verbose debug printing throughout the mimetic-operator/mass-matrix assembly code. |
| `dump` | `PetscInt` | DECK | `MHD_Config/dump`, default `0` | `ts_functions.c:17711,17715,17765` (`Monitor`'s `.vtr`/cell-solution dump path) | Enables saving solution snapshots in `.vtr`/cell-dump files via `Monitor`. |
| `monitor` | `PetscInt` | DECK | `MHD_Config/monitor`, default `0` | `mhd.c:729` (`if (user->monitor) { ... Monitor(...) }`) | **Recently added.** Gates the per-step `Monitor` diagnostics (step norms, max\|div B\|, toroidal currents) that `mhd_step` now calls directly, since PETSc's own `TSMonitorSet` registration (`mhd.c:227`) never fires (the driver advances via `TSStep`, not `TSSolve`). Off by default because `Monitor` performs a nested linear solve per step; the regression harness turns it on and uses its printed output as a per-step fingerprint. |
| `inertia` | `PetscInt` | DECK | `MHD_Config/inertia`, default `0` | `FormIFunction_Vperp_viscosity` in `ts_functions.c`, which picks one of two constant `MFD_ResidualTerms` from it | Enables the advective inertia term n_i (V.grad)V in the production momentum residual. Off by default: the term was historically multiplied by a literal `0.0`, and with the flag off it still is, verbatim, because deleting it moves results by 1-14 ULP. Read in the wrapper, not inside the shared body, so the default path compiles exactly as before. |
| `prestep` | `PetscInt` | NEVER WRITTEN | — | none in `src/mhd/*.c` besides `mfd_config.c:69` (`AppCtxView` print) | Documented as "activate the prestep to approximate the runaway current contribution," but no assignment exists anywhere (not in `kinetic.cpp`, not in any `.c` file). **NEVER WRITTEN**, confirmed. |
| `Ebc` | `PetscInt` | HARDCODED | `kinetic.cpp:235` sets it to `0` | `ts_functions.c:16489,17519` (`if (ictype==9 && user->Ebc)`) | Type of boundary condition for the E field; always `0` under the current deck wiring, so both gated branches are permanently dead in practice unless something else sets it (nothing does). |
| `savecoords` | `PetscInt` | DECK | `MHD_Config/savecoords`, default `0` | `mhd.c:219` (`if (user->savecoords) { SaveCoordinates(...); return(0); }`) | If set, `mhd_initialize` saves cell/face/edge-center coordinates to `.m` files and returns immediately — a diagnostics-only early-exit mode. |
| `savesol` | `PetscInt` | NEVER WRITTEN | — | `mhd.c:699` (`if (user->savesol) { SaveSolution(...); }`) | Documented as "save face-centered B values in `.m` files." **No assignment anywhere** (not DECK, not HARDCODED, not C-INTERNAL) — `grep -n "savesol" *.c` in `kinetic.cpp` returns nothing, and the only `.c` writes are absent; `mhd.c:699`'s guard is therefore permanently false. **NEVER WRITTEN**, confirmed. |
| `tempdump` | `PetscInt` | HARDCODED | `kinetic.cpp:236` sets it to `0` | `ts_functions.c:17604` (`Monitor`'s "save intermediate solution in .dat files" gate) | Always `0` under current wiring; the intermediate-`.dat`-file dump path in `Monitor` is dead unless something else sets this later (nothing does). |
| `dumpfreq` | `PetscInt` | DERIVED | `kinetic.cpp:237`: `ceil(ftime / (10.0*dt))` | `ts_functions.c:17605` (`step % user->dumpfreq == 0`) | Frequency (in steps) for the (currently dead, since `tempdump==0`) intermediate-dump path. |
| `testSpGD` | `PetscInt` | HARDCODED | `kinetic.cpp:238` sets it to `0` | **none** — `grep -n "testSpGD\b" *.c` (outside `mfd_config.c`'s print) returns zero hits | Documented as a flag for testing `S` vs. `S - dt*GD` Schur-complement variants; the comparison code it would gate does not appear to exist (or was removed) in the current `ts_functions.c`. Set but unused. |
| `testSpGDsamerhs` | `PetscInt` | HARDCODED | `kinetic.cpp:239` sets it to `0` | **none** — same as `testSpGD` | Same situation: set but unused. |
| `ic_binary_mode` | `char` | DERIVED | `kinetic.cpp:233`: `pin->GetOrAddInteger("MHD_Config","ic_binary_load",1)==1 ? 'l' : 'c'` | `ts_functions.c:20389` (`if (user->ic_binary_mode == 'l')`) | Selects binary-IC load (`'l'`) vs. compute (`'c'`) mode for `FormInitialSolution_LargeData`'s companion binary-vector path. |
| `input_folder` | `char[PETSC_MAX_PATH_LEN]` | DECK | `MHD_Config/input_folder`, default `"../../inputs/mhd"` | `mhd.c:108,122,131` (`PetscSNPrintf` building the `veceta_grid.../vecpsi_grid.../vecg_grid...` filenames) | Directory containing the grid-resolution-specific IC text files. |
| `ic_binary_path` | `char[PETSC_MAX_PATH_LEN]` | DECK | `MHD_Config/ic_binary_path`, default `""` | `ts_functions.c:20391-20392,20681` (binary `Vec` load/save viewer path) | Path to a PETSc binary-format solution vector, used by the `ic_binary_load`/large-data IC path. |

---

## 7. Counters / mutable run state

One-line meaning: quantities recomputed every step (or every `Monitor` call)
that track diagnostic scalars across the run.

| Field | Type | Producer | Value / source | Read sites | Meaning |
|---|---|---|---|---|---|
| `Iphi1` | `PetscScalar` | C-INTERNAL | `monitor_functions.c:2585` (zeroed), then `MPI_Reduce(&I1, &user->Iphi1, ...)` at `~2626` inside `ComputeCurrent` | `ts_functions.c:17798-17799` (`Monitor` prints it, appends to `TotalCurrent.txt`) | Toroidal current intensity inside the plasma (integral of `J` over the `dataC≈1.5±0.7` region). Only ever computed when `Monitor`→`ComputeCurrent` runs, i.e. gated by `user->monitor`. |
| `Iphi2` | `PetscScalar` | C-INTERNAL | `monitor_functions.c:2586`, reduced at `~2632` | `ts_functions.c:17799` | Current intensity outside the plasma (`dataC≈0` region, "wall"). Same gating as `Iphi1`. |
| `Iphi3` | `PetscScalar` | C-INTERNAL | `monitor_functions.c:2587`, reduced at `~2635` | `ts_functions.c:17800` | Current intensity inside the vacuum vessel (`dataC≈-1` region). Same gating. |
| `prev_current` | `double` | NEVER WRITTEN | — | **none** — zero hits for `user *-> *prev_current` in any `.c`, and zero hits in `kinetic.cpp` | Documented as "runaway current of the previous time step." No producer and no consumer anywhere in `mhd_core` or the kinetic coupling code. **NEVER WRITTEN and read-nowhere**, confirmed. |
| `present_current` | `double` | NEVER WRITTEN | — | **none** — same as above | No producer or consumer anywhere. **NEVER WRITTEN and read-nowhere**, confirmed. |

---

## 8. The C++ coupling window (`jre`)

One-line meaning: the single field through which the C++/Kokkos kinetic side
hands the MHD solver a runaway-current source term, without any other shared
mutable state crossing the C/C++ boundary through `User`.

| Field | Type | Producer | Value / source | Read sites | Meaning |
|---|---|---|---|---|---|
| `jre` | `view3d_t` | C-INTERNAL (struct itself) / cross-language wiring | `src/kinetic/kinetic.cpp:244`: `mhd_context->jre = wrap_view(Jre_mhd.view_host())` — wraps the host mirror of a Kokkos `DualView3` (`Jre_mhd`, shape `[NR][NZ][3]`) allocated in `kinetic.cpp:161` | `mhd.c:791-796` (`mhd_getF(..., fid_Jre, ...)` copies `jre.data` out via strided indexing), `ts_functions.c:2929,3603-3685` (`FormIFunction_BImplicit`-family functions read `jre.data[...]` as a source term for the runaway-current contribution to the B-field/E-field equations) | `view3d_t` is a plain `{data, dim0..2, stride0..2}` struct (no ownership) — `User` does not own this memory; it merely holds a pointer into the C++-side Kokkos view's host buffer. The C++ task `AdvanceBackgroundFields` (`src/tasks/BackgroundFields.cpp:4-10`) calls `Jre.sync_host()` before `mhd_step`, and `AccumulateView` (`src/tasks/ViewManipulation.cpp:28-46`) accumulates the kinetic particle current into `Jre_mhd`'s device view between MHD steps — so this is the one deliberate, explicit two-way handshake between the languages, as opposed to every other field which flows one direction (deck → `User` → `mhd_core`). No writes to `jre.data[...]` were found inside `src/mhd/*.c` itself (only reads) — the struct is written to from the C++ side (via the Kokkos view it aliases) and read from the C side. |

---

## Hazards

1. **Unterminated comment swallowing two struct fields.** At
   `mfd_config.h:168-170`:
   ```c
   PetscInt    oldstep;           /* Last step in previous simulation when a restart is used
   //PetscInt    dtchange;        /* Flag to indicate if a restart is used with a new time step : olditime...olddt...newitime...newdt...newftime */
   //PetscReal   olddt;           /* Length of old timestep when a restart is used with a new time step */
   ```
   The doc-comment on `oldstep` opens with `/*` on its own line and does
   **not** close there; the comment body runs through the next two
   `//`-prefixed lines and only closes at the first `*/` it encounters, which
   is the one at the end of the `dtchange` line (line 169). That leaves the
   *text* of the `olddt` line (line 170) outside any comment — but that line
   is itself `//`-prefixed, so it still doesn't declare a field. Net effect:
   `dtchange` and `olddt` are **not** struct members today (both lines are
   fully commented out, one via the runaway block comment, one via its own
   `//`), so the struct correctly has no fields for them. The hazard is
   purely for future edits: naively "closing" the `oldstep` comment on its
   own line, or uncommenting `dtchange`/`olddt` without noticing the
   overlapping `/* ... */`, will either (a) silently turn the `dtchange`/
   `olddt` line bodies into live (but bogus, comment-fragment) code, or (b)
   change what `oldstep`'s doc-comment covers, or (c) — if someone "fixes"
   the unterminated comment by adding a closing `*/` right after
   "restart is used" — expose the literal text `//PetscInt dtchange...` as
   real (broken) source. Any edit to lines 168-170 must be done by rewriting
   all three lines together, not by touching one in isolation.

2. **~27 unscoped macros aliasing `DMSTAG_*` constants to bare identifiers**
   (`mfd_config.h:68-94`), included transitively by every consumer of
   `mhd.h`:
   ```
   BACK_DOWN_LEFT, BACK_DOWN, BACK_DOWN_RIGHT, BACK_LEFT, BACK, BACK_RIGHT,
   BACK_UP_LEFT, BACK_UP, BACK_UP_RIGHT, DOWN_LEFT, DOWN, DOWN_RIGHT, LEFT,
   ELEMENT, RIGHT, UP_LEFT, UP, UP_RIGHT, FRONT_DOWN_LEFT, FRONT_DOWN,
   FRONT_DOWN_RIGHT, FRONT_LEFT, FRONT, FRONT_RIGHT, FRONT_UP_LEFT, FRONT_UP,
   FRONT_UP_RIGHT
   ```
   That is 27 macros, all `#define`d with no namespace/scoping (plain
   `#define NAME DMSTAG_NAME`, no `#undef` anywhere in the header). Because
   they are C preprocessor macros, they apply to **every** translation unit
   that (transitively) includes `mfd_config.h` — and that is not limited to
   `.c` files. `grep -rln '#include.*mhd\.h'` shows these C++ files pull it
   in: `src/tasks/BackgroundFields.cpp`, `src/tasks/Interpolate.cpp`,
   `src/mhd/mhd.cpp`, `src/mhd/MHDDriver.h`, `src/kinetic/hybrid.cpp`,
   `src/kinetic/ConservationDriver.h`, `src/kinetic/AvalancheDriver.h`,
   `src/kinetic/SingleParticleDriver.h`, `src/kinetic/profile.cpp`,
   `src/kinetic/ProfileDriver.h`, `src/kinetic/kinetic.cpp`,
   `src/kinetic/HybridDriver.h`. Any of these files (or anything they
   `#include` afterward) that declares a local variable, enum member, or
   calls a library function/macro literally named `LEFT`, `RIGHT`, `UP`,
   `DOWN`, `BACK`, `FRONT`, or `ELEMENT` — all disturbingly common short
   identifiers — will silently collide with these macros and get rewritten
   to `DMSTAG_LEFT` etc. by the preprocessor, with no compiler warning if the
   substitution happens to still typecheck. A targeted grep for bare uses of
   those seven short names in the affected C++ files found no current
   collisions, but the risk is structural (the macros are global once
   `mhd.h` is included) and worth flagging before adding new code to any of
   the 12 files above, especially anything that pulls in another
   header/library that itself defines short enum-style names.

3. **The `calloc` zero-init at `mhd.c:64` is load-bearing across the whole
   struct.** As shown above, at least 8 fields (`isni_boundary`,
   `DiagMe1`/downstream `MeVec1`, `prestep`, `savesol`, `testSpGD`,
   `testSpGDsamerhs`, `prev_current`, `present_current`, `dataz`/`dataphi`/
   `datar`/`numz`/`numphi`/`numr`) are read somewhere in the codebase while
   relying on `calloc`'s zero-fill as their only defined value (`0`/`NULL`).
   Splitting `User` into multiple structs, reordering its members, or
   replacing the blanket `calloc(1, sizeof(User))` with any form of
   designated/partial initialization would change some of these from
   "reliably zero" to "uninitialized," which for the `PetscInt`/`double`
   flags would likely still read as garbage-but-plausible (dangerous), and
   for the pointer fields (`DiagMe1`, `dataz`, etc.) would turn a currently
   inert `NULL`-guarded dead path into an actual wild-pointer dereference the
   first time one of those guarded branches executes.

---

## Verification summary

- **Total struct fields:** 99 (parsed directly from `mfd_config.h:97-200`
  with a comment-aware scanner that correctly attributes `oldstep`,
  `input_folder`, `ic_binary_path` despite the unterminated-comment hazard
  above). Per-section row counts: 18 in §1 (Physical/normalization), 17 in
  §2 (Geometry/mesh), 7 in §3 (Time control), 24 in §4 (Live PETSc objects),
  12 in §5 (Input data arrays), 15 in §6 (Flags), 5 in §7 (Counters/mutable
  run state), 1 in §8 (`jre`). `18+17+7+24+12+15+5+1 = 99`, and every one of
  the 99 names parsed from the header (`density` ... `jre`) was
  cross-checked to appear as exactly one table row above.
- **NEVER WRITTEN fields found:** `mi`(*) `density`(*), `savesol`, `prestep`,
  `prev_current`, `present_current`, `dataz`, `dataphi`, `datar`, `numz`,
  `numphi`, `numr`, `isni_boundary`, `DiagMe1`, `testSpGD`,
  `testSpGDsamerhs`. (*`mi`/`density` **are** written, by `kinetic.cpp`, but
  are read nowhere in `mhd_core` — listed here only because they're
  otherwise inert; they are not NEVER-WRITTEN in the strict sense and are
  correctly marked DERIVED/HARDCODED in §1, not flagged as NEVER WRITTEN in
  the tables.) The task's predicted NEVER-WRITTEN/read-nowhere list —
  `savesol`, `prestep`, `prev_current`, `present_current`, and
  `dataz`/`dataphi`/`datar` with `numz`/`numphi`/`numr` — is **fully
  confirmed** by grep: none of those 10 fields has any assignment site in
  `src/mhd/*.c` or `src/kinetic/kinetic.cpp`. Two additional fields not
  called out in the prompt's prediction were also found to be effectively
  dead: `isE_boundary`/`isB_boundary` (their would-be producer functions
  write to a same-named *local variable*, not `user->isE_boundary`/
  `user->isB_boundary`, at every call site) and `isni_boundary` (no producer
  at all), plus `DiagMe1` (documented as populated but no code path assigns
  it) and `testSpGD`/`testSpGDsamerhs` (set to `0` by `kinetic.cpp` but never
  read anywhere in `mhd_core`).
- **Read-nowhere-in-mhd_core fields (written but never consumed by the C
  solver):** `mi`, `density` (both only appear in `AppCtxView`'s print).
- **Contradictions with the stated VERIFIED FACTS:** none found. `dataC`'s
  tag→material mapping, layout, and duplicate-source relationship with
  `inputs/AxisSymmetricGeometry.dat` all check out as described. `datapsi`/
  `datag` are confirmed read only on the `ictype` 9/15 path. `numC` is
  confirmed to now hold a real, checked value. `monitor` is confirmed as a
  recently-added DECK flag (default 0) gating `Monitor` from `mhd_step`.

## Files touched

Only this file was created/edited:
`/workspace/docs/mhd/reference/user-struct.md`. No file under `src/` was
modified.
