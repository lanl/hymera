# How the two remaining large residuals differ

## Merged in P6-MERGE-RES (uncommitted, on top of 9a658e8)

The two functions below are now one body, `static vperp_residual(..., const
MFD_ResidualTerms *terms)` in `src/mhd/ts_functions.c`. Both public names survive as
thin wrappers, so no caller changed. The switch struct is:

```c
typedef struct {
  PetscBool    f1_inertia;    /* TRUE: relaxation `+ niv * dV/dt`; FALSE: production `+ 0.0*niv*(V.gradV ...)` */
  PetscBool    f3_resistive;  /* TRUE: `tau - curl2(B)` (betaf/betae2 term); FALSE: ideal `tau` */
  PetscBool    f3_jre;        /* TRUE: add the edge-averaged runaway current user->jre */
  const char * label;         /* PETSc log-event name */
} MFD_ResidualTerms;
```

| Wrapper | f1_inertia | f3_resistive | f3_jre |
|---|---|---|---|
| `FormIFunction_Vperp_viscosity` (production) | FALSE | TRUE | TRUE |
| `FormIFunction_newequilibrium_Vperp` (relaxation) | TRUE | FALSE | FALSE |

Owner decisions encoded here: the relaxation is **ideal** (no resistive term), for
consistency with previous runs. `jre` is its own switch, separate from
`f3_resistive`. The production `0.0 *` inertia term is kept verbatim.

How it is written:
- Every switch selects whole statements, never a `0/1` factor or an inline ternary,
  so each variant keeps its original expression tree and its FMA contraction.
- f1: `if (f1_inertia) { relaxation r,z statements } else { production r,z statements }`.
- f3: in each of the 5 inner-edge statements, a 4-way `if` on
  (`f3_resistive`, `f3_jre`). (TRUE, TRUE) is production's statement byte for byte.
  (FALSE, FALSE) is the relaxation's `F = tau;` byte for byte. The two mixed
  combinations are the natural sub-expressions. No wrapper uses them, so
  **they are not covered by the t1 regression baseline**.
- `vperp_residual` is `static inline __attribute__((always_inline))`. As a plain
  out-of-line static taking runtime `terms`, GCC shared code across the branches
  and 32 entries of the t1 fixture changed in the last ULP (30 production, 2
  relaxation). Once inlined into each wrapper with constant `terms`, both are
  bit-identical to the pre-merge functions, with the same fmadd/fmsub counts (61
  and 52).
- The production debug-print blocks (`if (user->debug)`) now also run in the
  relaxation. They only print and never write F. Those after the f3 statements
  still print the resistive curl pieces even when `f3_resistive` is FALSE.
- The relaxation's unused `curlBxB` declarations (A3) were dropped.
- `#line` directives keep every PETSc macro line of the merged body at production's
  original `__LINE__`, and every later function at its old numbering.
  `check_line_numbers.py` passes with no mismatches: the relaxation's macro
  lines were deleted, not moved, so none of them survive to be compared.

Verification: `t1_residuals` (both functions, nonzero `jre`) is bit-identical to
the baseline. With the relaxation wrapper's `f3_jre` flipped to TRUE, t1 fails
(118547 relaxation entries differ), which shows the gate is sensitive to this
switch.

The rest of this document is the pre-merge map. Its line numbers refer to the
file before the merge.


Task P6-MAP. A read-only, term-by-term comparison of the two residuals left in
`src/mhd/ts_functions.c`:

| Function | Lines | Registered at | Role |
|---|---|---|---|
| `FormIFunction_Vperp_viscosity` (called **VV** below) | 229-1058 | `src/mhd/mhd.c:316` | production residual, every time step |
| `FormIFunction_newequilibrium_Vperp` (called **NE** below) | 1526-2196 | `src/mhd/ts_functions.c:5595`, inside `FormInitialSolution_psi` | initial-condition relaxation (`ictype` 9/15, `itime == 0`, no binary IC load; see 5557-5566) |

All line numbers are for the current working tree and are **file lines**. The file
contains many `#line` directives (for example `#line 3002` at line 289), so line
numbers in compiler diagnostics, `objdump -l` and debuggers do *not* match file
lines. Do not cross-reference them.

## Summary

After comments and whitespace are removed, a token-level diff (Python `difflib` over
C tokens) gives 30 hunks. 14 of them are `if (user -> debug) { ... }` print blocks
that only VV has. With those removed, 16 hunks remain. They group into **11
differences**. There are **6 in category A** (cosmetic), **1 in B** (auxiliary
input: `jre`), **4 in C** (term differences: f1 r-row, f1 z-row, f3 resistive term,
f3 `jre` source), **0 in D** (coefficients), **0 in E** (boundaries or regions)
and **0 in F** (unclassified).

One further effect, f2, shows no token difference. It follows from the f3 terms
through `ApplyDerivedDivergence`, so it is listed as C5 (derived) and adds no
parameter.

Every material predicate, every boundary-index branch, every coefficient call
outside the f3 resistive term, the f4 and f5 rows, the source-potential block and
the full set of auxiliary-field constructions are token-identical.

The merge itself is mechanical. Three boolean switches plus the log label rebuild
both functions token-for-token, and this was checked by script (see Verification).
It does still touch two physics points: the f1 inertia model and the resistive
term in Ohm's law. VV's header comments (`ts_functions.h:32,35`) document both
choices as intended, so the merge does not need to *decide* them, only keep both.
The open questions are smaller: a dead `0.0 *` term and one stale log message.

## Table

| # | Category | Row/stratum | Vperp_viscosity (VV) | newequilibrium_Vperp (NE) | Lines (VV / NE) |
|---|---|---|---|---|---|
| A1 | A | function name | `FormIFunction_Vperp_viscosity` | `FormIFunction_newequilibrium_Vperp` | 229 / 1526 |
| A2 | A | PETSc log-event label | `"FormIFunction_Vperp_viscosity"` | `"FormIFunction_newequilibrium_Vperp"` | 236 / 1533 |
| A3 | A | declarations | none | `Vec curlBxB, curlBxBLocal` and `PetscScalar ****arrcurlBxB`, declared but never used | — / 1541, 1571 |
| A4 | A | debug output | 14 extra `if (user->debug)` print blocks (only `PetscPrintf`; they never write F) | none (both share the blocks at 477/1769, 487/1779, 1044/2182) | 575-586, 589-591, 594-596, 624-635, 641-643, 649-651, 655-657, 662-664, 669-671, 783-785, 798-800, 813-835, 848-850, 862-882 / — |
| A5 | A | comments | header says f1 = -(J×B + Re⁻¹∇²V), f2 includes -div(curl₂B), f3 includes curl₂B | header says f1 = nᵢ dV/dt - J×B - Re⁻¹∇²V, f2 = -div grad EP, f3 = tau - grad EP | 427-440, 521-527, 770, 892 / 1720-1732, 1813-1819, 2004, 2030 |
| A6 | A | `jre` indexing inside VV | the phi-z edge with `phibtype==0` uses `... + (ez)*jre.stride1]`; the `phibtype!=0` copy writes `+ 0 * jre.stride2` | n/a | 794-795 vs 844-845 (both in VV; same integer index) |
| B1 | B | auxiliary input | `view3d_t jre = user->jre;` | not read | 249 / — |
| C1 | C | f1, plasma vertex (`ivVrmphimzm[0]`, r) | `-(curlB×B)_r + 0.0*niv*(V·GradV1 - V_phi²/r) - Re⁻¹ LapV_r` | `-(curlB×B)_r + niv*dV_r/dt - Re⁻¹ LapV_r` | 535 / 1827 |
| C2 | C | f1, plasma vertex (`ivVrmphimzm[2]`, z) | `-(curlB×B)_z + 0.0*niv*(V·GradV3) - Re⁻¹ LapV_z` | `-(curlB×B)_z + niv*dV_z/dt - Re⁻¹ LapV_z` | 537 / 1829 |
| C3 | C | f3, inner edges (5 statements) | `tau - curl₂(B)` | `tau` | 772-775, 789-792, 803-806, 839-842, 853-856 / 2006, 2010, 2013, 2017, 2020 |
| C4 | C | f3, inner edges (same 5 statements) | `+ avg(jre)` | none | 776-781, 793-796, 807-810, 843-846, 857-860 / — |
| C5 | C (derived) | f2, inner vertices | includes `div(-curl₂B + avg(jre))` | does not | no token diff: 997-1026 / 2135-2164, via `mimetic_operators.c:186,352-353` |

## The C items

### C1, C2: f1 momentum rows on plasma vertices

Both functions use the same vertex predicate: `er > 0 && ez > 0` and all four
adjacent cells satisfy `fabs(dataC-1.5) < 0.7`. Both share the Lorentz term built
from `arrcurlBv`×`arrBv` and the viscous term `-(1.0 / user->Re) * arrLapV`. The φ
component (`ivVrmphimzm[1]`) is the identical constraint `V·B = 0` (539 / 1831).

VV, at 535 (r) and 537 (z):
```c
+ 0.0 * arrniv[..][ivVrmphimzm[0]] * (arrX[..][ivVrmphimzm[0]] * arrGradV1[..][ivVrmphimzm[0]]
    + arrX[..][ivVrmphimzm[1]] * arrGradV1[..][ivVrmphimzm[1]] + arrX[..][ivVrmphimzm[2]] * arrGradV1[..][ivVrmphimzm[2]]
    - arrX[..][ivVrmphimzm[1]] * arrX[..][ivVrmphimzm[1]] / arrCoorda[..][icrmphimzm[0]] )
/* z row: same shape with arrGradV3 and no V_phi^2/r term */
```
NE, at 1827 (r) and 1829 (z):
```c
+ arrniv[..][ivVrmphimzm[0]] * arrXdot[..][ivVrmphimzm[0]]      /* z row: [2] */
```

What the code shows: VV's f1 has no time derivative. On plasma vertices it is the
algebraic balance `J×B + Re⁻¹∇²V = 0`. Its advective term
`n (V·∇)V` (with the `-V_φ²/r` part on the r row) is multiplied by the literal
`0.0`. NE uses `n dV/dt - J×B - Re⁻¹∇²V`, so it has inertia and no advection.
`ts_functions.h:32` and `:35` describe exactly these two forms, so the split looks
deliberate.

Side effect: in VV, `niv`, `GradV1` and `GradV3` are consumed *only* through the
`0.0 *` term. In NE, `GradV1` and `GradV3` are never read. `GradV2` is built and
destroyed in both and read by neither (398-404 / 1693-1699). Both functions build
exactly the same auxiliary fields.

What a merge must parameterize: which of the two sub-expressions is added between
the Lorentz and viscous terms, on two statements (r and z). VV's expression must
be kept literally, `0.0 *` included. Without `-ffast-math` the compiler cannot fold
`0.0 * x`. The term contributes a signed zero, or NaN when `x` is Inf or NaN, so
deleting it is not guaranteed bit-exact.

### C3: f3 resistive term (derived curl of B)

These are the five inner-edge assignments in the second cell loop: the r-z edge,
plus the φ-z and r-φ edges in each `phibtype` branch. The guard conditions
(`!(er == 0 || ez == 0)` and the rest) are identical in both functions. Only the
right-hand side differs.

VV, for example at 772-775:
```c
arrF[..][ivErmzm] = arrX[..][ivErmzm] - ((arrX[..][ivBrm] * betaf(.., LEFT, user) / surface(.., LEFT, user) - ...
      ) * rmzmedgelength / betae2(er, ephi, ez, BACK_LEFT, user))
```
NE, at 2006:
```c
arrF[..][ivErmzm] = arrX[..][ivErmzm];
```

What the code shows: VV computes `tau - M_e⁻¹ Curlᵀ M_f B`. The face weight is
`betaf` and the edge weight is `betae2`, which `alphaec2` builds from `user->eta0`
and `user->etawall` (`mass_matrix_coefficients.c:312-314, 652-655`). Its header
comment calls this `derived_mimetic_curl2(B)` (434-435). NE has no such term, so
its Ohm's law is `tau = grad EP - R_ve(V×B)` with no resistive part. This is not a
coefficient substitution (category D): the term is absent entirely. `betae2`
appears nowhere else in either function except VV's debug prints.

What a merge must parameterize: whether the five statements subtract the curl
expression. The five expressions themselves are identical between the functions'
shared structure, because only VV has them.

### C4: f3 runaway-electron current source

These are the same five statements, as a trailing addend in VV:

| Edge | VV lines | Expression |
|---|---|---|
| r-z edge (`BACK_LEFT`, φ-directed) | 776-781 | `+ .25 * (jre[er,ez,1] + jre[er-1,ez,1] + jre[er,ez-1,1] + jre[er-1,ez-1,1])` |
| φ-z edge (`BACK_DOWN`, r-directed) | 793-796, 843-846 | `+ .5 * (jre[er,ez,0] + jre[er,ez-1,0])` |
| r-φ edge (`DOWN_LEFT`, z-directed) | 807-810, 857-860 | `+ .5 * (jre[er-1,ez,2] + jre[er,ez,2])` |

NE has nothing here. What the code shows: the cell-centred, axisymmetric
`(NR, NZ, 3)` kinetic current is averaged onto each edge type, taking the
component along that edge. The scaling (`eta_a3VaB0 / dtCD`) is applied on the C++
side (`src/tasks/PredictorCorrector.cpp:52-57`), which matches the comment at 770,
`+ (eta/(V_A*B_0)) j_RE`.

What a merge must parameterize: whether the addend is present. Per the owner
requirement and open question Q2, this must be an explicit switch and must be ON
for production.

### C5: f2 (derived; no token difference)

Both functions then form `Fcopy = F - tau` on inner edges (997-1021 / 2135-2159)
and call `ApplyDerivedDivergence(ts, Fcopy, F, user)` (1026 / 2164). That call
overwrites the EP rows on inner vertices (`mimetic_operators.c:352-353`). VV's
edge residual carries `-curl₂B + jre`, so VV's f2 also contains
`div(-curl₂B + avg jre)`, as its header says (432). NE's f2 does not. No
separate parameter is needed, because C3 and C4 carry it.

## Proposed parameters

This is the minimal set for one body, with the debug-free VV body as the base.

1. **`f1_inertia` (bool), on 2 statements (535, 537).** It selects between
   expression P, VV's `+ 0.0 * arrniv[..][ivVrmphimzm[0]] * (V·GradV1 - V_φ²/r)`
   on the r row and `+ 0.0 * arrniv[..][ivVrmphimzm[0]] * (V·GradV3)` on the z
   row, and expression Q, NE's `+ arrniv[..][ivVrmphimzm[0]] * arrXdot[..][ivVrmphimzm[c]]`
   with c = 0 or 2. VV = false, NE = true.
2. **`f3_resistive` (bool), on 5 statements.** It selects between P, VV's
   `arrX[tau] - ((Σ ±arrX[B]·betaf/surface) * edgelength / betae2(edge))`, and Q,
   NE's `arrX[tau]`. VV = true, NE = false.
3. **`f3_jre` (bool), on the same 5 statements plus the `view3d_t jre = user->jre`
   local.** It selects between P, VV's `+ .25*(4-point avg of jre comp 1)`,
   `+ .5*(2-point avg comp 0)` or `+ .5*(2-point avg comp 2)`, and Q, NE's empty
   addend. VV = true, NE = false.
4. **`label` (const char\*).** It selects between `"FormIFunction_Vperp_viscosity"`
   and `"FormIFunction_newequilibrium_Vperp"`, as `initialize_ep` already does
   (1089, 1512-1519).

Switches 2 and 3 are always set together in both existing functions, so strictly
the minimum is 1, 2+3 combined, and 4. Keep 2 and 3 separate anyway: Q2 requires
`jre` to be its own explicit parameter.

Not parameters (cosmetic): A3, the unused `curlBxB` declarations, can be dropped.
A4, the debug blocks, only print. The merged body can keep VV's blocks, ideally
gating the ones that print curl terms (813-835, 862-882) on `f3_resistive`. Nothing
in F depends on them.

Implementation caution, not a classification: the build is `-O2 -g -DNDEBUG` on
aarch64 GCC, with the default `-ffp-contract=fast`. Both functions' object code
contains fused multiply-add instructions (61 and 52 `fmadd`/`fmsub`-family
instructions). An inline ternary inside the expression, such as
`a + (flag ? niv*xdot : ...)`, can change which operations get fused, so
`niv*xdot` might no longer fuse into the surrounding add. Write each variant as a
complete statement under an `if`, so each branch has the original expression
tree. Then confirm with `tests/regression/t1_residuals.c`, which covers both
functions with a nonzero `jre`.

## Questions for the owner

1. **The `0.0 *` advective term in production f1 (C1/C2).** Is it a deliberate
   switch-off of `n (V·∇)V`, or a temporary disable? If deliberate, may it be
   deleted, together with `GradV1`, `GradV3` and `niv` in VV? Doing so is not
   bit-exact in signed-zero and non-finite cases, but it saves three vertex fields
   per residual evaluation. (`GradV2` is dead in both functions regardless.)
2. **No resistive term in relaxation f3 (C3).** NE's Ohm's law has no
   `curl₂B / betae2` term, yet after the relaxation solve `FormInitialSolution_psi`
   logs `"Saving solution after relaxation with resitivity %le"` with
   `user->etaplasma` (5847). Is ideal Ohm's law intended during relaxation, which
   would make the message stale, or was the resistive term meant to be there?
3. **`jre` in relaxation (C4).** Is its absence in NE intended? In every current
   call path it makes no practical difference: `Jre_mhd` is zero-filled
   (`src/kinetic/kinetic.cpp:161-163`) and wired to `user->jre` (244) just before
   `mhd_initialize` (246), which runs the relaxation. Turning `jre` on for NE would
   therefore change F only by signed zeros. The answer decides whether the merged
   relaxation call passes `f3_jre = true` (uniform physics) or `false` (bit-exact
   with today).

## Verification

`/tmp/verify.py` (not committed) does the following:
- tokenizes both bodies with comments and `#line` lines stripped;
- removes every balanced `if ( user -> debug ) { ... }` block, which takes VV from
  17 blocks to 3 and NE stays at 3;
- applies the parameters to VV's token stream: renames under `label`, inserts the
  A3 declarations, deletes the `jre` local (`f3_jre`), replaces the two f1
  sub-expressions (`f1_inertia`), and in each of the 5 f3 statements deletes the
  `- ( ( ... betae2 ... ) )` span (`f3_resistive`) and the following
  `+ .25/.5 * ( ...jre... )` span (`f3_jre`), asserting each span's shape.

Result: `f3 statements rewritten: 5` and
`IDENTICAL to newequilibrium_Vperp (debug-stripped): True`. With the VV settings
the body is VV's own token stream, debug blocks aside.

Coverage of the 16 non-debug hunks:
- hunks 1-2 are `label`;
- hunks 3 and 5 are A3;
- hunk 4 is `f3_jre`;
- hunks 6-8 and 9-11 are `f1_inertia` (r and z);
- hunks 21, 23, 25, 27 and 29 each contain a `f3_resistive` span immediately
  followed by a `f3_jre` span. `difflib` reports each pair as one hunk, but each
  splits cleanly at the `) ) + .25|.5 *` boundary.

No hunk needs a parameter beyond these, and none is covered twice.
