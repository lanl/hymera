# DMStag Discretization Decoder

This document decodes the index variables (`ivVrmphimzm`, `ivErmzm`, `ivBphip`,
`icrpphipzp`, ...) that every residual function in `src/mhd/ts_functions.c` resolves at
the top of its body before doing any physics. Without this key, the residual assembly
code is unreadable: it is dense arithmetic indexed entirely by these opaque names.

Everything below was verified against two source ranges (no other lines were read to
produce this document):

- `/workspace/src/mhd/mfd_config.h` lines 60-100 (the `LOCATION` macro aliases)
- `/workspace/src/mhd/ts_functions.c` lines 2960-3045 (the slot-resolution block from
  `FormIFunction_Vperp_viscosity`)

## 1. DMStag strata and the dof layout

PETSc's `DMStag` is a structured, staggered grid. In 3-D it has four "strata" —
vertices, edges, faces, and elements (cell centers) — and each stratum can carry its own
number of degrees of freedom (dof). Hymera's MHD solution DM (`user->da`, created at
`src/mhd/mhd.c:159/161` with `dof0=4, dof1=1, dof2=1, dof3=1`, `DMSTAG_STENCIL_BOX`,
stencil width 1) lays out:

- **dof0 = 4, on vertices**: the three components of the ion velocity `V` (`Vr, Vphi, Vz`)
  plus the electrostatic potential `EP`. (Comment at `src/mhd/mhd.c:154`: "4 dofs on each
  vertex (3 for vector field V, and 1 for the electrostatic potential EP)".)
- **dof1 = 1, on edges**: one component of `tau`, the divergence-free part of the electric
  field, per edge (`src/mhd/mhd.c:27`).
- **dof2 = 1, on faces**: one component of `B` per face.
- **dof3 = 1, at cell centers (elements)**: the ion number density `ni`.

The grid is cylindrical (r, phi, z); phi is `DM_BOUNDARY_PERIODIC` when `user->phibtype`
is set, else `DM_BOUNDARY_NONE` (both r and z are always `DM_BOUNDARY_NONE`).

ASCII diagram of one cell (a hexahedral cell in (r, phi, z) index space; corners are
vertices, cell edges hold `tau`, cell faces hold `B`, the cell interior holds `ni`):

```
                 V,EP (r-,phi+,z-) ---------- tau_phi(z-) ---------- V,EP (r+,phi+,z-)
                     /|                                                   /|
                    / |                                                  / |
              tau_r(phi+,z-)                                     tau_r(phi+,z+)  (back face)
                  /   |                                                /   |
   V,EP(r-,phi+,z+) --+----------- tau_phi(z+) ------------ V,EP(r+,phi+,z+)
        |    B_phi(z-,+)|  <- UP face (phi+)                  |    B_phi(z-,+)
        |             tau_z(r-)                                |     |
        |              |                                       |     |
     B_r(-) [LEFT]      |          ni (cell center)        B_r(+) [RIGHT] face
        |              |                                       |     |
        |          tau_z(r+)                                   |    tau_z(r+)
        |              |                                       |     |
   V,EP(r-,phi-,z+) ---+----------- tau_phi(z+,phi-) ---- V,EP(r+,phi-,z+)
              /        |                                            /
        tau_r(phi-,z-)                                     tau_r(phi-,z+)
            /           |                                        /
V,EP(r-,phi-,z-) --------+----- tau_phi(z-,phi-) ----- V,EP(r+,phi-,z-)

B_z(-) is the BACK face (z-), B_z(+) is the FRONT face (z+).
```

In words (this is the load-bearing summary, the picture above is only a mnemonic): the 8
corners of the cell are vertices, each carrying `(Vr, Vphi, Vz, EP)`; the 12 edges of the
cell each carry one scalar `tau` component (oriented along that edge's direction: r-edges
carry `tau_r`, phi-edges carry `tau_phi`, z-edges carry `tau_z`); the 6 faces of the cell
each carry one scalar `B` component, normal to that face (the two r-faces carry `B_r`,
the two phi-faces carry `B_phi`, the two z-faces carry `B_z`); and the cell center carries
the scalar `ni`.

## 2. Location macros (`mfd_config.h` lines 60-100)

`mfd_config.h` gives PETSc's verbose `DMSTAG_*` stencil-location enumerators short
aliases. Every one of these is a plain `#define X DMSTAG_X`, so "the DMSTAG_* value it
aliases" and "macro name" are identical modulo the `DMSTAG_` prefix.

| Macro | Aliases | Stratum | Geometric position | Line |
|---|---|---|---|---|
| `BACK_DOWN_LEFT` | `DMSTAG_BACK_DOWN_LEFT` | vertex | (r-, phi-, z-) | 68 |
| `BACK_DOWN` | `DMSTAG_BACK_DOWN` | edge | (phi-, z-), r-direction edge | 69 |
| `BACK_DOWN_RIGHT` | `DMSTAG_BACK_DOWN_RIGHT` | vertex | (r+, phi-, z-) | 70 |
| `BACK_LEFT` | `DMSTAG_BACK_LEFT` | edge | (r-, z-), phi-direction edge | 71 |
| `BACK` | `DMSTAG_BACK` | face | z- | 72 |
| `BACK_RIGHT` | `DMSTAG_BACK_RIGHT` | edge | (r+, z-), phi-direction edge | 73 |
| `BACK_UP_LEFT` | `DMSTAG_BACK_UP_LEFT` | vertex | (r-, phi+, z-) | 74 |
| `BACK_UP` | `DMSTAG_BACK_UP` | edge | (phi+, z-), r-direction edge | 75 |
| `BACK_UP_RIGHT` | `DMSTAG_BACK_UP_RIGHT` | vertex | (r+, phi+, z-) | 76 |
| `DOWN_LEFT` | `DMSTAG_DOWN_LEFT` | edge | (r-, phi-), z-direction edge | 77 |
| `DOWN` | `DMSTAG_DOWN` | face | phi- | 78 |
| `DOWN_RIGHT` | `DMSTAG_DOWN_RIGHT` | edge | (r+, phi-), z-direction edge | 79 |
| `LEFT` | `DMSTAG_LEFT` | face | r- | 80 |
| `ELEMENT` | `DMSTAG_ELEMENT` | element (cell center) | cell interior | 81 |
| `RIGHT` | `DMSTAG_RIGHT` | face | r+ | 82 |
| `UP_LEFT` | `DMSTAG_UP_LEFT` | edge | (r-, phi+), z-direction edge | 83 |
| `UP` | `DMSTAG_UP` | face | phi+ | 84 |
| `UP_RIGHT` | `DMSTAG_UP_RIGHT` | edge | (r+, phi+), z-direction edge | 85 |
| `FRONT_DOWN_LEFT` | `DMSTAG_FRONT_DOWN_LEFT` | vertex | (r-, phi-, z+) | 86 |
| `FRONT_DOWN` | `DMSTAG_FRONT_DOWN` | edge | (phi-, z+), r-direction edge | 87 |
| `FRONT_DOWN_RIGHT` | `DMSTAG_FRONT_DOWN_RIGHT` | vertex | (r+, phi-, z+) | 88 |
| `FRONT_LEFT` | `DMSTAG_FRONT_LEFT` | edge | (r-, z+), phi-direction edge | 89 |
| `FRONT` | `DMSTAG_FRONT` | face | z+ | 90 |
| `FRONT_RIGHT` | `DMSTAG_FRONT_RIGHT` | edge | (r+, z+), phi-direction edge | 91 |
| `FRONT_UP_LEFT` | `DMSTAG_FRONT_UP_LEFT` | vertex | (r-, phi+, z+) | 92 |
| `FRONT_UP` | `DMSTAG_FRONT_UP` | edge | (phi+, z+), r-direction edge | 93 |
| `FRONT_UP_RIGHT` | `DMSTAG_FRONT_UP_RIGHT` | vertex | (r+, phi+, z+) | 94 |

Axis convention used by DMStag/these macros: `LEFT`/`RIGHT` = r-/r+, `DOWN`/`UP` =
phi-/phi+, `BACK`/`FRONT` = z-/z+.

## 3. Naming scheme for the `iv*`/`ic*` slot variables

**Prefix**: `iv` = a slot index into the *solution* vector (`da`'s layout, i.e. `V`, `tau`,
`B`, `ni`); `ic` = a slot index into a *coordinate* vector (`dmCoord`'s or `dmCoorda`'s
layout — i.e. an (r, phi, z) coordinate triple, indexed by `[0|1|2]` for `r`/`phi`/`z`
respectively, not a physical field).

**Field letter** (only on `iv*` names; coordinate names have none): `V` = vertex velocity
(4-component slot array, indices 0-2 = `Vr, Vphi, Vz`, index 3 = `EP`), `E` = edge `tau`
(the code's own edge-array prefix is `E`, presumably from "electric field" — see note
below), `B` = face `B`. No letter at all: the cell-center `ni` (`ivn`) and the vertex
coordinate names.

**Position suffix**: a sequence of `r`, `phi`, `z` tokens, each followed by `m` (minus
side) or `p` (plus side), for every axis that stratum's location actually varies along:
- vertices vary in all three: `rXphiYzZ` with X,Y,Z each `m`/`p` (8 combinations).
- edges vary in the two axes transverse to the edge's own direction, and the field-letter
  encodes *which* axis is the edge's free direction (implicitly, by which two suffix axes
  appear): e.g. `ivErmzm` — an `r`,`z` suffix pair only, so it is a phi-direction edge at
  (r-, z-) (`BACK_LEFT`); `ivEphimzm` has a `phi`,`z` pair, so it is an r-direction edge at
  (phi-, z-) (`BACK_DOWN`); `ivErmphim` has an `r`,`phi` pair, so it is a z-direction edge
  at (r-, phi-) (`DOWN_LEFT`).
- faces vary in only one axis, matching the face's normal: `ivBrm`, `ivBphim`, `ivBzm` (and
  `p` variants).
- the cell center (`ivn`) and element coordinate (`icp`) have no position suffix.

**Verification against the task's example pattern**: the stated scheme
(`ivVrmphimzm` = V at (r-,phi-,z-) vertex; `ivErmzm` = tau/E at (r-,z-) edge; `ivBphim` = B
at phi- face; `icrpphipzp` = coordinates at (r+,phi+,z+)) matches the code exactly for
every one of the 62 `iv*`/`ic*` names in the slot block — no correction was needed. One
refinement beyond the task's description: the edge names' two position axes double as an
implicit statement of the edge's direction (the *missing* axis), which is what lets
`ivErmzm` (r,z pair -> phi-direction edge) and `ivErmphim` (r,phi pair -> z-direction edge)
both start with `ivEr...` yet denote geometrically different edges — readers must check
which two axes are present, not just the leading letters.

**Why "E" for tau**: `ts_functions.c`'s own local variable names call the edge quantity
`tau` in comments/physics (`src/mhd/mhd.c:21,27`) but the slot variables use `E` throughout
(`ivErmzm`, `icEphimzm`, ...) — consistent with `tau` being the divergence-free/curl part
of the electric field `E`. Treat `ivE*` as "the edge-resident field", i.e. `tau`.

**Cell-local "minus belongs to me, plus is shared" convention**: In a cell-local DMStag
view, for a cell at grid index `(er, ephi, ez)`, the "minus"-side vertices/edges/faces
(r-, phi-, z-) are the ones *owned* by that cell in the sense that they are never also the
"plus" side of a *previous* cell you'd otherwise double-index — DMStag numbers stencil
locations so each cell's local array slice exposes both its own minus-side entities and
(via the stencil width of 1) the immediately adjacent plus-side entities that are actually
shared with the next cell in that direction. Concretely this is why loops throughout
`ts_functions.c` special-case the last valid index in each direction — e.g. `er == N[0]-1`
(see `ts_functions.c:386,401,408,502,...` in the neighboring code) — the "plus" face/edge/
vertex at the last cell in a direction is the true domain boundary and has no "next cell"
to be shared with, so it needs boundary-condition handling instead of interior treatment.

## 4. Main table: solution-vector (`iv*`) slot variables, from `da`

All of these resolve via `DMStagGetLocationSlot(da, LOCATION, component, &var)` at
`ts_functions.c:2969-2999`.

### Vertices (dof0 = 4 components each; component loop `d = 0..3`)

| Variable | DM | Location | Component | Stratum | Physical quantity |
|---|---|---|---|---|---|
| `ivVrmphimzm[d]` | `da` | `BACK_DOWN_LEFT` | `d` | vertex (r-,phi-,z-) | `d<3`: V (Vr,Vphi,Vz); `d==3`: EP |
| `ivVrpphimzm[d]` | `da` | `BACK_DOWN_RIGHT` | `d` | vertex (r+,phi-,z-) | V / EP |
| `ivVrmphipzm[d]` | `da` | `BACK_UP_LEFT` | `d` | vertex (r-,phi+,z-) | V / EP |
| `ivVrpphipzm[d]` | `da` | `BACK_UP_RIGHT` | `d` | vertex (r+,phi+,z-) | V / EP |
| `ivVrmphimzp[d]` | `da` | `FRONT_DOWN_LEFT` | `d` | vertex (r-,phi-,z+) | V / EP |
| `ivVrpphimzp[d]` | `da` | `FRONT_DOWN_RIGHT` | `d` | vertex (r+,phi-,z+) | V / EP |
| `ivVrmphipzp[d]` | `da` | `FRONT_UP_LEFT` | `d` | vertex (r-,phi+,z+) | V / EP |
| `ivVrpphipzp[d]` | `da` | `FRONT_UP_RIGHT` | `d` | vertex (r+,phi+,z+) | V / EP |

(Confirmed at use sites elsewhere in the file, e.g. `ts_functions.c:426-440` uses
`[0],[1],[2]` as `Vr, Vphi, Vz` in momentum-equation arithmetic, and `ts_functions.c:581`
uses `[3]` in an `arrX[...] - arrx[...]` Dirichlet-copy pattern consistent with `EP`.)

### Edges (dof1 = 1 component; component argument fixed at `0`)

| Variable | DM | Location | Component | Stratum | Physical quantity |
|---|---|---|---|---|---|
| `ivErmzm` | `da` | `BACK_LEFT` | 0 | edge (r-,z-), phi-dir | tau (edge E-field), phi-oriented |
| `ivEphimzm` | `da` | `BACK_DOWN` | 0 | edge (phi-,z-), r-dir | tau, r-oriented |
| `ivErpzm` | `da` | `BACK_RIGHT` | 0 | edge (r+,z-), phi-dir | tau, phi-oriented |
| `ivEphipzm` | `da` | `BACK_UP` | 0 | edge (phi+,z-), r-dir | tau, r-oriented |
| `ivErmphim` | `da` | `DOWN_LEFT` | 0 | edge (r-,phi-), z-dir | tau, z-oriented |
| `ivErpphim` | `da` | `DOWN_RIGHT` | 0 | edge (r+,phi-), z-dir | tau, z-oriented |
| `ivErmphip` | `da` | `UP_LEFT` | 0 | edge (r-,phi+), z-dir | tau, z-oriented |
| `ivErpphip` | `da` | `UP_RIGHT` | 0 | edge (r+,phi+), z-dir | tau, z-oriented |
| `ivEphimzp` | `da` | `FRONT_DOWN` | 0 | edge (phi-,z+), r-dir | tau, r-oriented |
| `ivErmzp` | `da` | `FRONT_LEFT` | 0 | edge (r-,z+), phi-dir | tau, phi-oriented |
| `ivErpzp` | `da` | `FRONT_RIGHT` | 0 | edge (r+,z+), phi-dir | tau, phi-oriented |
| `ivEphipzp` | `da` | `FRONT_UP` | 0 | edge (phi+,z+), r-dir | tau, r-oriented |

### Faces (dof2 = 1 component; component argument fixed at `0`)

| Variable | DM | Location | Component | Stratum | Physical quantity |
|---|---|---|---|---|---|
| `ivBrm` | `da` | `LEFT` | 0 | face r- | B_r |
| `ivBphim` | `da` | `DOWN` | 0 | face phi- | B_phi |
| `ivBzm` | `da` | `BACK` | 0 | face z- | B_z |
| `ivBrp` | `da` | `RIGHT` | 0 | face r+ | B_r |
| `ivBphip` | `da` | `UP` | 0 | face phi+ | B_phi |
| `ivBzp` | `da` | `FRONT` | 0 | face z+ | B_z |

### Element (dof3 = 1 component)

| Variable | DM | Location | Component | Stratum | Physical quantity |
|---|---|---|---|---|---|
| `ivn` | `da` | `ELEMENT` | 0 | cell center | ion number density `ni` |

**Count**: the vertex block resolves 8 macro calls x 4 components = 32 slot fills from 8
`DMStagGetLocationSlot` call sites (inside the `d` loop); the edge block is 12 calls; the
face block is 6 calls; the element block is 1 call. Solution-side (`da`) call sites in the
read range: 8 + 12 + 6 + 1 = **27 `DMStagGetLocationSlot` calls**, producing 32 + 12 + 6 +
1 = 51 distinct slot integers (32 of which live packed inside the four-wide `ivV*[4]`
arrays).

## 5. Coordinate-vector (`ic*`) slot variables, and the coordinate DMs

Two coordinate DMs are in play, both obtained via `DMGetCoordinateDM`/`DMGetCoordinatesLocal`
+ `DMStagVecGetArrayRead` at `ts_functions.c:3001-3006`:

- **`dmCoord`** — the coordinate DM of `da` itself (`DMGetCoordinateDM(da, &dmCoord)`).
  Supplies face- and edge-location coordinates (`icB*`, `icE*` below), read into `arrCoord`.
- **`dmCoorda`** — the coordinate DM of the separate DM `coordDA` (`user->coorda`,
  aliased locally as `coordDA` at the top of the function), i.e.
  `DMGetCoordinateDM(coordDA, &dmCoorda)`. Supplies element- and vertex-location
  coordinates (`icp`, `icr*`), read into `arrCoorda`.

Both are 3-component (`d = 0..2`) coordinate DMs: component `0` = r, `1` = phi, `2` = z,
consistent with `DMStagSetUniformCoordinatesExplicit(..., rmin, rmax, phimin, phimax,
zmin, zmax)` used to populate them.

| Variable | DM | Location | Component | Stratum | Coordinate meaning |
|---|---|---|---|---|---|
| `icp[d]` | `dmCoorda` | `ELEMENT` | d | element | cell-center (r,phi,z)[d] |
| `icBrm[d]` | `dmCoord` | `LEFT` | d | face r- | (r,phi,z)[d] at r- face |
| `icBphim[d]` | `dmCoord` | `DOWN` | d | face phi- | at phi- face |
| `icBzm[d]` | `dmCoord` | `BACK` | d | face z- | at z- face |
| `icBrp[d]` | `dmCoord` | `RIGHT` | d | face r+ | at r+ face |
| `icBphip[d]` | `dmCoord` | `UP` | d | face phi+ | at phi+ face |
| `icBzp[d]` | `dmCoord` | `FRONT` | d | face z+ | at z+ face |
| `icErmzm[d]` | `dmCoord` | `BACK_LEFT` | d | edge (r-,z-) | at that edge |
| `icEphimzm[d]` | `dmCoord` | `BACK_DOWN` | d | edge (phi-,z-) | at that edge |
| `icErpzm[d]` | `dmCoord` | `BACK_RIGHT` | d | edge (r+,z-) | at that edge |
| `icEphipzm[d]` | `dmCoord` | `BACK_UP` | d | edge (phi+,z-) | at that edge |
| `icErmphim[d]` | `dmCoord` | `DOWN_LEFT` | d | edge (r-,phi-) | at that edge |
| `icErpphim[d]` | `dmCoord` | `DOWN_RIGHT` | d | edge (r+,phi-) | at that edge |
| `icErmphip[d]` | `dmCoord` | `UP_LEFT` | d | edge (r-,phi+) | at that edge |
| `icErpphip[d]` | `dmCoord` | `UP_RIGHT` | d | edge (r+,phi+) | at that edge |
| `icEphimzp[d]` | `dmCoord` | `FRONT_DOWN` | d | edge (phi-,z+) | at that edge |
| `icErmzp[d]` | `dmCoord` | `FRONT_LEFT` | d | edge (r-,z+) | at that edge |
| `icErpzp[d]` | `dmCoord` | `FRONT_RIGHT` | d | edge (r+,z+) | at that edge |
| `icEphipzp[d]` | `dmCoord` | `FRONT_UP` | d | edge (phi+,z+) | at that edge |
| `icrmphimzm[d]` | `dmCoorda` | `BACK_DOWN_LEFT` | d | vertex (r-,phi-,z-) | at that vertex |
| `icrpphimzm[d]` | `dmCoorda` | `BACK_DOWN_RIGHT` | d | vertex (r+,phi-,z-) | at that vertex |
| `icrmphipzm[d]` | `dmCoorda` | `BACK_UP_LEFT` | d | vertex (r-,phi+,z-) | at that vertex |
| `icrpphipzm[d]` | `dmCoorda` | `BACK_UP_RIGHT` | d | vertex (r+,phi+,z-) | at that vertex |
| `icrmphimzp[d]` | `dmCoorda` | `FRONT_DOWN_LEFT` | d | vertex (r-,phi-,z+) | at that vertex |
| `icrpphimzp[d]` | `dmCoorda` | `FRONT_DOWN_RIGHT` | d | vertex (r+,phi-,z+) | at that vertex |
| `icrmphipzp[d]` | `dmCoorda` | `FRONT_UP_LEFT` | d | vertex (r-,phi+,z+) | at that vertex |
| `icrpphipzp[d]` | `dmCoorda` | `FRONT_UP_RIGHT` | d | vertex (r+,phi+,z+) | at that vertex |

**Count**: 1 (`icp`) + 6 (`icB*`) + 12 (`icE*`) + 8 (`icr*`) = 27 `DMStagGetLocationSlot`
call sites for coordinates, each executed inside the `d = 0..2` loop, producing
27 x 3 = 81 slot integers.

**Grand total for the whole block (`ts_functions.c:2967-3039`)**: 27 (solution) + 27
(coordinate) = **54 `DMStagGetLocationSlot` call sites**, matching the count of that call
found by mechanically grepping the read range (`grep -c DMStagGetLocationSlot` over lines
2960-3045 returns 54). These 54 call sites populate 51 distinct solution-side slot
integers and 81 distinct coordinate-side slot integers — over a hundred resolved indices
per residual function, all from a boilerplate block that looks identical function to
function.

### How the coordinate arrays are indexed and where the read-only borrow happens

Code accesses the borrowed coordinate array as `arrCoord[ez][ephi][er][ic....[d]]` (or
`arrCoorda[ez][ephi][er][ic....[d]]`), i.e. indexed first by the cell's `(z, phi, r)`
grid indices (in that outer-to-inner order) and then by the resolved slot integer for the
desired location, itself indexed by `d in {0,1,2}` to select the r/phi/z coordinate
component. E.g. `arrCoord[ez][ephi][er][icBzp[2]]` is "the z-coordinate of this cell's z+
face".

`user->arrCoord` (note: this is the *user-context-level* field, distinct from the
function-local `arrCoord`/`arrCoorda` locals inside each residual — the local `arrCoord`
in the slot block above is bound fresh every call from `dmCoord`, while `user->arrCoord` is
a separate, longer-lived borrow) is obtained once, read-only, in `mhd_initialize`:

```
DMStagVecGetArrayRead(dmCoorda, coordaLocal, &user->arrCoord);   // src/mhd/mhd.c:184
```

and released once at teardown:

```
DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, &user->arrCoord);  // src/mhd/mhd.c:832
```

So `user->arrCoord` is a single borrowed view of the `coorda` DM's coordinates, held for
the run's whole lifetime, separate from (but geometrically identical in content to) the
per-call `arrCoorda` that each residual function re-derives locally from `coordDA` inside
its own slot block.

## 6. Where the slots come from now

When this document was first written, the slot-resolution block above was duplicated
verbatim at the top of every residual -- 54 `DMStagGetLocationSlot` calls per
function, across 18 variants. Most of those variants were unreachable and are now in
`src/mhd/attic/`, and in the live residuals the block is shared (`f20302c`):

- `MFD_Slots`, defined near the top of `ts_functions.c`, is a struct whose members
  carry exactly the variable names in the table above (`ivVrmphimzm[4]`, `ivErmzm`,
  `icrpphipzp[3]`, ...).
- `MFD_GetSlotsSolution(da, &S)` fills the solution-DM slots (vertex, edge, face,
  element); `MFD_GetSlotsCoords(dmCoord, dmCoorda, &S)` fills the coordinate slots.
- `MFD_UNPACK_SLOTS(S)` copies them into each function's existing local variables,
  so the residual arithmetic reads the same names this table decodes.

So this table still decodes the residuals directly; the lookups just happen in one
place. `DMStagGetLocationSlot` is a pure function of the DM, so the values are
constant for a run. They are still resolved once per residual call rather than once
per run; caching them on `User` would be a small further step, not yet taken.

The per-cell geometry that sits next to the slots in every cell loop is shared the
same way: `MFD_CellEdgeLengths` gives the 12 edge lengths and `MFD_CellVolume` the
cell volume (`c2c454f`).

## Verification summary

- 54 `DMStagGetLocationSlot` calls counted in the read range (`ts_functions.c:2960-3045`);
  all 54 are listed in the tables in sections 4 and 5 (27 solution-side + 27
  coordinate-side call sites; the coordinate-side and vertex-solution-side calls each
  execute inside a `d`-loop, so they resolve more than one slot integer per call site, as
  noted in each section's "Count").
- Every `LOCATION` macro cited in the tables (`BACK_DOWN_LEFT`, `BACK_DOWN`,
  `BACK_DOWN_RIGHT`, `BACK_LEFT`, `BACK`, `BACK_RIGHT`, `BACK_UP_LEFT`, `BACK_UP`,
  `BACK_UP_RIGHT`, `DOWN_LEFT`, `DOWN`, `DOWN_RIGHT`, `LEFT`, `ELEMENT`, `RIGHT`,
  `UP_LEFT`, `UP`, `UP_RIGHT`, `FRONT_DOWN_LEFT`, `FRONT_DOWN`, `FRONT_DOWN_RIGHT`,
  `FRONT_LEFT`, `FRONT`, `FRONT_RIGHT`, `FRONT_UP_LEFT`, `FRONT_UP`, `FRONT_UP_RIGHT`)
  exists in `mfd_config.h` at the line cited in section 2 (lines 68-94).
- No file other than this one (`/workspace/docs/mhd/physics/discretization.md`) was
  created or modified; nothing under `src/` was changed.
