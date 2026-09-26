# The material region model

`user->dataC` is the field that tells every mass-matrix coefficient function
which material a cell is made of. It is not a resistivity field, despite
being loaded from a file named `veceta`, and it is not a level-set function,
despite the header comment in `src/mhd/mfd_config.h:155` calling it one. It
is a five-valued integer tag, stored as `PetscReal`, compared against its
five values with a floating-point tolerance.

## Where it comes from

`ReadInitialData` (`src/mhd/ts_functions.c:18461`) is a generic ASCII-array
reader: it reads a leading count, then that many `%lf` tokens, into a freshly
`PetscMalloc1`'d buffer. `dataC` is one of three arrays it is used to fill;
the call site is in `src/mhd/mhd.c`, in a block starting at line 96 with a
comment that documents the three files' layouts:

```
src/mhd/mhd.c:96-103
  /* Grid data from inputs/mhd. The element count of each file is a fixed
   * function of the grid, so check it rather than trusting the count in the
   * file header:
   *   veceta -> dataC    cell-centred material tags   Nr     * Nphi     * Nz
   *   vecpsi -> datapsi  poloidal flux, r-z edge      (Nr+1) * Nphi     * (Nz+1)
   *   vecg   -> datag    R*B_phi, phi-face            Nr     * (Nphi+1) * Nz
```

The `dataC` load itself is `src/mhd/mhd.c:108-109`:

```c
PetscSNPrintf(filename, sizeof(filename), "%s/veceta_grid%.3Dx%.2Dx%.3D.txt", user->input_folder, user->Nr, user->Nphi, user->Nz);
PetscCall(ReadInitialData( & (user->dataC), & (user->numC), filename));
```

immediately followed (`mhd.c:110-115`) by a `PetscCheck` that `numC == Nr *
Nphi * Nz`, i.e. the file's own element count is trusted, but its *shape* is
verified against the grid the run was configured with — a mismatch aborts
with a message rather than reading out of bounds, per the comment above it.
For the production grid the file is `inputs/mhd/veceta_grid100x02x200.txt`
(confirmed present in-tree), and the only file matching that naming
convention in `inputs/mhd/`.

Layout: `dataC[er + ephi*Nr + ez*Nphi*Nr]`, i.e. `er` varies fastest, then
`ephi`, then `ez` — the same indexing expression appears at every consumption
site below.

Read the file directly and it contains exactly five distinct values:
`-2, -1, 0, 1, 2` (verified below, FACT confirmed). Not a continuum, not a
signed distance — five discrete tags.

## The tag-to-material map

The map from integer tag to physical material and resistivity parameter is
established independently in every `alpha*` coefficient function that
switches on `dataC`, of which `alphaec2` (`src/mhd/mass_matrix_coefficients.c:2372`)
is representative. Its material branch is at lines 2432-2477:

```
src/mhd/mass_matrix_coefficients.c:2434  fabs(dataC[...] - 1.0) < 1e-12  -> user->etaplasma   (tag 1, plasma)
src/mhd/mass_matrix_coefficients.c:2436  fabs(dataC[...])       < 1e-12  -> user->etawall     (tag 0, blanket wall)
src/mhd/mass_matrix_coefficients.c:2438  fabs(dataC[...] - 2.0) < 1e-12  -> user->etasepwal   (tag 2, separatrix-wall)
src/mhd/mass_matrix_coefficients.c:2440  fabs(dataC[...] + 1.0) < 1e-12  -> user->etaVV       (tag -1, vacuum vessel)
src/mhd/mass_matrix_coefficients.c:2442  else                             -> user->etaout      (tag -2, exterior)
```

i.e.

| tag | material | `User` field |
|---|---|---|
| `1` | plasma | `etaplasma` |
| `2` | separatrix-wall | `etasepwal` |
| `0` | blanket wall | `etawall` |
| `-1` | vacuum vessel | `etaVV` |
| `-2` | exterior | `etaout` |

Its sibling `alphaec` (`src/mhd/mass_matrix_coefficients.c:2658`) reproduces
the identical branch structure at lines 2718-2762, against `user->mu0`
instead of `user->eta0` as the numerator constant but with the same tag
order and the same five destination fields (`etaplasma`, `etawall`,
`etasepwal`, `etaVV`, `etaout`). The map is not centralized: every
`alpha*`/`beta*` function that is material-aware reimplements this same
five-branch `if`/`else if` chain over `dataC`, each with its own copy of the
`1e-12` tolerance and `+1.0`/`-1.0`/`-2.0` literals.

Note the comparisons are written as `fabs(dataC[...] - value) < 1e-12`,
floating-point closeness, not integer equality (`dataC` is a `PetscReal*`
even though it only ever holds five integer-valued doubles). This is
distinct from, and much tighter than, the plasma/separatrix predicate
discussed next.

## The plasma/separatrix predicate — a material test, not a geometric one

Distinct from the five-way exact-tag dispatch above, 112 sites in
`src/mhd/ts_functions.c` use a different, looser test:

```c
fabs(user->dataC[...] - 1.5) < 0.7
```

Count verified with `grep -oE 'dataC\[.*\] - 1\.5\) < 0\.7' src/mhd/ts_functions.c
| wc -l` = 98, plus the complementary `>= 0.7` form (the same predicate,
negated) at 14 more sites, for 112 total — matching
`grep -oE 'dataC\[.*\] - 1\.5\)( < | >= )0\.7' src/mhd/ts_functions.c | wc -l` = 112.

Algebraically, `|x - 1.5| < 0.7` is `x` in the open interval `(0.8, 2.2)`.
The only tag values `dataC` ever holds are `{-2, -1, 0, 1, 2}`, and the only
members of `{-2,-1,0,1,2}` inside `(0.8, 2.2)` are `1` and `2` — plasma and
separatrix-wall. So the predicate is exactly `tag in {1, 2}`, i.e. "this cell
is plasma or separatrix-wall," restated with an interval test instead of two
equality tests.

**This is a material predicate, not a geometric level set.** The `1.5` is
the midpoint between adjacent tag integers `1` and `2`; the `0.7` is a
tag-matching half-width chosen to catch both of those integers while
excluding `0` and `-1` (distance `1.5` and `2.5` respectively) — it is not a
physical length, a grid spacing, or a distance-to-boundary tolerance. Reading
it as "cells within `0.7` of some geometric reference `1.5`" is the natural
and wrong interpretation, and is the most likely way a future maintainer
introduces a bug here — e.g. by "fixing" the tolerance to track a mesh
refinement, when the tolerance's only job is to separate two adjacent
integers and is correct at any resolution. A single representative site,
`src/mhd/ts_functions.c:421`:

```c
if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7) {
```

## Verified against the data on disk

Reading `inputs/mhd/veceta_grid100x02x200.txt` directly (100x2x200 = 40000
cells, header count 40000, confirmed matching) gives exactly five distinct
values, with these counts:

Both φ-planes combined (80000 samples over the two planes, since the array
also carries an `ephi` axis of extent `Nphi = 2`):

| tag | count (both planes) |
|---|---|
| `-2` (exterior) | 11696 |
| `-1` (vacuum vessel) | 6604 |
| `0` (blanket wall) | 7960 |
| `1` (plasma) | 10460 |
| `2` (separatrix-wall) | 3280 |

Per φ-plane (20000 cells each; the two planes are identical, confirmed —
axisymmetry means the same tag map is duplicated at both `ephi=0` and
`ephi=1`):

| tag | material | count |
|---|---|---|
| `-2` | exterior | 5848 |
| `1` | plasma | 5230 |
| `0` | blanket wall | 3980 |
| `-1` | vacuum vessel | 3302 |
| `2` | separatrix-wall | 1640 |

Sum: 5848+5230+3980+3302+1640 = 20000, the full per-plane cell count.

The plasma/separatrix predicate `|dataC - 1.5| < 0.7` selects **6870** cells
per plane — that is `5230 (tag 1) + 1640 (tag 2) = 6870`, confirming the
predicate is exactly the union of the plasma and separatrix-wall tag classes
and nothing else.

## One dataset, two formats, two languages, no cross-check

`inputs/mhd/veceta_grid100x02x200.txt` (read by the C side, per above) and
`inputs/AxisSymmetricGeometry.dat` (read by the C++ side) are **bit-identical**
copies of the same material map, verified cell-by-cell: comparing
`dataC[er + 0*Nr + ez*Nphi*Nr]` (the `ephi=0` plane of the 3-D array) against
the `.dat` file's row `er`, column `ez` gives **zero mismatches over all
20000 cells** (and likewise zero mismatches against the `ephi=1` plane,
confirming both φ-planes of the 3-D array are identical, as expected for an
axisymmetric problem).

Precisely on the indexing: `AxisSymmetricGeometry.dat` is a 100-row,
200-column whitespace-separated ASCII matrix (confirmed: `wc -l` reports 100
lines, and the first line has 200 tokens). Row index corresponds to `er`,
column index corresponds to `ez`. Relative to the usual `(row, col) = (z, r)`
convention used when plotting a poloidal cross-section, this file's
`(row, col) = (er, ez)` layout is therefore the transpose of that plotting
convention — stated here as the precise index correspondence rather than
relying on the word "transposed" alone.

The `.dat` file is read independently by the C++ kinetic side at
`src/kinetic/kinetic.cpp:318-323`:

```c++
std::ifstream ifs(configurationdomain_file);
auto indicator_h = Kokkos::create_mirror_view(indicator);
for (int i = 0; i < NR; ++i) {
  for (int j = 0; j < NZ; ++j) {
    ifs >> indicator_h(i,j);
  }
}
```

This loop has no check that `ifs` opened successfully and no check that the
stream still has data at each `>>` — a missing file, a truncated file, or a
file with the wrong dimensions is read silently as zeros with no diagnostic.
Compare `ReadInitialData` on the C side, which at least checks the element
count against the expected `Nr*Nphi*Nz` (`src/mhd/mhd.c:110-115`) even though
it too has weak error handling (`i != *num` after the read loop only
`printf`s a warning, per `src/mhd/ts_functions.c:18485-18487`).

So: one physical dataset, maintained as two on-disk files in two different
formats, read by two different languages through two different (and
differently fragile) parsers, with nothing in the build or at runtime that
verifies the two agree. `tests/regression/viz/plot_materials.py` (below) is
the only place in the repository that currently checks this.

## Why `mass_matrix_coefficients.c` dominates the `dataC` references

```
$ grep -c 'dataC' src/mhd/*.c
src/mhd/main.c:0
src/mhd/mfd_config.c:1
src/mhd/default_petsc_options.c:0
src/mhd/mhd.c:2
src/mhd/monitor_functions.c:4
src/mhd/mass_matrix_coefficients.c:139
src/mhd/geometry.c:1
src/mhd/mimetic_operators.c:21
src/mhd/ts_functions.c:113
```

(`mfd_config.c:1` is just the struct dump's `"dataC = %s\n"` line, not a use
of the array.)

`mass_matrix_coefficients.c` is the diagonal mass-matrix coefficient file —
every `alpha*`/`beta*` function there computes a coefficient for one
mimetic-operator entry, and per the file's own physical role (see
`docs/mhd/README.md`'s file table: "diagonal mass-matrix coefficients ...
material-weighted") *every* such coefficient is material-weighted, i.e. it
divides a cell volume or edge/face length by a resistivity that depends on
which of the five materials the cell belongs to. Because the file defines
many near-duplicate coefficient functions (`alphaec2`, `alphaec`,
`alphaecnores`, `alphaecnomp`, `alphaecperp`, `alphaecphi`, `alphaecperp2`,
`alphaecphi2`, `alphaecphi_isolcell`, ...), each with its own copy of the
same five-branch tag dispatch, the same handful of `dataC[...]` index
expressions is repeated function-by-function throughout the file, producing
139 occurrences from what is conceptually one lookup table.

`ts_functions.c` is next (113, including the 112 plasma/separatrix
predicate sites plus one `>= 0.7` variant covered in the same grep) because
the residual-assembly functions (`FormIFunction_*`) branch on material to
decide which physics term applies at a given edge (Ohmic term inside plasma
vs. none in the wall, etc.) — the same material information, consumed at the
point of use rather than at the point of coefficient construction.

`mimetic_operators.c` (21) uses `dataC` for isolated special-cased operator
variants; `monitor_functions.c` (4), `mhd.c` (2), and `geometry.c` (1) touch
it only incidentally (diagnostics, the initial load, and one boundary
helper).

## Recommendation (not implemented here)

The 112 plasma/separatrix predicate sites in `ts_functions.c`, and the
repeated five-branch tag dispatch scattered across
`mass_matrix_coefficients.c`, are candidates for a single named helper, e.g.:

```c
PetscBool MFD_IsPlasmaRegion(User *user, PetscInt er, PetscInt ephi, PetscInt ez);
```

replacing each `fabs(user->dataC[er + ephi*N[0] + ez*N[1]*N[0]] - 1.5) < 0.7`
(and its `>= 0.7` negation) with a call to this helper. This is a mechanical,
provably behaviour-preserving change: the helper's body would be exactly the
existing expression, so every call site produces the identical boolean for
identical inputs — the substitution is not attempted here, only recommended,
since the fact-finding task this document answers does not extend to
modifying `src/`.

## Figure

`tests/regression/viz/plot_materials.py` renders the tag map (colored by the
five materials) with the plasma/separatrix predicate region shaded and the
three hardcoded "isolated cell" indices (see below) marked, plus it performs
the bit-identical cross-check against `AxisSymmetricGeometry.dat` described
above and prints connected-component analysis of the wall region. The
generated figure lives at
[`tests/regression/viz/out/materials.png`](../../../tests/regression/viz/out/materials.png)
(relative to this document).

Regenerate with:

```sh
python3 tests/regression/viz/plot_materials.py --nr 100 --nphi 2 --nz 200 -o tests/regression/viz/out
```

(This document does not attempt to view the image; see the script's own
console output for the numeric cross-check and connected-component report.)

## Consequences for changing the grid

The material data's filename encodes the grid resolution
(`veceta_grid<NR>x<Nphi>x<NZ>.txt`, checked in `src/mhd/mhd.c:96-115` as
described above), so any change to `Nr`, `Nphi`, or `Nz` requires a
regenerated file at the new resolution — the production tree only ships
`veceta_grid100x02x200.txt`. `tests/regression/gen_grid_data.py` can generate
tag data at an arbitrary resolution, either uniformly plasma-tagged (for
manufactured-solution convergence tests, where a single material keeps mass-
matrix coefficients smooth) or by nearest-neighbour downsampling of the
production 100x200 map (to keep all five material code paths exercised, at
the cost of not being a physics-equivalent coarsening).

Regenerating the tag field at a new resolution is not sufficient to make the
solver behave identically at that resolution, though, because
`alphaecphi_isolcell` (`src/mhd/mass_matrix_coefficients.c:3361`) hardcodes
three *absolute* cell indices, not a material tag or a coordinate:

```
src/mhd/mass_matrix_coefficients.c:3426:  if((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176)){
src/mhd/mass_matrix_coefficients.c:3451:  if((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176)){
src/mhd/mass_matrix_coefficients.c:3469:  if((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176)){
```

(all three verified at those exact line numbers). These three `(er, ez)`
pairs identify specific wall cells (tag `0`) that get a separate, lower
toroidal resistivity (`user->etawallphi_isol_cell` instead of
`user->etawallphi`) — evidently cells the code owner identified by hand as
needing special treatment on the production 100x200 mesh. On any other
`(Nr, Nz)` mesh, indices `(28, 183)`, `(13, 176)`, and `(47, 176)` address
different physical cells, or fall out of range entirely — there is no
mechanism that re-locates them relative to the new grid. Changing the mesh
resolution therefore silently changes which cells (if any) get the isolated
low-resistivity treatment, with no error or warning.

(Note: `alphaecphi_isolcell` is the function that contains these three
hardcoded-index literals. A closely related function, `betaephi_isolcell`
(`src/mhd/mass_matrix_coefficients.c:1326`), is the edge-length-weighted
helper that *calls* `alphaecphi_isolcell` to assemble the isolated-cell
resistivity onto φ-edges — it does not itself contain the `(er, ez)`
literals. Some comments in the regression tooling attribute the hardcoding
to `betaephi_isolcell` by association; the literals themselves live in
`alphaecphi_isolcell` at the three lines above.)
