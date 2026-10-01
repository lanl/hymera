/* T0 regression harness: the discrete mass-matrix coefficients.
 *
 * Evaluates every live alpha-, beta-, rese and condu function in
 * src/mhd/mass_matrix_coefficients.c over every grid index, stencil location and
 * material tag, and writes each result as an exact hexadecimal double. Two runs
 * that agree line for line agree bit for bit.
 *
 * This is the gate for deduplicating that file. Its 13 beta* functions are
 * near-identical clones -- betae vs betae2 differs only in which alpha* it sums --
 * so a refactor that collapses them can be proved bit-exact here in seconds,
 * with no input files, no solver and no time step.
 *
 * Coverage is chosen to reach every branch the coefficients take:
 *   - a small grid, so that boundary handling is a large fraction of all cells
 *     (er = 0 and Nr-1, the phi seam, ez = 0 and Nz-1, all corners);
 *   - phi both periodic and non-periodic (user->phibtype), which selects
 *     different branch ladders;
 *   - ictype 9 (material-dependent resistivity) and a non-9 value (uniform eta);
 *   - a tag field containing all five materials, laid out so that every cell has
 *     neighbours of differing tag;
 *   - the three hardcoded "isolated" cells (28,183), (13,176), (47,176) inside
 *     alphaecphi_isolcell, so the grid is sized to contain them.
 *
 * Usage:
 *   t0_coefficients > out.txt          then diff against a recorded baseline
 */

#include <petscdmstag.h>
#include <stdio.h>

#include "mfd_config.h"
#include "mass_matrix_coefficients.h"

/* Must cover the isolcell literals (er up to 47, ez up to 183) while staying
 * small enough to run in seconds. phi needs at least 2 cells for DMStag. */
#define NR   50
#define NPHI 2
#define NZ   186

/* Locations accepted by each function family. Each beta* SETERRQs on any other
 * location, so passing only valid ones keeps the harness from tripping that. */
static const DMStagStencilLocation EDGES[] = {
  DMSTAG_BACK_DOWN, DMSTAG_BACK_LEFT, DMSTAG_BACK_RIGHT, DMSTAG_BACK_UP,
  DMSTAG_DOWN_LEFT, DMSTAG_DOWN_RIGHT, DMSTAG_UP_LEFT, DMSTAG_UP_RIGHT,
  DMSTAG_FRONT_DOWN, DMSTAG_FRONT_LEFT, DMSTAG_FRONT_RIGHT, DMSTAG_FRONT_UP,
};
static const DMStagStencilLocation FACES[] = {
  DMSTAG_LEFT, DMSTAG_RIGHT, DMSTAG_DOWN, DMSTAG_UP, DMSTAG_BACK, DMSTAG_FRONT,
};
static const DMStagStencilLocation VERTICES[] = {
  DMSTAG_BACK_DOWN_LEFT, DMSTAG_BACK_DOWN_RIGHT, DMSTAG_BACK_UP_LEFT,
  DMSTAG_BACK_UP_RIGHT, DMSTAG_FRONT_DOWN_LEFT, DMSTAG_FRONT_DOWN_RIGHT,
  DMSTAG_FRONT_UP_LEFT, DMSTAG_FRONT_UP_RIGHT,
};
#define LEN(a) ((int)(sizeof(a) / sizeof((a)[0])))

typedef PetscScalar (*loc_fn)(PetscInt, PetscInt, PetscInt, DMStagStencilLocation, void *);
typedef PetscScalar (*cell_fn)(PetscInt, PetscInt, PetscInt, void *);

/* The live functions -- reachable from the residual, the operators or the
 * diagnostics. Dead ones are excluded: they are candidates for quarantine, and
 * the harness must not keep them alive by referencing them. */
static const struct { const char *name; loc_fn fn; const DMStagStencilLocation *locs; int nloc; }
LOC_FNS[] = {
  {"betaf",             betaf,             FACES,    LEN(FACES)},
  {"betae",             betae,             EDGES,    LEN(EDGES)},
  {"betaenores",        betaenores,        EDGES,    LEN(EDGES)},
  {"betaenomp",         betaenomp,         EDGES,    LEN(EDGES)},
  {"betae2",            betae2,            EDGES,    LEN(EDGES)},
  {"betaephi_isolcell", betaephi_isolcell, EDGES,    LEN(EDGES)},
  {"betaeperp2",        betaeperp2,        EDGES,    LEN(EDGES)},
  {"betavnomp",         betavnomp,         VERTICES, LEN(VERTICES)},
  {"rese",              rese,              EDGES,    LEN(EDGES)},
  {"condu",             condu,             EDGES,    LEN(EDGES)},
};

/* Cell functions. Only those reachable from a live beta* are listed. phi runs
 * over -1..NPHI because the beta* ladders call them with ephi-1 and ephi+1 at
 * the periodic seam, and the alpha* functions special-case exactly those values. */
static const struct { const char *name; cell_fn fn; } CELL_FNS[] = {
  {"alphaec",             alphaec},
  {"alphaec2",            alphaec2},
  {"alphaecnores",        alphaecnores},
  {"alphaecnomp",         alphaecnomp},
  {"alphaecperp2",        alphaecperp2},
  {"alphaecphi",          alphaecphi},
  {"alphaecphi_isolcell", alphaecphi_isolcell},
  {"alphafc",             alphafc},
  {"alphavcnomp",         alphavcnomp},
  {"resec",               resec},
  {"conduc",              conduc},
};

/* A tag field with all five materials and maximal tag variation between
 * neighbours. The three isolcell cells must be blanket wall (tag 0) for their
 * special branch to be reached, matching the production geometry. */
static double tag_at(int er, int ez) {
  if ((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176))
    return 0.0;
  static const double tags[5] = {1.0, 2.0, 0.0, -1.0, -2.0};
  return tags[(er + 2 * ez) % 5];
}

static PetscErrorCode setup(User *u, PetscInt phibtype, PetscInt ictype) {
  PetscFunctionBeginUser;

  u->Nr = NR; u->Nphi = NPHI; u->Nz = NZ;
  u->phibtype = phibtype;
  u->ictype = ictype;

  /* Physical constants in the ranges the production run prints from AppCtxView,
   * deliberately all distinct so that a substitution of one for another in a
   * refactor changes the output. */
  u->L0 = 2.0;
  u->mu0 = 1.25663706212e-06;
  u->eta0 = 29.1599;
  u->eta = 1.0;
  u->etaplasma = 2.46016e-05;
  u->etasepwal = 3.17e-05;
  u->etawall = 0.044;
  u->etawallperp = 0.051;
  u->etawallphi = 0.037;
  u->etawallphi_isol_cell = 0.029;
  u->etaVV = 1.30288e-06;
  u->etaout = 0.00130288;
  u->rmin = 3.05; u->rmax = 9.95;
  u->zmin = -5.95; u->zmax = 5.95;
  u->phimin = 0.0; u->phimax = 2.0 * PETSC_PI;
  u->dphi = (u->phimax - u->phimin) / NPHI;
  u->debug = 0;

  PetscCall(PetscMalloc1(NR * NPHI * NZ, &u->dataC));
  for (int ez = 0; ez < NZ; ++ez)
    for (int ephi = 0; ephi < NPHI; ++ephi)
      for (int er = 0; er < NR; ++er)
        u->dataC[er + ephi * NR + ez * NPHI * NR] = tag_at(er, ez);
  u->numC = NR * NPHI * NZ;

  /* Mirror mhd_initialize: the coordinate DM, uniform coordinates, and the
   * borrowed read-only coordinate array the coefficients index into. */
  const DMBoundaryType bphi = phibtype ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
  PetscCall(DMStagCreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, bphi, DM_BOUNDARY_NONE,
                           NR, NPHI, NZ, 1, 1, 1, 4, 1, 1, 1, DMSTAG_STENCIL_BOX, 1,
                           NULL, NULL, NULL, &u->coorda));
  PetscCall(DMSetUp(u->coorda));
  PetscCall(DMStagSetUniformCoordinatesExplicit(u->coorda, u->rmin / u->L0, u->rmax / u->L0,
                                                u->phimin, u->phimax,
                                                u->zmin / u->L0, u->zmax / u->L0));
  DM  dmc;
  Vec cl;
  PetscCall(DMGetCoordinateDM(u->coorda, &dmc));
  PetscCall(DMGetCoordinatesLocal(u->coorda, &cl));
  PetscCall(DMStagVecGetArrayRead(dmc, cl, &u->arrCoord));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode teardown(User *u) {
  PetscFunctionBeginUser;
  DM  dmc;
  Vec cl;
  PetscCall(DMGetCoordinateDM(u->coorda, &dmc));
  PetscCall(DMGetCoordinatesLocal(u->coorda, &cl));
  PetscCall(DMStagVecRestoreArrayRead(dmc, cl, &u->arrCoord));
  PetscCall(DMDestroy(&u->coorda));
  PetscCall(PetscFree(u->dataC));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  const PetscInt phibtypes[] = {0, 1};
  const PetscInt ictypes[] = {9, 1};
  long count = 0;

  for (int b = 0; b < 2; ++b) {
    for (int c = 0; c < 2; ++c) {
      User u = {0};
      PetscCall(setup(&u, phibtypes[b], ictypes[c]));
      printf("# phibtype=%d ictype=%d\n", (int)phibtypes[b], (int)ictypes[c]);

      for (int f = 0; f < LEN(LOC_FNS); ++f) {
        for (int ez = 0; ez < NZ; ++ez)
          for (int ephi = 0; ephi < NPHI; ++ephi)
            for (int er = 0; er < NR; ++er)
              for (int l = 0; l < LOC_FNS[f].nloc; ++l) {
                const double v = LOC_FNS[f].fn(er, ephi, ez, LOC_FNS[f].locs[l], &u);
                printf("%s %d %d %d %d %a\n", LOC_FNS[f].name, er, ephi, ez,
                       (int)LOC_FNS[f].locs[l], v);
                ++count;
              }
      }
      for (int f = 0; f < LEN(CELL_FNS); ++f) {
        for (int ez = 0; ez < NZ; ++ez)
          for (int ephi = -1; ephi <= NPHI; ++ephi)
            for (int er = 0; er < NR; ++er) {
              const double v = CELL_FNS[f].fn(er, ephi, ez, &u);
              printf("%s %d %d %d %a\n", CELL_FNS[f].name, er, ephi, ez, v);
              ++count;
            }
      }
      PetscCall(teardown(&u));
    }
  }

  fprintf(stderr, "t0_coefficients: %ld evaluations\n", count);
  PetscCall(PetscFinalize());
  return 0;
}
