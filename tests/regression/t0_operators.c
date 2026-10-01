/* T0 regression harness: the discrete mimetic operators.
 *
 * Builds the solver's DMStag (the production DOF layout: 4 per vertex, 1 per edge,
 * face and element) on a small grid, fills an input vector with a fixed
 * deterministic pattern, applies each live Vec-to-Vec operator in
 * mimetic_operators.c and geometry.c, and prints every entry of every result as
 * an exact hexadecimal double. Two runs that agree line for line agree bit for
 * bit.
 *
 * Sibling of t0_coefficients.c, one level up: those test the mass-matrix
 * coefficients directly, this tests the operators that consume them. It is the
 * gate for deduplicating the curl family (FormDerivedCurl, FormDerivedCurlnores,
 * FormDerivedCurlnomp are token-identical except for the edge coefficient) and
 * any later operator refactor.
 *
 * The input pattern deliberately has no symmetry and no zeros, so that an
 * operator reading the wrong slot, the wrong neighbour or the wrong component
 * produces a different result rather than coincidentally the same one.
 *
 * Usage:   t0_operators > out.txt
 */

#include <petscdmstag.h>
#include <petscts.h>
#include <stdio.h>

#include "mfd_config.h"
#include "mimetic_operators.h"
#include "geometry.h"

/* Must contain the three hardcoded isolated cells used by the coefficients. */
#define NR   50
#define NPHI 2
#define NZ   186

typedef PetscErrorCode (*vec_op)(TS, Vec, Vec, void *);

static const struct { const char *name; vec_op fn; } OPS[] = {
  {"FormPrimaryCurl",              FormPrimaryCurl},
  {"FormDerivedCurl",              FormDerivedCurl},
  {"FormDerivedCurlnores",         FormDerivedCurlnores},
  {"FormDerivedCurlnomp",          FormDerivedCurlnomp},
  {"ApplyVectorLaplacian",         ApplyVectorLaplacian},
  {"FormDiscreteGradientEP_noMat", FormDiscreteGradientEP_noMat},
  {"FormElectricField",            FormElectricField},
  {"CellToVertexProjectionScalar", CellToVertexProjectionScalar},
  {"VertexToEdgeReconstruction",   VertexToEdgeReconstruction},
  {"VertexToFaceReconstruction",   VertexToFaceReconstruction},
  {"FaceToVertexProjection",       FaceToVertexProjection},
  {"EdgeToVertexProjection",       EdgeToVertexProjection},
  {"CellToFaceProjection",         CellToFaceProjection},
};
#define NOPS ((int)(sizeof(OPS) / sizeof(OPS[0])))

static double tag_at(int er, int ez) {
  if ((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176))
    return 0.0;
  static const double tags[5] = {1.0, 2.0, 0.0, -1.0, -2.0};
  return tags[(er + 2 * ez) % 5];
}

/* A smooth-ish but asymmetric, nowhere-zero pattern. Irrational-looking
 * coefficients avoid accidental cancellation between neighbouring entries. */
static double pattern(PetscInt i) {
  const double x = (double)i;
  return 1.0 + 0.37 * sin(0.013 * x) + 0.21 * cos(0.0071 * x * x) + 1e-3 * x;
}

static PetscErrorCode setup(User *u, TS *ts, PetscInt phibtype) {
  PetscFunctionBeginUser;

  u->Nr = NR; u->Nphi = NPHI; u->Nz = NZ;
  u->phibtype = phibtype;
  u->ictype = 9;
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
  u->Re = 200.0;
  u->rmin = 3.05; u->rmax = 9.95;
  u->zmin = -5.95; u->zmax = 5.95;
  u->phimin = 0.0; u->phimax = 2.0 * PETSC_PI;
  u->dphi = (u->phimax - u->phimin) / NPHI;
  u->dr = (u->rmax - u->rmin) / NR;
  u->dz = (u->zmax - u->zmin) / NZ;
  u->debug = 0;

  PetscCall(PetscMalloc1(NR * NPHI * NZ, &u->dataC));
  for (int ez = 0; ez < NZ; ++ez)
    for (int ephi = 0; ephi < NPHI; ++ephi)
      for (int er = 0; er < NR; ++er)
        u->dataC[er + ephi * NR + ez * NPHI * NR] = tag_at(er, ez);
  u->numC = NR * NPHI * NZ;

  /* Mirror mhd_initialize: the solution DM, the coordinate DM, uniform
   * coordinates on both, and the borrowed coordinate array. */
  const DMBoundaryType bphi = phibtype ? DM_BOUNDARY_PERIODIC : DM_BOUNDARY_NONE;
  PetscCall(DMStagCreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, bphi, DM_BOUNDARY_NONE,
                           NR, NPHI, NZ, 1, 1, 1, 4, 1, 1, 1, DMSTAG_STENCIL_BOX, 1,
                           NULL, NULL, NULL, &u->da));
  PetscCall(DMSetUp(u->da));
  PetscCall(DMStagCreate3d(PETSC_COMM_SELF, DM_BOUNDARY_NONE, bphi, DM_BOUNDARY_NONE,
                           NR, NPHI, NZ, 1, 1, 1, 4, 1, 1, 1, DMSTAG_STENCIL_BOX, 1,
                           NULL, NULL, NULL, &u->coorda));
  PetscCall(DMSetUp(u->coorda));
  PetscCall(DMStagSetUniformCoordinatesExplicit(u->da, u->rmin / u->L0, u->rmax / u->L0,
                                                u->phimin, u->phimax,
                                                u->zmin / u->L0, u->zmax / u->L0));
  PetscCall(DMStagSetUniformCoordinatesExplicit(u->coorda, u->rmin / u->L0, u->rmax / u->L0,
                                                u->phimin, u->phimax,
                                                u->zmin / u->L0, u->zmax / u->L0));
  DM  dmc;
  Vec cl;
  PetscCall(DMGetCoordinateDM(u->coorda, &dmc));
  PetscCall(DMGetCoordinatesLocal(u->coorda, &cl));
  PetscCall(DMStagVecGetArrayRead(dmc, cl, &u->arrCoord));

  PetscCall(TSCreate(PETSC_COMM_SELF, ts));
  PetscCall(TSSetDM(*ts, u->da));
  u->ts = *ts;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode teardown(User *u, TS *ts) {
  PetscFunctionBeginUser;
  DM  dmc;
  Vec cl;
  PetscCall(DMGetCoordinateDM(u->coorda, &dmc));
  PetscCall(DMGetCoordinatesLocal(u->coorda, &cl));
  PetscCall(DMStagVecRestoreArrayRead(dmc, cl, &u->arrCoord));
  PetscCall(TSDestroy(ts));
  PetscCall(DMDestroy(&u->coorda));
  PetscCall(DMDestroy(&u->da));
  PetscCall(PetscFree(u->dataC));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  long count = 0;
  for (PetscInt b = 0; b <= 1; ++b) {
    User u = {0};
    TS   ts;
    PetscCall(setup(&u, &ts, b));
    printf("# phibtype=%d\n", (int)b);

    Vec X;
    PetscCall(DMCreateGlobalVector(u.da, &X));
    PetscInt n;
    PetscCall(VecGetSize(X, &n));
    {
      PetscScalar *x;
      PetscCall(VecGetArray(X, &x));
      for (PetscInt i = 0; i < n; ++i) x[i] = pattern(i);
      PetscCall(VecRestoreArray(X, &x));
    }

    for (int k = 0; k < NOPS; ++k) {
      Vec F;
      PetscCall(DMCreateGlobalVector(u.da, &F));
      PetscCall(VecZeroEntries(F));
      PetscCall(OPS[k].fn(ts, X, F, &u));
      const PetscScalar *f;
      PetscCall(VecGetArrayRead(F, &f));
      for (PetscInt i = 0; i < n; ++i) {
        printf("%s %d %a\n", OPS[k].name, (int)i, f[i]);
        ++count;
      }
      PetscCall(VecRestoreArrayRead(F, &f));
      PetscCall(VecDestroy(&F));
    }
    PetscCall(VecDestroy(&X));
    PetscCall(teardown(&u, &ts));
  }

  fprintf(stderr, "t0_operators: %ld values\n", count);
  PetscCall(PetscFinalize());
  return 0;
}
