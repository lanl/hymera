/* T1 regression harness: the live implicit residuals.
 *
 * Evaluates each residual the solver registers with TSSetIFunction exactly once,
 * on a fixed state, and prints every entry of F as an exact hexadecimal double.
 * Two runs that agree line for line agree bit for bit.
 *
 * This is the gate for refactoring ts_functions.c. Its residuals are ~1000-line
 * functions that share a long skeleton -- the same ~70 DMStagGetLocationSlot
 * lookups, the same per-cell edge-length block, the same auxiliary-field
 * construction -- so collapsing that skeleton must leave F unchanged to the last
 * bit. A full time step cannot provide that check: the residual's output is
 * consumed by an inexact Newton-Krylov solve, which can amplify a one-ULP change
 * into a change at solver tolerance. A single evaluation has no such
 * amplification.
 *
 * Residuals covered, with where each is registered:
 *   FormIFunction_Vperp_viscosity        mhd.c         production time step
 *   FormIFunction_newequilibrium_Vperp   ts_functions  initial-condition relaxation
 *   FormIFunction_InitializeEP           ts_functions  initial-condition EP solve
 *   FormIFunction_InitializeEP_halo      ts_functions  initial-condition EP solve
 *
 * State: X and Xdot are filled with fixed asymmetric, nowhere-zero patterns, so
 * that reading a wrong slot, neighbour or component changes the result. The
 * runaway-electron current user->jre, which only the production residual reads,
 * is nonzero for the same reason. The material tags contain all five materials
 * and the three hardcoded isolated cells, matching t0_coefficients.
 *
 * ictype is 9, the production value. The residuals call FormExactSolution for
 * their reference state; for ictype 9 with Ebc == 0 that only assigns fixed
 * values, so no nested solve or input file is involved.
 *
 * Usage:   t1_residuals > out.txt
 */

#include <petscdmstag.h>
#include <petscts.h>
#include <math.h>
#include <stdio.h>

#include "mfd_config.h"
#include "ts_functions.h"

#define NR   50
#define NPHI 2
#define NZ   186

typedef PetscErrorCode (*ifunction)(TS, PetscReal, Vec, Vec, Vec, void *);

static const struct { const char *name; ifunction fn; } RESIDUALS[] = {
  {"FormIFunction_Vperp_viscosity",      FormIFunction_Vperp_viscosity},
  {"FormIFunction_newequilibrium_Vperp", FormIFunction_newequilibrium_Vperp},
  {"FormIFunction_InitializeEP",         FormIFunction_InitializeEP},
  {"FormIFunction_InitializeEP_halo",    FormIFunction_InitializeEP_halo},
};
#define NRES ((int)(sizeof(RESIDUALS) / sizeof(RESIDUALS[0])))

/* Material tags. A tokamak-like layout -- nested shells of plasma, separatrix,
 * blanket wall, vacuum vessel and exterior -- so that large connected plasma
 * regions exist and the residuals' interior plasma branches (which require all
 * four cells around a vertex to be plasma) are exercised. A band of rapidly
 * alternating tags crosses the middle so that every material also borders every
 * other, exercising the interface branches. The three hardcoded isolated cells
 * are blanket wall, as in the production geometry.
 *
 * An earlier version used only the alternating pattern; that left no vertex with
 * four plasma neighbours, so whole branches were never evaluated and a
 * deliberate perturbation in one of them went undetected. */
static double tag_at(int er, int ez) {
  if ((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176))
    return 0.0;
  static const double band[5] = {1.0, 2.0, 0.0, -1.0, -2.0};
  if (ez >= 40 && ez <= 60 && er >= 5 && er <= 45) return band[(er + 2 * ez) % 5];
  const double r = hypot((er - 24.5) / 25.0, (ez - 92.5) / 92.0);
  if (r < 0.45) return 1.0;
  if (r < 0.60) return 2.0;
  if (r < 0.80) return 0.0;
  if (r < 0.92) return -1.0;
  return -2.0;
}

static double pattern(PetscInt i, double phase) {
  const double x = (double)i;
  return 1.0 + 0.37 * sin(0.013 * x + phase) + 0.21 * cos(0.0071 * x * x) + 1e-3 * x;
}

static PetscErrorCode setup(User *u, TS *ts, PetscInt phibtype, double *jre) {
  PetscFunctionBeginUser;

  u->Nr = NR; u->Nphi = NPHI; u->Nz = NZ;
  u->phibtype = phibtype;
  u->ictype = 9;
  u->Ebc = 0;
  u->itime = 0.0;
  u->L0 = 2.0;
  u->B0 = 5.3;
  u->V_A = 1.16024e+07;
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
  u->dampV = 0.17;
  u->dt = 100.0;
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

  /* jre: an NR x NZ x 3 host view, laid out as the C++ side's Kokkos DualView
   * (LayoutRight: last index fastest). */
  for (int i = 0; i < NR * NZ * 3; ++i) jre[i] = 0.3 * pattern(i, 0.7);
  u->jre.data = jre;
  u->jre.dim0 = NR;  u->jre.dim1 = NZ;  u->jre.dim2 = 3;
  u->jre.stride0 = NZ * 3;  u->jre.stride1 = 3;  u->jre.stride2 = 1;

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
  PetscCall(TSSetTimeStep(*ts, u->dt));
  u->ts = *ts;

  /* X0 is the previous-step state some residuals difference against. */
  PetscCall(DMCreateGlobalVector(u->da, &u->X0));
  {
    PetscInt     n;
    PetscScalar *x;
    PetscCall(VecGetSize(u->X0, &n));
    PetscCall(VecGetArray(u->X0, &x));
    for (PetscInt i = 0; i < n; ++i) x[i] = pattern(i, 2.1);
    PetscCall(VecRestoreArray(u->X0, &x));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode teardown(User *u, TS *ts) {
  PetscFunctionBeginUser;
  DM  dmc;
  Vec cl;
  PetscCall(DMGetCoordinateDM(u->coorda, &dmc));
  PetscCall(DMGetCoordinatesLocal(u->coorda, &cl));
  PetscCall(DMStagVecRestoreArrayRead(dmc, cl, &u->arrCoord));
  PetscCall(VecDestroy(&u->X0));
  PetscCall(TSDestroy(ts));
  PetscCall(DMDestroy(&u->coorda));
  PetscCall(DMDestroy(&u->da));
  PetscCall(PetscFree(u->dataC));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode fill(Vec v, double phase) {
  PetscInt     n;
  PetscScalar *x;
  PetscFunctionBeginUser;
  PetscCall(VecGetSize(v, &n));
  PetscCall(VecGetArray(v, &x));
  for (PetscInt i = 0; i < n; ++i) x[i] = pattern(i, phase);
  PetscCall(VecRestoreArray(v, &x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  static double jre[NR * NZ * 3];
  long count = 0;

  for (PetscInt b = 0; b <= 1; ++b) {
    User u = {0};
    TS   ts;
    PetscCall(setup(&u, &ts, b, jre));
    printf("# phibtype=%d\n", (int)b);

    Vec X, Xdot;
    PetscCall(DMCreateGlobalVector(u.da, &X));
    PetscCall(DMCreateGlobalVector(u.da, &Xdot));
    PetscCall(fill(X, 0.0));
    PetscCall(fill(Xdot, 1.3));

    for (int k = 0; k < NRES; ++k) {
      Vec F;
      PetscCall(DMCreateGlobalVector(u.da, &F));
      PetscCall(VecZeroEntries(F));
      PetscCall(RESIDUALS[k].fn(ts, 0.0, X, Xdot, F, &u));
      PetscInt           n;
      const PetscScalar *f;
      PetscCall(VecGetSize(F, &n));
      PetscCall(VecGetArrayRead(F, &f));
      for (PetscInt i = 0; i < n; ++i) {
        printf("%s %d %a\n", RESIDUALS[k].name, (int)i, f[i]);
        ++count;
      }
      PetscCall(VecRestoreArrayRead(F, &f));
      PetscCall(VecDestroy(&F));
    }
    PetscCall(VecDestroy(&X));
    PetscCall(VecDestroy(&Xdot));
    PetscCall(teardown(&u, &ts));
  }

  fprintf(stderr, "t1_residuals: %ld values\n", count);
  PetscCall(PetscFinalize());
  return 0;
}
