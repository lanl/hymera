//========================================================================================
// (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
// for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================

#include "mfd_config.h"
#include "ts_functions.h"
#include "monitor_functions.h"
#include "geometry.h"
#include "mass_matrix_coefficients.h"
#include "mimetic_operators.h"

/* DMStag location-slot indices shared by the FormIFunction_* residuals.
   Member names match the residuals' local variables one-for-one. */
typedef struct {
  /* solution DM (da): vertex (4 components), edge, face, element */
  PetscInt ivVrmphimzm[4], ivVrpphimzm[4], ivVrmphipzm[4], ivVrpphipzm[4];
  PetscInt ivVrmphimzp[4], ivVrpphimzp[4], ivVrmphipzp[4], ivVrpphipzp[4];
  PetscInt ivErmzm, ivEphimzm, ivErpzm, ivEphipzm;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  PetscInt ivEphimzp, ivErmzp, ivErpzp, ivEphipzp;
  PetscInt ivBrm, ivBphim, ivBzm, ivBrp, ivBphip, ivBzp;
  PetscInt ivn;
  /* coordinate DMs (dmCoord: faces/edges; dmCoorda: element/vertices), 3 components */
  PetscInt icp[3];
  PetscInt icBrm[3], icBphim[3], icBzm[3], icBrp[3], icBphip[3], icBzp[3];
  PetscInt icErmzm[3], icEphimzm[3], icErpzm[3], icEphipzm[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icEphimzp[3], icErmzp[3], icErpzp[3], icEphipzp[3];
  PetscInt icrmphimzm[3], icrpphimzm[3], icrmphipzm[3], icrpphipzm[3];
  PetscInt icrmphimzp[3], icrpphimzp[3], icrmphipzp[3], icrpphipzp[3];
} MFD_Slots;

static PetscErrorCode MFD_GetSlotsSolution(DM da, MFD_Slots * s) {
  PetscInt d;

  PetscFunctionBeginUser;
  for (d = 0; d < 4; ++d) {
    /* Vertex locations */
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & s->ivVrmphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & s->ivVrpphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & s->ivVrmphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & s->ivVrpphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & s->ivVrmphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & s->ivVrpphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & s->ivVrmphipzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & s->ivVrpphipzp[d]));
  }
  /* Edge locations */
  PetscCall(DMStagGetLocationSlot(da, BACK_LEFT, 0, & s->ivErmzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_DOWN, 0, & s->ivEphimzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_RIGHT, 0, & s->ivErpzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_UP, 0, & s->ivEphipzm));
  PetscCall(DMStagGetLocationSlot(da, DOWN_LEFT, 0, & s->ivErmphim));
  PetscCall(DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & s->ivErpphim));
  PetscCall(DMStagGetLocationSlot(da, UP_LEFT, 0, & s->ivErmphip));
  PetscCall(DMStagGetLocationSlot(da, UP_RIGHT, 0, & s->ivErpphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN, 0, & s->ivEphimzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_LEFT, 0, & s->ivErmzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & s->ivErpzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_UP, 0, & s->ivEphipzp));
  /* Face locations */
  PetscCall(DMStagGetLocationSlot(da, LEFT, 0, & s->ivBrm));
  PetscCall(DMStagGetLocationSlot(da, DOWN, 0, & s->ivBphim));
  PetscCall(DMStagGetLocationSlot(da, BACK, 0, & s->ivBzm));
  PetscCall(DMStagGetLocationSlot(da, RIGHT, 0, & s->ivBrp));
  PetscCall(DMStagGetLocationSlot(da, UP, 0, & s->ivBphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT, 0, & s->ivBzp));
  /* Cell locations */
  PetscCall(DMStagGetLocationSlot(da, ELEMENT, 0, & s->ivn));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode MFD_GetSlotsCoords(DM dmCoord, DM dmCoorda, MFD_Slots * s) {
  PetscInt d;

  PetscFunctionBeginUser;
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & s->icp[d]));
    /* Face coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, LEFT, d, & s->icBrm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN, d, & s->icBphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK, d, & s->icBzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, RIGHT, d, & s->icBrp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP, d, & s->icBphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT, d, & s->icBzp[d]));
    /* Edge coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & s->icErmzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & s->icEphimzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & s->icErpzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_UP, d, & s->icEphipzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & s->icErmphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & s->icErpphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & s->icErmphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & s->icErpphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & s->icEphimzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & s->icErmzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & s->icErpzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & s->icEphipzp[d]));
    /* Vertex coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & s->icrmphimzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & s->icrpphimzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & s->icrmphipzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & s->icrpphipzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & s->icrmphimzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & s->icrpphimzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & s->icrmphipzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & s->icrpphipzp[d]));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Copy every slot index from an MFD_Slots value into the identically named
   locals of the calling residual. Uses the caller's loop index `d`. */
#define MFD_UNPACK_SLOTS(S) do { \
    for (d = 0; d < 4; ++d) { \
      ivVrmphimzm[d] = (S).ivVrmphimzm[d]; ivVrpphimzm[d] = (S).ivVrpphimzm[d]; \
      ivVrmphipzm[d] = (S).ivVrmphipzm[d]; ivVrpphipzm[d] = (S).ivVrpphipzm[d]; \
      ivVrmphimzp[d] = (S).ivVrmphimzp[d]; ivVrpphimzp[d] = (S).ivVrpphimzp[d]; \
      ivVrmphipzp[d] = (S).ivVrmphipzp[d]; ivVrpphipzp[d] = (S).ivVrpphipzp[d]; \
    } \
    ivErmzm = (S).ivErmzm; ivEphimzm = (S).ivEphimzm; ivErpzm = (S).ivErpzm; ivEphipzm = (S).ivEphipzm; \
    ivErmphim = (S).ivErmphim; ivErpphim = (S).ivErpphim; ivErmphip = (S).ivErmphip; ivErpphip = (S).ivErpphip; \
    ivEphimzp = (S).ivEphimzp; ivErmzp = (S).ivErmzp; ivErpzp = (S).ivErpzp; ivEphipzp = (S).ivEphipzp; \
    ivBrm = (S).ivBrm; ivBphim = (S).ivBphim; ivBzm = (S).ivBzm; \
    ivBrp = (S).ivBrp; ivBphip = (S).ivBphip; ivBzp = (S).ivBzp; \
    ivn = (S).ivn; \
    for (d = 0; d < 3; ++d) { \
      icp[d] = (S).icp[d]; \
      icBrm[d] = (S).icBrm[d]; icBphim[d] = (S).icBphim[d]; icBzm[d] = (S).icBzm[d]; \
      icBrp[d] = (S).icBrp[d]; icBphip[d] = (S).icBphip[d]; icBzp[d] = (S).icBzp[d]; \
      icErmzm[d] = (S).icErmzm[d]; icEphimzm[d] = (S).icEphimzm[d]; icErpzm[d] = (S).icErpzm[d]; icEphipzm[d] = (S).icEphipzm[d]; \
      icErmphim[d] = (S).icErmphim[d]; icErpphim[d] = (S).icErpphim[d]; icErmphip[d] = (S).icErmphip[d]; icErpphip[d] = (S).icErpphip[d]; \
      icEphimzp[d] = (S).icEphimzp[d]; icErmzp[d] = (S).icErmzp[d]; icErpzp[d] = (S).icErpzp[d]; icEphipzp[d] = (S).icEphipzp[d]; \
      icrmphimzm[d] = (S).icrmphimzm[d]; icrpphimzm[d] = (S).icrpphimzm[d]; icrmphipzm[d] = (S).icrmphipzm[d]; icrpphipzm[d] = (S).icrpphipzm[d]; \
      icrmphimzp[d] = (S).icrmphimzp[d]; icrpphimzp[d] = (S).icrpphimzp[d]; icrmphipzp[d] = (S).icrmphipzp[d]; icrpphipzp[d] = (S).icrpphipzp[d]; \
    } \
  } while (0)

/* Per-cell volume, INT_c(r dphi dr dz). Expressions are verbatim from the
   residuals' cell loops; do not rearrange them (results must stay bit-identical). */
static inline PetscScalar MFD_CellVolume(PetscScalar **** arrCoord, PetscInt er, PetscInt ephi, PetscInt ez,
  const PetscInt N[3], PetscReal dphi,
  const PetscInt icBrm[3], const PetscInt icBphim[3], const PetscInt icBzm[3],
  const PetscInt icBrp[3], const PetscInt icBphip[3], const PetscInt icBzp[3]) {
  PetscScalar cellvolume;

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }
  return cellvolume;
}

/* The 12 edge lengths of one cell, from its vertex coordinates. Expressions are
   verbatim from the residuals' cell loops; do not rearrange them. */
static inline void MFD_CellEdgeLengths(PetscScalar **** arrCoorda, PetscInt er, PetscInt ephi, PetscInt ez,
  const PetscInt N[3], PetscReal dphi,
  const PetscInt icrmphimzm[3], const PetscInt icrpphimzm[3], const PetscInt icrmphipzm[3], const PetscInt icrpphipzm[3],
  const PetscInt icrmphimzp[3], const PetscInt icrpphimzp[3], const PetscInt icrmphipzp[3], const PetscInt icrpphipzp[3],
  PetscScalar * rmzmedgelength, PetscScalar * rmphimedgelength, PetscScalar * rmzpedgelength, PetscScalar * rmphipedgelength,
  PetscScalar * phimzmedgelength, PetscScalar * phimzpedgelength, PetscScalar * rpphimedgelength, PetscScalar * rpzmedgelength,
  PetscScalar * phipzmedgelength, PetscScalar * rpzpedgelength, PetscScalar * rpphipedgelength, PetscScalar * phipzpedgelength) {
  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    * rmzmedgelength = dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
  } else {
    * rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
  }

  * rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    * rmzpedgelength = dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
  } else {
    * rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
  }

  * rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

  * phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

  * phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

  * rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    * rpzmedgelength = dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
  } else {
    * rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
  }

  * phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    * rpzpedgelength = dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
  } else {
    * rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
  }

  * rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

  * phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */
}
#line 20

PetscErrorCode FormIJacobian_BImplicit(TS ts, PetscReal t, Vec X, Vec Xdot, PetscReal a, Mat J, Mat Jpre, void * ptr) {
  PetscFunctionBeginUser;
  PetscFunctionReturn(PETSC_SUCCESS);
}


#line 973

#line 1934

/* Which terms the shared V_perp residual body includes. FormIFunction_Vperp_viscosity
 * (production time step) and FormIFunction_newequilibrium_Vperp (initial-condition
 * relaxation) were 830- and 671-line copies of vperp_residual differing only in:
 *
 *   f1_inertia    r/z momentum rows on plasma vertices:
 *                   TRUE : - (curl B x B) + n_i dV/dt - Re^-1 Lap V            (relaxation)
 *                   FALSE: - (curl B x B) + 0.0 * n_i (V.grad V ...) - Re^-1 Lap V
 *                          (production; the `0.0 *` advective term is kept verbatim)
 *   f3_resistive  inner-edge Ohm's-law rows subtract derived_mimetic_curl2(B)
 *                   (the betaf/betae2 term). FALSE = ideal Ohm's law, which the
 *                   relaxation uses by owner decision (consistency with previous runs).
 *   f3_jre        inner-edge Ohm's-law rows add the edge-averaged runaway current
 *                   user->jre (open question Q2). Kept separate from f3_resistive.
 *   label         PETSc log-event name.
 *
 *                                       f1_inertia  f3_resistive  f3_jre
 *   FormIFunction_Vperp_viscosity        FALSE       TRUE          TRUE
 *   FormIFunction_newequilibrium_Vperp   TRUE        FALSE         FALSE
 *
 * Every switch selects between complete, verbatim statements: the build contracts
 * to FMA (-ffp-contract=fast), so changing an expression's shape could change bits.
 * Only the two combinations above are covered by tests/regression/t1_residuals.c;
 * the f3 branches for (resistive, no jre) and (jre, not resistive) are untested.
 *
 * vperp_residual is always_inline so each wrapper gets a body specialised on its
 * constant `terms`, with the unused branches folded away before optimisation. As
 * a plain out-of-line static (runtime `terms`) GCC merged code across the
 * branches and the residual changed in the last ULP at 32 entries of the t1
 * fixture; inlined, both wrappers are bit-identical to the pre-merge functions
 * (same fmadd/fmsub counts, 61 and 52). Do not drop the attribute. */
typedef struct {
  PetscBool    f1_inertia;
  PetscBool    f3_resistive;
  PetscBool    f3_jre;
  PetscBool    f1_advection;  /* production only: include n_i (V.grad)V for real */
  const char * label;
} MFD_ResidualTerms;

#line 2909

static inline __attribute__((always_inline)) PetscErrorCode vperp_residual(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void * ptr, const MFD_ResidualTerms * terms) {
  PetscFunctionBeginUser;

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister(terms->label,classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  PetscScalar dt, cellvolume;
  Vec fLocal, xLocal, bcLocal, xdotLocal, pLocal;
  Vec VxBe, VxBeLocal, VxB, Vf, VfLocal, nif, nifLocal, niv = NULL, nivLocal = NULL, Bv, BvLocal, curlBv, curlBvLocal, GradEP, GradEPLocal, F1 = NULL, F2 = NULL, F3 = NULL, GradV1 = NULL, GradV1Local = NULL, GradV2 = NULL, GradV2Local = NULL, GradV3 = NULL, GradV3Local = NULL, Fcopy, FcopyLocal, LapV, LapVLocal ;
  Vec x, potential;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  view3d_t jre = user->jre;

  PetscInt icp[3] PETSC_UNUSED;
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3] PETSC_UNUSED, icErmzp[3] PETSC_UNUSED, icErpzm[3] PETSC_UNUSED, icErpzp[3] PETSC_UNUSED;

  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3] PETSC_UNUSED, icErpphim[3] PETSC_UNUSED, icErmphip[3] PETSC_UNUSED, icErpphip[3] PETSC_UNUSED;
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  PetscInt ivn;

  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, ** ** arrP, ** ** arrx, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength, ** ** arrXdot, ** ** arrBv, ** ** arrcurlBv, ** ** arrnif, ** ** arrniv = NULL, ** ** arrVf, ** ** arrVxBe, ** ** arrGradEP, ** ** arrFcopy, ** ** arrGradV3 = NULL, ** ** arrGradV2 = NULL, ** ** arrGradV1 = NULL, ** ** arrLapV;

  PetscInt steps=0;

  PetscCall(TSGetStepNumber(ts,&steps));
  PetscCall(VecZeroEntries(F));
  PetscCall(TSGetDM(ts, & da));

  PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
  PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));

  MFD_Slots S;
  PetscCall(MFD_GetSlotsSolution(da, & S));
#line 3002
  PetscCall(DMGetCoordinateDM(da, & dmCoord));
  PetscCall(DMGetCoordinatesLocal(da, & coordLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord));
  PetscCall(DMGetCoordinateDM(coordDA, & dmCoorda));
  PetscCall(DMGetCoordinatesLocal(coordDA, & coordaLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda));
  PetscCall(MFD_GetSlotsCoords(dmCoord, dmCoorda, & S));
  MFD_UNPACK_SLOTS(S);
#line 3041
  /* Compute the source term potential for time-dependent manufactured solution */
  PetscCall(DMCreateGlobalVector(da, & potential));
  FormSourceTermPotential(ts, t, potential, user);
  PetscCall(DMGetLocalVector(da, & pLocal));
  PetscCall(DMGlobalToLocalBegin(da, potential, INSERT_VALUES, pLocal));
  PetscCall(DMGlobalToLocalEnd(da, potential, INSERT_VALUES, pLocal));
  PetscCall(DMStagVecGetArray(da, pLocal, & arrP));

  /* Compute the exact solution to set boundary conditions */
  PetscCall(DMCreateGlobalVector(da, & x));
  FormExactSolution(t, ts, & x, user);
  PetscCall(DMGetLocalVector(da, & bcLocal));
  PetscCall(DMGlobalToLocalBegin(da, x, INSERT_VALUES, bcLocal));
  PetscCall(DMGlobalToLocalEnd(da, x, INSERT_VALUES, bcLocal));
  PetscCall(DMStagVecGetArrayRead(da, bcLocal, & arrx));
  PetscCall(TSGetTimeStep(ts, & dt));
  {
    /* Compute the gradient of EP */
    PetscCall(VecDuplicate(X, & GradEP));
    PetscCall(VecCopy(X, GradEP));
    FormDiscreteGradientEP_noMat(ts, X, GradEP, user);
    PetscCall(DMGetLocalVector(da, & GradEPLocal));
    PetscCall(DMGlobalToLocalBegin(da, GradEP, INSERT_VALUES, GradEPLocal));
    PetscCall(DMGlobalToLocalEnd(da, GradEP, INSERT_VALUES, GradEPLocal));
    PetscCall(DMStagVecGetArrayRead(da, GradEPLocal, & arrGradEP));
  }
  /* niv feeds n_i dV/dt (relaxation) and the advective inertia term (production);
     GradV1/2/3 feed only the latter. Production builds both whether or not
     user->inertia is set: with the flag off the term is still evaluated, multiplied
     by 0.0, because removing it changes the result in the last bit (see the comment
     at the f1 statements). Relaxation never reads GradV, so it skips that work. */
  const PetscBool build_niv   = PETSC_TRUE;
  const PetscBool build_gradV = (PetscBool) (!terms->f1_inertia);
  if (build_gradV) {
#line 3068
    /* Compute the gradient of V */
    PetscCall(VecDuplicate(X, & F1));
    PetscCall(VecCopy(X, F1));
    PetscCall(VecDuplicate(X, & F2));
    PetscCall(VecCopy(X, F2));
    PetscCall(VecDuplicate(X, & F3));
    PetscCall(VecCopy(X, F3));
    FormDiscreteGradientVectorField(ts, X, F1, F2, F3, user);
  }
  {
    /* Compute the vector laplacian of V */
    PetscCall(VecDuplicate(X, & LapV));
    PetscCall(VecZeroEntries(LapV));
    ApplyVectorLaplacian(ts, X, LapV, user);
    PetscCall(DMGetLocalVector(da, & LapVLocal));
    PetscCall(DMGlobalToLocalBegin(da, LapV, INSERT_VALUES, LapVLocal));
    PetscCall(DMGlobalToLocalEnd(da, LapV, INSERT_VALUES, LapVLocal));
    PetscCall(DMStagVecGetArrayRead(da, LapVLocal, & arrLapV));
  }
  /* Compute the projection vectors */
  /* P_{c->v}(ni) */
  if (build_niv) {
#line 3089
  PetscCall(DMCreateGlobalVector(da, & niv));
  CellToVertexProjectionScalar(ts, X, niv, user);
  PetscCall(DMGetLocalVector(da, & nivLocal));
  PetscCall(DMGlobalToLocalBegin(da, niv, INSERT_VALUES, nivLocal));
  PetscCall(DMGlobalToLocalEnd(da, niv, INSERT_VALUES, nivLocal));
  PetscCall(DMStagVecGetArrayRead(da, nivLocal, & arrniv));
  }
#line 3095
  /* P_{c->f}(ni) */
  PetscCall(DMCreateGlobalVector(da, & nif));
  CellToFaceProjection(ts, X, nif, user);
  PetscCall(DMGetLocalVector(da, & nifLocal));
  PetscCall(DMGlobalToLocalBegin(da, nif, INSERT_VALUES, nifLocal));
  PetscCall(DMGlobalToLocalEnd(da, nif, INSERT_VALUES, nifLocal));
  PetscCall(DMStagVecGetArrayRead(da, nifLocal, & arrnif));
  /* P_{f->v}(B) */
  PetscCall(DMCreateGlobalVector(da, & Bv));
  FaceToVertexProjection(ts, X, Bv, user);
  PetscCall(DMGetLocalVector(da, & BvLocal));
  PetscCall(DMGlobalToLocalBegin(da, Bv, INSERT_VALUES, BvLocal));
  PetscCall(DMGlobalToLocalEnd(da, Bv, INSERT_VALUES, BvLocal));
  PetscCall(DMStagVecGetArrayRead(da, BvLocal, & arrBv));
  /* P_{e->v}(der_curl_no_mp(B)) */
  PetscCall(DMCreateGlobalVector(da, & curlBv));
  PetscCall(VecZeroEntries(curlBv));
  Vec curlB;
  PetscCall(VecDuplicate(X, & curlB));
  PetscCall(VecCopy(X, curlB));
  FormDerivedCurlnomp(ts, X, curlB, user); //This updates only the E field part in curlB by computing the derived mimetic curl operator applied to B field of X that does not include material properties
  EdgeToVertexProjection(ts, curlB, curlBv, user);
  //PetscBarrier((PetscObject) curlB);
  PetscCall(VecDestroy( & curlB));
  PetscCall(DMGetLocalVector(da, & curlBvLocal));
  PetscCall(DMGlobalToLocalBegin(da, curlBv, INSERT_VALUES, curlBvLocal));
  PetscCall(DMGlobalToLocalEnd(da, curlBv, INSERT_VALUES, curlBvLocal));
  PetscCall(DMStagVecGetArrayRead(da, curlBvLocal, & arrcurlBv));
  /* P_{e->v}(prim_grad(V)) */
  if (build_gradV) {
#line 3124
  PetscCall(DMCreateGlobalVector(da, & GradV1));
  EdgeToVertexProjection(ts, F1, GradV1, user);
  PetscCall(VecDestroy( & F1));
  PetscCall(DMGetLocalVector(da, & GradV1Local));
  PetscCall(DMGlobalToLocalBegin(da, GradV1, INSERT_VALUES, GradV1Local));
  PetscCall(DMGlobalToLocalEnd(da, GradV1, INSERT_VALUES, GradV1Local));
  PetscCall(DMStagVecGetArrayRead(da, GradV1Local, & arrGradV1));

  PetscCall(DMCreateGlobalVector(da, & GradV3));
  EdgeToVertexProjection(ts, F3, GradV3, user);
  PetscCall(VecDestroy( & F3));
  PetscCall(DMGetLocalVector(da, & GradV3Local));
  PetscCall(DMGlobalToLocalBegin(da, GradV3, INSERT_VALUES, GradV3Local));
  PetscCall(DMGlobalToLocalEnd(da, GradV3, INSERT_VALUES, GradV3Local));
  PetscCall(DMStagVecGetArrayRead(da, GradV3Local, & arrGradV3));

  PetscCall(DMCreateGlobalVector(da, & GradV2));
  EdgeToVertexProjection(ts, F2, GradV2, user);
  PetscCall(VecDestroy( & F2));
  PetscCall(DMGetLocalVector(da, & GradV2Local));
  PetscCall(DMGlobalToLocalBegin(da, GradV2, INSERT_VALUES, GradV2Local));
  PetscCall(DMGlobalToLocalEnd(da, GradV2, INSERT_VALUES, GradV2Local));
  PetscCall(DMStagVecGetArrayRead(da, GradV2Local, & arrGradV2));
  }
#line 3147
  /* Compute the reconstruction vectors */
  /* R_{v->f}(V) */
  PetscCall(DMCreateGlobalVector(da, & Vf));
  VertexToFaceReconstruction(ts, X, Vf, user);
  PetscCall(DMGetLocalVector(da, & VfLocal));
  PetscCall(DMGlobalToLocalBegin(da, Vf, INSERT_VALUES, VfLocal));
  PetscCall(DMGlobalToLocalEnd(da, Vf, INSERT_VALUES, VfLocal));
  PetscCall(DMStagVecGetArrayRead(da, VfLocal, & arrVf));
  /* R_{v->e}(VxP_{f->v}(B)) */
  PetscCall(DMCreateGlobalVector(da, & VxB));
  PetscCall(VecZeroEntries(VxB));
  VertexCrossProduct(ts, X, Bv, VxB, user);
  PetscCall(DMCreateGlobalVector(da, & VxBe));
  VertexToEdgeReconstruction(ts, VxB, VxBe, user);
  PetscCall(VecDestroy( & VxB));
  PetscCall(DMGetLocalVector(da, & VxBeLocal));
  PetscCall(DMGlobalToLocalBegin(da, VxBe, INSERT_VALUES, VxBeLocal));
  PetscCall(DMGlobalToLocalEnd(da, VxBe, INSERT_VALUES, VxBeLocal));
  PetscCall(DMStagVecGetArrayRead(da, VxBeLocal, & arrVxBe));

#line 3177

  /* Compute function over the locally owned part of the grid */
  /* f1(V,EP,tau,B,ni) . e_r = - (P_{e->v}(der_mim_curl_no_mp(B)) x P_{f->v}(B) + Re^{-1} (\nabla^2 V)) . e_r ; on plasma vertices
     f1(V,EP,tau,B,ni) . e_z = - (P_{e->v}(der_mim_curl_no_mp(B)) x P_{f->v}(B) + Re^{-1} (\nabla^2 V)) . e_z ; on plasma vertices
     f1(V,EP,tau,B,ni) . e_phi = V. P_{f->v}(B) ; on plasma vertices
     f1(V,EP,tau,B,ni) = V; on other vertices
     f2(V,EP,tau,B,ni) = - derived_mimetic_div(primary_mimetic_grad(EP)) - (1/(L0*V_A)) * derived_mimetic_div(derived_mimetic_curl(B)) = - derived_mimetic_div(primary_mimetic_grad(EP)) - derived_mimetic_div(derived_mimetic_curl_2(B)); on all vertices
     f2(V,EP,tau,B,ni) += derived_mimetic_div(R_ve(VxP_fv(B))); on plasma vertices
     f3(V,EP,tau,B,ni) = tau - primary_mimetic_grad(EP) - (1/(L0*V_A)) * derived_mimetic_curl(B) = tau - primary_mimetic_grad(EP) - (1/(L0*V_A)) * primary_mimetic_curl^T (beta_f B)/beta_e ≡ tau - primary_mimetic_grad(EP) - (1/(L0*V_A)) * M_e^{-1} Curl^T M_f B ; on all edges
     f3(V,EP,tau,B,ni) = tau - primary_mimetic_grad(EP) - derived_mimetic_curl2(B) ; on all edges
     f3(V,EP,tau,B,ni) += R_ve(VxP_fv(B)) ; on inner plasma edges
     f4(V,EP,tau,B,ni) = dB/dt + primary_mimetic_curl(tau); on all faces
     f5(V,EP,tau,B,ni) = dni/dt ; in all cells
     f5(V,EP,tau,B,ni) += primary_mimetic_divergence(P_{c->f}(ni)*R_{v->f}(V)); in plasma cells
     */
  PetscCall(DMGetLocalVector(da, & fLocal));
  PetscCall(DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal));
  PetscCall(DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal));
  PetscCall(DMStagVecGetArray(da, fLocal, & arrF));

  PetscCall(DMGetLocalVector(da, & xLocal));
  PetscCall(DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMStagVecGetArrayRead(da, xLocal, & arrX));

  PetscCall(DMGetLocalVector(da, & xdotLocal));
  PetscCall(DMGlobalToLocalBegin(da, Xdot, INSERT_VALUES, xdotLocal));
  PetscCall(DMGlobalToLocalEnd(da, Xdot, INSERT_VALUES, xdotLocal));
  PetscCall(DMStagVecGetArrayRead(da, xdotLocal, & arrXdot));

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        cellvolume = MFD_CellVolume(arrCoord, er, ephi, ez, N, user -> dphi, icBrm, icBphim, icBzm, icBrp, icBphip, icBzp);
#line 3217

        MFD_CellEdgeLengths(arrCoorda, er, ephi, ez, N, user -> dphi,
          icrmphimzm, icrpphimzm, icrmphipzm, icrpphipzm, icrmphimzp, icrpphimzp, icrmphipzp, icrpphipzp,
          & rmzmedgelength, & rmphimedgelength, & rmzpedgelength, & rmphipedgelength, & phimzmedgelength, & phimzpedgelength,
          & rpphimedgelength, & rpzmedgelength, & phipzmedgelength, & rpzpedgelength, & rpphipedgelength, & phipzpedgelength);
#line 3257

        /* Set boundary conditions for tau field */
        /* f3(V,EP,tau,B,ni) = (tau - Eboundarycondition) */
        if (er == 0 || ez == 0) {
          arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivErmzm] - arrx[ez][ephi][er][ivErmzm]);
        }
        if (!(user -> phibtype)) {
          if (er == 0 || ephi == 0) {
            arrF[ez][ephi][er][ivErmphim] = (arrX[ez][ephi][er][ivErmphim] - arrx[ez][ephi][er][ivErmphim]);
            if (user -> debug) {
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Ermphim,%d,%d,%d) = %g\n", er, ephi, ez, (double) arrF[ez][ephi][er][ivErmphim]));
            }
          }
          if (ephi == 0 || ez == 0) {
            arrF[ez][ephi][er][ivEphimzm] = (arrX[ez][ephi][er][ivEphimzm] - arrx[ez][ephi][er][ivEphimzm]);
          }
        } else {
          if (er == 0) {
            arrF[ez][ephi][er][ivErmphim] = (arrX[ez][ephi][er][ivErmphim] - arrx[ez][ephi][er][ivErmphim]);
            if (user -> debug) {
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Ermphim,%d,%d,%d) = %g\n", er, ephi, ez, (double) arrF[ez][ephi][er][ivErmphim]));
            }
          }
          if (ez == 0) {
            arrF[ez][ephi][er][ivEphimzm] = (arrX[ez][ephi][er][ivEphimzm] - arrx[ez][ephi][er][ivEphimzm]);
          }
        }
        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivErpzm] = (arrX[ez][ephi][er][ivErpzm] - arrx[ez][ephi][er][ivErpzm]);
          arrF[ez][ephi][er][ivErpphim] = (arrX[ez][ephi][er][ivErpphim] - arrx[ez][ephi][er][ivErpphim]);
        }
        if (!(user -> phibtype)) {
          if (ephi == N[1] - 1) {
            arrF[ez][ephi][er][ivEphipzm] = (arrX[ez][ephi][er][ivEphipzm] - arrx[ez][ephi][er][ivEphipzm]);
            arrF[ez][ephi][er][ivErmphip] = (arrX[ez][ephi][er][ivErmphip] - arrx[ez][ephi][er][ivErmphip]);
          }
        }
        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivErmzp] = (arrX[ez][ephi][er][ivErmzp] - arrx[ez][ephi][er][ivErmzp]);
          arrF[ez][ephi][er][ivEphimzp] = (arrX[ez][ephi][er][ivEphimzp] - arrx[ez][ephi][er][ivEphimzp]);
        }
        if (!(user -> phibtype)) {
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrF[ez][ephi][er][ivErpphip] = (arrX[ez][ephi][er][ivErpphip] - arrx[ez][ephi][er][ivErpphip]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrF[ez][ephi][er][ivEphipzp] = (arrX[ez][ephi][er][ivEphipzp] - arrx[ez][ephi][er][ivEphipzp]);
          }
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrF[ez][ephi][er][ivErpzp] = (arrX[ez][ephi][er][ivErpzp] - arrx[ez][ephi][er][ivErpzp]);
        }

        /* f1(V,EP,tau,B,ni) . e_r = - (P_{e->v}(der_mim_curl_no_mp(B)) x P_{f->v}(B) + Re^{-1}(\nabla^2 V)) . e_r ; on plasma vertices
           f1(V,EP,tau,B,ni) . e_z = - (P_{e->v}(der_mim_curl_no_mp(B)) x P_{f->v}(B) + Re^{-1}(\nabla^2 V)) . e_z ; on plasma vertices
           f1(V,EP,tau,B,ni) . e_phi = V . P_{f->v}(B) ; on plasma vertices
           f1(V,EP,tau,B,ni) = V; on other vertices
           f5(V,EP,tau,B,ni) = dni/dt ; in all cells
           f5(V,EP,tau,B,ni) += primary_mimetic_divergence(P_{c->f}(ni)*R_{v->f}(V)); in plasma cells
           */
        arrF[ez][ephi][er][ivn] = arrXdot[ez][ephi][er][ivn];

        if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7) {
          arrF[ez][ephi][er][ivn] += (-surface(er, ephi, ez, LEFT, user) / cellvolume) * arrnif[ez][ephi][er][ivBrm] * arrVf[ez][ephi][er][ivBrm] + (surface(er, ephi, ez, RIGHT, user) / cellvolume) * arrnif[ez][ephi][er][ivBrp] * arrVf[ez][ephi][er][ivBrp] + (-surface(er, ephi, ez, DOWN, user) / cellvolume) * arrnif[ez][ephi][er][ivBphim] * arrVf[ez][ephi][er][ivBphim] + (surface(er, ephi, ez, UP, user) / cellvolume) * arrnif[ez][ephi][er][ivBphip] * arrVf[ez][ephi][er][ivBphip] + (-surface(er, ephi, ez, BACK, user) / cellvolume) * arrnif[ez][ephi][er][ivBzm] * arrVf[ez][ephi][er][ivBzm] + (surface(er, ephi, ez, FRONT, user) / cellvolume) * arrnif[ez][ephi][er][ivBzp] * arrVf[ez][ephi][er][ivBzp];
        }

        if (er > 0 && ez > 0 && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7 && fabs(user -> dataC[er + ephi * N[0] + (ez - 1) * N[1] * N[0]] - 1.5) < 0.7 && fabs(user -> dataC[er - 1 + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7 && fabs(user -> dataC[er - 1 + ephi * N[0] + (ez - 1) * N[1] * N[0]] - 1.5) < 0.7) {
          /* f1_inertia selects whole statements (never a 0/1 factor) so each variant keeps
             its original expression tree, and hence its FMA contraction, bit for bit. */
          if (terms->f1_inertia) {
            arrF[ez][ephi][er][ivVrmphimzm[0]] = - (arrcurlBv[ez][ephi][er][ivVrmphimzm[1]] * arrBv[ez][ephi][er][ivVrmphimzm[2]] - arrcurlBv[ez][ephi][er][ivVrmphimzm[2]] * arrBv[ez][ephi][er][ivVrmphimzm[1]]) + arrniv[ez][ephi][er][ivVrmphimzm[0]] * arrXdot[ez][ephi][er][ivVrmphimzm[0]] - (1.0 / user->Re) * arrLapV[ez][ephi][er][ivVrmphimzm[0]];

            arrF[ez][ephi][er][ivVrmphimzm[2]] = - (arrcurlBv[ez][ephi][er][ivVrmphimzm[0]] * arrBv[ez][ephi][er][ivVrmphimzm[1]] - arrcurlBv[ez][ephi][er][ivVrmphimzm[1]] * arrBv[ez][ephi][er][ivVrmphimzm[0]]) + arrniv[ez][ephi][er][ivVrmphimzm[0]] * arrXdot[ez][ephi][er][ivVrmphimzm[2]] - (1.0 / user->Re) * arrLapV[ez][ephi][er][ivVrmphimzm[2]];
          } else {
            /* Advective inertia n_i (V.grad)V, enabled by MHD_Config/inertia (default 0).
             * The default branch keeps the original statement, term multiplied by 0.0,
             * verbatim: deleting the term instead lets GCC contract the remaining
             * products into FMAs differently, which moves 30 entries of the residual by
             * 1-14 ULP and so breaks bit-for-bit agreement with earlier runs. Do not
             * simplify it without re-baselining. */
            if (terms->f1_advection) {
              arrF[ez][ephi][er][ivVrmphimzm[0]] = - (arrcurlBv[ez][ephi][er][ivVrmphimzm[1]] * arrBv[ez][ephi][er][ivVrmphimzm[2]] - arrcurlBv[ez][ephi][er][ivVrmphimzm[2]] * arrBv[ez][ephi][er][ivVrmphimzm[1]]) + arrniv[ez][ephi][er][ivVrmphimzm[0]] * (arrX[ez][ephi][er][ivVrmphimzm[0]] * arrGradV1[ez][ephi][er][ivVrmphimzm[0]] + arrX[ez][ephi][er][ivVrmphimzm[1]] * arrGradV1[ez][ephi][er][ivVrmphimzm[1]] + arrX[ez][ephi][er][ivVrmphimzm[2]] * arrGradV1[ez][ephi][er][ivVrmphimzm[2]] - arrX[ez][ephi][er][ivVrmphimzm[1]] * arrX[ez][ephi][er][ivVrmphimzm[1]] / arrCoorda[ez][ephi][er][icrmphimzm[0]] ) - (1.0 / user->Re) * arrLapV[ez][ephi][er][ivVrmphimzm[0]];

              arrF[ez][ephi][er][ivVrmphimzm[2]] = - (arrcurlBv[ez][ephi][er][ivVrmphimzm[0]] * arrBv[ez][ephi][er][ivVrmphimzm[1]] - arrcurlBv[ez][ephi][er][ivVrmphimzm[1]] * arrBv[ez][ephi][er][ivVrmphimzm[0]]) + arrniv[ez][ephi][er][ivVrmphimzm[0]] * (arrX[ez][ephi][er][ivVrmphimzm[0]] * arrGradV3[ez][ephi][er][ivVrmphimzm[0]] + arrX[ez][ephi][er][ivVrmphimzm[1]] * arrGradV3[ez][ephi][er][ivVrmphimzm[1]] + arrX[ez][ephi][er][ivVrmphimzm[2]] * arrGradV3[ez][ephi][er][ivVrmphimzm[2]]) - (1.0 / user->Re) * arrLapV[ez][ephi][er][ivVrmphimzm[2]];
            } else {
              arrF[ez][ephi][er][ivVrmphimzm[0]] = - (arrcurlBv[ez][ephi][er][ivVrmphimzm[1]] * arrBv[ez][ephi][er][ivVrmphimzm[2]] - arrcurlBv[ez][ephi][er][ivVrmphimzm[2]] * arrBv[ez][ephi][er][ivVrmphimzm[1]]) + 0.0 * arrniv[ez][ephi][er][ivVrmphimzm[0]] * (arrX[ez][ephi][er][ivVrmphimzm[0]] * arrGradV1[ez][ephi][er][ivVrmphimzm[0]] + arrX[ez][ephi][er][ivVrmphimzm[1]] * arrGradV1[ez][ephi][er][ivVrmphimzm[1]] + arrX[ez][ephi][er][ivVrmphimzm[2]] * arrGradV1[ez][ephi][er][ivVrmphimzm[2]] - arrX[ez][ephi][er][ivVrmphimzm[1]] * arrX[ez][ephi][er][ivVrmphimzm[1]] / arrCoorda[ez][ephi][er][icrmphimzm[0]] ) - (1.0 / user->Re) * arrLapV[ez][ephi][er][ivVrmphimzm[0]];

              arrF[ez][ephi][er][ivVrmphimzm[2]] = - (arrcurlBv[ez][ephi][er][ivVrmphimzm[0]] * arrBv[ez][ephi][er][ivVrmphimzm[1]] - arrcurlBv[ez][ephi][er][ivVrmphimzm[1]] * arrBv[ez][ephi][er][ivVrmphimzm[0]]) + 0.0 * arrniv[ez][ephi][er][ivVrmphimzm[0]] * (arrX[ez][ephi][er][ivVrmphimzm[0]] * arrGradV3[ez][ephi][er][ivVrmphimzm[0]] + arrX[ez][ephi][er][ivVrmphimzm[1]] * arrGradV3[ez][ephi][er][ivVrmphimzm[1]] + arrX[ez][ephi][er][ivVrmphimzm[2]] * arrGradV3[ez][ephi][er][ivVrmphimzm[2]]) - (1.0 / user->Re) * arrLapV[ez][ephi][er][ivVrmphimzm[2]];
            }
          }
#line 3327

          arrF[ez][ephi][er][ivVrmphimzm[1]] = (arrX[ez][ephi][er][ivVrmphimzm[0]] * arrBv[ez][ephi][er][ivVrmphimzm[0]] + arrX[ez][ephi][er][ivVrmphimzm[1]] * arrBv[ez][ephi][er][ivVrmphimzm[1]] + arrX[ez][ephi][er][ivVrmphimzm[2]] * arrBv[ez][ephi][er][ivVrmphimzm[2]]) ;

        } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) >= 0.7) {
          arrF[ez][ephi][er][ivVrmphimzm[0]] = arrX[ez][ephi][er][ivVrmphimzm[0]] - arrx[ez][ephi][er][ivVrmphimzm[0]];
          arrF[ez][ephi][er][ivVrmphimzm[1]] = arrX[ez][ephi][er][ivVrmphimzm[1]] - arrx[ez][ephi][er][ivVrmphimzm[1]];
          arrF[ez][ephi][er][ivVrmphimzm[2]] = arrX[ez][ephi][er][ivVrmphimzm[2]] - arrx[ez][ephi][er][ivVrmphimzm[2]];
          arrF[ez][ephi][er][ivVrmphimzp[0]] = arrX[ez][ephi][er][ivVrmphimzp[0]] - arrx[ez][ephi][er][ivVrmphimzp[0]];
          arrF[ez][ephi][er][ivVrmphimzp[1]] = arrX[ez][ephi][er][ivVrmphimzp[1]] - arrx[ez][ephi][er][ivVrmphimzp[1]];
          arrF[ez][ephi][er][ivVrmphimzp[2]] = arrX[ez][ephi][er][ivVrmphimzp[2]] - arrx[ez][ephi][er][ivVrmphimzp[2]];
          arrF[ez][ephi][er][ivVrmphipzm[0]] = arrX[ez][ephi][er][ivVrmphipzm[0]] - arrx[ez][ephi][er][ivVrmphipzm[0]];
          arrF[ez][ephi][er][ivVrmphipzm[1]] = arrX[ez][ephi][er][ivVrmphipzm[1]] - arrx[ez][ephi][er][ivVrmphipzm[1]];
          arrF[ez][ephi][er][ivVrmphipzm[2]] = arrX[ez][ephi][er][ivVrmphipzm[2]] - arrx[ez][ephi][er][ivVrmphipzm[2]];
          arrF[ez][ephi][er][ivVrmphipzp[0]] = arrX[ez][ephi][er][ivVrmphipzp[0]] - arrx[ez][ephi][er][ivVrmphipzp[0]];
          arrF[ez][ephi][er][ivVrmphipzp[1]] = arrX[ez][ephi][er][ivVrmphipzp[1]] - arrx[ez][ephi][er][ivVrmphipzp[1]];
          arrF[ez][ephi][er][ivVrmphipzp[2]] = arrX[ez][ephi][er][ivVrmphipzp[2]] - arrx[ez][ephi][er][ivVrmphipzp[2]];
          arrF[ez][ephi][er][ivVrpphimzm[0]] = arrX[ez][ephi][er][ivVrpphimzm[0]] - arrx[ez][ephi][er][ivVrpphimzm[0]];
          arrF[ez][ephi][er][ivVrpphimzm[1]] = arrX[ez][ephi][er][ivVrpphimzm[1]] - arrx[ez][ephi][er][ivVrpphimzm[1]];
          arrF[ez][ephi][er][ivVrpphimzm[2]] = arrX[ez][ephi][er][ivVrpphimzm[2]] - arrx[ez][ephi][er][ivVrpphimzm[2]];
          arrF[ez][ephi][er][ivVrpphimzp[0]] = arrX[ez][ephi][er][ivVrpphimzp[0]] - arrx[ez][ephi][er][ivVrpphimzp[0]];
          arrF[ez][ephi][er][ivVrpphimzp[1]] = arrX[ez][ephi][er][ivVrpphimzp[1]] - arrx[ez][ephi][er][ivVrpphimzp[1]];
          arrF[ez][ephi][er][ivVrpphimzp[2]] = arrX[ez][ephi][er][ivVrpphimzp[2]] - arrx[ez][ephi][er][ivVrpphimzp[2]];
          arrF[ez][ephi][er][ivVrpphipzm[0]] = arrX[ez][ephi][er][ivVrpphipzm[0]] - arrx[ez][ephi][er][ivVrpphipzm[0]];
          arrF[ez][ephi][er][ivVrpphipzm[1]] = arrX[ez][ephi][er][ivVrpphipzm[1]] - arrx[ez][ephi][er][ivVrpphipzm[1]];
          arrF[ez][ephi][er][ivVrpphipzm[2]] = arrX[ez][ephi][er][ivVrpphipzm[2]] - arrx[ez][ephi][er][ivVrpphipzm[2]];
          arrF[ez][ephi][er][ivVrpphipzp[0]] = arrX[ez][ephi][er][ivVrpphipzp[0]] - arrx[ez][ephi][er][ivVrpphipzp[0]];
          arrF[ez][ephi][er][ivVrpphipzp[1]] = arrX[ez][ephi][er][ivVrpphipzp[1]] - arrx[ez][ephi][er][ivVrpphipzp[1]];
          arrF[ez][ephi][er][ivVrpphipzp[2]] = arrX[ez][ephi][er][ivVrpphipzp[2]] - arrx[ez][ephi][er][ivVrpphipzp[2]];
        } else {
          arrF[ez][ephi][er][ivVrmphimzm[0]] = arrX[ez][ephi][er][ivVrmphimzm[0]] - arrx[ez][ephi][er][ivVrmphimzm[0]];
          arrF[ez][ephi][er][ivVrmphimzm[1]] = arrX[ez][ephi][er][ivVrmphimzm[1]] - arrx[ez][ephi][er][ivVrmphimzm[1]];
          arrF[ez][ephi][er][ivVrmphimzm[2]] = arrX[ez][ephi][er][ivVrmphimzm[2]] - arrx[ez][ephi][er][ivVrmphimzm[2]];
        }

        /* f4(V,EP,tau,B,ni) = dB/dt + primary_mimetic_curl(tau) */
        arrF[ez][ephi][er][ivBrm] = arrXdot[ez][ephi][er][ivBrm] + (rmzmedgelength * arrX[ez][ephi][er][ivErmzm] - rmzpedgelength * arrX[ez][ephi][er][ivErmzp] + rmphipedgelength * arrX[ez][ephi][er][ivErmphip] - rmphimedgelength * arrX[ez][ephi][er][ivErmphim]) / surface(er, ephi, ez, LEFT, user); /* Left face */
        /* DEBUG PRINT*/
        if (user -> debug) {
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Brm) = %g\n", (double) arrF[ez][ephi][er][ivBrm]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "1st term in F(Brm) = %g\n", (double) arrX[ez][ephi][er][ivErmzm]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 1st term in F(Brm) = %g\n", (double) rmzmedgelength));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "2nd term in F(Brm) = %g\n", (double) arrX[ez][ephi][er][ivErmzp]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 2nd term in F(Brm) = %g\n", (double) rmzpedgelength));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "3rd term in F(Brm) = %g\n", (double) arrX[ez][ephi][er][ivErmphip]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 3rd term in F(Brm) = %g\n", (double) rmphipedgelength));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "4th term in F(Brm) = %g\n", (double) arrX[ez][ephi][er][ivErmphim]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 4th term in F(Brm) = %g\n", (double) rmphimedgelength));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "surface in F(Brm) = %g\n", (double) surface(er, ephi, ez, LEFT, user)));
        }
        arrF[ez][ephi][er][ivBphim] = arrXdot[ez][ephi][er][ivBphim] + (-phimzmedgelength * arrX[ez][ephi][er][ivEphimzm] + phimzpedgelength * arrX[ez][ephi][er][ivEphimzp] - rpphimedgelength * arrX[ez][ephi][er][ivErpphim] + rmphimedgelength * arrX[ez][ephi][er][ivErmphim]) / surface(er, ephi, ez, DOWN, user); /* Down face */
        /* DEBUG PRINT*/
        if (user -> debug) {
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Bphim) = %g\n", (double) arrF[ez][ephi][er][ivBphim]));
        }
        arrF[ez][ephi][er][ivBzm] = arrXdot[ez][ephi][er][ivBzm] + (phimzmedgelength * arrX[ez][ephi][er][ivEphimzm] - phipzmedgelength * arrX[ez][ephi][er][ivEphipzm] + rpzmedgelength * arrX[ez][ephi][er][ivErpzm] - rmzmedgelength * arrX[ez][ephi][er][ivErmzm]) / surface(er, ephi, ez, BACK, user); /* Back face */
        /* DEBUG PRINT*/
        if (user -> debug) {
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Bzm) = %g\n", (double) arrF[ez][ephi][er][ivBzm]));
        }

        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivBrp] = arrXdot[ez][ephi][er][ivBrp] + (rpzmedgelength * arrX[ez][ephi][er][ivErpzm] - rpzpedgelength * arrX[ez][ephi][er][ivErpzp] + rpphipedgelength * arrX[ez][ephi][er][ivErpphip] - rpphimedgelength * arrX[ez][ephi][er][ivErpphim]) / surface(er, ephi, ez, RIGHT, user); /* Right face */
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivBphip] = arrXdot[ez][ephi][er][ivBphip] + (-phipzmedgelength * arrX[ez][ephi][er][ivEphipzm] + phipzpedgelength * arrX[ez][ephi][er][ivEphipzp] - rpphipedgelength * arrX[ez][ephi][er][ivErpphip] + rmphipedgelength * arrX[ez][ephi][er][ivErmphip]) / surface(er, ephi, ez, UP, user); /* Up face */
        }

        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivBzp] = arrXdot[ez][ephi][er][ivBzp] + (phimzpedgelength * arrX[ez][ephi][er][ivEphimzp] - phipzpedgelength * arrX[ez][ephi][er][ivEphipzp] + rpzpedgelength * arrX[ez][ephi][er][ivErpzp] - rmzpedgelength * arrX[ez][ephi][er][ivErmzp]) / surface(er, ephi, ez, FRONT, user); /* Front face */
        }

        /* Adding source term to balance the PDE: f4(B,E) = f4(B,E) - primary_mimetic_curl(potential) */
        if(user->ictype == 11 && ephi == 0){
          arrP[ez][ephi][er][ivEphimzm] = - condu(er, ephi, ez, BACK_DOWN, user) * PetscSinReal((condu(er, ephi, ez, BACK_DOWN, user) / user->mu0) * t) * (arrCoord[ez][ephi][er][icEphimzm[2]] ) / user->mu0;
          arrP[ez][ephi][er][ivEphimzp] = - condu(er, ephi, ez, FRONT_DOWN, user) * PetscSinReal((condu(er, ephi, ez, FRONT_DOWN, user) / user->mu0) * t) * (arrCoord[ez][ephi][er][icEphimzp[2]] ) / user->mu0;
        }
        if(user->ictype == 11 && ephi == N[1]-1){
          arrP[ez][ephi][er][ivEphipzm] = - condu(er, ephi, ez, BACK_UP, user) * PetscSinReal((condu(er, ephi, ez, BACK_UP, user) / user->mu0) * t) * (arrCoord[ez][ephi][er][icEphipzm[2]] - arrCoord[ez][ephi][er][icEphipzm[0]] * PETSC_PI * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icEphipzm[0]]))) / user->mu0;
          if(ez == N[2]-1){
            arrP[ez][ephi][er][ivEphipzp] = - condu(er, ephi, ez, FRONT_UP, user) * PetscSinReal((condu(er, ephi, ez, FRONT_UP, user) / user->mu0) * t) * (arrCoord[ez][ephi][er][icEphipzp[2]] - arrCoord[ez][ephi][er][icEphipzp[0]] * PETSC_PI * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icEphipzp[0]]))) / user->mu0;
          }
        }

        arrF[ez][ephi][er][ivBrm] -= (rmzmedgelength * arrP[ez][ephi][er][ivErmzm] - rmzpedgelength * arrP[ez][ephi][er][ivErmzp] + rmphipedgelength * arrP[ez][ephi][er][ivErmphip] - rmphimedgelength * arrP[ez][ephi][er][ivErmphim]) / surface(er, ephi, ez, LEFT, user); /* Left face */
        /* DEBUG PRINT*/
        if (user -> debug) {
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Source(Brm) = %g\n", (double) arrF[ez][ephi][er][ivBrm]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "1st term in Source(Brm) = %g\n", (double) arrP[ez][ephi][er][ivErmzm]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 1st term in Source(Brm) = %g\n", (double) rmzmedgelength));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "2nd term in Source(Brm) = %g\n", (double) arrP[ez][ephi][er][ivErmzp]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 2nd term in Source(Brm) = %g\n", (double) rmzpedgelength));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "3rd term in Source(Brm) = %g\n", (double) arrP[ez][ephi][er][ivErmphip]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 3rd term in Source(Brm) = %g\n", (double) rmphipedgelength));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "4th term in Source(Brm) = %g\n", (double) arrP[ez][ephi][er][ivErmphim]));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 4th term in Source(Brm) = %g\n", (double) rmphimedgelength));
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "surface in Source(Brm) = %g\n", (double) surface(er, ephi, ez, LEFT, user)));
        }

        if (!(user -> phibtype)) {
          if (ephi == N[1] - 1) {
            arrF[ez][ephi][er][ivBphip] -= (-phipzmedgelength * arrP[ez][ephi][er][ivEphipzm] + phipzpedgelength * arrP[ez][ephi][er][ivEphipzp] - rpphipedgelength * arrP[ez][ephi][er][ivErpphip] + rmphipedgelength * arrP[ez][ephi][er][ivErmphip]) / surface(er, ephi, ez, UP, user); /* Up face */
            /* DEBUG PRINT*/
            if (user -> debug) {
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Source(Bphip) = %g\n", (double) arrF[ez][ephi][er][ivBphip]));
            }
          }
        }

        arrF[ez][ephi][er][ivBphim] -= (-phimzmedgelength * arrP[ez][ephi][er][ivEphimzm] + phimzpedgelength * arrP[ez][ephi][er][ivEphimzp] - rpphimedgelength * arrP[ez][ephi][er][ivErpphim] + rmphimedgelength * arrP[ez][ephi][er][ivErmphim]) / surface(er, ephi, ez, DOWN, user); /* Face down */
        /* DEBUG PRINT*/
        if (user -> debug) {
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Source(Bphim) = %g\n", (double) arrF[ez][ephi][er][ivBphim]));
        }

        arrF[ez][ephi][er][ivBzm] -= (phimzmedgelength * arrP[ez][ephi][er][ivEphimzm] - phipzmedgelength * arrP[ez][ephi][er][ivEphipzm] + rpzmedgelength * arrP[ez][ephi][er][ivErpzm] - rmzmedgelength * arrP[ez][ephi][er][ivErmzm]) / surface(er, ephi, ez, BACK, user); /* Back face */
        /* DEBUG PRINT*/
        if (user -> debug) {
          PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Source(Bzm) = %g\n", (double) arrF[ez][ephi][er][ivBzm]));
        }

        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivBrp] -= (rpzmedgelength * arrP[ez][ephi][er][ivErpzm] - rpzpedgelength * arrP[ez][ephi][er][ivErpzp] + rpphipedgelength * arrP[ez][ephi][er][ivErpphip] - rpphimedgelength * arrP[ez][ephi][er][ivErpphim]) / surface(er, ephi, ez, RIGHT, user); /* Right face */
          /* DEBUG PRINT*/
          if (user -> debug) {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Source(Brp) = %g\n", (double) arrF[ez][ephi][er][ivBrp]));
          }
        }
        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivBzp] -= (phimzpedgelength * arrP[ez][ephi][er][ivEphimzp] - phipzpedgelength * arrP[ez][ephi][er][ivEphipzp] + rpzpedgelength * arrP[ez][ephi][er][ivErpzp] - rmzpedgelength * arrP[ez][ephi][er][ivErmzp]) / surface(er, ephi, ez, FRONT, user); /* Front face */
          /* DEBUG PRINT*/
          if (user -> debug) {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Source(Bzp) = %g\n", (double) arrF[ez][ephi][er][ivBzp]));
          }
        }

        /* Set boundary conditions for EP field */
        /* f2(V,EP,tau,B,ni) = (EP - EPboundarycondition) */
        if (er == 0 || ez == 0 || (ephi == 0 && !(user -> phibtype))) {
          arrF[ez][ephi][er][ivVrmphimzm[3]] = arrX[ez][ephi][er][ivVrmphimzm[3]] - arrx[ez][ephi][er][ivVrmphimzm[3]];
        }
        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivVrpphimzm[3]] = arrX[ez][ephi][er][ivVrpphimzm[3]] - arrx[ez][ephi][er][ivVrpphimzm[3]];
        }
        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivVrmphimzp[3]] = arrX[ez][ephi][er][ivVrmphimzp[3]] - arrx[ez][ephi][er][ivVrmphimzp[3]];
        }
        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivVrmphipzm[3]] = arrX[ez][ephi][er][ivVrmphipzm[3]] - arrx[ez][ephi][er][ivVrmphipzm[3]];
        }
        if (ez == N[2] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivVrmphipzp[3]] = arrX[ez][ephi][er][ivVrmphipzp[3]] - arrx[ez][ephi][er][ivVrmphipzp[3]];
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrF[ez][ephi][er][ivVrpphimzp[3]] = arrX[ez][ephi][er][ivVrpphimzp[3]] - arrx[ez][ephi][er][ivVrpphimzp[3]];
        }
        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivVrpphipzm[3]] = arrX[ez][ephi][er][ivVrpphipzm[3]] - arrx[ez][ephi][er][ivVrpphipzm[3]];
        }
        if (er == N[0] - 1 && ez == N[2] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivVrpphipzp[3]] = arrX[ez][ephi][er][ivVrpphipzp[3]] - arrx[ez][ephi][er][ivVrpphipzp[3]];
        }
      }
    }
  }
  /* End of triple for loop */

#line 3499

  PetscCall(DMStagVecRestoreArrayRead(da, LapVLocal, & arrLapV));
  PetscCall(DMRestoreLocalVector(da, & LapVLocal));
  PetscCall(VecDestroy(& LapV));

  if (build_gradV) {
#line 3504
  PetscCall(DMStagVecRestoreArrayRead(da, GradV1Local, & arrGradV1));
  PetscCall(DMRestoreLocalVector(da, & GradV1Local));
  PetscCall(VecDestroy( & GradV1));

  PetscCall(DMStagVecRestoreArrayRead(da, GradV2Local, & arrGradV2));
  PetscCall(DMRestoreLocalVector(da, & GradV2Local));
  PetscCall(VecDestroy( & GradV2));

  PetscCall(DMStagVecRestoreArrayRead(da, GradV3Local, & arrGradV3));
  PetscCall(DMRestoreLocalVector(da, & GradV3Local));
  PetscCall(VecDestroy( & GradV3));
  }
#line 3515

  PetscCall(DMStagVecRestoreArray(da, pLocal, & arrP));
  PetscCall(DMRestoreLocalVector(da, & pLocal));
  //PetscBarrier((PetscObject) potential);
  PetscCall(VecDestroy( & potential));

  PetscCall(DMStagVecRestoreArrayRead(da, nifLocal, & arrnif));
  PetscCall(DMRestoreLocalVector(da, & nifLocal));
  PetscCall(VecDestroy( & nif));

  if (build_niv) {
#line 3525
  PetscCall(DMStagVecRestoreArrayRead(da, nivLocal, & arrniv));
  PetscCall(DMRestoreLocalVector(da, & nivLocal));
  PetscCall(VecDestroy( & niv));
  }
#line 3528

  PetscCall(DMStagVecRestoreArrayRead(da, BvLocal, & arrBv));
  PetscCall(DMRestoreLocalVector(da, & BvLocal));
  PetscCall(VecDestroy( & Bv));

  PetscCall(DMStagVecRestoreArrayRead(da, VfLocal, & arrVf));
  PetscCall(DMRestoreLocalVector(da, & VfLocal));
  PetscCall(VecDestroy( & Vf));

  PetscCall(DMStagVecRestoreArrayRead(da, curlBvLocal, & arrcurlBv));
  PetscCall(DMRestoreLocalVector(da, & curlBvLocal));
  PetscCall(VecDestroy( & curlBv));

  if (user -> phibtype) {
    PetscCall(DMStagVecRestoreArray(da, fLocal, & arrF));
    PetscCall(DMLocalToGlobal(da, fLocal, INSERT_VALUES, F));
    /*DMRestoreLocalVector(da,&fLocal);*/

    /*DMGetLocalVector(da,&fLocal);*/
    PetscCall(DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(da, fLocal, & arrF));
  }

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        MFD_CellEdgeLengths(arrCoorda, er, ephi, ez, N, user -> dphi,
          icrmphimzm, icrpphimzm, icrmphipzm, icrpphipzm, icrmphimzp, icrpphimzp, icrmphipzp, icrpphipzp,
          & rmzmedgelength, & rmphimedgelength, & rmzpedgelength, & rmphipedgelength, & phimzmedgelength, & phimzpedgelength,
          & rpphimedgelength, & rpzmedgelength, & phipzmedgelength, & rpzpedgelength, & rpphipedgelength, & phipzpedgelength);
#line 3595


        /* f3(V,EP,tau,B,ni) = tau - derived_mimetic_curl2(B) + (eta/(V_A*B_0)) j_RE; for all inner edges */
        if (!(er == 0 || ez == 0)) {
          if (terms->f3_resistive && terms->f3_jre) {
            arrF[ez][ephi][er][ivErmzm] = arrX[ez][ephi][er][ivErmzm] - ((arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
                  arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
                  arrX[ez - 1][ephi][er][ivBrm] * betaf(er, ephi, ez - 1, LEFT, user) / surface(er, ephi, ez - 1, LEFT, user) +
                  arrX[ez][ephi][er - 1][ivBzm] * betaf(er - 1, ephi, ez, BACK, user) / surface(er - 1, ephi, ez, BACK, user)) * rmzmedgelength / betae2(er, ephi, ez, BACK_LEFT, user))
              + .25 * (
                  jre.data[(er    ) * jre.stride0 +  (ez    ) * jre.stride1 + 1 * jre.stride2] +
                  jre.data[(er - 1) * jre.stride0 +  (ez    ) * jre.stride1 + 1 * jre.stride2] +
                  jre.data[(er    ) * jre.stride0 +  (ez - 1) * jre.stride1 + 1 * jre.stride2] +
                  jre.data[(er - 1) * jre.stride0 +  (ez - 1) * jre.stride1 + 1 * jre.stride2]
                );
          } else if (terms->f3_resistive) { /* not used by any wrapper; NOT covered by the t1 baseline */
            arrF[ez][ephi][er][ivErmzm] = arrX[ez][ephi][er][ivErmzm] - ((arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
                  arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
                  arrX[ez - 1][ephi][er][ivBrm] * betaf(er, ephi, ez - 1, LEFT, user) / surface(er, ephi, ez - 1, LEFT, user) +
                  arrX[ez][ephi][er - 1][ivBzm] * betaf(er - 1, ephi, ez, BACK, user) / surface(er - 1, ephi, ez, BACK, user)) * rmzmedgelength / betae2(er, ephi, ez, BACK_LEFT, user));
          } else if (terms->f3_jre) { /* not used by any wrapper; NOT covered by the t1 baseline */
            arrF[ez][ephi][er][ivErmzm] = arrX[ez][ephi][er][ivErmzm]
              + .25 * (
                  jre.data[(er    ) * jre.stride0 +  (ez    ) * jre.stride1 + 1 * jre.stride2] +
                  jre.data[(er - 1) * jre.stride0 +  (ez    ) * jre.stride1 + 1 * jre.stride2] +
                  jre.data[(er    ) * jre.stride0 +  (ez - 1) * jre.stride1 + 1 * jre.stride2] +
                  jre.data[(er - 1) * jre.stride0 +  (ez - 1) * jre.stride1 + 1 * jre.stride2]
                );
          } else {
            arrF[ez][ephi][er][ivErmzm] = arrX[ez][ephi][er][ivErmzm];
          }
#line 3609
          /* DEBUG PRINT*/
          if (user -> debug) {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Ermzm) = %E\n", (double) arrF[ez][ephi][er][ivErmzm]));
          }
        }
        if (!(user -> phibtype)) {
          if (!(ephi == 0 || ez == 0)) {
            if (terms->f3_resistive && terms->f3_jre) {
              arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm] - ((-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                    arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
                    arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
                    arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae2(er, ephi, ez, BACK_DOWN, user))
                + .5 * (
                    jre.data[er * jre.stride0 + (ez    ) * jre.stride1] +
                    jre.data[er * jre.stride0 + (ez - 1) * jre.stride1]
                  );
            } else if (terms->f3_resistive) { /* not used by any wrapper; NOT covered by the t1 baseline */
              arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm] - ((-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                    arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
                    arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
                    arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae2(er, ephi, ez, BACK_DOWN, user));
            } else if (terms->f3_jre) { /* not used by any wrapper; NOT covered by the t1 baseline */
              arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm]
                + .5 * (
                    jre.data[er * jre.stride0 + (ez    ) * jre.stride1] +
                    jre.data[er * jre.stride0 + (ez - 1) * jre.stride1]
                  );
            } else {
              arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm];
            }
#line 3624
            /* DEBUG PRINT*/
            if (user -> debug) {
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Ephimzm) = %g\n", (double) arrF[ez][ephi][er][ivEphimzm]));
            }
          }
          if (!(er == 0 || ephi == 0)) {
            if (terms->f3_resistive && terms->f3_jre) {
              arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim] - ((-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
                    arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                    arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
                    arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user))
                + .5 * (
                    jre.data[(er - 1) * jre.stride0 + ez * jre.stride1 + 2 * jre.stride2] +
                    jre.data[(er    ) * jre.stride0 + ez * jre.stride1 + 2 * jre.stride2]
                  );
            } else if (terms->f3_resistive) { /* not used by any wrapper; NOT covered by the t1 baseline */
              arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim] - ((-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
                    arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                    arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
                    arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user));
            } else if (terms->f3_jre) { /* not used by any wrapper; NOT covered by the t1 baseline */
              arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim]
                + .5 * (
                    jre.data[(er - 1) * jre.stride0 + ez * jre.stride1 + 2 * jre.stride2] +
                    jre.data[(er    ) * jre.stride0 + ez * jre.stride1 + 2 * jre.stride2]
                  );
            } else {
              arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim];
            }
#line 3638

            /* DEBUG PRINT*/
            if (user -> debug) {
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Ermphim,%d,%d,%d) = %g\n", er, ephi, ez, (double) arrF[ez][ephi][er][ivErmphim]));

              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "1st term in F(Ermphim) = %E\n", (double)(-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user)));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B in 1st term of F(Ermphim) = %E\n", (double)(arrX[ez][ephi][er][ivBrm])));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "betaf in 1st term of F(Ermphim) = %E\n", (double)(betaf(er, ephi, ez, LEFT, user))));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "surface in 1st term of F(Ermphim) = %E\n", (double)(surface(er, ephi, ez, LEFT, user))));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 1st term of F(Ermphim) = %E\n", (double)(rmphimedgelength)));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "betae in 1st term in F(Ermphim) = %E\n", (double)(betae2(er, ephi, ez, DOWN_LEFT, user))));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "2nd term in F(Ermphim) = %E\n", (double)(+arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user)));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B in 2nd term of F(Ermphim) = %E\n", (double) arrX[ez][ephi][er][ivBphim]));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "beta/surf in 2nd term of F(Ermphim) = %E\n", (double)(betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user))));

              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "3rd term in F(Ermphim) = %E\n", (double)(+arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user)));

              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B in 3rd term of F(Ermphim) = %E\n", (double) arrX[ez][ephi - 1][er][ivBrm]));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B(Ermphim,%d,%d,%d) = %E\n", er, ephi - 1, ez, (double) arrX[ez][ephi - 1][er][ivBrm]));

              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "4th term in F(Ermphim) = %E\n", (double)(-arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user)));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B in 4th term of F(Ermphim) = %E\n", (double)(-arrX[ez][ephi][er - 1][ivBphim])));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "beta/surf in 4th term of F(Ermphim) = %E\n", (double)(betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user))));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Last term in F(Ermphim) = %E\n", (double)(-arrX[ez][ephi][er][ivErmphim])));
            }
          }
        } else {
          if (!(ez == 0)) {
            if (terms->f3_resistive && terms->f3_jre) {
              arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm] - ((-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                    arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
                    arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
                    arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae2(er, ephi, ez, BACK_DOWN, user))
                + .5 * (
                    jre.data[er * jre.stride0 + (ez    ) * jre.stride1 + 0 * jre.stride2] +
                    jre.data[er * jre.stride0 + (ez - 1) * jre.stride1 + 0 * jre.stride2]
                  );
            } else if (terms->f3_resistive) { /* not used by any wrapper; NOT covered by the t1 baseline */
              arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm] - ((-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                    arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
                    arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
                    arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae2(er, ephi, ez, BACK_DOWN, user));
            } else if (terms->f3_jre) { /* not used by any wrapper; NOT covered by the t1 baseline */
              arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm]
                + .5 * (
                    jre.data[er * jre.stride0 + (ez    ) * jre.stride1 + 0 * jre.stride2] +
                    jre.data[er * jre.stride0 + (ez - 1) * jre.stride1 + 0 * jre.stride2]
                  );
            } else {
              arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm];
            }
#line 3674
            /* DEBUG PRINT*/
            if (user -> debug) {
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Ephimzm) = %g\n", (double) arrF[ez][ephi][er][ivEphimzm]));
            }
          }
          if (!(er == 0)) {
            if (terms->f3_resistive && terms->f3_jre) {
              arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim] - ((-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
                    arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                    arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
                    arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user))
                + .5 * (
                    jre.data[(er - 1) * jre.stride0 + ez * jre.stride1 + 2 * jre.stride2] +
                    jre.data[(er    ) * jre.stride0 + ez * jre.stride1 + 2 * jre.stride2]
                  );
            } else if (terms->f3_resistive) { /* not used by any wrapper; NOT covered by the t1 baseline */
              arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim] - ((-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
                    arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                    arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
                    arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user));
            } else if (terms->f3_jre) { /* not used by any wrapper; NOT covered by the t1 baseline */
              arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim]
                + .5 * (
                    jre.data[(er - 1) * jre.stride0 + ez * jre.stride1 + 2 * jre.stride2] +
                    jre.data[(er    ) * jre.stride0 + ez * jre.stride1 + 2 * jre.stride2]
                  );
            } else {
              arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim];
            }
#line 3688
            /* DEBUG PRINT*/
            if (user -> debug) {
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Ermphim,%d,%d,%d) = %g\n", er, ephi, ez, (double) arrF[ez][ephi][er][ivErmphim]));

              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "1st term in F(Ermphim) = %E\n", (double)(-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user)));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B in 1st term of F(Ermphim) = %E\n", (double)(arrX[ez][ephi][er][ivBrm])));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "betaf in 1st term of F(Ermphim) = %E\n", (double)(betaf(er, ephi, ez, LEFT, user))));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "surface in 1st term of F(Ermphim) = %E\n", (double)(surface(er, ephi, ez, LEFT, user))));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "edge length in 1st term of F(Ermphim) = %E\n", (double)(rmphimedgelength)));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "betae in 1st term in F(Ermphim) = %E\n", (double)(betae2(er, ephi, ez, DOWN_LEFT, user))));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "2nd term in F(Ermphim) = %E\n", (double)(+arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user)));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B in 2nd term of F(Ermphim) = %E\n", (double) arrX[ez][ephi][er][ivBphim]));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "beta/surf in 2nd term of F(Ermphim) = %E\n", (double)(betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user))));

              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "3rd term in F(Ermphim) = %E\n", (double)(+arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user)));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B in 3rd term of F(Ermphim) = %E\n", (double) arrX[ez][ephi - 1][er][ivBrm]));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B(Ermphim,%d,%d,%d) = %E\n", er, ephi - 1, ez, (double) arrX[ez][ephi - 1][er][ivBrm]));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "4th term in F(Ermphim) = %E\n", (double)(-arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae2(er, ephi, ez, DOWN_LEFT, user)));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "B in 4th term of F(Ermphim) = %E\n", (double)(-arrX[ez][ephi][er - 1][ivBphim])));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "beta/surf in 4th term of F(Ermphim) = %E\n", (double)(betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user))));
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Last term in F(Ermphim) = %E\n", (double)(-arrX[ez][ephi][er][ivErmphim])));
            }
          }
        }
      }
    }
  }

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        //F3(V,EP,tau,B,ni) += R_ve(VxP_fv(B)) ; on inner plasma edges
        if (!(er == 0 || ez == 0)) {
          if ((fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7) && (fabs(user -> dataC[er - 1 + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7) && (fabs(user -> dataC[er + ephi * N[0] + (ez - 1) * N[1] * N[0]] - 1.5) < 0.7) && (fabs(user -> dataC[er - 1 + ephi * N[0] + (ez - 1) * N[1] * N[0]] - 1.5) < 0.7)) {
            arrF[ez][ephi][er][ivErmzm] += arrVxBe[ez][ephi][er][ivErmzm];
          }
        }
        if (!(user -> phibtype)) {
          if (!(ephi == 0 || ez == 0)) {
            if ((fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7) && (fabs(user -> dataC[er + ephi * N[0] + (ez - 1) * N[1] * N[0]] - 1.5) < 0.7)) {
              arrF[ez][ephi][er][ivEphimzm] += arrVxBe[ez][ephi][er][ivEphimzm];
            }
          }
          if (!(er == 0 || ephi == 0)) {
            if ((fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7) && (fabs(user -> dataC[er - 1 + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7)) {
              arrF[ez][ephi][er][ivErmphim] += arrVxBe[ez][ephi][er][ivErmphim];
            }
          }
        } else {
          if (!(ez == 0)) {
            if ((fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7) && (fabs(user -> dataC[er + ephi * N[0] + (ez - 1) * N[1] * N[0]] - 1.5) < 0.7)) {
              arrF[ez][ephi][er][ivEphimzm] += arrVxBe[ez][ephi][er][ivEphimzm];
            }
          }
          if (!(er == 0)) {
            if ((fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7) && (fabs(user -> dataC[er - 1 + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7)) {
              arrF[ez][ephi][er][ivErmphim] += arrVxBe[ez][ephi][er][ivErmphim];
            }
          }
        }
      }
    }
  }

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        //F3(V,EP,tau,B,ni) -= prim_grad(EP) ; on inner edges
        if (!(er == 0 || ez == 0)) {
          arrF[ez][ephi][er][ivErmzm] -= arrGradEP[ez][ephi][er][ivErmzm];
        }
        if (!(user -> phibtype)) {
          if (!(ephi == 0 || ez == 0)) {
            arrF[ez][ephi][er][ivEphimzm] -= arrGradEP[ez][ephi][er][ivEphimzm];
          }
          if (!(er == 0 || ephi == 0)) {
            arrF[ez][ephi][er][ivErmphim] -= arrGradEP[ez][ephi][er][ivErmphim];
          }
        } else {
          if (!(ez == 0)) {
            arrF[ez][ephi][er][ivEphimzm] -= arrGradEP[ez][ephi][er][ivEphimzm];
          }
          if (!(er == 0)) {
            arrF[ez][ephi][er][ivErmphim] -= arrGradEP[ez][ephi][er][ivErmphim];
          }
        }
        /* Set boundary conditions for velocity */
        /* f1(V,EP,tau,B,ni) = dV/dt */
        if (er == 0 || ez == 0) {
          arrF[ez][ephi][er][ivVrmphimzm[0]] = arrXdot[ez][ephi][er][ivVrmphimzm[0]];
          arrF[ez][ephi][er][ivVrmphimzm[1]] = arrXdot[ez][ephi][er][ivVrmphimzm[1]];
          arrF[ez][ephi][er][ivVrmphimzm[2]] = arrXdot[ez][ephi][er][ivVrmphimzm[2]];
          arrF[ez][ephi][er][ivVrmphipzm[0]] = arrXdot[ez][ephi][er][ivVrmphipzm[0]];
          arrF[ez][ephi][er][ivVrmphipzm[1]] = arrXdot[ez][ephi][er][ivVrmphipzm[1]];
          arrF[ez][ephi][er][ivVrmphipzm[2]] = arrXdot[ez][ephi][er][ivVrmphipzm[2]];
        }
        if (er == 0 || ez == N[2]-1) {
          arrF[ez][ephi][er][ivVrmphimzp[0]] = arrXdot[ez][ephi][er][ivVrmphimzp[0]];
          arrF[ez][ephi][er][ivVrmphimzp[1]] = arrXdot[ez][ephi][er][ivVrmphimzp[1]];
          arrF[ez][ephi][er][ivVrmphimzp[2]] = arrXdot[ez][ephi][er][ivVrmphimzp[2]];
          arrF[ez][ephi][er][ivVrmphipzp[0]] = arrXdot[ez][ephi][er][ivVrmphipzp[0]];
          arrF[ez][ephi][er][ivVrmphipzp[1]] = arrXdot[ez][ephi][er][ivVrmphipzp[1]];
          arrF[ez][ephi][er][ivVrmphipzp[2]] = arrXdot[ez][ephi][er][ivVrmphipzp[2]];
        }
        if (ez == 0 || er == N[0]-1) {
          arrF[ez][ephi][er][ivVrpphimzm[0]] = arrXdot[ez][ephi][er][ivVrpphimzm[0]];
          arrF[ez][ephi][er][ivVrpphimzm[1]] = arrXdot[ez][ephi][er][ivVrpphimzm[1]];
          arrF[ez][ephi][er][ivVrpphimzm[2]] = arrXdot[ez][ephi][er][ivVrpphimzm[2]];
          arrF[ez][ephi][er][ivVrpphipzm[0]] = arrXdot[ez][ephi][er][ivVrpphipzm[0]];
          arrF[ez][ephi][er][ivVrpphipzm[1]] = arrXdot[ez][ephi][er][ivVrpphipzm[1]];
          arrF[ez][ephi][er][ivVrpphipzm[2]] = arrXdot[ez][ephi][er][ivVrpphipzm[2]];
        }
        if (ez == N[2]-1 || er == N[0]-1) {
          arrF[ez][ephi][er][ivVrpphimzp[0]] = arrXdot[ez][ephi][er][ivVrpphimzp[0]];
          arrF[ez][ephi][er][ivVrpphimzp[1]] = arrXdot[ez][ephi][er][ivVrpphimzp[1]];
          arrF[ez][ephi][er][ivVrpphimzp[2]] = arrXdot[ez][ephi][er][ivVrpphimzp[2]];
          arrF[ez][ephi][er][ivVrpphipzp[0]] = arrXdot[ez][ephi][er][ivVrpphipzp[0]];
          arrF[ez][ephi][er][ivVrpphipzp[1]] = arrXdot[ez][ephi][er][ivVrpphipzp[1]];
          arrF[ez][ephi][er][ivVrpphipzp[2]] = arrXdot[ez][ephi][er][ivVrpphipzp[2]];
        }
      }
    }
  }

  PetscCall(DMStagVecRestoreArray(da, fLocal, & arrF));
  PetscCall(DMLocalToGlobal(da, fLocal, INSERT_VALUES, F));
  PetscCall(DMRestoreLocalVector(da,&fLocal));

  PetscCall(VecDuplicate(F, & Fcopy));
  PetscCall(VecCopy(F, Fcopy));
  //VecScale(Fcopy, -1.0);
  PetscCall(DMGetLocalVector(da, & FcopyLocal));
  PetscCall(DMGlobalToLocalBegin(da, Fcopy, INSERT_VALUES, FcopyLocal));
  PetscCall(DMGlobalToLocalEnd(da, Fcopy, INSERT_VALUES, FcopyLocal));
  PetscCall(DMStagVecGetArray(da, FcopyLocal, & arrFcopy));

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        //F3copy(V,EP,tau,B,ni) -= tau ; on inner edges
        if (!(er == 0 || ez == 0)) {
          arrFcopy[ez][ephi][er][ivErmzm] -= arrX[ez][ephi][er][ivErmzm];
        }
        if (!(user -> phibtype)) {
          if (!(ephi == 0 || ez == 0)) {
            arrFcopy[ez][ephi][er][ivEphimzm] -= arrX[ez][ephi][er][ivEphimzm];
          }
          if (!(er == 0 || ephi == 0)) {
            arrFcopy[ez][ephi][er][ivErmphim] -= arrX[ez][ephi][er][ivErmphim];
          }
        } else {
          if (!(ez == 0)) {
            arrFcopy[ez][ephi][er][ivEphimzm] -= arrX[ez][ephi][er][ivEphimzm];
          }
          if (!(er == 0)) {
            arrFcopy[ez][ephi][er][ivErmphim] -= arrX[ez][ephi][er][ivErmphim];
          }
        }
      }
    }
  }
  PetscCall(DMStagVecRestoreArray(da, FcopyLocal, & arrFcopy));
  PetscCall(DMLocalToGlobal(da, FcopyLocal, INSERT_VALUES, Fcopy));
  PetscCall(DMRestoreLocalVector(da, & FcopyLocal));

  ApplyDerivedDivergence(ts, Fcopy, F, user);

  /* Restore vectors */
  PetscCall(DMStagVecRestoreArrayRead(da, GradEPLocal, & arrGradEP));
  PetscCall(DMRestoreLocalVector(da, & GradEPLocal));
  PetscCall(DMStagVecRestoreArrayRead(da, VxBeLocal, & arrVxBe));
  PetscCall(DMRestoreLocalVector(da, & VxBeLocal));


  PetscCall(DMStagVecRestoreArrayRead(da, xLocal, & arrX));
  PetscCall(DMRestoreLocalVector(da, & xLocal));
  PetscCall(DMStagVecRestoreArrayRead(da, xdotLocal, & arrXdot));
  PetscCall(DMRestoreLocalVector(da, & xdotLocal));
  PetscCall(DMStagVecRestoreArrayRead(da, bcLocal, & arrx));
  PetscCall(DMRestoreLocalVector(da, & bcLocal));

  PetscCall(DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord));
  PetscCall(DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda));
  if (user -> debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F = \n"));
    PetscCall(VecView(F, PETSC_VIEWER_STDOUT_WORLD));
  }

  PetscCall(VecDestroy( & GradEP));
  PetscCall(VecDestroy( & VxBe));
  PetscCall(VecDestroy( & Fcopy));
  PetscCall(VecDestroy( & x));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));


  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Production residual, every time step (registered in mhd.c). */
PetscErrorCode FormIFunction_Vperp_viscosity(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void * ptr) {
  /* MHD_Config/inertia picks between two separately compiled specializations.
   * The flag must not reach vperp_residual as a runtime value: a runtime branch
   * there lets GCC share code across the alternatives, which changed 30 entries
   * of the default (inertia off) residual by 1-14 ULP. With the choice made here
   * and each call inlined on a constant, the default path compiles exactly as it
   * did before the flag existed. */
  static const MFD_ResidualTerms terms = {
    .f1_inertia = PETSC_FALSE, .f3_resistive = PETSC_TRUE, .f3_jre = PETSC_TRUE,
    .f1_advection = PETSC_FALSE, .label = "FormIFunction_Vperp_viscosity"};
  static const MFD_ResidualTerms terms_advection = {
    .f1_inertia = PETSC_FALSE, .f3_resistive = PETSC_TRUE, .f3_jre = PETSC_TRUE,
    .f1_advection = PETSC_TRUE, .label = "FormIFunction_Vperp_viscosity"};
  User * user = (User *) ptr;
  if (user->inertia) return vperp_residual(ts, t, X, Xdot, F, ptr, &terms_advection);
  return vperp_residual(ts, t, X, Xdot, F, ptr, &terms);
}

#line 4866

#line 5847

#line 6822

#line 7516

#line 8192

#line 9003

#line 9400

/* Residual for the initial-condition solve of the electrostatic potential EP and
 * the tau field, holding B, V and n_i fixed (their rows are dX/dt).
 *
 * FormIFunction_InitializeEP and FormIFunction_InitializeEP_halo were 423- and
 * 427-line copies of this body, differing only in the edge coefficient that
 * divides the derived curl in the tau (Ohm's law) rows:
 *
 *                       r-z edges (BACK_LEFT)    phi-z / r-phi edges
 *   InitializeEP        betae2                   betae2
 *   InitializeEP_halo   betaephi_isolcell        betaeperp2
 *
 * i.e. the halo variant gives the wall anisotropic resistivity -- toroidal on
 * r-z edges, poloidal on the others -- with the three hardcoded isolated cells
 * (open question Q1). `label` keeps each variant's PETSc log-event name. */
#line 9401
static PetscErrorCode initialize_ep(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void * ptr, mfd_edge_coeff beta_rz, mfd_edge_coeff beta_other, const char * label) {
  PetscFunctionBeginUser;

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister(label,classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));


  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  PetscScalar dt;
  Vec fLocal, xLocal, bcLocal, xdotLocal, pLocal;
  Vec GradEP, GradEPLocal, Fcopy, FcopyLocal;
  Vec x, potential;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icp[3] PETSC_UNUSED;
  PetscInt icBrp[3] PETSC_UNUSED, icBphip[3] PETSC_UNUSED, icBzp[3] PETSC_UNUSED, icBrm[3] PETSC_UNUSED, icBphim[3] PETSC_UNUSED, icBzm[3] PETSC_UNUSED;
  PetscInt icErmzm[3] PETSC_UNUSED, icErmzp[3] PETSC_UNUSED, icErpzm[3] PETSC_UNUSED, icErpzp[3] PETSC_UNUSED;

  PetscInt icEphimzm[3] PETSC_UNUSED, icEphipzm[3] PETSC_UNUSED, icEphimzp[3] PETSC_UNUSED, icEphipzp[3] PETSC_UNUSED;
  PetscInt icErmphim[3] PETSC_UNUSED, icErpphim[3] PETSC_UNUSED, icErmphip[3] PETSC_UNUSED, icErpphip[3] PETSC_UNUSED;
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  PetscInt ivn;

  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, ** ** arrP, ** ** arrx, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength, ** ** arrXdot, ** ** arrGradEP, ** ** arrFcopy;

  PetscCall(VecZeroEntries(F));
  PetscCall(TSGetDM(ts, & da));

  PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
  PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));

  MFD_Slots S;
  PetscCall(MFD_GetSlotsSolution(da, & S));
#line 9489
  PetscCall(DMGetCoordinateDM(da, & dmCoord));
  PetscCall(DMGetCoordinatesLocal(da, & coordLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord));
  PetscCall(DMGetCoordinateDM(coordDA, & dmCoorda));
  PetscCall(DMGetCoordinatesLocal(coordDA, & coordaLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda));
  PetscCall(MFD_GetSlotsCoords(dmCoord, dmCoorda, & S));
  MFD_UNPACK_SLOTS(S);
#line 9528
  /* Compute the source term potential for time-dependent manufactured solution */
  PetscCall(DMCreateGlobalVector(da, & potential));
  FormSourceTermPotential(ts, t, potential, user);
  PetscCall(DMGetLocalVector(da, & pLocal));
  PetscCall(DMGlobalToLocalBegin(da, potential, INSERT_VALUES, pLocal));
  PetscCall(DMGlobalToLocalEnd(da, potential, INSERT_VALUES, pLocal));
  PetscCall(DMStagVecGetArrayRead(da, pLocal, & arrP));

  /* Compute the exact solution to set boundary conditions */
  PetscCall(DMCreateGlobalVector(da, & x));
  FormExactSolution(t, ts, & x, user);
  PetscCall(DMGetLocalVector(da, & bcLocal));
  PetscCall(DMGlobalToLocalBegin(da, x, INSERT_VALUES, bcLocal));
  PetscCall(DMGlobalToLocalEnd(da, x, INSERT_VALUES, bcLocal));
  PetscCall(DMStagVecGetArrayRead(da, bcLocal, & arrx));
  PetscCall(TSGetTimeStep(ts, & dt));
  {
    /* Compute the gradient of EP */
    PetscCall(VecDuplicate(X, & GradEP));
    PetscCall(VecCopy(X, GradEP));
    FormDiscreteGradientEP_noMat(ts, X, GradEP, user);
    PetscCall(DMGetLocalVector(da, & GradEPLocal));
    PetscCall(DMGlobalToLocalBegin(da, GradEP, INSERT_VALUES, GradEPLocal));
    PetscCall(DMGlobalToLocalEnd(da, GradEP, INSERT_VALUES, GradEPLocal));
    PetscCall(DMStagVecGetArrayRead(da, GradEPLocal, & arrGradEP));
  }


  /* Compute function over the locally owned part of the grid */
  /* f1(V,EP,tau,B,ni) = dV/dt; on all vertices
     f2(V,EP,tau,B,ni) = -derived_mimetic_div(primary_mimetic_grad(EP)) - derived_mimetic_div(derived_mimetic_curl2(B)); on all vertices
     f3(V,EP,tau,B,ni) = tau - primary_mimetic_grad(EP) - (1/(L0*V_A)) * derived_mimetic_curl(B) = tau - primary_mimetic_grad(EP) - derived_mimetic_curl2(B); on all edges
     f4(V,EP,tau,B,ni) = dB/dt; on all faces
     f5(V,EP,tau,B,ni) = dni/dt; in all cells */
  PetscCall(DMGetLocalVector(da, & fLocal));
  PetscCall(DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal));
  PetscCall(DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal));
  PetscCall(DMStagVecGetArray(da, fLocal, & arrF));

  PetscCall(DMGetLocalVector(da, & xLocal));
  PetscCall(DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMStagVecGetArrayRead(da, xLocal, & arrX));

  PetscCall(DMGetLocalVector(da, & xdotLocal));
  PetscCall(DMGlobalToLocalBegin(da, Xdot, INSERT_VALUES, xdotLocal));
  PetscCall(DMGlobalToLocalEnd(da, Xdot, INSERT_VALUES, xdotLocal));
  PetscCall(DMStagVecGetArrayRead(da, xdotLocal, & arrXdot));

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

#line 9587

        MFD_CellEdgeLengths(arrCoorda, er, ephi, ez, N, user -> dphi,
          icrmphimzm, icrpphimzm, icrmphipzm, icrpphipzm, icrmphimzp, icrpphimzp, icrmphipzp, icrpphipzp,
          & rmzmedgelength, & rmphimedgelength, & rmzpedgelength, & rmphipedgelength, & phimzmedgelength, & phimzpedgelength,
          & rpphimedgelength, & rpzmedgelength, & phipzmedgelength, & rpzpedgelength, & rpphipedgelength, & phipzpedgelength);
#line 9627

        /* Set boundary conditions for tau field */
        /* f3(V,EP,tau,B,ni) = tau */
        if (er == 0 || ez == 0) {
          arrF[ez][ephi][er][ivErmzm] = arrX[ez][ephi][er][ivErmzm];
        }
        if (!(user -> phibtype)) {
          if (er == 0 || ephi == 0) {
            arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim];
            if (user -> debug) {
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Ermphim,%d,%d,%d) = %g\n", er, ephi, ez, (double) arrF[ez][ephi][er][ivErmphim]));
            }
          }
          if (ephi == 0 || ez == 0) {
            arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm];
          }
        } else {
          if (er == 0) {
            arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim];
            if (user -> debug) {
              PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F(Ermphim,%d,%d,%d) = %g\n", er, ephi, ez, (double) arrF[ez][ephi][er][ivErmphim]));
            }
          }
          if (ez == 0) {
            arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm];
          }
        }
        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivErpzm] = arrX[ez][ephi][er][ivErpzm];
          arrF[ez][ephi][er][ivErpphim] = arrX[ez][ephi][er][ivErpphim];
        }
        if (!(user -> phibtype)) {
          if (ephi == N[1] - 1) {
            arrF[ez][ephi][er][ivEphipzm] = arrX[ez][ephi][er][ivEphipzm];
            arrF[ez][ephi][er][ivErmphip] = arrX[ez][ephi][er][ivErmphip];
          }
        }
        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivErmzp] = arrX[ez][ephi][er][ivErmzp];
          arrF[ez][ephi][er][ivEphimzp] = arrX[ez][ephi][er][ivEphimzp];
        }
        if (!(user -> phibtype)) {
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrF[ez][ephi][er][ivErpphip] = arrX[ez][ephi][er][ivErpphip];
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrF[ez][ephi][er][ivEphipzp] = arrX[ez][ephi][er][ivEphipzp];
          }
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrF[ez][ephi][er][ivErpzp] = arrX[ez][ephi][er][ivErpzp];
        }

        /* f1(V,EP,tau,B,ni) = dV/dt; on all vertices
           f5(V,EP,tau,B,ni) = dni/dt ; in all cells
           */
        arrF[ez][ephi][er][ivn] = arrXdot[ez][ephi][er][ivn];

        arrF[ez][ephi][er][ivVrmphimzm[0]] = arrXdot[ez][ephi][er][ivVrmphimzm[0]];
        arrF[ez][ephi][er][ivVrmphimzm[1]] = arrXdot[ez][ephi][er][ivVrmphimzm[1]];
        arrF[ez][ephi][er][ivVrmphimzm[2]] = arrXdot[ez][ephi][er][ivVrmphimzm[2]];
        arrF[ez][ephi][er][ivVrmphimzp[0]] = arrXdot[ez][ephi][er][ivVrmphimzp[0]];
        arrF[ez][ephi][er][ivVrmphimzp[1]] = arrXdot[ez][ephi][er][ivVrmphimzp[1]];
        arrF[ez][ephi][er][ivVrmphimzp[2]] = arrXdot[ez][ephi][er][ivVrmphimzp[2]];
        arrF[ez][ephi][er][ivVrmphipzm[0]] = arrXdot[ez][ephi][er][ivVrmphipzm[0]];
        arrF[ez][ephi][er][ivVrmphipzm[1]] = arrXdot[ez][ephi][er][ivVrmphipzm[1]];
        arrF[ez][ephi][er][ivVrmphipzm[2]] = arrXdot[ez][ephi][er][ivVrmphipzm[2]];
        arrF[ez][ephi][er][ivVrmphipzp[0]] = arrXdot[ez][ephi][er][ivVrmphipzp[0]];
        arrF[ez][ephi][er][ivVrmphipzp[1]] = arrXdot[ez][ephi][er][ivVrmphipzp[1]];
        arrF[ez][ephi][er][ivVrmphipzp[2]] = arrXdot[ez][ephi][er][ivVrmphipzp[2]];
        arrF[ez][ephi][er][ivVrpphimzm[0]] = arrXdot[ez][ephi][er][ivVrpphimzm[0]];
        arrF[ez][ephi][er][ivVrpphimzm[1]] = arrXdot[ez][ephi][er][ivVrpphimzm[1]];
        arrF[ez][ephi][er][ivVrpphimzm[2]] = arrXdot[ez][ephi][er][ivVrpphimzm[2]];
        arrF[ez][ephi][er][ivVrpphimzp[0]] = arrXdot[ez][ephi][er][ivVrpphimzp[0]];
        arrF[ez][ephi][er][ivVrpphimzp[1]] = arrXdot[ez][ephi][er][ivVrpphimzp[1]];
        arrF[ez][ephi][er][ivVrpphimzp[2]] = arrXdot[ez][ephi][er][ivVrpphimzp[2]];
        arrF[ez][ephi][er][ivVrpphipzm[0]] = arrXdot[ez][ephi][er][ivVrpphipzm[0]];
        arrF[ez][ephi][er][ivVrpphipzm[1]] = arrXdot[ez][ephi][er][ivVrpphipzm[1]];
        arrF[ez][ephi][er][ivVrpphipzm[2]] = arrXdot[ez][ephi][er][ivVrpphipzm[2]];
        arrF[ez][ephi][er][ivVrpphipzp[0]] = arrXdot[ez][ephi][er][ivVrpphipzp[0]];
        arrF[ez][ephi][er][ivVrpphipzp[1]] = arrXdot[ez][ephi][er][ivVrpphipzp[1]];
        arrF[ez][ephi][er][ivVrpphipzp[2]] = arrXdot[ez][ephi][er][ivVrpphipzp[2]];

        /* f4(V,EP,tau,B,ni) = dB/dt */
        arrF[ez][ephi][er][ivBrm] = arrXdot[ez][ephi][er][ivBrm]; /* Left face */

        arrF[ez][ephi][er][ivBphim] = arrXdot[ez][ephi][er][ivBphim]; /* Down face */

        arrF[ez][ephi][er][ivBzm] = arrXdot[ez][ephi][er][ivBzm]; /* Back face */

        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivBrp] = arrXdot[ez][ephi][er][ivBrp]; /* Right face */
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivBphip] = arrXdot[ez][ephi][er][ivBphip]; /* Up face */
        }

        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivBzp] = arrXdot[ez][ephi][er][ivBzp]; /* Front face */
        }

        /* Set boundary conditions for EP field */
        /* f2(V,EP,tau,B,ni) = (EP - EPboundarycondition) */
        if (er == 0 || ez == 0 || (ephi == 0 && !(user -> phibtype))) {
          arrF[ez][ephi][er][ivVrmphimzm[3]] = arrX[ez][ephi][er][ivVrmphimzm[3]];
        }
        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivVrpphimzm[3]] = arrX[ez][ephi][er][ivVrpphimzm[3]];
        }
        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivVrmphimzp[3]] = arrX[ez][ephi][er][ivVrmphimzp[3]];
        }
        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivVrmphipzm[3]] = arrX[ez][ephi][er][ivVrmphipzm[3]];
        }
        if (ez == N[2] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivVrmphipzp[3]] = arrX[ez][ephi][er][ivVrmphipzp[3]];
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrF[ez][ephi][er][ivVrpphimzp[3]] = arrX[ez][ephi][er][ivVrpphimzp[3]];
        }
        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivVrpphipzm[3]] = arrX[ez][ephi][er][ivVrpphipzm[3]];
        }
        if (er == N[0] - 1 && ez == N[2] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivVrpphipzp[3]] = arrX[ez][ephi][er][ivVrpphipzp[3]];
        }
      }
    }
  }
  /* End of triple for loop */

  PetscCall(DMStagVecRestoreArray(da, pLocal, & arrP));
  PetscCall(DMRestoreLocalVector(da, & pLocal));
  //PetscBarrier((PetscObject) potential);
  PetscCall(VecDestroy( & potential));

  if (user -> phibtype) {
    PetscCall(DMStagVecRestoreArray(da, fLocal, & arrF));
    PetscCall(DMLocalToGlobal(da, fLocal, INSERT_VALUES, F));
    /*DMRestoreLocalVector(da,&fLocal);*/

    /*DMGetLocalVector(da,&fLocal);*/
    PetscCall(DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal));
    PetscCall(DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal));
    PetscCall(DMStagVecGetArray(da, fLocal, & arrF));
  }

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        MFD_CellEdgeLengths(arrCoorda, er, ephi, ez, N, user -> dphi,
          icrmphimzm, icrpphimzm, icrmphipzm, icrpphipzm, icrmphimzp, icrpphimzp, icrmphipzp, icrpphipzp,
          & rmzmedgelength, & rmphimedgelength, & rmzpedgelength, & rmphipedgelength, & phimzmedgelength, & phimzpedgelength,
          & rpphimedgelength, & rpzmedgelength, & phipzmedgelength, & rpzpedgelength, & rpphipedgelength, & phipzpedgelength);
#line 9819

        /* f3(V,EP,tau,B,ni) = tau - derived_mimetic_curl2(B) = tau - primary_mimetic_curl^T (beta_f B)/beta_e2; for all inner edges */
        if (!(er == 0 || ez == 0)) {
          arrF[ez][ephi][er][ivErmzm] = arrX[ez][ephi][er][ivErmzm] - ((arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
                arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
                arrX[ez - 1][ephi][er][ivBrm] * betaf(er, ephi, ez - 1, LEFT, user) / surface(er, ephi, ez - 1, LEFT, user) +
                arrX[ez][ephi][er - 1][ivBzm] * betaf(er - 1, ephi, ez, BACK, user) / surface(er - 1, ephi, ez, BACK, user)) * rmzmedgelength / beta_rz(er, ephi, ez, BACK_LEFT, user)) ;
        }
        if (!(user -> phibtype)) {
          if (!(ephi == 0 || ez == 0)) {
            arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm] - ((-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                  arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
                  arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
                  arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / beta_other(er, ephi, ez, BACK_DOWN, user)) ;
          }
          if (!(er == 0 || ephi == 0)) {
            arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim] - ((-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
                  arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                  arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
                  arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / beta_other(er, ephi, ez, DOWN_LEFT, user)) ;
          }
        } else {
          if (!(ez == 0)) {
            arrF[ez][ephi][er][ivEphimzm] = arrX[ez][ephi][er][ivEphimzm] - ((-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                  arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
                  arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
                  arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / beta_other(er, ephi, ez, BACK_DOWN, user)) ;
          }
          if (!(er == 0)) {
            arrF[ez][ephi][er][ivErmphim] = arrX[ez][ephi][er][ivErmphim] - ((-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
                  arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
                  arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
                  arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / beta_other(er, ephi, ez, DOWN_LEFT, user)) ;
          }
        }
      }
    }
  }
  //PetscBarrier((PetscObject) F);

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        //F3(V,EP,tau,B,ni) -= prim_grad(EP) ; on inner edges
        if (!(er == 0 || ez == 0)) {
          arrF[ez][ephi][er][ivErmzm] -= arrGradEP[ez][ephi][er][ivErmzm];
        }
        if (!(user -> phibtype)) {
          if (!(ephi == 0 || ez == 0)) {
            arrF[ez][ephi][er][ivEphimzm] -= arrGradEP[ez][ephi][er][ivEphimzm];
          }
          if (!(er == 0 || ephi == 0)) {
            arrF[ez][ephi][er][ivErmphim] -= arrGradEP[ez][ephi][er][ivErmphim];
          }
        } else {
          if (!(ez == 0)) {
            arrF[ez][ephi][er][ivEphimzm] -= arrGradEP[ez][ephi][er][ivEphimzm];
          }
          if (!(er == 0)) {
            arrF[ez][ephi][er][ivErmphim] -= arrGradEP[ez][ephi][er][ivErmphim];
          }
        }
      }
    }
  }

  PetscCall(DMStagVecRestoreArray(da, fLocal, & arrF));
  PetscCall(DMLocalToGlobal(da, fLocal, INSERT_VALUES, F));

  PetscCall(VecDuplicate(F, & Fcopy));
  PetscCall(VecCopy(F, Fcopy));
  //VecScale(Fcopy, -1.0);

  PetscCall(DMGetLocalVector(da, & FcopyLocal));
  PetscCall(DMGlobalToLocalBegin(da, Fcopy, INSERT_VALUES, FcopyLocal));
  PetscCall(DMGlobalToLocalEnd(da, Fcopy, INSERT_VALUES, FcopyLocal));
  PetscCall(DMStagVecGetArray(da, FcopyLocal, & arrFcopy));

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        //F3copy(V,EP,tau,B,ni) -= tau ; on inner edges
        if (!(er == 0 || ez == 0)) {
          arrFcopy[ez][ephi][er][ivErmzm] -= arrX[ez][ephi][er][ivErmzm];
        }
        if (!(user -> phibtype)) {
          if (!(ephi == 0 || ez == 0)) {
            arrFcopy[ez][ephi][er][ivEphimzm] -= arrX[ez][ephi][er][ivEphimzm];
          }
          if (!(er == 0 || ephi == 0)) {
            arrFcopy[ez][ephi][er][ivErmphim] -= arrX[ez][ephi][er][ivErmphim];
          }
        } else {
          if (!(ez == 0)) {
            arrFcopy[ez][ephi][er][ivEphimzm] -= arrX[ez][ephi][er][ivEphimzm];
          }
          if (!(er == 0)) {
            arrFcopy[ez][ephi][er][ivErmphim] -= arrX[ez][ephi][er][ivErmphim];
          }
        }
      }
    }
  }
  PetscCall(DMStagVecRestoreArray(da, FcopyLocal, & arrFcopy));
  PetscCall(DMLocalToGlobal(da, FcopyLocal, INSERT_VALUES, Fcopy));
  PetscCall(DMRestoreLocalVector(da, & FcopyLocal));

  ApplyDerivedDivergence(ts, Fcopy, F, user);

  /* Restore vectors */
  PetscCall(DMStagVecRestoreArrayRead(da, GradEPLocal, & arrGradEP));
  PetscCall(DMRestoreLocalVector(da, & GradEPLocal));
  PetscCall(VecDestroy( & GradEP));

  PetscCall(DMRestoreLocalVector(da, & fLocal));
  PetscCall(DMStagVecRestoreArrayRead(da, xLocal, & arrX));
  PetscCall(DMRestoreLocalVector(da, & xLocal));
  PetscCall(DMStagVecRestoreArrayRead(da, xdotLocal, & arrXdot));
  PetscCall(DMRestoreLocalVector(da, & xdotLocal));
  PetscCall(DMStagVecRestoreArrayRead(da, bcLocal, & arrx));
  PetscCall(DMRestoreLocalVector(da, & bcLocal));

  PetscCall(DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord));
  PetscCall(DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda));
  if (user -> debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F = \n"));
    PetscCall(VecView(F, PETSC_VIEWER_STDOUT_WORLD));
  }

  //PetscBarrier((PetscObject) x);
  PetscCall(VecDestroy( & x));
  PetscCall(VecDestroy( & Fcopy));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));


  PetscFunctionReturn(PETSC_SUCCESS);
}
PetscErrorCode FormIFunction_InitializeEP(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void * ptr) {
  return initialize_ep(ts, t, X, Xdot, F, ptr, betae2, betae2, "FormIFunction_InitializeEP");
}
#line 1510

PetscErrorCode FormIFunction_InitializeEP_halo(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void * ptr) {
  return initialize_ep(ts, t, X, Xdot, F, ptr, betaephi_isolcell, betaeperp2, "FormIFunction_InitializeEP_halo");
}
#line 1938

#line 11303

#line 12103

/* Initial-condition relaxation residual (ictype 9/15, see FormInitialSolution_psi):
 * f1 has inertia n_i dV/dt and no advection; Ohm's law is IDEAL (no resistive
 * curl_2(B) term -- owner decision, kept for consistency with previous runs);
 * no runaway-current source. */
PetscErrorCode FormIFunction_newequilibrium_Vperp(TS ts, PetscReal t, Vec X, Vec Xdot, Vec F, void * ptr) {
  static const MFD_ResidualTerms terms = {
    .f1_inertia = PETSC_TRUE, .f3_resistive = PETSC_FALSE, .f3_jre = PETSC_FALSE,
    .label = "FormIFunction_newequilibrium_Vperp"};
  return vperp_residual(ts, t, X, Xdot, F, ptr, &terms);
}

#line 13805

PetscErrorCode FormRHSFunction_BImplicit(TS ts, PetscReal t, Vec X, Vec F, void * ptr) {
  PetscFunctionBeginUser;

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("FormRHSFunction_BImplicit",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec fLocal, xLocal;

#line 13823
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];

  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  PetscInt ivn;

  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[3], ivVrmphimzp[3], ivVrmphipzm[3], ivVrmphipzp[3];
  PetscInt ivVrpphimzm[3], ivVrpphipzm[3], ivVrpphipzp[3], ivVrpphimzp[3];

  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;

  PetscCall(VecZeroEntries(F));
  PetscCall(TSGetDM(ts, & da));

  PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
  PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));

  for (d = 0; d < 3; ++d) {
    /* Vertex locations */
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]));
  }
  /* Edge locations */
  PetscCall(DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm));
  PetscCall(DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim));
  PetscCall(DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim));
  PetscCall(DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip));
  PetscCall(DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp));
  /* Face locations */
  PetscCall(DMStagGetLocationSlot(da, LEFT, 0, & ivBrm));
  PetscCall(DMStagGetLocationSlot(da, DOWN, 0, & ivBphim));
  PetscCall(DMStagGetLocationSlot(da, BACK, 0, & ivBzm));
  PetscCall(DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp));
  PetscCall(DMStagGetLocationSlot(da, UP, 0, & ivBphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT, 0, & ivBzp));
  /* Cell locations */
  PetscCall(DMStagGetLocationSlot(da, ELEMENT, 0, & ivn));

  PetscCall(DMGetCoordinateDM(da, & dmCoord));
  PetscCall(DMGetCoordinatesLocal(da, & coordLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord));
  PetscCall(DMGetCoordinateDM(coordDA, & dmCoorda));
  PetscCall(DMGetCoordinatesLocal(coordDA, & coordaLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda));
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]));
    /* Face coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]));
    /* Edge coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]));
    /* Vertex coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]));
  }


#line 13961

  /*Compute function over the locally owned part of the grid */
  /* f1(V,E,B,ni) = 0 ; on all vertices
     f2(V,E,B,ni) = 0; on all edges
     f3(V,E,B,ni) = 0; on all faces
     f4(V,E,B,ni) = - primary_mimetic_div(P_cf(ni)R_vf(V)); in plasma cells
     f4(V,E,B,ni) = 0; in other cells */
  PetscCall(DMGetLocalVector(da, & fLocal));
  PetscCall(DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal));
  PetscCall(DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal));
  PetscCall(DMStagVecGetArray(da, fLocal, & arrF));

  PetscCall(DMGetLocalVector(da, & xLocal));
  PetscCall(DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMStagVecGetArrayRead(da, xLocal, & arrX));

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

#line 13988

        MFD_CellEdgeLengths(arrCoorda, er, ephi, ez, N, user -> dphi,
          icrmphimzm, icrpphimzm, icrmphipzm, icrpphipzm, icrmphimzp, icrpphimzp, icrmphipzp, icrpphipzp,
          & rmzmedgelength, & rmphimedgelength, & rmzpedgelength, & rmphipedgelength, & phimzmedgelength, & phimzpedgelength,
          & rpphimedgelength, & rpzmedgelength, & phipzmedgelength, & rpzpedgelength, & rpphipedgelength, & phipzpedgelength);
#line 14028

        /* f1(V,E,B,ni) = 0; on all vertices

           f4(V,E,B,ni) = - prim_div(P_{c->f}(ni)R_{v->f}(V)); in plasma cells
           f4(V,E,B,ni) = 0; in wall cells
           */
        arrF[ez][ephi][er][ivn] = 0.0;
#line 14038

        arrF[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
        arrF[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
        arrF[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

        arrF[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
        arrF[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
        arrF[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

        arrF[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
        arrF[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
        arrF[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

        arrF[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
        arrF[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
        arrF[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

        arrF[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
        arrF[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
        arrF[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

        arrF[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
        arrF[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
        arrF[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

        arrF[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
        arrF[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
        arrF[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

        arrF[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
        arrF[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
        arrF[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

        if (user -> ictype == 12 && er!=0 && ez!=0) {
          /* curlB x B = - (1/r) * (cos(2*r_max/r)*(exp(- z / z_max)/z_max + 2*r_max*sin(2*r_max/r)/r) + (r-z)*(2*r-z)) e_r +
             (- cos(2*r_max/r) + (2r-z) * exp(-z/z_max)/r^2) e_phi +
             (-z+r + exp(-z/z_max)/r * (exp(-z/z_max)/(r*z_max) + sin(2*r_max/r)* 2r_max/r^2) )e_z */
          if (er > 0 && ez > 0 && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7 && fabs(user -> dataC[er + ephi * N[0] + (ez - 1) * N[1] * N[0]] - 1.5) < 0.7 && fabs(user -> dataC[er - 1 + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7 && fabs(user -> dataC[er - 1 + ephi * N[0] + (ez - 1) * N[1] * N[0]] - 1.5) < 0.7) {
            arrF[ez][ephi][er][ivVrmphimzm[0]] = -(1.0/ arrCoord[ez][ephi][er][icrmphimzm[0]]) * (PetscCosReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icrmphimzm[0]]) * (PetscExpReal(- arrCoord[ez][ephi][er][icrmphimzm[2]] / user->zmax) / user->zmax + PetscSinReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icrmphimzm[0]]) * 2.0 * user->rmax / arrCoord[ez][ephi][er][icrmphimzm[0]]) + (arrCoord[ez][ephi][er][icrmphimzm[0]] - arrCoord[ez][ephi][er][icrmphimzm[2]]) * (2.0 * arrCoord[ez][ephi][er][icrmphimzm[0]] - arrCoord[ez][ephi][er][icrmphimzm[2]]) ) ;
            arrF[ez][ephi][er][ivVrmphimzm[1]] = - PetscCosReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icrmphimzm[0]]) + PetscExpReal(- arrCoord[ez][ephi][er][icrmphimzm[2]] / user->zmax) * (2.0 * arrCoord[ez][ephi][er][icrmphimzm[0]] - arrCoord[ez][ephi][er][icrmphimzm[2]]) / PetscSqr(arrCoord[ez][ephi][er][icrmphimzm[0]]);
            arrF[ez][ephi][er][ivVrmphimzm[2]] = arrCoord[ez][ephi][er][icrmphimzm[0]] - arrCoord[ez][ephi][er][icrmphimzm[2]] + ( PetscExpReal(- arrCoord[ez][ephi][er][icrmphimzm[2]] / user->zmax) / (user->zmax * arrCoord[ez][ephi][er][icrmphimzm[0]]) + PetscSinReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icrmphimzm[0]]) * 2.0 * user->rmax / PetscSqr(arrCoord[ez][ephi][er][icrmphimzm[0]]) ) * PetscExpReal(- arrCoord[ez][ephi][er][icrmphimzm[2]] / user->zmax) / arrCoord[ez][ephi][er][icrmphimzm[0]];
          }
        }

        /* f3(B,E) = 0 */
        arrF[ez][ephi][er][ivBrm] = 0.0;
        arrF[ez][ephi][er][ivBphim] = 0.0;
        arrF[ez][ephi][er][ivBzm] = 0.0;
        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivBrp] = 0.0;
        }
        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivBphip] = 0.0;
        }
        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivBzp] = 0.0;
        }

        /* f2(B,E) = 0; for edges outside of the plasma region and the wall/plasma interface */
        if (!(er == 0 || ez == 0)) {
          arrF[ez][ephi][er][ivErmzm] = 0.0;
        }
        if (!(user -> phibtype)) {
          if (!(ephi == 0 || ez == 0)) {
            arrF[ez][ephi][er][ivEphimzm] = 0.0;
          }
          if (!(er == 0 || ephi == 0)) {
            arrF[ez][ephi][er][ivErmphim] = 0.0;
          }
        } else {
          if (!(ez == 0)) {
            arrF[ez][ephi][er][ivEphimzm] = 0.0;
          }
          if (!(er == 0)) {
            arrF[ez][ephi][er][ivErmphim] = 0.0;
          }
        }

      }
    }
  }
  /* End of triple for loop */

  /* Restore vectors */
#line 14136

  PetscCall(DMStagVecRestoreArray(da, fLocal, & arrF));
  PetscCall(DMLocalToGlobal(da, fLocal, INSERT_VALUES, F));
  PetscCall(DMRestoreLocalVector(da, & fLocal));
  PetscCall(DMStagVecRestoreArrayRead(da, xLocal, & arrX));
  PetscCall(DMRestoreLocalVector(da, & xLocal));

  PetscCall(DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord));
  PetscCall(DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda));
  if (user -> debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "F = \n"));
    PetscCall(VecView(F, PETSC_VIEWER_STDOUT_WORLD));
  }

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormInitialSolution(TS ts, Vec X, void * ptr) {
  PetscFunctionBeginUser;

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("FormInitialSolution",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  User * user = (User * ) ptr;
  DM da;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec xLocal;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];

  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];
  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];

  PetscInt ivn;

  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  DM dmCoord;
  PetscScalar ** ** arrCoord, ** ** arrX, time = 0.0;
  PetscInt countz = 0, countphi = 0, countr = 0;

  PetscCall(TSGetDM(ts, & da));
  PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
  PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
  for (d = 0; d < 4; ++d) {
    /* Vertex locations */
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]));
  }
  /* Edge locations */
  PetscCall(DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm));
  PetscCall(DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim));
  PetscCall(DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim));
  PetscCall(DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip));
  PetscCall(DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp));
  /* Face locations */
  PetscCall(DMStagGetLocationSlot(da, LEFT, 0, & ivBrm));
  PetscCall(DMStagGetLocationSlot(da, DOWN, 0, & ivBphim));
  PetscCall(DMStagGetLocationSlot(da, BACK, 0, & ivBzm));
  PetscCall(DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp));
  PetscCall(DMStagGetLocationSlot(da, UP, 0, & ivBphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT, 0, & ivBzp));
  /* Cell locations */
  PetscCall(DMStagGetLocationSlot(da, ELEMENT, 0, & ivn));

  PetscCall(DMGetCoordinateDM(da, & dmCoord));
  PetscCall(DMGetCoordinatesLocal(da, & coordLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord));
  for (d = 0; d < 3; ++d) {
    /* Face coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]));
    /* Edge coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]));
  }

#line 14270
  if (1) {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       Cancel all components except B
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Vec tau, EP, ni, Vi_perp;
    PetscCall(VecGetSubVector( X, user -> isEP, & EP));
    PetscCall(VecScale(EP, 0.0));
    PetscCall(VecRestoreSubVector( X, user -> isEP, & EP));
    PetscCall(VecGetSubVector( X, user -> istau, & tau));
    PetscCall(VecScale(tau, 0.0));
    PetscCall(VecRestoreSubVector( X, user -> istau, & tau));
    PetscCall(VecGetSubVector( X, user -> isni, & ni));
    PetscCall(VecScale(ni, 0.0));
    PetscCall(VecRestoreSubVector( X, user -> isni, & ni));
    PetscCall(VecGetSubVector( X, user -> isV, & Vi_perp));
    PetscCall(VecScale(Vi_perp, 0.0));
    PetscCall(VecRestoreSubVector( X, user -> isV, & Vi_perp));
  }
  /* Compute function over the locally owned part of the grid */
  PetscCall(DMGetLocalVector(da, & xLocal));
  PetscCall(DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMStagVecGetArray(da, xLocal, & arrX));
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        if (user -> ictype == 1) {
          /* B(r, phi, z, t) = (cos(pi*t)/r) * ((z-0.5)^2 + sin(phi)) e_r */
          arrX[ez][ephi][er][ivBrm] = PetscCosReal(PETSC_PI * time) * (PetscSinReal(arrCoord[ez][ephi][er][icBrm[1]]) + PetscSqr(arrCoord[ez][ephi][er][icBrm[2]] - 0.5)) / arrCoord[ez][ephi][er][icBrm[0]];
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = PetscCosReal(PETSC_PI * time) * (PetscSinReal(arrCoord[ez][ephi][er][icBrp[1]]) + PetscSqr(arrCoord[ez][ephi][er][icBrp[2]] - 0.5)) / arrCoord[ez][ephi][er][icBrp[0]];
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t) = (cos(pi*t)*(2*z-1)/r) e_phi - (cos(pi*t)*cos(phi)/r^2) e_z */
          arrX[ez][ephi][er][ivErmzm] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErmzm[2]] - 1.0) / arrCoord[ez][ephi][er][icErmzm[0]];
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = -PetscCosReal(arrCoord[ez][ephi][er][icErmphim[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErmphim[0]]);

          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErpzm[2]] - 1.0) / arrCoord[ez][ephi][er][icErpzm[0]];
            arrX[ez][ephi][er][ivErpphim] = -PetscCosReal(arrCoord[ez][ephi][er][icErpphim[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = -PetscCosReal(arrCoord[ez][ephi][er][icErmphip[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErmzp[2]] - 1.0) / arrCoord[ez][ephi][er][icErmzp[0]];
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = -PetscCosReal(arrCoord[ez][ephi][er][icErpphip[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErpzp[2]] - 1.0) / arrCoord[ez][ephi][er][icErpzp[0]];
          }
        } else if (user -> ictype == 2) {
          /* B(r, phi, z, t=0) = sin(pi*r) e_phi */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icBphim[0]]);
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icBphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t=0) = (pi*cos(pi*r) + sin(pi*r)/r) e_z */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphim[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphim[0]]) / arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphim[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphim[0]]) / arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphip[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphip[0]]) / arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphip[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphip[0]]) / arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        } else if (user -> ictype == 3) {
          /* B(r, phi, z, t=0) = cos(phi) e_r - sin(phi) e_phi */
          arrX[ez][ephi][er][ivBrm] = PetscCosReal(arrCoord[ez][ephi][er][icBrm[1]]);
          arrX[ez][ephi][er][ivBphim] = -PetscSinReal(arrCoord[ez][ephi][er][icBphim[1]]);
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = PetscCosReal(arrCoord[ez][ephi][er][icBrp[1]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = -PetscSinReal(arrCoord[ez][ephi][er][icBphip[1]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        } else if (user -> ictype == 4) {
          /* B(r, phi, z, t=0) = e_z */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = 1.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 1.0;
          }
          /* E(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        } else if (user -> ictype == 5) {
          /* B(r, phi, z, t=0) = (1/r) e_r */
          arrX[ez][ephi][er][ivBrm] = 1.0 / arrCoord[ez][ephi][er][icBrm[0]];
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 1.0 / arrCoord[ez][ephi][er][icBrp[0]];
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        } else if (user -> ictype == 6) {
          /* B(r, phi, z, t=0) = e_phi */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 1.0;
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 1.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t=0) = (1/r) e_z */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = 1.0 / arrCoord[ez][ephi][er][icErmphim[0]];
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 1.0 / arrCoord[ez][ephi][er][icErpphim[0]];
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 1.0 / arrCoord[ez][ephi][er][icErmphip[0]];
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 1.0 / arrCoord[ez][ephi][er][icErpphip[0]];
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        } else if (user -> ictype == 7) {
          /* B(r, phi, z, t=0) = (1 + sin(phi)) e_z */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = PetscSinReal(arrCoord[ez][ephi][er][icBzm[1]]) + 1.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = PetscSinReal(arrCoord[ez][ephi][er][icBzp[1]]) + 1.0;
          }
          /* E(r, phi, z, t=0) = (cos(phi)/r) e_r */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = PetscCosReal(arrCoord[ez][ephi][er][icEphimzm[1]]) / arrCoord[ez][ephi][er][icEphimzm[0]];
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = PetscCosReal(arrCoord[ez][ephi][er][icEphipzm[1]]) / arrCoord[ez][ephi][er][icEphipzm[0]];
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = PetscCosReal(arrCoord[ez][ephi][er][icEphimzp[1]]) / arrCoord[ez][ephi][er][icEphimzp[0]];
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = PetscCosReal(arrCoord[ez][ephi][er][icEphipzp[1]]) / arrCoord[ez][ephi][er][icEphipzp[0]];
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        } else if (user -> ictype == 8) {
          /* B(r, phi, z, t=0) = exp(sin(phi)*r) e_z */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icBzm[1]]) * arrCoord[ez][ephi][er][icBzm[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icBzp[1]]) * arrCoord[ez][ephi][er][icBzp[0]]);
          }
          /* E(r, phi, z, t=0) = exp(sin(phi)*r) * (cos(phi) e_r - sin(phi) e_phi) */
          arrX[ez][ephi][er][ivErmzm] = -PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icErmzm[1]]) * arrCoord[ez][ephi][er][icErmzm[0]]) * PetscSinReal(arrCoord[ez][ephi][er][icErmzm[1]]);
          arrX[ez][ephi][er][ivEphimzm] = PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icEphimzm[1]]) * arrCoord[ez][ephi][er][icEphimzm[0]]) * PetscCosReal(arrCoord[ez][ephi][er][icEphimzm[1]]);
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = -PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icErpzm[1]]) * arrCoord[ez][ephi][er][icErpzm[0]]) * PetscSinReal(arrCoord[ez][ephi][er][icErpzm[1]]);
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icEphipzm[1]]) * arrCoord[ez][ephi][er][icEphipzm[0]]) * PetscCosReal(arrCoord[ez][ephi][er][icEphipzm[1]]);
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = -PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icErmzp[1]]) * arrCoord[ez][ephi][er][icErmzp[0]]) * PetscSinReal(arrCoord[ez][ephi][er][icErmzp[1]]);
            arrX[ez][ephi][er][ivEphimzp] = PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icEphimzp[1]]) * arrCoord[ez][ephi][er][icEphimzp[0]]) * PetscCosReal(arrCoord[ez][ephi][er][icEphimzp[1]]);
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icEphipzp[1]]) * arrCoord[ez][ephi][er][icEphipzp[0]]) * PetscCosReal(arrCoord[ez][ephi][er][icEphipzp[1]]);
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = -PetscExpReal(PetscSinReal(arrCoord[ez][ephi][er][icErpzp[1]]) * arrCoord[ez][ephi][er][icErpzp[0]]) * PetscSinReal(arrCoord[ez][ephi][er][icErpzp[1]]);
          }
        } else if (user -> ictype == 10) {
          /* B(r, phi, z, t=0) = e_phi + sqrt(2*ln(2*r_max/r)) e_z */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 1.0;
          arrX[ez][ephi][er][ivBzm] = PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzm[0]]));
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 1.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzp[0]]));
          }
          /* tau(r, phi, z, t=0) = 0 */
          /*
             arrX[ez][ephi][er][ivErmzm] = 0.0;
             arrX[ez][ephi][er][ivEphimzm] = 0.0;
             arrX[ez][ephi][er][ivErmphim] = 0.0;
             if (er == N[0] - 1) {
             arrX[ez][ephi][er][ivErpzm] = 0.0;
             arrX[ez][ephi][er][ivErpphim] = 0.0;
             }
             if (ephi == N[1] - 1) {
             arrX[ez][ephi][er][ivEphipzm] = 0.0;
             arrX[ez][ephi][er][ivErmphip] = 0.0;
             }
             if (ez == N[2] - 1) {
             arrX[ez][ephi][er][ivErmzp] = 0.0;
             arrX[ez][ephi][er][ivEphimzp] = 0.0;
             }
             if (er == N[0] - 1 && ephi == N[1] - 1) {
             arrX[ez][ephi][er][ivErpphip] = 0.0;
             }
             if (ephi == N[1] - 1 && ez == N[2] - 1) {
             arrX[ez][ephi][er][ivEphipzp] = 0.0;
             }
             if (er == N[0] - 1 && ez == N[2] - 1) {
             arrX[ez][ephi][er][ivErpzp] = 0.0;
             }*/

          /* tau(r, phi, z, t=0) = eta / (mu0 * V_A) / r * [e_z + sqrt(2*ln(2*r_max/r))^{-1} e_phi] */
          arrX[ez][ephi][er][ivErmzm] = condu(er, ephi, ez, BACK_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzm[0]])));
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = condu(er, ephi, ez, DOWN_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = condu(er, ephi, ez, BACK_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzm[0]])));
            arrX[ez][ephi][er][ivErpphim] = condu(er, ephi, ez, DOWN_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = condu(er, ephi, ez, UP_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = condu(er, ephi, ez, FRONT_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzp[0]])));
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = condu(er, ephi, ez, UP_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = condu(er, ephi, ez, FRONT_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzp[0]])));
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        } else if (user -> ictype == 9 && user -> itime == 0.0) {
          if(1){
            /* B(r, phi, z, t=0) Loaded from ReadInitialData output */
            arrX[ez][ephi][er][ivBrm] = user -> datar[er + ephi * (N[0] + 1) + ez * N[1] * (N[0] + 1)] / user->B0;
            countr = countr + 1;
            //PetscPrintf(PETSC_COMM_SELF,"arrX[%d][%d][%d][%d] = %15.10e\n",ez,ephi,er,ivBrm,arrX[ez][ephi][er][ivBrm]);

            arrX[ez][ephi][er][ivBphim] = user -> dataphi[er + ephi * N[0] + ez * N[0] * (N[1] + 1)] / user->B0;
            countphi = countphi + 1;
            //PetscPrintf(PETSC_COMM_SELF,"arrX[%d][%d][%d][%d] = %15.10e\n",ez,ephi,er,ivBphim,arrX[ez][ephi][er][ivBphim]);

            arrX[ez][ephi][er][ivBzm] = user -> dataz[er + ephi * N[0] + ez * N[0] * N[1]] / user->B0;
            countz = countz + 1;
            //PetscPrintf(PETSC_COMM_SELF,"arrX[%d][%d][%d][%d] = %15.10e\n",ez,ephi,er,ivBzm,arrX[ez][ephi][er][ivBzm]);

            if (er == N[0] - 1) {
              arrX[ez][ephi][er][ivBrp] = user -> datar[er + 1 + ephi * (N[0] + 1) + ez * N[1] * (N[0] + 1)] / user->B0;
              countr = countr + 1;
              //PetscPrintf(PETSC_COMM_SELF,"arrX[%d][%d][%d][%d] = %15.10e\n",ez,ephi,er,ivBrp,arrX[ez][ephi][er][ivBrp]);
            }
            if (ephi == N[1] - 1) {
              arrX[ez][ephi][er][ivBphip] = user -> dataphi[er + (ephi + 1) * N[0] + ez * N[0] * (N[1] + 1)] / user->B0;
              countphi = countphi + 1;
              //PetscPrintf(PETSC_COMM_SELF,"arrX[%d][%d][%d][%d] = %15.10e\n",ez,ephi,er,ivBphip,arrX[ez][ephi][er][ivBphip]);
            }
            if (ez == N[2] - 1) {
              arrX[ez][ephi][er][ivBzp] = user -> dataz[er + ephi * N[0] + (ez + 1) * N[0] * N[1]] / user->B0;
              countz = countz + 1;
              //PetscPrintf(PETSC_COMM_SELF,"arrX[%d][%d][%d][%d] = %15.10e\n",ez,ephi,er,ivBzp,arrX[ez][ephi][er][ivBzp]);
            }

            if (countr > user -> numr) printf("countr = [%d] exceeds numr = [%d]\n", countr, user -> numr);
            if (countphi > user -> numphi) printf("countphi = [%d] exceeds numphi = [%d]\n", countphi, user -> numphi);
            if (countz > user -> numz) printf("countz = [%d] exceeds numz = [%d]\n", countz, user -> numz);
          }
          /* tau(r, phi, z, t=0) initialized with zero values. The tau field will be computed later from B field with the derived curl operator */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        } else if (user -> ictype == 11) {
          /* B(r, phi, z, t=0) = (2+cos(eta/mu0 * t)) * [e_phi + sqrt(2*ln(2*r_max/r)) e_z] */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 2.0 + PetscCosReal((user->eta / user->mu0) * time);
          arrX[ez][ephi][er][ivBzm] = (2.0 + PetscCosReal((user->eta / user->mu0) * time)) * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzm[0]]));
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 2.0 + PetscCosReal((user->eta / user->mu0) * time);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = (2.0 + PetscCosReal((user->eta / user->mu0) * time)) * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzp[0]]));
          }
          /* tau(r, phi, z, t=0) = eta / (mu0 * V_A) * (cos(eta/mu0 * t)+2) / r * [e_z + sqrt(2*ln(2*r_max/r))^{-1} e_phi] */
          arrX[ez][ephi][er][ivErmzm] = condu(er, ephi, ez, BACK_LEFT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, BACK_LEFT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzm[0]])));
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = condu(er, ephi, ez, DOWN_LEFT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, DOWN_LEFT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = condu(er, ephi, ez, BACK_RIGHT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, BACK_RIGHT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzm[0]])));
            arrX[ez][ephi][er][ivErpphim] = condu(er, ephi, ez, DOWN_RIGHT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, DOWN_RIGHT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = condu(er, ephi, ez, UP_LEFT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, UP_LEFT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = condu(er, ephi, ez, FRONT_LEFT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, FRONT_LEFT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzp[0]])));
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = condu(er, ephi, ez, UP_RIGHT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, UP_RIGHT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = condu(er, ephi, ez, FRONT_RIGHT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, FRONT_RIGHT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzp[0]])));
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        }
        else if (user -> ictype == 12) {
          /* B(r, phi, z, t=0) = exp(- z / z_max) / r e_r + (r - z) e_phi + cos(2*r_max/r) e_z */
          arrX[ez][ephi][er][ivBrm] = PetscExpReal(- arrCoord[ez][ephi][er][icBrm[2]] / user->zmax) / arrCoord[ez][ephi][er][icBrm[0]];
          arrX[ez][ephi][er][ivBphim] = arrCoord[ez][ephi][er][icBphim[0]] - arrCoord[ez][ephi][er][icBphim[2]];
          arrX[ez][ephi][er][ivBzm] = PetscCosReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzm[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = PetscExpReal(- arrCoord[ez][ephi][er][icBrp[2]] / user->zmax) / arrCoord[ez][ephi][er][icBrp[0]];
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = arrCoord[ez][ephi][er][icBphip[0]] - arrCoord[ez][ephi][er][icBphip[2]];
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = PetscCosReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzp[0]]);
          }
          /* tau(r, phi, z, t=0) = eta / (mu0 * V_A) / r * [e_z + sqrt(2*ln(2*r_max/r))^{-1} e_phi] */
          arrX[ez][ephi][er][ivErmzm] = condu(er, ephi, ez, BACK_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzm[0]])));
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = condu(er, ephi, ez, DOWN_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = condu(er, ephi, ez, BACK_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzm[0]])));
            arrX[ez][ephi][er][ivErpphim] = condu(er, ephi, ez, DOWN_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = condu(er, ephi, ez, UP_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = condu(er, ephi, ez, FRONT_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzp[0]])));
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = condu(er, ephi, ez, UP_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = condu(er, ephi, ez, FRONT_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzp[0]])));
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        }

        else if (user -> ictype == 13) {
          /* B(r, phi, z, t=0) = e_phi + sqrt(2*ln(2*r_max/r)) e_z */
          if(er!=0) {
            arrX[ez][ephi][er][ivBrm] = 0.0;
          }
          arrX[ez][ephi][er][ivBphim] = 1.0;
          if(ez!=0) {
            arrX[ez][ephi][er][ivBzm] = PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzm[0]]));
          }
          if(er==0) {
            arrX[ez][ephi][er][ivBrm] = NAN;
          }
          if(ez==0) {
            arrX[ez][ephi][er][ivBzm] = NAN;
          }
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = NAN;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 1.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = NAN;
          }
          /* tau(r, phi, z, t=0) = 0 */
          /*
             arrX[ez][ephi][er][ivErmzm] = 0.0;
             arrX[ez][ephi][er][ivEphimzm] = 0.0;
             arrX[ez][ephi][er][ivErmphim] = 0.0;
             if (er == N[0] - 1) {
             arrX[ez][ephi][er][ivErpzm] = 0.0;
             arrX[ez][ephi][er][ivErpphim] = 0.0;
             }
             if (ephi == N[1] - 1) {
             arrX[ez][ephi][er][ivEphipzm] = 0.0;
             arrX[ez][ephi][er][ivErmphip] = 0.0;
             }
             if (ez == N[2] - 1) {
             arrX[ez][ephi][er][ivErmzp] = 0.0;
             arrX[ez][ephi][er][ivEphimzp] = 0.0;
             }
             if (er == N[0] - 1 && ephi == N[1] - 1) {
             arrX[ez][ephi][er][ivErpphip] = 0.0;
             }
             if (ephi == N[1] - 1 && ez == N[2] - 1) {
             arrX[ez][ephi][er][ivEphipzp] = 0.0;
             }
             if (er == N[0] - 1 && ez == N[2] - 1) {
             arrX[ez][ephi][er][ivErpzp] = 0.0;
             }*/

          /* tau(r, phi, z, t=0) = eta / (mu0 * V_A) / r * [e_z + sqrt(2*ln(2*r_max/r))^{-1} e_phi] */
          arrX[ez][ephi][er][ivErmzm] = condu(er, ephi, ez, BACK_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzm[0]])));
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = condu(er, ephi, ez, DOWN_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = condu(er, ephi, ez, BACK_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzm[0]])));
            arrX[ez][ephi][er][ivErpphim] = condu(er, ephi, ez, DOWN_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = condu(er, ephi, ez, UP_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = condu(er, ephi, ez, FRONT_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzp[0]])));
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = condu(er, ephi, ez, UP_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = condu(er, ephi, ez, FRONT_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzp[0]])));
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        }
        else if (user -> ictype > 13) PetscPrintf(PETSC_COMM_WORLD, "Test case not set\n");

      }
    }
  }

  /* Restore vectors */
  PetscCall(DMStagVecRestoreArray(da, xLocal, & arrX));
  if (user -> itime == 0.0) {
    PetscCall(DMLocalToGlobal(da, xLocal, INSERT_VALUES, X));
    PetscCall(DMRestoreLocalVector(da, & xLocal));
  } else {
    char filename[PETSC_MAX_PATH_LEN];
    PetscViewer viewerX;

    /* Read X in binary file */
    PetscCall(PetscSNPrintf(filename, sizeof(filename), "%s/X_ic%.2d_grid%.2dx%.2dx%.2d_step%.3d_time%5.7f.dat", user->input_folder, user -> ictype, user -> Nr, user -> Nphi, user -> Nz, (int)(user -> oldstep), (double) user -> itime));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Reading X vector from file %s ...\n", filename));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, & viewerX));
    PetscCall(VecLoad(X, viewerX));
    /* Destroy the viewer */
    PetscCall(PetscViewerDestroy( & viewerX));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Reading from file %s is over.\n", filename));
  }

#line 15192

  if ((user -> ictype == 9 || user -> ictype == 15) && user -> itime == 0.0) {
    SNES dummysnes;
    TS dummyts;
    KSP dummyKSP;
    PC dummypc;
    Mat J;
    PetscInt steps;
    PetscReal ftime;
    TSConvergedReason reason;

    PetscCall(TSCreate(PETSC_COMM_WORLD, & dummyts));
    PetscCall(TSSetDM(dummyts, user->coorda));

    PetscCall(TSSetProblemType(dummyts, TS_LINEAR));

    PetscCall(TSGetSNES(dummyts, & dummysnes));
    //SNESSetType(dummysnes, SNESKSPONLY);
    PetscCall(TSSetType(dummyts, TSBEULER));
    PetscCall(DMCreateMatrix(da, & J));
    /* Use coloring to compute finite difference J efficiently */
    PetscCall(SNESSetJacobian(dummysnes, J, J, SNESComputeJacobianDefaultColor, PETSC_NULLPTR));
    PetscCall(TSSetIFunction(dummyts, NULL, FormIFunction_InitializeEP, user));

    PetscCall(SNESSetUseMatrixFree(dummysnes,PETSC_TRUE,PETSC_FALSE));
    PetscCall(SNESSetOptionsPrefix(dummysnes, "dummySNES_"));
    PetscCall(SNESSetFromOptions(dummysnes));

    PetscCall(TSSetTime(dummyts, 0.0));
    PetscCall(TSSetMaxTime(dummyts, 1e-2));
    PetscCall(TSSetExactFinalTime(dummyts, TS_EXACTFINALTIME_STEPOVER));
    PetscCall(TSSetTimeStep(dummyts, 1e-2));
    PetscCall(TSSetSolution(dummyts, X));

    //TSGetSNES(dummyts, & dummysnes);
    //SNESGetKSP(dummysnes, & dummyKSP);
    TSGetKSP(dummyts, & dummyKSP);
    PetscCall(KSPGetPC(dummyKSP, & dummypc));
    //PCSetType(dummypc, PCNONE);

    PetscCall(PCSetType(dummypc,PCBJACOBI));
    //PCASMSetOverlap(dummypc,7);
    //PCFactorSetMatSolverType(dummypc,MATSOLVERMUMPS);
    //PCFactorSetUpMatSolverType(dummypc);


    PetscCall(KSPSetOptionsPrefix(dummyKSP, "dummyKSP_"));

    /* Logic below modifies the PC directly, so this is the last chance to change the solver from the command line */
    PetscCall(KSPSetFromOptions(dummyKSP));

    {/* first level -> split ni from the rest : {ni}, {V Phi tau B}
        second level -> split tau from {V Phi B} : {ni}, {{tau},{V Phi B}}
        third level -> split V from {Phi B} : {ni}, {{tau},{{Phi B}, {V}}}
        fourth level -> split Phi from B : {ni}, {{tau},{{{Phi}, {B}}, {V}}}
        */
      IS            is[2];
      DMStagStencil stencil0[1], stencil1[10];
      PC            pc_notc, pc_noe;

      const char *name[2] = {"ni", "TEBV"};

      //PetscCall(KSPGetPC(dummyKSP,&dummypc));
      //PetscCall(PCSetType(dummypc,PCFIELDSPLIT));

      // First split is cells
      stencil0[0].loc = DMSTAG_ELEMENT;
      stencil0[0].c = 0;

      // Second split is the rest
      for (PetscInt c=0; c<4; ++c) {
        stencil1[c].loc = DMSTAG_BACK_DOWN_LEFT;
        stencil1[c].c = c;
      }
      stencil1[4].loc = DMSTAG_LEFT;
      stencil1[4].c = 0;
      stencil1[5].loc = DMSTAG_BACK;
      stencil1[5].c = 0;
      stencil1[6].loc = DMSTAG_DOWN;
      stencil1[6].c = 0;
      stencil1[7].loc = DMSTAG_BACK_DOWN;
      stencil1[7].c = 0;
      stencil1[8].loc = DMSTAG_BACK_LEFT;
      stencil1[8].c = 0;
      stencil1[9].loc = DMSTAG_DOWN_LEFT;
      stencil1[9].c = 0;

      PetscCall(DMStagCreateISFromStencils(da,1,stencil0,&is[0]));
      PetscCall(DMStagCreateISFromStencils(da,10,stencil1,&is[1]));

      for (PetscInt i=0; i<2; ++i) {
        PetscCall(PCFieldSplitSetIS(dummypc,name[i],is[i]));
      }

      for (PetscInt i=0; i<2; ++i) {
        PetscCall(ISDestroy(&is[i]));
      }

      /* If the fieldsplit PC wasn't overridden, further split the second split */
      {
        PCType pc_type;
        PetscBool is_fieldsplit;

        PetscCall(KSPGetPC(dummyKSP, &dummypc));
        PetscCall(PCGetType(dummypc,&pc_type));
        PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
        if (is_fieldsplit) {
          DM            dm_notc;
          KSP           *sub_ksp;

          PetscInt      n_splits;
          DMStagStencil stencil_notc_edges[3], stencil_notc_notedges[7];
          IS            is_notc[2];
          const char    *name_notc[2] = {"tau","EBV"};

          PetscCall(PCSetUp(dummypc)); // Set up the Fieldsplit PC
          PetscCall(PCFieldSplitGetSubKSP(dummypc,&n_splits,&sub_ksp));
          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
          PetscCall(KSPGetPC(sub_ksp[1],&pc_notc));
          PetscCall(PetscFree(sub_ksp));

          PetscCall(DMStagCreateCompatibleDMStag(da,4,1,1,0,&dm_notc));

          // First split within notc is edges
          stencil_notc_edges[0].loc = DMSTAG_BACK_DOWN;
          stencil_notc_edges[0].c = 0;
          stencil_notc_edges[1].loc = DMSTAG_BACK_LEFT;
          stencil_notc_edges[1].c = 0;
          stencil_notc_edges[2].loc = DMSTAG_DOWN_LEFT;
          stencil_notc_edges[2].c = 0;

          // Second split within notc is faces and vertices
          for (PetscInt c=0; c<3; ++c) {
            stencil_notc_notedges[c].loc = DMSTAG_BACK_DOWN_LEFT;
            stencil_notc_notedges[c].c = c;
          }
          stencil_notc_notedges[3].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_notc_notedges[3].c = 3;
          stencil_notc_notedges[4].loc = DMSTAG_LEFT;
          stencil_notc_notedges[4].c = 0;
          stencil_notc_notedges[5].loc = DMSTAG_BACK;
          stencil_notc_notedges[5].c = 0;
          stencil_notc_notedges[6].loc = DMSTAG_DOWN;
          stencil_notc_notedges[6].c = 0;

          PetscCall(DMStagCreateISFromStencils(dm_notc,3,stencil_notc_edges,&is_notc[0]));
          PetscCall(DMStagCreateISFromStencils(dm_notc,7,stencil_notc_notedges,&is_notc[1]));

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(PCFieldSplitSetIS(pc_notc,name_notc[i],is_notc[i]));
          }

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(ISDestroy(&is_notc[i]));
          }
          PetscCall(DMDestroy(&dm_notc));
        }
      }

      /* If the fieldsplit PC wasn't overridden, further split the second split of the second level */
      {
        PCType pc_type;
        PetscBool is_fieldsplit;

        PetscCall(PCGetType(pc_notc,&pc_type));
        PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
        if (is_fieldsplit) {
          DM            dm_noe;
          KSP           *sub_ksp;

          PetscInt      n_splits;
          DMStagStencil stencil_noe_velocity[3], stencil_noe_notvelocity[4];
          IS            is_noe[2];
          const char    *name_noe[2] = {"EB", "V"};

          PetscCall(PCSetUp(pc_notc)); // Set up the Fieldsplit PC
          PetscCall(PCFieldSplitGetSubKSP(pc_notc,&n_splits,&sub_ksp));
          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
          PetscCall(KSPGetPC(sub_ksp[1],&pc_noe));
          PetscCall(PetscFree(sub_ksp));

          PetscCall(DMStagCreateCompatibleDMStag(da,4,0,1,0,&dm_noe));

          // First split within notv is faces and 4th dofs on vertices
          stencil_noe_notvelocity[0].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_noe_notvelocity[0].c = 3;

          stencil_noe_notvelocity[1].loc = DMSTAG_LEFT;
          stencil_noe_notvelocity[1].c = 0;
          stencil_noe_notvelocity[2].loc = DMSTAG_BACK;
          stencil_noe_notvelocity[2].c = 0;
          stencil_noe_notvelocity[3].loc = DMSTAG_DOWN;
          stencil_noe_notvelocity[3].c = 0;

          // Second split within notv is the first 3 dofs on vertices
          for (PetscInt c=0; c<3; ++c) {
            stencil_noe_velocity[c].loc = DMSTAG_BACK_DOWN_LEFT;
            stencil_noe_velocity[c].c = c;
          }

          PetscCall(DMStagCreateISFromStencils(dm_noe,4,stencil_noe_notvelocity,&is_noe[0]));
          PetscCall(DMStagCreateISFromStencils(dm_noe,3,stencil_noe_velocity,&is_noe[1]));

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(PCFieldSplitSetIS(pc_noe,name_noe[i],is_noe[i]));
          }

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(ISDestroy(&is_noe[i]));
          }
          PetscCall(DMDestroy(&dm_noe));
        }
      }

      /* If the fieldsplit PC wasn't overridden, further split the first split of the third level */
      {
        PCType pc_type;
        PetscBool is_fieldsplit;

        PetscCall(PCGetType(pc_noe,&pc_type));
        PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
        if (is_fieldsplit) {
          DM            dm_notv;
          KSP           *sub_ksp;
          PC            pc_noe_2;
          PetscInt      n_splits;
          DMStagStencil stencil_notv_faces[3], stencil_notv_notfaces[3];
          IS            is_notv[2];
          const char    *name_notv[2] = {"EP", "B"};

          PetscCall(PCSetUp(pc_noe)); // Set up the Fieldsplit PC
          PetscCall(PCFieldSplitGetSubKSP(pc_noe,&n_splits,&sub_ksp));
          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
          PetscCall(KSPGetPC(sub_ksp[0],&pc_noe_2));
          PetscCall(PetscFree(sub_ksp));

          PetscCall(DMStagCreateCompatibleDMStag(da,1,0,1,0,&dm_notv));

          // First split within notv is the 4th dofs on vertices
          stencil_notv_notfaces[0].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_notv_notfaces[0].c = 0;

          // Second split within notv is faces
          stencil_notv_faces[0].loc = DMSTAG_LEFT;
          stencil_notv_faces[0].c = 0;
          stencil_notv_faces[1].loc = DMSTAG_BACK;
          stencil_notv_faces[1].c = 0;
          stencil_notv_faces[2].loc = DMSTAG_DOWN;
          stencil_notv_faces[2].c = 0;

          PetscCall(DMStagCreateISFromStencils(dm_notv,1,stencil_notv_notfaces,&is_notv[0]));
          PetscCall(DMStagCreateISFromStencils(dm_notv,3,stencil_notv_faces,&is_notv[1]));

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(PCFieldSplitSetIS(pc_noe_2,name_notv[i],is_notv[i]));
          }

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(ISDestroy(&is_notv[i]));
          }
          PetscCall(DMDestroy(&dm_notv));
        }
      }
    }

    KSPSetTolerances(dummyKSP,1e-8,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

    PetscCall(TSSolve(dummyts, X));

    // if(0){
    // KSPConvergedReasonView(dummyKSP, PETSC_VIEWER_DEFAULT);
    // KSPMonitorSet(dummyKSP, (PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*))KSPMonitorResidual, NULL, NULL);
    // KSPView(dummyKSP, PETSC_VIEWER_STDOUT_WORLD);
    // }

    PetscCall(TSGetSolveTime(dummyts, & ftime));
    PetscCall(TSGetStepNumber(dummyts, & steps));
    PetscCall(TSGetConvergedReason(dummyts, & reason));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Setting Initial conditions: %s at time %g after %d steps\n", TSConvergedReasons[reason], (double) ftime, steps));
    //DumpSolution(ts, 50, X, user);
    PetscCall(MatDestroy( & J));
    PetscCall(TSDestroy( & dummyts));
  }

  PetscCall(DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord));
  if (user -> debug) {
    /*This print is just for debugging*/
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Initial solution vector\n"));
    PetscCall(VecView(X, PETSC_VIEWER_STDOUT_WORLD));
  }

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormExactSolution(PetscReal time, TS ts, Vec * X, void * ptr) {
  PetscFunctionBeginUser;

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("FormExactSolution",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  User * user = (User * ) ptr;
  DM da;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec xLocal;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];

  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];
  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];

  PetscInt ivn;

  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  DM dmCoord;
  PetscScalar ** ** arrCoord, ** ** arrX, dt;
  PetscInt countz = 0, countphi = 0, countr = 0;

  PetscCall(TSGetDM(ts, & da));
  PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
  PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));

  for (d = 0; d < 4; ++d) {
    /* Vertex locations */
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]));
  }
  /* Edge locations */
  PetscCall(DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm));
  PetscCall(DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim));
  PetscCall(DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim));
  PetscCall(DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip));
  PetscCall(DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp));
  /* Face locations */
  PetscCall(DMStagGetLocationSlot(da, LEFT, 0, & ivBrm));
  PetscCall(DMStagGetLocationSlot(da, DOWN, 0, & ivBphim));
  PetscCall(DMStagGetLocationSlot(da, BACK, 0, & ivBzm));
  PetscCall(DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp));
  PetscCall(DMStagGetLocationSlot(da, UP, 0, & ivBphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT, 0, & ivBzp));
  /* Cell locations */
  PetscCall(DMStagGetLocationSlot(da, ELEMENT, 0, & ivn));

  PetscCall(DMGetCoordinateDM(da, & dmCoord));
  PetscCall(DMGetCoordinatesLocal(da, & coordLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord));
  for (d = 0; d < 3; ++d) {
    /* Face coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]));
    /* Edge coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]));
  }
  /* Compute function over the locally owned part of the grid */
  PetscCall(DMGetLocalVector(da, & xLocal));
  PetscCall(DMStagVecGetArray(da, xLocal, & arrX));
  PetscCall(TSGetTimeStep(ts, & dt));
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        if (user -> ictype == 1) {
          /* B(r, phi, z, t) = (cos(pi*t)/r) * ((z-0.5)^2 + sin(phi)) e_r */
          arrX[ez][ephi][er][ivBrm] = PetscCosReal(PETSC_PI * time) * (PetscSinReal(arrCoord[ez][ephi][er][icBrm[1]]) + PetscSqr(arrCoord[ez][ephi][er][icBrm[2]] - 0.5)) / arrCoord[ez][ephi][er][icBrm[0]];
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = PetscCosReal(PETSC_PI * time) * (PetscSinReal(arrCoord[ez][ephi][er][icBrp[1]]) + PetscSqr(arrCoord[ez][ephi][er][icBrp[2]] - 0.5)) / arrCoord[ez][ephi][er][icBrp[0]];
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t) = (cos(pi*t)*(2*z-1)/r) e_phi - (cos(pi*t)*cos(phi)/r^2) e_z */
          arrX[ez][ephi][er][ivErmzm] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErmzm[2]] - 1.0) / arrCoord[ez][ephi][er][icErmzm[0]];
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = -PetscCosReal(arrCoord[ez][ephi][er][icErmphim[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErpzm[2]] - 1.0) / arrCoord[ez][ephi][er][icErpzm[0]];
            arrX[ez][ephi][er][ivErpphim] = -PetscCosReal(arrCoord[ez][ephi][er][icErpphim[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = -PetscCosReal(arrCoord[ez][ephi][er][icErmphip[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErmzp[2]] - 1.0) / arrCoord[ez][ephi][er][icErmzp[0]];
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = -PetscCosReal(arrCoord[ez][ephi][er][icErpphip[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErpzp[2]] - 1.0) / arrCoord[ez][ephi][er][icErpzp[0]];
          }
        }
        if (user -> ictype == 2) {
          /* B(r, phi, z, t) = sin(pi*r) exp(-t) e_phi */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icBphim[0]]) * PetscExpReal(-time);
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icBphip[0]]) * PetscExpReal(-time);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t) = (pi*cos(pi*r) + sin(pi*r)/r) exp(-t) e_z */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphim[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphim[0]]) / arrCoord[ez][ephi][er][icErmphim[0]]) * PetscExpReal(-time);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphim[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphim[0]]) / arrCoord[ez][ephi][er][icErpphim[0]]) * PetscExpReal(-time);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphip[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphip[0]]) / arrCoord[ez][ephi][er][icErmphip[0]]) * PetscExpReal(-time);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphip[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphip[0]]) / arrCoord[ez][ephi][er][icErpphip[0]]) * PetscExpReal(-time);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        }
        if (user -> ictype == 3) {
          /* B(r, phi, z, t) = cos(phi) e_r - sin(phi) e_phi */
          arrX[ez][ephi][er][ivBrm] = PetscCosReal(arrCoord[ez][ephi][er][icBrm[1]]);
          arrX[ez][ephi][er][ivBphim] = -PetscSinReal(arrCoord[ez][ephi][er][icBphim[1]]);
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = PetscCosReal(arrCoord[ez][ephi][er][icBrp[1]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = -PetscSinReal(arrCoord[ez][ephi][er][icBphip[1]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t) = 0 */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        }
        if (user -> ictype == 4) {
          /* B(r, phi, z, t) = e_z */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = 1.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 1.0;
          }
          /* E(r, phi, z, t) = 0 */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        }
        if (user -> ictype == 5) {
          /* B(r, phi, z, t) = (1/r) e_r */
          arrX[ez][ephi][er][ivBrm] = 1.0 / arrCoord[ez][ephi][er][icBrm[0]];
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 1.0 / arrCoord[ez][ephi][er][icBrp[0]];
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t) = 0 */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        }
        if (user -> ictype == 6) {
          /* B(r, phi, z, t<=t_1) = (1 - t/r^2) e_phi */
          /* B(r, phi, z, t_1<t<=t_2) = (1 - 2dt/r^2 + time/r^2 - 3(time^2 - dt^2)/(2r^4) ) e_phi */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          /*arrX[ez][ephi][er][ivBphim] = 1.0 - time / PetscPowReal(arrCoord[ez][ephi][er][icBphim[0]],2) ;*/
          arrX[ez][ephi][er][ivBphim] = 1.0 - (time > 0.0) * (2.0 * dt - time) / PetscPowReal(arrCoord[ez][ephi][er][icBphim[0]], 2) - (time > 0.0) * 3.0 * (PetscSqr(time) - PetscSqr(dt)) / (2.0 * PetscPowReal(arrCoord[ez][ephi][er][icBphim[0]], 4));
          if (user -> debug) {
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "X_exact(Bphim) = %g\n", (double) arrX[ez][ephi][er][ivBphim]));
          }
          arrX[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            /*arrX[ez][ephi][er][ivBphip] = 1.0 - time / PetscPowReal(arrCoord[ez][ephi][er][icBphip[0]],2) ;*/
            arrX[ez][ephi][er][ivBphip] = 1.0 - (time > 0.0) * (2.0 * dt - time) / PetscPowReal(arrCoord[ez][ephi][er][icBphip[0]], 2) - (time > 0.0) * 3.0 * (PetscSqr(time) - PetscSqr(dt)) / (2.0 * PetscPowReal(arrCoord[ez][ephi][er][icBphip[0]], 4));
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t<=t_1) = (1/r + t/r^3) e_z */
          /* E(r, phi, z, t<=t_2) = (1/r + 2dt/r^3 - time/r^3 + 9(time^2 - dt^2)/(2r^5) ) e_z */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          /*arrX[ez][ephi][er][ivErmphim] = 1.0 / arrCoord[ez][ephi][er][icErmphim[0]] + time / PetscPowReal(arrCoord[ez][ephi][er][icErmphim[0]],3) ; */
          arrX[ez][ephi][er][ivErmphim] = (1.0 + (time > 0.0) * (2.0 * dt - time) / PetscPowReal(arrCoord[ez][ephi][er][icErmphim[0]], 2) + (time > 0.0) * 9.0 * (PetscSqr(time) - PetscSqr(dt)) / (2.0 * PetscPowReal(arrCoord[ez][ephi][er][icErmphim[0]], 4))) / arrCoord[ez][ephi][er][icErmphim[0]];
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            /*arrX[ez][ephi][er][ivErpphim] = 1.0 / arrCoord[ez][ephi][er][icErpphim[0]] + time / PetscPowReal(arrCoord[ez][ephi][er][icErpphim[0]],3);*/
            arrX[ez][ephi][er][ivErpphim] = (1.0 + (time > 0.0) * (2.0 * dt - time) / PetscPowReal(arrCoord[ez][ephi][er][icErpphim[0]], 2) + (time > 0.0) * 9.0 * (PetscSqr(time) - PetscSqr(dt)) / (2.0 * PetscPowReal(arrCoord[ez][ephi][er][icErpphim[0]], 4))) / arrCoord[ez][ephi][er][icErpphim[0]];
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            /*arrX[ez][ephi][er][ivErmphip] = 1.0 / arrCoord[ez][ephi][er][icErmphip[0]] + time / PetscPowReal(arrCoord[ez][ephi][er][icErmphip[0]],3);*/
            arrX[ez][ephi][er][ivErmphip] = (1.0 + (time > 0.0) * (2.0 * dt - time) / PetscPowReal(arrCoord[ez][ephi][er][icErmphip[0]], 2) + (time > 0.0) * 9.0 * (PetscSqr(time) - PetscSqr(dt)) / (2.0 * PetscPowReal(arrCoord[ez][ephi][er][icErmphip[0]], 4))) / arrCoord[ez][ephi][er][icErmphip[0]];
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            /*arrX[ez][ephi][er][ivErpphip] = 1.0 / arrCoord[ez][ephi][er][icErpphip[0]] + time / PetscPowReal(arrCoord[ez][ephi][er][icErpphip[0]],3);*/
            arrX[ez][ephi][er][ivErpphip] = (1.0 + (time > 0.0) * (2.0 * dt - time) / PetscPowReal(arrCoord[ez][ephi][er][icErpphip[0]], 2) + (time > 0.0) * 9.0 * (PetscSqr(time) - PetscSqr(dt)) / (2.0 * PetscPowReal(arrCoord[ez][ephi][er][icErpphip[0]], 4))) / arrCoord[ez][ephi][er][icErpphip[0]];
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        }
        if (user -> ictype == 7) {
          /* B(r, phi, z, t) = (cos(pi*t) + sin(phi)) e_z */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = PetscSinReal(arrCoord[ez][ephi][er][icBzm[1]]) + PetscCosReal(PETSC_PI * time);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = PetscSinReal(arrCoord[ez][ephi][er][icBzp[1]]) + PetscCosReal(PETSC_PI * time);
          }
          /* E(r, phi, z, t) = (cos(phi)/r) e_r */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = PetscCosReal(arrCoord[ez][ephi][er][icEphimzm[1]]) / arrCoord[ez][ephi][er][icEphimzm[0]];
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = PetscCosReal(arrCoord[ez][ephi][er][icEphipzm[1]]) / arrCoord[ez][ephi][er][icEphipzm[0]];
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = PetscCosReal(arrCoord[ez][ephi][er][icEphimzp[1]]) / arrCoord[ez][ephi][er][icEphimzp[0]];
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = PetscCosReal(arrCoord[ez][ephi][er][icEphipzp[1]]) / arrCoord[ez][ephi][er][icEphipzp[0]];
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }
        }
        if (user -> ictype == 8) {
          /* B(r, phi, z, t) = exp(sin(phi)*r+t) e_z */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 0.0;
          arrX[ez][ephi][er][ivBzm] = PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icBzm[1]]) * arrCoord[ez][ephi][er][icBzm[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icBzp[1]]) * arrCoord[ez][ephi][er][icBzp[0]]);
          }
          /* E(r, phi, z, t) = exp(sin(phi)*r+t) * (cos(phi) e_r - sin(phi) e_phi) */
          arrX[ez][ephi][er][ivErmzm] = -PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icErmzm[1]]) * arrCoord[ez][ephi][er][icErmzm[0]]) * PetscSinReal(arrCoord[ez][ephi][er][icErmzm[1]]);
          arrX[ez][ephi][er][ivEphimzm] = PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icEphimzm[1]]) * arrCoord[ez][ephi][er][icEphimzm[0]]) * PetscCosReal(arrCoord[ez][ephi][er][icEphimzm[1]]);
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = -PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icErpzm[1]]) * arrCoord[ez][ephi][er][icErpzm[0]]) * PetscSinReal(arrCoord[ez][ephi][er][icErpzm[1]]);
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icEphipzm[1]]) * arrCoord[ez][ephi][er][icEphipzm[0]]) * PetscCosReal(arrCoord[ez][ephi][er][icEphipzm[1]]);
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = -PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icErmzp[1]]) * arrCoord[ez][ephi][er][icErmzp[0]]) * PetscSinReal(arrCoord[ez][ephi][er][icErmzp[1]]);
            arrX[ez][ephi][er][ivEphimzp] = PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icEphimzp[1]]) * arrCoord[ez][ephi][er][icEphimzp[0]]) * PetscCosReal(arrCoord[ez][ephi][er][icEphimzp[1]]);
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icEphipzp[1]]) * arrCoord[ez][ephi][er][icEphipzp[0]]) * PetscCosReal(arrCoord[ez][ephi][er][icEphipzp[1]]);
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = -PetscExpReal(time + PetscSinReal(arrCoord[ez][ephi][er][icErpzp[1]]) * arrCoord[ez][ephi][er][icErpzp[0]]) * PetscSinReal(arrCoord[ez][ephi][er][icErpzp[1]]);
          }
        }

        if (user -> ictype == 9) {
          /* B(r, phi, z, t=0) Loaded from ReadInitialData output */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          countr = countr + 1;

          arrX[ez][ephi][er][ivBphim] = 0.0;
          countphi = countphi + 1;

          arrX[ez][ephi][er][ivBzm] = 0.0;
          countz = countz + 1;

          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
            countr = countr + 1;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 0.0;
            countphi = countphi + 1;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = 0.0;
            countz = countz + 1;
          }

          //if (countr > user -> numr) printf("countr = [%d] exceeds numr = [%d]\n", countr, user -> numr);
          //if (countphi > user -> numphi) printf("countphi = [%d] exceeds numphi = [%d]\n", countphi, user -> numphi);
          //if (countz > user -> numz) printf("countz = [%d] exceeds numz = [%d]\n", countz, user -> numz);
          /* tau(r, phi, z, t=0) initialized with zero values. The E field will be computed later from B field with the derived curl operator */
          arrX[ez][ephi][er][ivErmzm] = 0.0;
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        }

        if (user -> ictype == 10) {
          /* B(r, phi, z, t=0) = e_phi + sqrt(2*ln(2*r_max/r)) e_z */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 1.0;
          arrX[ez][ephi][er][ivBzm] = PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzm[0]]));
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 1.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzp[0]]));
          }
          /* tau(r, phi, z, t=0) initialized with zero values. The tau field can be computed later from B field with the derived curl operator */
          /*arrX[ez][ephi][er][ivErmzm] = 0.0;
            arrX[ez][ephi][er][ivEphimzm] = 0.0;
            arrX[ez][ephi][er][ivErmphim] = 0.0;
            if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
            }
            if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = 0.0;
            }
            if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
            }
            if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
            }
            if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
            }
            if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
            }*/
          /* tau(r, phi, z, t=0) = eta / (mu0 * V_A) / r * [e_z + sqrt(2*ln(2*r_max/r))^{-1} e_phi] */
          arrX[ez][ephi][er][ivErmzm] = condu(er, ephi, ez, BACK_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzm[0]])));
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = condu(er, ephi, ez, DOWN_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = condu(er, ephi, ez, BACK_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzm[0]])));
            arrX[ez][ephi][er][ivErpphim] = condu(er, ephi, ez, DOWN_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = condu(er, ephi, ez, UP_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = condu(er, ephi, ez, FRONT_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzp[0]])));
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = condu(er, ephi, ez, UP_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = condu(er, ephi, ez, FRONT_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzp[0]])));
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        }

        if (user -> ictype == 11) {
          /* B(r, phi, z, t=0) = [2+cos(eta/mu0 * t)] * [e_phi + sqrt(2*ln(2*r_max/r)) e_z] */
          arrX[ez][ephi][er][ivBrm] = 0.0;
          arrX[ez][ephi][er][ivBphim] = 2.0 + PetscCosReal((user->eta / user->mu0) * time);
          arrX[ez][ephi][er][ivBzm] = (2.0 + PetscCosReal((user->eta / user->mu0) * time)) * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzm[0]]));
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 2.0 + PetscCosReal((user->eta / user->mu0) * time);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = (2.0 + PetscCosReal((user->eta / user->mu0) * time)) * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzp[0]]));
          }
          /* tau(r, phi, z, t=0) = eta / (mu0 * V_A) * (cos(eta/mu0 * t) + 2)/r * [e_z + sqrt(2*ln(2*r_max/r))^{-1} e_phi] */
          if(er == 0 || ez == 0){
            arrX[ez][ephi][er][ivErmzm] = 0.0;
          }
          else{
            arrX[ez][ephi][er][ivErmzm] = condu(er, ephi, ez, BACK_LEFT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, BACK_LEFT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzm[0]])));
          }
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          if((ephi == 0 && !user->phibtype) || er == 0){
            arrX[ez][ephi][er][ivErmphim] = 0.0;
          }
          else{
            arrX[ez][ephi][er][ivErmphim] = condu(er, ephi, ez, DOWN_LEFT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, DOWN_LEFT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphim[0]]);
          }
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = 0.0;
            arrX[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = condu(er, ephi, ez, UP_LEFT, user) * (2.0 + PetscCosReal(condu(er, ephi, ez, UP_LEFT, user) / user->mu0 * time)) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = 0.0;
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = 0.0;
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        }

        if (user -> ictype == 12) {
          /* B(r, phi, z, t=0) = exp(- z / z_max) / r e_r + (r - z) e_phi + cos(2*r_max/r) e_z */
          arrX[ez][ephi][er][ivBrm] = PetscExpReal(- arrCoord[ez][ephi][er][icBrm[2]] / user->zmax) / arrCoord[ez][ephi][er][icBrm[0]];
          arrX[ez][ephi][er][ivBphim] = arrCoord[ez][ephi][er][icBphim[0]] - arrCoord[ez][ephi][er][icBphim[2]];
          arrX[ez][ephi][er][ivBzm] = PetscCosReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzm[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = PetscExpReal(- arrCoord[ez][ephi][er][icBrp[2]] / user->zmax) / arrCoord[ez][ephi][er][icBrp[0]];
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = arrCoord[ez][ephi][er][icBphip[0]] - arrCoord[ez][ephi][er][icBphip[2]];
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = PetscCosReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzp[0]]);
          }
          /* tau(r, phi, z, t=0) = eta / (mu0 * V_A) / r * [e_z + sqrt(2*ln(2*r_max/r))^{-1} e_phi] */
          arrX[ez][ephi][er][ivErmzm] = condu(er, ephi, ez, BACK_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzm[0]])));
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = condu(er, ephi, ez, DOWN_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = condu(er, ephi, ez, BACK_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzm[0]])));
            arrX[ez][ephi][er][ivErpphim] = condu(er, ephi, ez, DOWN_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = condu(er, ephi, ez, UP_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = condu(er, ephi, ez, FRONT_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzp[0]])));
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = condu(er, ephi, ez, UP_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = condu(er, ephi, ez, FRONT_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzp[0]])));
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        }

        if (user -> ictype == 13) {
          /* B(r, phi, z, t=0) = e_phi + sqrt(2*ln(2*r_max/r)) e_z */
          if(er!=0) {
            arrX[ez][ephi][er][ivBrm] = 0.0;
          }
          arrX[ez][ephi][er][ivBphim] = 1.0;
          if(ez!=0) {
            arrX[ez][ephi][er][ivBzm] = PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icBzm[0]]));
          }
          if(er==0) {
            arrX[ez][ephi][er][ivBrm] = NAN;
          }
          if(ez==0) {
            arrX[ez][ephi][er][ivBzm] = NAN;
          }
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivBrp] = NAN;
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivBphip] = 1.0;
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivBzp] = NAN;
          }
          /* tau(r, phi, z, t=0) = 0 */
          /*
             arrX[ez][ephi][er][ivErmzm] = 0.0;
             arrX[ez][ephi][er][ivEphimzm] = 0.0;
             arrX[ez][ephi][er][ivErmphim] = 0.0;
             if (er == N[0] - 1) {
             arrX[ez][ephi][er][ivErpzm] = 0.0;
             arrX[ez][ephi][er][ivErpphim] = 0.0;
             }
             if (ephi == N[1] - 1) {
             arrX[ez][ephi][er][ivEphipzm] = 0.0;
             arrX[ez][ephi][er][ivErmphip] = 0.0;
             }
             if (ez == N[2] - 1) {
             arrX[ez][ephi][er][ivErmzp] = 0.0;
             arrX[ez][ephi][er][ivEphimzp] = 0.0;
             }
             if (er == N[0] - 1 && ephi == N[1] - 1) {
             arrX[ez][ephi][er][ivErpphip] = 0.0;
             }
             if (ephi == N[1] - 1 && ez == N[2] - 1) {
             arrX[ez][ephi][er][ivEphipzp] = 0.0;
             }
             if (er == N[0] - 1 && ez == N[2] - 1) {
             arrX[ez][ephi][er][ivErpzp] = 0.0;
             }*/

          /* tau(r, phi, z, t=0) = eta / (mu0 * V_A) / r * [e_z + sqrt(2*ln(2*r_max/r))^{-1} e_phi] */
          arrX[ez][ephi][er][ivErmzm] = condu(er, ephi, ez, BACK_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzm[0]])));
          arrX[ez][ephi][er][ivEphimzm] = 0.0;
          arrX[ez][ephi][er][ivErmphim] = condu(er, ephi, ez, DOWN_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrX[ez][ephi][er][ivErpzm] = condu(er, ephi, ez, BACK_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzm[0]])));
            arrX[ez][ephi][er][ivErpphim] = condu(er, ephi, ez, DOWN_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivEphipzm] = 0.0;
            arrX[ez][ephi][er][ivErmphip] = condu(er, ephi, ez, UP_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErmzp] = condu(er, ephi, ez, FRONT_LEFT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzp[0]])));
            arrX[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrX[ez][ephi][er][ivErpphip] = condu(er, ephi, ez, UP_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrX[ez][ephi][er][ivErpzp] = condu(er, ephi, ez, FRONT_RIGHT, user) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzp[0]])));
          }

          /* n_i(r, phi, z, t=0) = ni_0 */
          arrX[ez][ephi][er][ivn] = 1.0;

          /* Vi_perp(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrmphipzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphimzp[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzm[2]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[0]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[1]] = 0.0;
          arrX[ez][ephi][er][ivVrpphipzp[2]] = 0.0;

          /* EP(r, phi, z, t=0) = 0 */
          arrX[ez][ephi][er][ivVrmphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrmphipzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphimzp[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzm[3]] = 0.0;

          arrX[ez][ephi][er][ivVrpphipzp[3]] = 0.0;
        }
      }
    }
  }
  /* Restore vectors */
  PetscCall(DMStagVecRestoreArray(da, xLocal, & arrX));
  PetscCall(DMLocalToGlobal(da, xLocal, INSERT_VALUES, * X));
  PetscCall(DMRestoreLocalVector(da, & xLocal));
  PetscCall(DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord));

  if (user -> ictype == 9 && user -> Ebc) {
    Vec Xcopy;
    PetscCall(VecDuplicate( * X, & Xcopy));
    PetscCall(VecCopy( * X, Xcopy));
    PetscCall(VecScale(Xcopy, 1.0 / 11000000.0));
    FormDerivedCurl(ts, Xcopy, * X, user); //This updates only the E field part in X by computing the derived mimetic curl operator applied to B field of Xcopy
                                           //PetscBarrier((PetscObject) Xcopy);
    PetscCall(VecDestroy( & Xcopy));
  }

#line 16521

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

#line 17555

PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal time, Vec X, void * ptr) {
  PetscFunctionBeginUser;

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("Monitor",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  User * user = (User * ) ptr;
  Mat D;
  Vec exactX, div, div2, curl, postprocX, B;
  PetscReal norm2, normmax, dt;
  PetscInt nr, nphi, nz;
  DM da, newda;


  PetscCall(TSGetTimeStep(ts, & dt));
  /* Adjust the viscosity coefficient and step size over time*/
  //TSSetTimeStep(ts,2*dt);
  //user->Re = 1.0 * PetscExpReal(4*(time + dt - user->itime)/10.0);
  //user->Re = 1.0 + (9.0/20.0) * (time + dt - user->itime);
  //user->Re = 1.0/(0.001 + 0.999 * PetscExpReal(-(time + dt - user->itime)/0.53));
  //user->Re = 1.0/(0.1 + 0.9 * PetscExpReal(-(time + dt - user->itime)/10.0));
  //if(step>100){
  //TSSetTimeStep(ts,200.0);
  //}
  //else{
  //TSSetTimeStep(ts,user->Re/5.0);
  //}
  //user->Re = 1.0/(0.01 + 0.09 * PetscExpReal(-(time + dt - user->itime)/100.0));

  // diffX := X - oldX = X^{n+1} - X^{n} for the corrector
  {
    Vec diffX;
    PetscReal norm_2, norm_max;
    PetscCall(VecDuplicate(X, & diffX));
    PetscCall(VecZeroEntries(diffX));
    PetscCall(VecWAXPY(diffX, -1.0, user -> X0, X));
    PetscCall(VecNorm(diffX, NORM_2, & norm_2));
    PetscCall(VecNorm(diffX, NORM_MAX, & norm_max));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Timestep %3D (CORRECTOR): step size = %g, time = %g, 2-norm of X^{n+1} - X^{n} = %g, max norm of X^{n+1} - X^{n} = %g\n", (int)(step + user -> oldstep), (double) dt, (double) time, (double) norm_2, (double) norm_max));
    PetscCall(VecDestroy(& diffX));
  }
  /* Compute the exact solution */
  PetscCall(TSGetDM(ts, & da));
  if (user -> ictype != 9){
    PetscCall(DMCreateGlobalVector(da, & exactX));
    FormExactSolution(time, ts, & exactX, user);
  }

  /* Save intermediate and final solutions in .dat files */
  if (user -> tempdump) {
    if ((time >= user -> ftime) || (step > 0 && (step % user -> dumpfreq == 0)))
    {
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Saving intermediate solution at time %g, step %d\n", (double) time, step));
      SaveIntermediateSolution(ts, (int)(step + user -> oldstep), time, X, user);
      char filename[PETSC_MAX_PATH_LEN];
      snprintf(filename, 50, "./snapshot_%d.h5", (int)(step + user -> oldstep));
      /* RuKS_H5Write(user->runaway_solver, filename); */
    }
  }
  //if (user->step_gl > -1) PushParticles(ts, user->oldX, X, user, 0);
  //user->step_gl++;
  /* Save solution in .vtr files */
  // PetscPrintf(PETSC_COMM_WORLD,"About to check for dump, dump = %g", (int) user -> dump);

  if (user -> dump) {
    // PetscPrintf(PETSC_COMM_WORLD,"dump flag checked");

    DumpSolution_Cell(ts, (int)(step + user -> oldstep), X, user);
    if (1 && user -> ictype == 9 && user -> oldstep == 0) {
      DumpLevelSet(ts, user);
    }
  }

#line 17653
  /* Save SNES residual */
#line 17670
  /* Save curl(B)xB and curl(B) */
#line 17697
  /* Save curl(tau) */
#line 17706

  /* Compute the 2-norm and max-norm of the error */
  if (user -> ictype != 9){
    PetscCall(VecAXPY(exactX, -1.0, X));
    PetscCall(VecNorm(exactX, NORM_2, & norm2));
    /* norm_2 = PetscSqrtReal(appctx->h)*norm_2; Scale the 2-norm by the grid spacing */
    norm2 = norm2 * PetscSqrtReal(user -> dr * user -> dphi * user -> dz); /* Scale the 2-norm by the grid spacing */
    PetscCall(VecNorm(exactX, NORM_MAX, & normmax));
  }
  /* Save Error in .vtr files */
  if (user -> dump && user -> ictype != 9 && user -> ictype != 10) {
    DumpError(ts, (int)(step + user -> oldstep), exactX, user);
  }
  /* Display information at each time step */
  if (user -> ictype != 9 && user -> ictype != 10) PetscPrintf(PETSC_COMM_WORLD, "Timestep %3D: step size = %g, time = %g, 2-norm error = %g, max norm error = %g\n", (int)(step + user -> oldstep), (double) dt, (double) time, (double) norm2, (double) normmax);

  /* Print debugging information if desired */
  if (user -> debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error vector\n"));
    PetscCall(VecView(exactX, PETSC_VIEWER_STDOUT_WORLD));
  }

  /* Back substitution : B^{n+1} := B^{n} - dt * primary_curl(E^{n+1}) */
  PetscCall(DMCreateGlobalVector(da, & curl));
  PetscCall(DMCreateGlobalVector(da, & postprocX));

  if (step != 0) {
    FormPrimaryCurl(ts, X, curl, user);
    VecAXPBYPCZ(postprocX, 1.0, -user -> dt, 0.0, user -> X0, curl); /* postprocX = oldX - dt * curl */
  }
  PetscCall(VecCopy(X, user -> X0));


  /* Check that the solution is divergence-free */
  /* Create a new DM to store cell-centered divergence values */
  PetscCall(DMStagGetNumRanks(da, & nr, & nphi, & nz));
  if (user -> phibtype) {
    PetscCall(DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE, user -> Nr, user -> Nphi, user -> Nz, nr, nphi, nz, 0, 0, 0, 1, DMSTAG_STENCIL_BOX, 1, NULL, NULL, NULL, & newda));
  } else {
    PetscCall(DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, user -> Nr, user -> Nphi, user -> Nz, nr, nphi, nz, 0, 0, 0, 1, DMSTAG_STENCIL_BOX, 1, NULL, NULL, NULL, & newda));
  }
  PetscCall(DMSetFromOptions(newda));
  PetscCall(DMSetUp(newda));
  PetscCall(DMStagSetUniformCoordinatesExplicit(newda, user->rmin/user->L0, user->rmax/user->L0, user->phimin, user->phimax, user->zmin/user->L0, user->zmax/user->L0));

  PetscCall(DMCreateGlobalVector(newda, & div));
  PetscCall(DMCreateMatrix(user -> coorda, & D));
  FormDiscreteDivergence(ts, newda, D, X, div, user);
  PetscCall(VecScale(div, 1.0/user->L0));
  if (user -> debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Divergence operator matrix: \n"));
    PetscCall(MatView(D, PETSC_VIEWER_STDOUT_WORLD));

    /* MatMult(D,X,Y); */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "X vector: \n"));
    PetscCall(VecView(X, PETSC_VIEWER_STDOUT_WORLD));

    /* PetscPrintf(PETSC_COMM_WORLD,"Divergence of B vector: \n");
       VecView(Y,PETSC_VIEWER_STDOUT_WORLD); */

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Divergence of B vector: \n"));
    PetscCall(VecView(div, PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(VecNorm(div, NORM_MAX, & normmax));
  if (user -> ictype != 9 && user -> ictype != 10){
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "max norm error of divergence of B = %g\n", (double) normmax));

  }
  else {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Timestep %3D: step size = %g, time = %g, max norm error of divergence of B = %g\n", (int)(step + user -> oldstep), (double) dt, (double) time, (double) normmax));
  }

  if (user -> dump) {
    DumpDivergence(ts, newda, (int)(step + user -> oldstep), div, user);
  }

  if (user -> tstype > 1) {
    PetscCall(VecGetSubVector(X, user -> isB, & B));
    PetscCall(VecNorm(B, NORM_MAX, & norm2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "max_norm(div B) / max_norm(B) = %g\n", (double) normmax / norm2));
    PetscCall(VecRestoreSubVector(X, user -> isB, & B));
  }

  /* Check that the post-processed solution is divergence-free */
  PetscCall(DMCreateGlobalVector(newda, & div2));
  if (step != 0) {
    FormDiscreteDivergence(ts, newda, D, postprocX, div2, user);
    PetscCall(VecNorm(div2, NORM_MAX, & normmax));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "max norm error of divergence of post-processed B = %g\n", (double) normmax));
    if (user -> tstype > 1) {
      PetscCall(VecGetSubVector(postprocX, user -> isB, & B));
      PetscCall(VecNorm(B, NORM_MAX, & norm2));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "max_norm(div postprocB) / max_norm(postprocB) = %g\n", (double) normmax / norm2));
      PetscCall(VecRestoreSubVector(postprocX, user -> isB, & B));
    }
  }

  ComputeCurrent(ts, X, user);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Current intensity inside plasma = %g\n", (double)(user -> Iphi1)));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Current intensity outside plasma = %g\n", (double)(user -> Iphi2)));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Current intensity inside vacuum vessel = %g\n", (double)(user -> Iphi3)));

  FILE *fp;
  fp = fopen("TotalCurrent.txt", "a");
  PetscCall(PetscFPrintf(PETSC_COMM_WORLD, fp, " %.8f  %.8f \n", (double) time, (double) (user -> Iphi1) ));
  fclose(fp);

  /* DESTROY VECTORS */
  if (user -> ictype != 9){
    PetscCall(VecDestroy( & exactX));
  }
  PetscCall(MatDestroy( & D));
  PetscCall(VecDestroy( & div));
  PetscCall(VecDestroy( & div2));
  PetscCall(VecDestroy( & curl));
  PetscCall(VecDestroy( & postprocX));
  PetscCall(DMDestroy( & newda));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormDummyIJacobian4(TS ts, Vec X, Vec Xdot, PetscReal a, Mat J, Mat Jpre, void * ptr) {
  PetscFunctionBeginUser;
  User * user = (User * ) ptr;
  DM da;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  PetscInt N[3], er, ephi, ez;

  PetscCall(MatZeroEntries(Jpre));
  PetscCall(TSGetDM(ts, & da));
  PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
  PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col[1];
        PetscScalar valJ[1];
        PetscInt nEntries = 1;

        /* Electrostatic potential part: Phi */

        nEntries = 1;
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_DOWN_LEFT;
        row.c = 3;

        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = BACK_DOWN_LEFT;
        col[0].c = 3;
        valJ[0] = 1.0;

        PetscCall(DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES));

        if (er == N[0] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_RIGHT;
          row.c = 3;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN_RIGHT;
          col[0].c = 3;
          valJ[0] = 1.0;

          PetscCall(DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES));
        }
        if (ephi == N[1] - 1 && !user -> phibtype) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_LEFT;
          row.c = 3;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_LEFT;
          col[0].c = 3;
          valJ[0] = 1.0;

          PetscCall(DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES));
        }
        if (ez == N[2] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_LEFT;
          row.c = 3;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN_LEFT;
          col[0].c = 3;
          valJ[0] = 1.0;

          PetscCall(DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES));
        }
        if (er == N[0] - 1 && ephi == N[1] - 1 && !user -> phibtype) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_RIGHT;
          row.c = 3;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_RIGHT;
          col[0].c = 3;
          valJ[0] = 1.0;

          PetscCall(DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES));
        }
        if (ez == N[2] - 1 && ephi == N[1] - 1 && !user -> phibtype) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_LEFT;
          row.c = 3;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_LEFT;
          col[0].c = 3;
          valJ[0] = 1.0;

          PetscCall(DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES));
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_RIGHT;
          row.c = 3;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN_RIGHT;
          col[0].c = 3;
          valJ[0] = 1.0;

          PetscCall(DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES));
        }
        if (er == N[0] - 1 && ez == N[2] - 1 && ephi == N[1] - 1 && !user -> phibtype) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_RIGHT;
          row.c = 3;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_RIGHT;
          col[0].c = 3;
          valJ[0] = 1.0;

          PetscCall(DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES));
        }
      }
    }
  }

  PetscCall(MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY));

  if (J != Jpre) {
    PetscCall(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY));
  }
  if (user -> debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Jpre:\n"));
    PetscCall(MatView(Jpre, PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode SampleShellPCSetUp(PC pc){
  PetscFunctionBeginUser;
  User * user;
  SNES snes;
  Mat J;
  KSP KSP_ETBN;
  PC PC_ETBN, PC_V;
  PetscErrorCode ierr = 0;

  PCShellGetContext(pc,&user);
  Mat myS;

  /* Creating matrices */
  ierr = MatCreate(PETSC_COMM_WORLD,&myS); CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,& user->DiagBlock_V); CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,& user->OffDiagBlock_L); CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,& user->OffDiagBlock_U); CHKERRQ(ierr);

  const IS islist2[5] = {user->isV, user->isEP, user->istau, user->isB, user->isni};
  IS isALL;

  ierr = ISConcatenate(PETSC_COMM_WORLD,5,islist2,&isALL); CHKERRQ(ierr);

  //PetscPrintf(PETSC_COMM_WORLD, "The global indices for all fields are as follows.\n");
  //ISView(isALL,PETSC_VIEWER_STDOUT_SELF);
  //return 0;

  ierr = ISDifference(isALL, user->isV, & user->isALL_V); CHKERRQ(ierr);

  //PetscPrintf(PETSC_COMM_WORLD, "The global indices for all fields except V are as follows.\n");
  //ISView(isALL_V,PETSC_VIEWER_STDOUT_SELF);
  //return 0;

  ierr = TSGetSNES(user->ts, & snes); CHKERRQ(ierr);
  SNESGetJacobian(snes,&J,NULL,NULL,NULL);

  ierr = MatGetSchurComplement (J, user->isALL_V, user->isALL_V, user->isV, user->isV, MAT_INITIAL_MATRIX, & user->DiagBlock_V, MAT_SCHUR_COMPLEMENT_AINV_DIAG, MAT_INITIAL_MATRIX, &myS); CHKERRQ(ierr);

  ierr = MatDestroy(&myS); CHKERRQ(ierr);
  ierr = ISDestroy( & isALL); CHKERRQ(ierr);

  PetscCall(MatCreateSubMatrix(J,user->isALL_V,user->isV,MAT_INITIAL_MATRIX,& user->OffDiagBlock_U));
  PetscCall(MatCreateSubMatrix(J,user->isV,user->isALL_V,MAT_INITIAL_MATRIX,& user->OffDiagBlock_L));

  /* Creating and setting KSP Solvers*/
  ierr = MatSchurComplementGetKSP(user->DiagBlock_V,&KSP_ETBN);CHKERRQ(ierr);
  //KSPSetType(KSP_ETBN,KSPPREONLY);
  KSPAppendOptionsPrefix(KSP_ETBN, "KSP_ETBN");
  ierr = KSPGetPC(KSP_ETBN,&PC_ETBN);CHKERRQ(ierr);
  PetscCall(PCSetType(PC_ETBN,PCASM));
  PCASMSetOverlap(PC_ETBN,7);
#line 18072

  //PCSetType(PC_ETBN,PCLU);
  //PCFactorSetMatSolverType(PC_ETBN,MATSOLVERSUPERLU_DIST);

  ierr = KSPCreate(PETSC_COMM_WORLD, & user->KSP_V);CHKERRQ(ierr);
  ierr = KSPSetOperators(user->KSP_V, user->DiagBlock_V, user->DiagBlock_V);CHKERRQ(ierr);
  ierr = KSPGetPC(user->KSP_V, & PC_V);CHKERRQ(ierr);
  PetscCall(PCSetType(PC_V, PCNONE));
  ierr = PCSetUp(PC_V);CHKERRQ(ierr);
  ierr = KSPSetUp(user->KSP_V);CHKERRQ(ierr);
  KSPAppendOptionsPrefix(user->KSP_V, "KSP_V");
  ierr = KSPSetTolerances(user->KSP_V, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 200); //use 200 outer iterations for the V solve
  CHKERRQ(ierr);

#line 18091

  PetscFunctionReturn(PETSC_SUCCESS);
}

#line 18151

PetscErrorCode SampleShellPCApply(PC pc,Vec x,Vec y){
  PetscFunctionBeginUser;
  User * user;
  Vec x_ETBN, x_V, y_ETBN, y_V;
  Vec z_ETBN, z_V;
  PetscErrorCode ierr = 0;
  KSP KSP_ETBN;

  ierr = PCShellGetContext(pc,&user);CHKERRQ(ierr);
  //VecZeroEntries(y);
  PetscCall(VecCopy(x, y));

  ierr = MatSchurComplementGetKSP(user->DiagBlock_V,&KSP_ETBN);CHKERRQ(ierr);

  PetscCall(VecGetSubVector( y, user -> isALL_V, & y_ETBN));
  PetscCall(VecGetSubVector( x, user -> isALL_V, & x_ETBN));
  PetscCall(VecDuplicate(x_ETBN, & z_ETBN));

  PetscCall(VecGetSubVector( y, user -> isV, & y_V));
  PetscCall(VecGetSubVector( x, user -> isV, & x_V));
  PetscCall(VecDuplicate(x_V, & z_V));

  // y_ETBN := J_{ETBN}^{-1} x_ETBN
  ierr = KSPSolve(KSP_ETBN, x_ETBN, y_ETBN);CHKERRQ(ierr);
  // z_V := x_V - J_{V,ETBN} y_ETBN
  ierr = VecScale(y_ETBN, -1.0);CHKERRQ(ierr);
  ierr = MatMultAdd(user->OffDiagBlock_L,y_ETBN,x_V,z_V);CHKERRQ(ierr);
  ierr = VecScale(y_ETBN, -1.0);CHKERRQ(ierr);
  // y_V := (J_V - J_{V,ETBN} * J_{ETBN}^{-1} J_{ETBN,V} )^{-1} z_V
  ierr = KSPSolve(user->KSP_V, z_V, y_V);CHKERRQ(ierr);
  // z_ETBN := J_{ETBN,V} y_V
  ierr = MatMult(user->OffDiagBlock_U, y_V, z_ETBN);CHKERRQ(ierr);
  // z_ETBN := J_{ETBN}^{-1} z_ETBN
  ierr = KSPSolve(KSP_ETBN, z_ETBN, z_ETBN);CHKERRQ(ierr);
  // y_ETBN := y_ETBN - z_ETBN
  ierr = VecAXPY(y_ETBN, -1.0, z_ETBN);CHKERRQ(ierr);

  PetscCall(VecRestoreSubVector( y, user -> isALL_V, & y_ETBN));
  PetscCall(VecRestoreSubVector( x, user -> isALL_V, & x_ETBN));

  PetscCall(VecRestoreSubVector( y, user -> isV, & y_V));
  PetscCall(VecRestoreSubVector( x, user -> isV, & x_V));

  PetscCall(VecDestroy(& z_V));
  PetscCall(VecDestroy(& z_ETBN));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode SampleShellPCDestroy(PC pc){
  PetscFunctionBeginUser;
  User * user;

  PCShellGetContext(pc,&user);
  /* Destroy matrices */
  PetscCall(MatDestroy( & user->DiagBlock_V));
  PetscCall(MatDestroy( & user->OffDiagBlock_L));
  PetscCall(MatDestroy( & user->OffDiagBlock_U));
  /* Destroy KSP Solvers */
  PetscCall(KSPDestroy( & user->KSP_V));
  /* Destroy Index Set */
  PetscCall(ISDestroy( & user->isALL_V));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#line 18273

#line 18304

#line 18317

#line 18403

#line 18445

#line 18460

PetscErrorCode ReadInitialData(PetscReal ** data, PetscInt * num,
    const char * file_name) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("ReadInitialData",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  PetscInt i;
  char filename[PETSC_MAX_PATH_LEN];

  snprintf(filename, sizeof(filename), "%s", file_name);
  FILE * pFile = fopen(filename, "r");
  if (pFile == NULL) SETERRQ(PETSC_COMM_SELF, 1, "Incorrect location of data");
  fscanf(pFile, "%d", num);

  PetscCall(PetscMalloc1( * num, data));

  i = 0;
  while (fscanf(pFile, "%lf", & (( * data)[i])) == 1) {
    //printf("%lf\n", (*data)[i]);
    i = i + 1;
  }
  fclose(pFile);

  if (i != * num) {
    printf("i = [%d] but num = [%d]\n", i, * num);
  }

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "============Finished reading in InitialData============\n"));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

#line 18549

#line 19664

#line 20206

PetscErrorCode FormInitialSolution_psi(TS ts, Vec X, void * ptr) {
  PetscFunctionBeginUser;

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("FormInitialSolution_psi",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  User * user = (User * ) ptr;
  DM da;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec xLocal, XcopyLocal, Xcopy;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];

  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];
  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];

  PetscInt ivn;

  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  DM dmCoord;
  PetscScalar ** ** arrCoord, ** ** arrX, ** ** arrXcopy;
  PetscInt countg = 0, countpsi = 0;

  PetscCall(VecZeroEntries(X));
  PetscCall(TSGetDM(ts, & da));
  PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
  PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
  for (d = 0; d < 4; ++d) {
    /* Vertex locations */
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]));
    PetscCall(DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]));
  }
  /* Edge locations */
  PetscCall(DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm));
  PetscCall(DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm));
  PetscCall(DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim));
  PetscCall(DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim));
  PetscCall(DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip));
  PetscCall(DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp));
  PetscCall(DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp));
  /* Face locations */
  PetscCall(DMStagGetLocationSlot(da, LEFT, 0, & ivBrm));
  PetscCall(DMStagGetLocationSlot(da, DOWN, 0, & ivBphim));
  PetscCall(DMStagGetLocationSlot(da, BACK, 0, & ivBzm));
  PetscCall(DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp));
  PetscCall(DMStagGetLocationSlot(da, UP, 0, & ivBphip));
  PetscCall(DMStagGetLocationSlot(da, FRONT, 0, & ivBzp));
  /* Cell locations */
  PetscCall(DMStagGetLocationSlot(da, ELEMENT, 0, & ivn));

  PetscCall(DMGetCoordinateDM(da, & dmCoord));
  PetscCall(DMGetCoordinatesLocal(da, & coordLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord));
  for (d = 0; d < 3; ++d) {
    /* Face coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]));
    /* Edge coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]));
  }

  PetscCall(VecDuplicate(X, & Xcopy));

  /* Compute function over the locally owned part of the grid */
  PetscCall(DMGetLocalVector(da, & XcopyLocal));
  PetscCall(DMGlobalToLocalBegin(da, Xcopy, INSERT_VALUES, XcopyLocal));
  PetscCall(DMGlobalToLocalEnd(da, Xcopy, INSERT_VALUES, XcopyLocal));
  PetscCall(DMStagVecGetArray(da, XcopyLocal , & arrXcopy));

  PetscCall(DMGetLocalVector(da, & xLocal));
  PetscCall(DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal));
  PetscCall(DMStagVecGetArray(da, xLocal, & arrX));
  //Set X's magnetic field with G(psi)/r e_phi
  //Set Xcopy's tau field with psi/r e_phi
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        /* G(psi(r, z, t=0)) Loaded from ReadInitialData output */
        arrX[ez][ephi][er][ivBphim] = user -> datag[er + ephi * N[0] + ez * N[0] * (N[1] + 1)] / (arrCoord[ez][ephi][er][icBphim[0]] * user->L0);
        countg = countg + 1;

        if (ephi == N[1] - 1) {
          arrX[ez][ephi][er][ivBphip] = user -> datag[er + (ephi + 1) * N[0] + ez * N[0] * (N[1] + 1)] / (arrCoord[ez][ephi][er][icBphip[0]] * user->L0);
          countg = countg + 1;
        }

        if (countg > user -> numg) printf("countg = [%d] exceeds numg = [%d]\n", countg, user -> numg);

        /* tau(r, phi, z, t=0) of Xcopy initialized with psi values. */
        arrXcopy[ez][ephi][er][ivErmzm] = user -> datapsi[er + ephi * (N[0] + 1) + ez * N[1] * (N[0] + 1)] / (arrCoord[ez][ephi][er][icErmzm[0]] * user->L0);
        countpsi = countpsi + 1;
        if (er == N[0] - 1) {
          arrXcopy[ez][ephi][er][ivErpzm] = user -> datapsi[er + 1 + ephi * (N[0] + 1) + ez * N[1] * (N[0] + 1)] / (arrCoord[ez][ephi][er][icErpzm[0]] * user->L0);
          countpsi = countpsi + 1;
        }
        if (ez == N[2] - 1) {
          arrXcopy[ez][ephi][er][ivErmzp] = user -> datapsi[er + ephi * (N[0] + 1) + (ez + 1) * N[1] * (N[0] + 1)] / (arrCoord[ez][ephi][er][icErmzp[0]] * user->L0);
          countpsi = countpsi + 1;
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrXcopy[ez][ephi][er][ivErpzp] = user -> datapsi[er + 1 + ephi * (N[0] + 1) + (ez + 1) * N[1] * (N[0] + 1)] / (arrCoord[ez][ephi][er][icErpzp[0]] * user->L0);
          countpsi = countpsi + 1;
        }
        if (countpsi > user -> numpsi) printf("countpsi = [%d] exceeds numpsi = [%d]\n", countpsi, user -> numpsi);

        /* n_i(r, phi, z, t=0) = ni_0 */
        arrX[ez][ephi][er][ivn] = 1.0;
      }
    }
  }

  /* Restore vectors */
  PetscCall(DMStagVecRestoreArray(da, xLocal, & arrX));
  PetscCall(DMStagVecRestoreArray(da, XcopyLocal , & arrXcopy));
  PetscCall(PetscFree(user->datag));
  PetscCall(PetscFree(user->datapsi));

  if (user -> itime == 0.0) {
    Vec curl;
    PetscCall(DMLocalToGlobal(da, xLocal, INSERT_VALUES, X));
    PetscCall(DMLocalToGlobal(da, XcopyLocal, INSERT_VALUES, Xcopy));
    PetscCall(DMRestoreLocalVector(da, & xLocal));
    PetscCall(DMRestoreLocalVector(da, & XcopyLocal));

    PetscCall(DMCreateGlobalVector(da, & curl));
    PetscCall(VecScale(curl,0.0));
    //Apply primary curl to Xcopy, and save in curl
    FormPrimaryCurl(ts, Xcopy, curl, user);
    //X := X - (1/L0)curl;
    PetscCall(VecAXPY(X,-1.0 / user->L0,curl));
    //Normalize wrt B_0
    PetscCall(VecGetSubVector( X, user -> isB, & curl));
    VecScale(curl, 1.0 / user->B0); // B := tilde{B} = (B_0^-1) B
    PetscCall(VecRestoreSubVector( X, user -> isB, & curl));
    //Destroy vectors
    PetscCall(VecDestroy(& Xcopy));
    PetscCall(VecDestroy(& curl));
  }


  // The following will replace the X with something from binary file, if no binary file is specified, the solution will come from efit and relaxation
  if(user->ic_binary_mode == 'l'){ // If loading binary
    PetscViewer viewer;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Reading X vector from file %s ...\n", user->ic_binary_path));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, user->ic_binary_path, FILE_MODE_READ, & viewer));
    PetscCall(stag_vec_io(user, viewer, X, PETSC_TRUE));
    PetscCall(PetscViewerDestroy(&viewer));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Reading from file %s is over.\n", user->ic_binary_path));
    /* Destroy the viewer */
  }  // If creating the binary
  else if ( (user -> ictype == 9 || user -> ictype == 15) && user -> itime == 0.0) {
    SNES dummysnes;
    TS dummyts;
    KSP dummyKSP;
    PC dummypc;
    Mat J;
    PetscInt steps;
    PetscReal ftime;
    TSConvergedReason reason;
    TSAdapt adapt;

    PetscCall(TSCreate(PETSC_COMM_WORLD, & dummyts));
    PetscCall(TSSetDM(dummyts, user->coorda));
    TSMonitorSet(dummyts, Monitor, user, NULL); /* Set optional user-defined monitoring routine */
    PetscCall(TSSetProblemType(dummyts, TS_NONLINEAR));

    PetscCall(TSGetSNES(dummyts, & dummysnes));
    //SNESSetType(dummysnes, SNESKSPONLY);
    TSGetAdapt(dummyts, & adapt);
    TSAdaptSetType(adapt, TSADAPTNONE);
    TSSetType(dummyts, TSARKIMEX); /* Additive Runge-Kutta IMEX method */
    TSARKIMEXSetFullyImplicit(dummyts, PETSC_TRUE);
    TSARKIMEXSetType(dummyts, TSARKIMEXL2);
    TSSetEquationType(dummyts,TS_EQ_IMPLICIT);

    PetscCall(DMCreateMatrix(da, & J));
    /* Use coloring to compute finite difference J efficiently */
    PetscCall(SNESSetJacobian(dummysnes, J, J, SNESComputeJacobianDefaultColor, PETSC_NULLPTR));
    PetscCall(TSSetIFunction(dummyts, NULL, FormIFunction_newequilibrium_Vperp, user));
    /* No RHS: FormRHSFunction_BImplicit is identically zero; see mhd.c. */

    PetscCall(SNESSetUseMatrixFree(dummysnes,PETSC_TRUE,PETSC_FALSE));
    PetscCall(SNESSetOptionsPrefix(dummysnes, "dummySNES_"));
    PetscCall(SNESSetFromOptions(dummysnes));

    //TSSetTime(dummyts, 0.0);
    PetscCall(TSSetTime(dummyts, user->itime));
    // TSSetDuration(dummyts,30,3000.0);
    PetscCall(TSSetMaxTime(dummyts, 3000.0));
    PetscCall(TSSetMaxSteps(dummyts, 30));
    PetscCall(TSSetExactFinalTime(dummyts, TS_EXACTFINALTIME_STEPOVER));
    PetscCall(TSSetTimeStep(dummyts, 100.0));
    PetscCall(TSSetSolution(dummyts, X));



    PetscCall(SNESGetKSP(dummysnes, & dummyKSP));
    //TSGetKSP(dummyts, & dummyKSP);
    PetscCall(KSPGetPC(dummyKSP, & dummypc));

    PetscCall(PCSetType(dummypc,PCBJACOBI));
    PetscCall(KSPSetOptionsPrefix(dummyKSP, "dummyKSP_"));

    /* Logic below modifies the PC directly, so this is the last chance to change the solver from the command line */
    PetscCall(KSPSetFromOptions(dummyKSP));

    PetscBool is_fieldsplit;
    {/* first level -> split ni from the rest : {ni}, {V Phi tau B}
        second level -> split tau from {V Phi B} : {ni}, {{tau},{V Phi B}}
        third level -> split V from {Phi B} : {ni}, {{tau},{{Phi B}, {V}}}
        fourth level -> split Phi from B : {ni}, {{tau},{{{Phi}, {B}}, {V}}}
        */
      IS            is[2];
      DMStagStencil stencil0[1], stencil1[10];
      PC            pc_notc, pc_noe;

      const char *name[2] = {"ni", "TEBV"};

      // First split is cells
      stencil0[0].loc = DMSTAG_ELEMENT;
      stencil0[0].c = 0;

      // Second split is the rest
      for (PetscInt c=0; c<4; ++c) {
        stencil1[c].loc = DMSTAG_BACK_DOWN_LEFT;
        stencil1[c].c = c;
      }
      stencil1[4].loc = DMSTAG_LEFT;
      stencil1[4].c = 0;
      stencil1[5].loc = DMSTAG_BACK;
      stencil1[5].c = 0;
      stencil1[6].loc = DMSTAG_DOWN;
      stencil1[6].c = 0;
      stencil1[7].loc = DMSTAG_BACK_DOWN;
      stencil1[7].c = 0;
      stencil1[8].loc = DMSTAG_BACK_LEFT;
      stencil1[8].c = 0;
      stencil1[9].loc = DMSTAG_DOWN_LEFT;
      stencil1[9].c = 0;

      PetscCall(DMStagCreateISFromStencils(da,1,stencil0,&is[0]));
      PetscCall(DMStagCreateISFromStencils(da,10,stencil1,&is[1]));

      for (PetscInt i=0; i<2; ++i) {
        PetscCall(PCFieldSplitSetIS(dummypc,name[i],is[i]));
      }

      for (PetscInt i=0; i<2; ++i) {
        PetscCall(ISDestroy(&is[i]));
      }

      /* If the fieldsplit PC wasn't overridden, further split the second split */
      {
        PCType pc_type;


        PetscCall(KSPGetPC(dummyKSP, &dummypc));
        PetscCall(PCGetType(dummypc,&pc_type));
        PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
        if (is_fieldsplit) {
          DM            dm_notc;
          KSP           *sub_ksp;

          PetscInt      n_splits;
          DMStagStencil stencil_notc_edges[3], stencil_notc_notedges[7];
          IS            is_notc[2];
          const char    *name_notc[2] = {"tau","EBV"};

          PetscCall(PCSetUp(dummypc)); // Set up the Fieldsplit PC
          PetscCall(PCFieldSplitGetSubKSP(dummypc,&n_splits,&sub_ksp));
          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
          PetscCall(KSPGetPC(sub_ksp[1],&pc_notc));
          PetscCall(PetscFree(sub_ksp));

          PetscCall(DMStagCreateCompatibleDMStag(da,4,1,1,0,&dm_notc));

          // First split within notc is edges
          stencil_notc_edges[0].loc = DMSTAG_BACK_DOWN;
          stencil_notc_edges[0].c = 0;
          stencil_notc_edges[1].loc = DMSTAG_BACK_LEFT;
          stencil_notc_edges[1].c = 0;
          stencil_notc_edges[2].loc = DMSTAG_DOWN_LEFT;
          stencil_notc_edges[2].c = 0;

          // Second split within notc is faces and vertices
          for (PetscInt c=0; c<3; ++c) {
            stencil_notc_notedges[c].loc = DMSTAG_BACK_DOWN_LEFT;
            stencil_notc_notedges[c].c = c;
          }
          stencil_notc_notedges[3].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_notc_notedges[3].c = 3;
          stencil_notc_notedges[4].loc = DMSTAG_LEFT;
          stencil_notc_notedges[4].c = 0;
          stencil_notc_notedges[5].loc = DMSTAG_BACK;
          stencil_notc_notedges[5].c = 0;
          stencil_notc_notedges[6].loc = DMSTAG_DOWN;
          stencil_notc_notedges[6].c = 0;

          PetscCall(DMStagCreateISFromStencils(dm_notc,3,stencil_notc_edges,&is_notc[0]));
          PetscCall(DMStagCreateISFromStencils(dm_notc,7,stencil_notc_notedges,&is_notc[1]));

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(PCFieldSplitSetIS(pc_notc,name_notc[i],is_notc[i]));
          }

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(ISDestroy(&is_notc[i]));
          }
          PetscCall(DMDestroy(&dm_notc));
        }
      }

      /* If the fieldsplit PC wasn't overridden, further split the second split of the second level */
      if (is_fieldsplit) {
        PCType pc_type;


        PetscCall(PCGetType(pc_notc,&pc_type));
        PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
        if (is_fieldsplit) {
          DM            dm_noe;
          KSP           *sub_ksp;

          PetscInt      n_splits;
          DMStagStencil stencil_noe_EP[1], stencil_noe_notEP[6];
          IS            is_noe[2];
          const char    *name_noe[2] = {"EP", "BV"};

          PetscCall(PCSetUp(pc_notc)); // Set up the Fieldsplit PC
          PetscCall(PCFieldSplitGetSubKSP(pc_notc,&n_splits,&sub_ksp));
          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
          PetscCall(KSPGetPC(sub_ksp[1],&pc_noe));
          PetscCall(PetscFree(sub_ksp));

          PetscCall(DMStagCreateCompatibleDMStag(da,4,0,1,0,&dm_noe));

          // First split within notv is 4th dofs on vertices
          stencil_noe_EP[0].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_noe_EP[0].c = 3;

          // Second split within notv is faces and the first 3 dofs on vertices
          for (PetscInt c=0; c<3; ++c) {
            stencil_noe_notEP[c].loc = DMSTAG_BACK_DOWN_LEFT;
            stencil_noe_notEP[c].c = c;
          }
          stencil_noe_notEP[3].loc = DMSTAG_LEFT;
          stencil_noe_notEP[3].c = 0;
          stencil_noe_notEP[4].loc = DMSTAG_BACK;
          stencil_noe_notEP[4].c = 0;
          stencil_noe_notEP[5].loc = DMSTAG_DOWN;
          stencil_noe_notEP[5].c = 0;

          PetscCall(DMStagCreateISFromStencils(dm_noe,1,stencil_noe_EP,&is_noe[0]));
          PetscCall(DMStagCreateISFromStencils(dm_noe,6,stencil_noe_notEP,&is_noe[1]));

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(PCFieldSplitSetIS(pc_noe,name_noe[i],is_noe[i]));
          }

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(ISDestroy(&is_noe[i]));
          }
          PetscCall(DMDestroy(&dm_noe));
        }
      }

      PC            pc_noe_2;

      /* If the fieldsplit PC wasn't overridden, further split the first split of the third level */
      if (is_fieldsplit) {
        PCType pc_type;


        PetscCall(PCGetType(pc_noe,&pc_type));
        PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
        if (is_fieldsplit) {
          DM            dm_notv;
          KSP           *sub_ksp;

          PetscInt      n_splits;
          DMStagStencil stencil_notv_faces[3], stencil_notv_notfaces[3];
          IS            is_notv[2];
          const char    *name_notv[2] = {"B", "V"};

          PetscCall(PCSetUp(pc_noe)); // Set up the Fieldsplit PC
          PetscCall(PCFieldSplitGetSubKSP(pc_noe,&n_splits,&sub_ksp));
          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
          PetscCall(KSPGetPC(sub_ksp[1],&pc_noe_2));
          PetscCall(PetscFree(sub_ksp));

          PetscCall(DMStagCreateCompatibleDMStag(da,3,0,1,0,&dm_notv));

          // First split within notv is faces
          stencil_notv_faces[0].loc = DMSTAG_LEFT;
          stencil_notv_faces[0].c = 0;
          stencil_notv_faces[1].loc = DMSTAG_BACK;
          stencil_notv_faces[1].c = 0;
          stencil_notv_faces[2].loc = DMSTAG_DOWN;
          stencil_notv_faces[2].c = 0;

          // Second split within notv is vertices
          for (PetscInt c=0; c<3; ++c) {
            stencil_notv_notfaces[c].loc = DMSTAG_BACK_DOWN_LEFT;
            stencil_notv_notfaces[c].c = c;
          }

          PetscCall(DMStagCreateISFromStencils(dm_notv,3,stencil_notv_faces,&is_notv[0]));
          PetscCall(DMStagCreateISFromStencils(dm_notv,3,stencil_notv_notfaces,&is_notv[1]));


          for (PetscInt i=0; i<2; ++i) {
            PetscCall(PCFieldSplitSetIS(pc_noe_2,name_notv[i],is_notv[i]));
          }

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(ISDestroy(&is_notv[i]));
          }
          PetscCall(DMDestroy(&dm_notv));
        }
      }
    }

    KSPSetTolerances(dummyKSP,1e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

    PetscCall(TSSolve(dummyts, X));
    PetscCall(TSGetSolveTime(dummyts, & ftime));
    PetscCall(TSGetStepNumber(dummyts, & steps));

    PetscMPIInt size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Saving relaxed solution (run etaplasma %le, not used by the ideal relaxation), MPI size = %d", user->etaplasma, size));
    PetscViewer viewer;
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, user->ic_binary_path, FILE_MODE_WRITE, & viewer));
    PetscCall(stag_vec_io(user, viewer, X, PETSC_FALSE));
    PetscCall(PetscViewerDestroy(&viewer));

    SaveIntermediateSolution(dummyts, steps, ftime, X, user);

    PetscCall(TSGetSolveTime(dummyts, & ftime));
    PetscCall(TSGetStepNumber(dummyts, & steps));
    PetscCall(TSGetConvergedReason(dummyts, & reason));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Setting Initial conditions for B and V: %s at time %g after %d steps\n", TSConvergedReasons[reason], (double) ftime, steps));
    //DumpSolution(ts, 50, X, user);
    SaveIntermediateSolution(dummyts, steps, ftime, X, user);
    PetscCall(MatDestroy( & J));
    PetscCall(TSDestroy( & dummyts));
  }

  if ( (user -> ictype == 9 || user -> ictype == 15) && (user -> itime == 0.0) && 1) {
    SNES dummysnes;
    TS dummyts;
    KSP dummyKSP;
    PC dummypc;
    Mat J;
    PetscInt steps;
    PetscReal ftime;
    TSConvergedReason reason;

    PetscCall(TSCreate(PETSC_COMM_WORLD, & dummyts));
    PetscCall(TSSetDM(dummyts, user->coorda));
    //TSMonitorSet(dummyts, Monitor, user, NULL); /* Set optional user-defined monitoring routine */
    PetscCall(TSSetProblemType(dummyts, TS_NONLINEAR));

    PetscCall(TSGetSNES(dummyts, & dummysnes));
    //SNESSetType(dummysnes, SNESKSPONLY);
    PetscCall(TSSetType(dummyts, TSBEULER));
    PetscCall(DMCreateMatrix(da, & J));
    /* Use coloring to compute finite difference J efficiently */
    PetscCall(SNESSetJacobian(dummysnes, J, J, SNESComputeJacobianDefaultColor, PETSC_NULLPTR));
    PetscCall(TSSetIFunction(dummyts, NULL, FormIFunction_InitializeEP_halo, user));

    PetscCall(SNESSetUseMatrixFree(dummysnes,PETSC_TRUE,PETSC_FALSE));
    PetscCall(SNESSetOptionsPrefix(dummysnes, "dummySNES_"));
    PetscCall(SNESSetFromOptions(dummysnes));

    PetscCall(TSSetTime(dummyts, 0.0));
    PetscCall(TSSetMaxTime(dummyts, 1e-2));
    PetscCall(TSSetExactFinalTime(dummyts, TS_EXACTFINALTIME_STEPOVER));
    PetscCall(TSSetTimeStep(dummyts, 1e-2));
    PetscCall(TSSetSolution(dummyts, X));

    PetscCall(SNESGetKSP(dummysnes, & dummyKSP));
    //TSGetKSP(dummyts, & dummyKSP);
    PetscCall(KSPGetPC(dummyKSP, & dummypc));

    PetscCall(PCSetType(dummypc,PCBJACOBI));
    PetscCall(KSPSetOptionsPrefix(dummyKSP, "dummyKSP_"));

    /* Logic below modifies the PC directly, so this is the last chance to change the solver from the command line */
    PetscCall(KSPSetFromOptions(dummyKSP));

    PetscBool is_fieldsplit;
    {/* first level -> split ni from the rest : {ni}, {V Phi tau B}
        second level -> split tau from {V Phi B} : {ni}, {{tau},{V Phi B}}
        third level -> split V from {Phi B} : {ni}, {{tau},{{Phi B}, {V}}}
        fourth level -> split Phi from B : {ni}, {{tau},{{{Phi}, {B}}, {V}}}
        */
      IS            is[2];
      DMStagStencil stencil0[1], stencil1[10];
      PC            pc_notc, pc_noe;

      const char *name[2] = {"ni", "TEBV"};

      // First split is cells
      stencil0[0].loc = DMSTAG_ELEMENT;
      stencil0[0].c = 0;

      // Second split is the rest
      for (PetscInt c=0; c<4; ++c) {
        stencil1[c].loc = DMSTAG_BACK_DOWN_LEFT;
        stencil1[c].c = c;
      }
      stencil1[4].loc = DMSTAG_LEFT;
      stencil1[4].c = 0;
      stencil1[5].loc = DMSTAG_BACK;
      stencil1[5].c = 0;
      stencil1[6].loc = DMSTAG_DOWN;
      stencil1[6].c = 0;
      stencil1[7].loc = DMSTAG_BACK_DOWN;
      stencil1[7].c = 0;
      stencil1[8].loc = DMSTAG_BACK_LEFT;
      stencil1[8].c = 0;
      stencil1[9].loc = DMSTAG_DOWN_LEFT;
      stencil1[9].c = 0;

      PetscCall(DMStagCreateISFromStencils(da,1,stencil0,&is[0]));
      PetscCall(DMStagCreateISFromStencils(da,10,stencil1,&is[1]));

      for (PetscInt i=0; i<2; ++i) {
        PetscCall(PCFieldSplitSetIS(dummypc,name[i],is[i]));
      }

      for (PetscInt i=0; i<2; ++i) {
        PetscCall(ISDestroy(&is[i]));
      }

      /* If the fieldsplit PC wasn't overridden, further split the second split */
      {
        PCType pc_type;


        PetscCall(KSPGetPC(dummyKSP, &dummypc));
        PetscCall(PCGetType(dummypc,&pc_type));
        PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
        if (is_fieldsplit) {
          DM            dm_notc;
          KSP           *sub_ksp;

          PetscInt      n_splits;
          DMStagStencil stencil_notc_edges[3], stencil_notc_notedges[7];
          IS            is_notc[2];
          const char    *name_notc[2] = {"tau","EBV"};

          PetscCall(PCSetUp(dummypc)); // Set up the Fieldsplit PC
          PetscCall(PCFieldSplitGetSubKSP(dummypc,&n_splits,&sub_ksp));
          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
          PetscCall(KSPGetPC(sub_ksp[1],&pc_notc));
          PetscCall(PetscFree(sub_ksp));

          PetscCall(DMStagCreateCompatibleDMStag(da,4,1,1,0,&dm_notc));

          // First split within notc is edges
          stencil_notc_edges[0].loc = DMSTAG_BACK_DOWN;
          stencil_notc_edges[0].c = 0;
          stencil_notc_edges[1].loc = DMSTAG_BACK_LEFT;
          stencil_notc_edges[1].c = 0;
          stencil_notc_edges[2].loc = DMSTAG_DOWN_LEFT;
          stencil_notc_edges[2].c = 0;

          // Second split within notc is faces and vertices
          for (PetscInt c=0; c<3; ++c) {
            stencil_notc_notedges[c].loc = DMSTAG_BACK_DOWN_LEFT;
            stencil_notc_notedges[c].c = c;
          }
          stencil_notc_notedges[3].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_notc_notedges[3].c = 3;
          stencil_notc_notedges[4].loc = DMSTAG_LEFT;
          stencil_notc_notedges[4].c = 0;
          stencil_notc_notedges[5].loc = DMSTAG_BACK;
          stencil_notc_notedges[5].c = 0;
          stencil_notc_notedges[6].loc = DMSTAG_DOWN;
          stencil_notc_notedges[6].c = 0;

          PetscCall(DMStagCreateISFromStencils(dm_notc,3,stencil_notc_edges,&is_notc[0]));
          PetscCall(DMStagCreateISFromStencils(dm_notc,7,stencil_notc_notedges,&is_notc[1]));

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(PCFieldSplitSetIS(pc_notc,name_notc[i],is_notc[i]));
          }

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(ISDestroy(&is_notc[i]));
          }
          PetscCall(DMDestroy(&dm_notc));
        }
      }

      /* If the fieldsplit PC wasn't overridden, further split the second split of the second level */
      if (is_fieldsplit) {
        PCType pc_type;


        PetscCall(PCGetType(pc_notc,&pc_type));
        PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
        if (is_fieldsplit) {
          DM            dm_noe;
          KSP           *sub_ksp;

          PetscInt      n_splits;
          DMStagStencil stencil_noe_EP[1], stencil_noe_notEP[6];
          IS            is_noe[2];
          const char    *name_noe[2] = {"EP", "BV"};

          PetscCall(PCSetUp(pc_notc)); // Set up the Fieldsplit PC
          PetscCall(PCFieldSplitGetSubKSP(pc_notc,&n_splits,&sub_ksp));
          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
          PetscCall(KSPGetPC(sub_ksp[1],&pc_noe));
          PetscCall(PetscFree(sub_ksp));

          PetscCall(DMStagCreateCompatibleDMStag(da,4,0,1,0,&dm_noe));

          // First split within notv is 4th dofs on vertices
          stencil_noe_EP[0].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_noe_EP[0].c = 3;

          // Second split within notv is faces and the first 3 dofs on vertices
          for (PetscInt c=0; c<3; ++c) {
            stencil_noe_notEP[c].loc = DMSTAG_BACK_DOWN_LEFT;
            stencil_noe_notEP[c].c = c;
          }
          stencil_noe_notEP[3].loc = DMSTAG_LEFT;
          stencil_noe_notEP[3].c = 0;
          stencil_noe_notEP[4].loc = DMSTAG_BACK;
          stencil_noe_notEP[4].c = 0;
          stencil_noe_notEP[5].loc = DMSTAG_DOWN;
          stencil_noe_notEP[5].c = 0;

          PetscCall(DMStagCreateISFromStencils(dm_noe,1,stencil_noe_EP,&is_noe[0]));
          PetscCall(DMStagCreateISFromStencils(dm_noe,6,stencil_noe_notEP,&is_noe[1]));

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(PCFieldSplitSetIS(pc_noe,name_noe[i],is_noe[i]));
          }

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(ISDestroy(&is_noe[i]));
          }
          PetscCall(DMDestroy(&dm_noe));
        }
      }

      PC            pc_noe_2;

      /* If the fieldsplit PC wasn't overridden, further split the first split of the third level */
      if (is_fieldsplit) {
        PCType pc_type;


        PetscCall(PCGetType(pc_noe,&pc_type));
        PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
        if (is_fieldsplit) {
          DM            dm_notv;
          KSP           *sub_ksp;

          PetscInt      n_splits;
          DMStagStencil stencil_notv_faces[3], stencil_notv_notfaces[3];
          IS            is_notv[2];
          const char    *name_notv[2] = {"B", "V"};

          PetscCall(PCSetUp(pc_noe)); // Set up the Fieldsplit PC
          PetscCall(PCFieldSplitGetSubKSP(pc_noe,&n_splits,&sub_ksp));
          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
          PetscCall(KSPGetPC(sub_ksp[1],&pc_noe_2));
          PetscCall(PetscFree(sub_ksp));

          PetscCall(DMStagCreateCompatibleDMStag(da,3,0,1,0,&dm_notv));

          // First split within notv is faces
          stencil_notv_faces[0].loc = DMSTAG_LEFT;
          stencil_notv_faces[0].c = 0;
          stencil_notv_faces[1].loc = DMSTAG_BACK;
          stencil_notv_faces[1].c = 0;
          stencil_notv_faces[2].loc = DMSTAG_DOWN;
          stencil_notv_faces[2].c = 0;

          // Second split within notv is vertices
          for (PetscInt c=0; c<3; ++c) {
            stencil_notv_notfaces[c].loc = DMSTAG_BACK_DOWN_LEFT;
            stencil_notv_notfaces[c].c = c;
          }

          PetscCall(DMStagCreateISFromStencils(dm_notv,3,stencil_notv_faces,&is_notv[0]));
          PetscCall(DMStagCreateISFromStencils(dm_notv,3,stencil_notv_notfaces,&is_notv[1]));


          for (PetscInt i=0; i<2; ++i) {
            PetscCall(PCFieldSplitSetIS(pc_noe_2,name_notv[i],is_notv[i]));
          }

          for (PetscInt i=0; i<2; ++i) {
            PetscCall(ISDestroy(&is_notv[i]));
          }
          PetscCall(DMDestroy(&dm_notv));
        }
      }
    }

    KSPSetTolerances(dummyKSP,1e-8,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

    PetscCall(TSSolve(dummyts, X));


#line 20977


    PetscCall(TSGetSolveTime(dummyts, & ftime));
    PetscCall(TSGetStepNumber(dummyts, & steps));
    PetscCall(TSGetConvergedReason(dummyts, & reason));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Setting Initial conditions for EP and tau: %s at time %g after %d steps\n", TSConvergedReasons[reason], (double) ftime, steps));
    //DumpSolution(ts, 50, X, user);
    PetscCall(MatDestroy( & J));
    PetscCall(TSDestroy( & dummyts));
  }

  PetscCall(DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord));
  if (user -> debug) {
    /*This print is just for debugging*/
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Initial solution vector\n"));
    PetscCall(VecView(X, PETSC_VIEWER_STDOUT_WORLD));
  }

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

#line 21819

/*
   Rank-independent restart I/O for the DMStag solution vector.

   A raw VecView/VecLoad on a DMStag global vector writes the data in PETSc
   parallel ordering (each rank's owned block concatenated in rank order),
   which depends on the number of MPI ranks and their domain decomposition.
   Restarting on a different rank count then scatters the entries into the
   wrong grid locations -> garbage state -> crash. DMStag provides no
   natural-ordering hook (unlike DMDA).

   Fix: split the stag vector into one single-component DMDA vector per
   (storage location, component). DMDA's VecView/VecLoad permute to a
   rank-independent natural ordering automatically, so the on-disk layout is
   the same regardless of decomposition. On load we push the DMDA values back
   into the stag global vector over the local owned corners (the split DMDA
   shares the stag ownership, so indices align 1:1 including the extra
   boundary points on the last rank).
*/

PetscErrorCode stag_vec_io(User *user, PetscViewer viewer,
                           Vec X, PetscBool load)
{
    DM       stagdm;
    PetscInt dof[4];

    PetscFunctionBeginUser;

    PetscCall(TSGetDM(user->ts, &stagdm));
    PetscCall(DMStagGetDOF(stagdm, &dof[0], &dof[1],
                                      &dof[2], &dof[3]));

    const DMStagStencilLocation locs[8] = {
        DMSTAG_BACK_DOWN_LEFT,
        DMSTAG_BACK_DOWN, DMSTAG_BACK_LEFT, DMSTAG_DOWN_LEFT,
        DMSTAG_LEFT, DMSTAG_DOWN, DMSTAG_BACK,
        DMSTAG_ELEMENT
    };

    const PetscInt ndof[8] = {
        dof[0],
        dof[1], dof[1], dof[1],
        dof[2], dof[2], dof[2],
        dof[3]
    };

    Vec Xloc = NULL;

    if (load) {
        /*
         * The file contains every stratum/component, so rebuild X.
         */
        PetscCall(VecZeroEntries(X));
        PetscCall(DMGetLocalVector(stagdm, &Xloc));
    }

    for (PetscInt s = 0; s < 8; ++s) {
        for (PetscInt c = 0; c < ndof[s]; ++c) {
            DM  da;
            Vec davec;
            char name[64];

            PetscCall(PetscSNPrintf(name, sizeof(name),
                                    "stag_%d_%d",
                                    (int)locs[s], (int)c));

            PetscCall(DMStagVecSplitToDMDA(stagdm, X, locs[s], c,
                                           &da, &davec));
            PetscCall(PetscObjectSetName((PetscObject)davec, name));

            if (load) {
                PetscInt slot;
                PetscInt xs, ys, zs, xm, ym, zm;
                PetscScalar ****stagarr;
                PetscScalar ***daarr;

                PetscCall(VecLoad(davec, viewer));

                /*
                 * Xloc is the ghosted DMStag vector.
                 */
                PetscCall(VecZeroEntries(Xloc));

                PetscCall(DMStagGetLocationSlot(stagdm, locs[s],
                                                c, &slot));

                PetscCall(DMStagVecGetArray(stagdm, Xloc, &stagarr));
                PetscCall(DMDAVecGetArrayRead(da, davec, &daarr));
                PetscCall(DMDAGetCorners(da, &xs, &ys, &zs,
                                             &xm, &ym, &zm));

                for (PetscInt k = zs; k < zs + zm; ++k)
                    for (PetscInt j = ys; j < ys + ym; ++j)
                        for (PetscInt i = xs; i < xs + xm; ++i)
                            stagarr[k][j][i][slot] = daarr[k][j][i];

                PetscCall(DMDAVecRestoreArrayRead(da, davec, &daarr));
                PetscCall(DMStagVecRestoreArray(stagdm, Xloc, &stagarr));

                /*
                 * Add only the populated owned entries into X.
                 * The other entries of Xloc are zero.
                 */
                PetscCall(DMLocalToGlobalBegin(stagdm, Xloc,
                                               ADD_VALUES, X));
                PetscCall(DMLocalToGlobalEnd(stagdm, Xloc,
                                             ADD_VALUES, X));
            } else {
                PetscCall(VecView(davec, viewer));
            }

            PetscCall(VecDestroy(&davec));
            PetscCall(DMDestroy(&da));
        }
    }

    if (load)
        PetscCall(DMRestoreLocalVector(stagdm, &Xloc));

    PetscFunctionReturn(PETSC_SUCCESS);
}
