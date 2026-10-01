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
//
// ATTIC -- NOT BUILT.
//
// This file is not built and is not part of any CMake target. It holds
// functions that had zero callers anywhere in src/, tests/, or the built
// objects (verified via nm) as of the commit that quarantined ts_functions.c
// (0823f34), re-verified for geometry.c in this commit. They were moved out
// of geometry.c verbatim (byte-identical bodies) to shrink that file while
// keeping the code in-tree and greppable for physics review.
//
// Retained for physics review only; scheduled for deletion in a later,
// separate commit. Known defects: see docs/mhd/physics/open-questions.md.
// Do NOT re-enable (re-add to a CMake target / re-wire callers) without
// first regenerating the regression baselines under tests/regression/.
//

#include <mfd_config.h>
#include <ts_functions.h>
#include <monitor_functions.h>
#include <geometry.h>
#include <mass_matrix_coefficients.h>
#include <mimetic_operators.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmpatch.h>
#include <petscsf.h>

#include "mhd.h"

#define PETSC_NULL_VEC PETSC_NULLPTR

PetscErrorCode EBoundaryAdjusters(TS ts, Mat LM, Mat Offset, void * ptr) {
  User * user = (User * ) ptr;
  DM da;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  PetscInt N[3], er, ephi, ez;

  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col[1];
        PetscScalar valJ[1];
        PetscInt nEntries = 1;

        /* B field part */

        /* Equation on left face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = LEFT;
        row.c = 0;

        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = LEFT;
        col[0].c = 0;
        valJ[0] = 1.0;

        DMStagMatSetValuesStencil(da, LM, 1, & row, nEntries, col, valJ, INSERT_VALUES);

        /* Equation on down face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN;
        row.c = 0;

        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = DOWN;
        col[0].c = 0;
        valJ[0] = 1.0;

        DMStagMatSetValuesStencil(da, LM, 1, & row, nEntries, col, valJ, INSERT_VALUES);

        /* Equation on back face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK;
        row.c = 0;

        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = BACK;
        col[0].c = 0;
        valJ[0] = 1.0;

        DMStagMatSetValuesStencil(da, LM, 1, & row, nEntries, col, valJ, INSERT_VALUES);

        if (er == N[0] - 1) {
          /* Equation on right boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, LM, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on up boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, LM, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on front boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, LM, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        /* E field part */

        if (er == 0 || ez == 0) {
          /* Equation on the back or left boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if ((user -> phibtype && (er == 0)) || (((er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on the left or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if ((user -> phibtype && (ez == 0)) || (((ephi == 0 || ez == 0)) && !(user -> phibtype))) {
          /* Equation on the left or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (er == N[0] - 1) {
          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP_RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && ez == N[2] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Equation on the front and right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Offset, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if ((user -> phibtype && !(er == 0)) || ((!(er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal down left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, LM, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if ((user -> phibtype && !(ez == 0)) || ((!(ez == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal back down edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, LM, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (!(er == 0 || ez == 0)) {
          /* Equation on internal back left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, LM, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

      }
    }
  }

  MatAssemblyBegin(Offset, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Offset, MAT_FINAL_ASSEMBLY);

  MatAssemblyBegin(LM, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(LM, MAT_FINAL_ASSEMBLY);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "LM Matrix:\n");
    MatView(LM, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "Offset Matrix:\n");
    MatView(Offset, PETSC_VIEWER_STDOUT_WORLD);
  }
  return (0);
}

PetscErrorCode ComputeIsBBoundary(TS ts, IS * isB_boundary, void * ptr) {

  User * user = (User * ) ptr;
  DM da;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  PetscInt N[3], er, ephi, ez, n = 1;

  PetscReal t;
  Vec dummyX;
  Mat Jpre;
  PC dummypc;
  KSP dummyksp, * subksp;
  PetscErrorCode ierr = 0;

  TSGetDM(ts, & da);
  DMCreateGlobalVector(da, & dummyX);
  FormInitialSolution(ts, dummyX, user);
  TSGetTime(ts, & t);
  DMCreateMatrix(da, & Jpre);
  MatZeroEntries(Jpre);

  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col[1];
        PetscScalar valJ[1];
        PetscInt nEntries;

        /* The edges are oriented in the directions of unit vectors: e_r, e_phi and e_z */

        /* B field part */

        if (er == 0) {
          /* Equation on left face */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ephi == 0 && !(user -> phibtype)) {
          /* Equation on down face */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ez == 0) {
          /* Equation on back face */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (er == N[0] - 1) {
          /* Equation on right boundary face */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on up boundary face */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on front boundary face */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

      }
    }
  }

  MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Jpre:\n");
    MatView(Jpre, PETSC_VIEWER_STDOUT_WORLD);
  }

  KSPCreate(PETSC_COMM_WORLD, & dummyksp);
  KSPSetFromOptions(dummyksp);
  KSPSetOperators(dummyksp, Jpre, Jpre);
  PetscBarrier((PetscObject)Jpre);
  MatDestroy( & Jpre);
  KSPGetPC(dummyksp, & dummypc);
  PCSetType(dummypc, PCFIELDSPLIT);
  PCFieldSplitSetDetectSaddlePoint(dummypc, PETSC_TRUE);
  PCSetUp(dummypc);
  PCFieldSplitSetSchurFactType(dummypc, PC_FIELDSPLIT_SCHUR_FACT_FULL);
  PCFieldSplitSetSchurPre(dummypc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
  KSPSetUp(dummyksp);
  KSPSetTolerances(dummyksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 outer iteration for the dummy solve
  PCFieldSplitGetSubKSP(dummypc, & n, & subksp);
  KSPSetTolerances(subksp[1], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 inner iteration maximum for the dummy solve
  KSPSolve(dummyksp, dummyX, dummyX);
  PetscBarrier((PetscObject)dummyX);
  VecDestroy( & dummyX);
  IS isdup;
  PCFieldSplitGetISByIndex(dummypc, 0, & isdup);
  ISDuplicate(isdup, isB_boundary);
  PetscBarrier((PetscObject)dummyksp);
  KSPDestroy( & dummyksp);

  return (0);
}

PetscErrorCode ComputeIsCBoundary(TS ts, IS * isC_boundary, void * ptr) {

  User * user = (User * ) ptr;
  DM coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  PetscInt N[3], er, ephi, ez, n = 1;

  PetscReal t;
  Vec dummyX;
  Mat Jpre;
  PC dummypc;
  KSP dummyksp, * subksp;
  PetscErrorCode ierr = 0;

  DMCreateGlobalVector(coordDA, & dummyX);
  VecSet(dummyX,1.0);
  DMCreateMatrix(coordDA, & Jpre);
  MatZeroEntries(Jpre);

  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col[1];
        PetscScalar valJ[1];
        PetscInt nEntries;

        /* The edges are oriented in the directions of unit vectors: e_r, e_phi and e_z */

        /* Element part */

        if (er == 0 || er == N[0]-1 || (ephi == 0 && !user->phibtype) || (ephi == N[1] - 1 && !user->phibtype) || ez == 0 || ez == N[2] - 1) {
          /* Equation on boundary cell */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = ELEMENT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = ELEMENT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(coordDA, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

      }
    }
  }

  MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Jpre:\n");
    MatView(Jpre, PETSC_VIEWER_STDOUT_WORLD);
  }

  KSPCreate(PETSC_COMM_WORLD, & dummyksp);
  KSPSetFromOptions(dummyksp);
  KSPSetOperators(dummyksp, Jpre, Jpre);
  PetscBarrier((PetscObject)Jpre);
  MatDestroy( & Jpre);
  KSPGetPC(dummyksp, & dummypc);
  PCSetType(dummypc, PCFIELDSPLIT);
  PCFieldSplitSetDetectSaddlePoint(dummypc, PETSC_TRUE);
  PCSetUp(dummypc);
  PCFieldSplitSetSchurFactType(dummypc, PC_FIELDSPLIT_SCHUR_FACT_FULL);
  PCFieldSplitSetSchurPre(dummypc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
  KSPSetUp(dummyksp);
  KSPSetTolerances(dummyksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 outer iteration for the dummy solve
  PCFieldSplitGetSubKSP(dummypc, & n, & subksp);
  KSPSetTolerances(subksp[1], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 inner iteration maximum for the dummy solve
  KSPSolve(dummyksp, dummyX, dummyX);
  PetscBarrier((PetscObject)dummyX);
  VecDestroy( & dummyX);
  IS isdup;
  PCFieldSplitGetISByIndex(dummypc, 0, & isdup);
  ISDuplicate(isdup, isC_boundary);
  PetscBarrier((PetscObject)dummyksp);
  KSPDestroy( & dummyksp);

  return (0);
}

PetscErrorCode ReadDataInVec(DM da, Vec F_r2, Vec F_phi2, Vec F_z2, Vec C2, void *ptr)
{
    User           *user = (User*)ptr;
    DM             coordDA = user->coorda;
    DM             daC, dmC;
    DM             dmFr, dmFphi,dmFz, daFr, daFphi,daFz;
    DM             dmEr, dmEphi,dmEz, daEr, daEphi,daEz;
    DM             dmV, daV;

    PetscInt       startr,startphi,startz,nr,nphi,nz;

    Vec            vecC, C, vecFr, vecFphi, vecFz, F_r, F_phi, F_z, vecEr, vecEphi, vecEz, E_r, E_r2, E_phi, E_phi2, E_z, E_z2, vecV, V, V2;

    Vec            C_rLocal, C_phiLocal, C_zLocal, F_rLocal, F_phiLocal, F_zLocal, V_rLocal, V_phiLocal, V_zLocal, xLocal;

    PetscInt       N[3],er,ephi,ez,d;

    PetscInt       icBrp[3],icBphip[3],icBzp[3],icBrm[3],icBphim[3],icBzm[3];
    PetscInt       icErmzm[3],icErmzp[3],icErpzm[3],icErpzp[3];
    PetscInt       icEphimzm[3],icEphipzm[3],icEphimzp[3],icEphipzp[3];
    PetscInt       icErmphim[3],icErpphim[3],icErmphip[3],icErpphip[3];
    PetscInt       icrmphimzm[3],icrpphimzm[3],icrmphipzm[3],icrpphipzm[3];
    PetscInt       icrmphimzp[3],icrpphimzp[3],icrmphipzp[3],icrpphipzp[3];
    PetscInt       icp[3];


    PetscInt        ivBrp,ivBphip,ivBzp,ivBrm,ivBphim,ivBzm;

    PetscInt          ivErmzm,ivErmzp,ivErpzm,ivErpzp;
    PetscInt        ivEphimzm,ivEphipzm,ivEphimzp,ivEphipzp;
    PetscInt        ivErmphim,ivErpphim,ivErmphip,ivErpphip;
    DM              dmCoord,dmCoorda;
    PetscScalar       ****arrX,****arrCr,****arrCphi,****arrCz,****arrFr,****arrFphi,****arrFz,****arrVr,****arrVphi,****arrVz;


    DMStagCreateCompatibleDMStag(da,0,0,1,0,&dmFr); /* 1 dof per face */
    DMSetUp(dmFr);
    DMStagSetUniformCoordinatesExplicit(dmFr,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmFr,&F_r2);

    DMStagCreateCompatibleDMStag(da,0,0,1,0,&dmFphi); /* 1 dof per face */
    DMSetUp(dmFphi);
    DMStagSetUniformCoordinatesExplicit(dmFphi,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmFphi,&F_phi2);

    DMStagCreateCompatibleDMStag(da,0,0,1,0,&dmFz); /* 1 dof per face */
    DMSetUp(dmFz);
    DMStagSetUniformCoordinatesExplicit(dmFz,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmFz,&F_z2);

    PetscPrintf(PETSC_COMM_WORLD,"Before creating compatible DMStag for E_r\n");

    DMStagCreateCompatibleDMStag(da,0,1,0,0,&dmEr); /* 1 dof1 per edge */
    DMSetUp(dmEr);
    DMStagSetUniformCoordinatesExplicit(dmEr,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmEr,&E_r2);

    DMStagCreateCompatibleDMStag(da,0,1,0,0,&dmEphi); /* 1 dof per edge */
    DMSetUp(dmEphi);
    DMStagSetUniformCoordinatesExplicit(dmEphi,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmEphi,&E_phi2);

    DMStagCreateCompatibleDMStag(da,0,1,0,0,&dmEz); /* 1 dof per edge */
    DMSetUp(dmEz);
    DMStagSetUniformCoordinatesExplicit(dmEz,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmEz,&E_z2);

    DMStagCreateCompatibleDMStag(da,0,0,0,1,&dmC); /* 1 dof per cell */
    DMSetUp(dmC);
    DMStagSetUniformCoordinatesExplicit(dmC,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmC,&C2);
    //STOPPED HERE


    PetscPrintf(PETSC_COMM_WORLD,"Before copying E_r coordinates \n");
    DMStagGetCorners(dmEr,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmEr,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK_DOWN;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK_DOWN;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK_DOWN;    from[2].c = 2;

                DMStagVecSetValuesStencil(dmEr,E_r2,3,from,valFrom,INSERT_VALUES);
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK_UP;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK_UP;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK_UP;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmEr,E_r2,1,from,valFrom,INSERT_VALUES);
                }
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT_DOWN;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT_DOWN;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT_DOWN;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmEr,E_r2,3,from,valFrom,INSERT_VALUES);
                }
                if(ephi == N[1]-1 && ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT_UP;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT_UP;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT_UP;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmEr,E_r2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(E_r2);
    VecAssemblyEnd(E_r2);

    DMStagVecSplitToDMDA(dmEr,E_r2,BACK_DOWN,-1,&daEr,&vecEr); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecEr,"rEdge_center_coordinates");


    PetscPrintf(PETSC_COMM_WORLD,"Before copying E_phi coordinates\n");
    DMStagGetCorners(dmEphi,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmEphi,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK_LEFT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK_LEFT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK_LEFT;    from[2].c = 2;

                DMStagVecSetValuesStencil(dmEphi,E_phi2,3,from,valFrom,INSERT_VALUES);
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK_RIGHT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK_RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK_RIGHT;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmEphi,E_phi2,3,from,valFrom,INSERT_VALUES);
                }
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT_LEFT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT_LEFT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT_LEFT;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmEphi,E_phi2,3,from,valFrom,INSERT_VALUES);
                }
                if(er == N[0]-1 && ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT_RIGHT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT_RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT_RIGHT;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmEphi,E_phi2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(E_phi2);
    VecAssemblyEnd(E_phi2);

    PetscPrintf(PETSC_COMM_WORLD,"Before calling DMStagVecSplitToDMDA for E_phi coordinates\n");
    DMStagVecSplitToDMDA(dmEphi,E_phi2,BACK_LEFT,-1,&daEphi,&vecEphi); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecEphi,"phiEdge_center_coordinates");


    PetscPrintf(PETSC_COMM_WORLD,"Before copying E_z coordinates\n");
    DMStagGetCorners(dmEz,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmEz,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN_LEFT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = DOWN_LEFT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = DOWN_LEFT;    from[2].c = 2;

                DMStagVecSetValuesStencil(dmEz,E_z2,3,from,valFrom,INSERT_VALUES);
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN_RIGHT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = DOWN_RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = DOWN_RIGHT;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmEz,E_z2,3,from,valFrom,INSERT_VALUES);
                }
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP_LEFT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = UP_LEFT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = UP_LEFT;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmEz,E_z2,3,from,valFrom,INSERT_VALUES);
                }
                if(er == N[0]-1 && ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP_RIGHT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = UP_RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = UP_RIGHT;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmEz,E_z2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(E_z2);
    VecAssemblyEnd(E_z2);

    PetscPrintf(PETSC_COMM_WORLD,"Before calling DMStagVecSplitToDMDA for E_z coordinates\n");
    DMStagVecSplitToDMDA(dmEz,E_z2,DOWN_LEFT,-1,&daEz,&vecEz); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecEz,"zEdge_center_coordinates");


    PetscPrintf(PETSC_COMM_WORLD,"Before copying F_r coordinates\n");
    DMStagGetCorners(dmFr,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmFr,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = LEFT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = LEFT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = LEFT;    from[2].c = 2;

                DMStagVecSetValuesStencil(dmFr,F_r2,3,from,valFrom,INSERT_VALUES);
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = RIGHT;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = RIGHT;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmFr,F_r2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(F_r2);
    VecAssemblyEnd(F_r2);

    DMStagVecSplitToDMDA(dmFr,F_r2,LEFT,-1,&daFr,&vecFr); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecFr,"rFace_center_coordinates");


    PetscPrintf(PETSC_COMM_WORLD,"Before copying F_phi coordinates\n");
    DMStagGetCorners(dmFphi,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmFphi,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = DOWN;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = DOWN;    from[2].c = 2;

                DMStagVecSetValuesStencil(dmFphi,F_phi2,3,from,valFrom,INSERT_VALUES);
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = UP;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = UP;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmFphi,F_phi2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(F_phi2);
    VecAssemblyEnd(F_phi2);

    DMStagVecSplitToDMDA(dmFphi,F_phi2,DOWN,-1,&daFphi,&vecFphi); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecFphi,"phiFace_center_coordinates");


    PetscPrintf(PETSC_COMM_WORLD,"Before copying F_z coordinates\n");
    DMStagGetCorners(dmFz,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmFz,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK;    from[2].c = 2;

                DMStagVecSetValuesStencil(dmFz,F_z2,3,from,valFrom,INSERT_VALUES);
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT;    from[2].c = 2;

                    DMStagVecSetValuesStencil(dmFz,F_z2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(F_z2);
    VecAssemblyEnd(F_z2);

    DMStagVecSplitToDMDA(dmFz,F_z2,BACK,-1,&daFz,&vecFz); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecFz,"zFace_center_coordinates");


    PetscPrintf(PETSC_COMM_WORLD,"Before copying cell coordinates\n");

    DMStagGetCorners(dmC,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = ELEMENT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = ELEMENT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = ELEMENT;    from[2].c = 2;

                DMStagVecSetValuesStencil(dmC,C2,3,from,valFrom,INSERT_VALUES);
            }
        }
    }
    VecAssemblyBegin(C2);
    VecAssemblyEnd(C2);

    DMStagVecSplitToDMDA(dmC,C2,ELEMENT,-1,&daC,&vecC); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecC,"Cell_center_coordinates");



    //PetscBarrier((PetscObject)viewerD);
    //PetscViewerDestroy(&viewerD);

    DMDestroy(&dmFr);
    DMDestroy(&dmFphi);
    DMDestroy(&dmFz);
    DMDestroy(&dmEr);
    DMDestroy(&dmEphi);
    DMDestroy(&dmEz);
    DMDestroy(&dmC);

    VecDestroy(&C2);
    VecDestroy(&E_r2);
    VecDestroy(&E_phi2);
    VecDestroy(&E_z2);
    VecDestroy(&F_r2);
    VecDestroy(&F_phi2);
    VecDestroy(&F_z2);

    return(0);
}

PetscErrorCode CellToVertexProjectionVector(TS ts, Vec C, Vec V, void *ptr)
{
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("CellToVertexProjectionVector",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;
    Vec            CLocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    VecZeroEntries(V);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
    DMGetLocalVector(da,&CLocal);
    DMGlobalToLocal(da,C,INSERT_VALUES,CLocal);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[24], to[3];
          PetscScalar valFrom[24], valTo[3];

          if(er>0 && ez>0 && (ephi>0 || user->phibtype)){
          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = ELEMENT;
          from[0].c = 0;
          from[1].i = er-1;
          from[1].j = ephi;
          from[1].k = ez;
          from[1].loc = ELEMENT;
          from[1].c = 0;
          from[2].i = er;
          from[2].j = ephi-1;
          from[2].k = ez;
          from[2].loc = ELEMENT;
          from[2].c = 0;
          from[3].i = er;
          from[3].j = ephi;
          from[3].k = ez-1;
          from[3].loc = ELEMENT;
          from[3].c = 0;
          from[4].i = er-1;
          from[4].j = ephi-1;
          from[4].k = ez;
          from[4].loc = ELEMENT;
          from[4].c = 0;
          from[5].i = er-1;
          from[5].j = ephi;
          from[5].k = ez-1;
          from[5].loc = ELEMENT;
          from[5].c = 0;
          from[6].i = er;
          from[6].j = ephi-1;
          from[6].k = ez-1;
          from[6].loc = ELEMENT;
          from[6].c = 0;
          from[7].i = er-1;
          from[7].j = ephi-1;
          from[7].k = ez-1;
          from[7].loc = ELEMENT;
          from[7].c = 0;
          from[8].i = er;
          from[8].j = ephi;
          from[8].k = ez;
          from[8].loc = ELEMENT;
          from[8].c = 1;
          from[9].i = er-1;
          from[9].j = ephi;
          from[9].k = ez;
          from[9].loc = ELEMENT;
          from[9].c = 1;
          from[10].i = er;
          from[10].j = ephi-1;
          from[10].k = ez;
          from[10].loc = ELEMENT;
          from[10].c = 1;
          from[11].i = er;
          from[11].j = ephi;
          from[11].k = ez-1;
          from[11].loc = ELEMENT;
          from[11].c = 1;
          from[12].i = er-1;
          from[12].j = ephi-1;
          from[12].k = ez;
          from[12].loc = ELEMENT;
          from[12].c = 1;
          from[13].i = er-1;
          from[13].j = ephi;
          from[13].k = ez-1;
          from[13].loc = ELEMENT;
          from[13].c = 1;
          from[14].i = er;
          from[14].j = ephi-1;
          from[14].k = ez-1;
          from[14].loc = ELEMENT;
          from[14].c = 1;
          from[15].i = er-1;
          from[15].j = ephi-1;
          from[15].k = ez-1;
          from[15].loc = ELEMENT;
          from[15].c = 1;
          from[16].i = er;
          from[16].j = ephi;
          from[16].k = ez;
          from[16].loc = ELEMENT;
          from[16].c = 2;
          from[17].i = er-1;
          from[17].j = ephi;
          from[17].k = ez;
          from[17].loc = ELEMENT;
          from[17].c = 2;
          from[18].i = er;
          from[18].j = ephi-1;
          from[18].k = ez;
          from[18].loc = ELEMENT;
          from[18].c = 2;
          from[19].i = er;
          from[19].j = ephi;
          from[19].k = ez-1;
          from[19].loc = ELEMENT;
          from[19].c = 2;
          from[20].i = er-1;
          from[20].j = ephi-1;
          from[20].k = ez;
          from[20].loc = ELEMENT;
          from[20].c = 2;
          from[21].i = er-1;
          from[21].j = ephi;
          from[21].k = ez-1;
          from[21].loc = ELEMENT;
          from[21].c = 2;
          from[22].i = er;
          from[22].j = ephi-1;
          from[22].k = ez-1;
          from[22].loc = ELEMENT;
          from[22].c = 2;
          from[23].i = er-1;
          from[23].j = ephi-1;
          from[23].k = ez-1;
          from[23].loc = ELEMENT;
          from[23].c = 2;
          DMStagVecGetValuesStencil(da, CLocal, 24, from, valFrom);
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 0;
          valTo[0] = 0.125 * (valFrom[0] + valFrom[1] + valFrom[2] + valFrom[3] + valFrom[4] + valFrom[5] + valFrom[6] + valFrom[7]);
          to[1].i = er;
          to[1].j = ephi;
          to[1].k = ez;
          to[1].loc = BACK_DOWN_LEFT;
          to[1].c = 1;
          valTo[1] = 0.125 * (valFrom[8] + valFrom[9] + valFrom[10] + valFrom[11] + valFrom[12] + valFrom[13] + valFrom[14] + valFrom[15]);
          to[2].i = er;
          to[2].j = ephi;
          to[2].k = ez;
          to[2].loc = BACK_DOWN_LEFT;
          to[2].c = 2;
          valTo[2] = 0.125 * (valFrom[16] + valFrom[17] + valFrom[18] + valFrom[19] + valFrom[20] + valFrom[21] + valFrom[22] + valFrom[23]);
          DMStagVecSetValuesStencil(da, V, 3, to, valTo, INSERT_VALUES);
          }
        }
      }
    }
    VecAssemblyBegin(V);
    VecAssemblyEnd(V);
    DMRestoreLocalVector(da,&CLocal);

    PetscLogEventEnd(USER_EVENT,0,0,0,0);

    return(0);
}

PetscErrorCode VertexToCellReconstruction(TS ts, Vec V, Vec C, void *ptr)
{
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("VertexToCellReconstruction",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;
    Vec            VLocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    VecZeroEntries(C);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
    DMGetLocalVector(da,&VLocal);
    DMGlobalToLocal(da,V,INSERT_VALUES,VLocal);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[24], to[3];
          PetscScalar valFrom[24], valTo[3];

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = FRONT_UP_RIGHT;
          from[0].c = 0;
          from[1].i = er;
          from[1].j = ephi;
          from[1].k = ez;
          from[1].loc = FRONT_UP_LEFT;
          from[1].c = 0;
          from[2].i = er;
          from[2].j = ephi;
          from[2].k = ez;
          from[2].loc = FRONT_DOWN_LEFT;
          from[2].c = 0;
          from[3].i = er;
          from[3].j = ephi;
          from[3].k = ez;
          from[3].loc = FRONT_DOWN_RIGHT;
          from[3].c = 0;
          from[4].i = er;
          from[4].j = ephi;
          from[4].k = ez;
          from[4].loc = BACK_UP_RIGHT;
          from[4].c = 0;
          from[5].i = er;
          from[5].j = ephi;
          from[5].k = ez;
          from[5].loc = BACK_UP_LEFT;
          from[5].c = 0;
          from[6].i = er;
          from[6].j = ephi;
          from[6].k = ez;
          from[6].loc = BACK_DOWN_LEFT;
          from[6].c = 0;
          from[7].i = er;
          from[7].j = ephi;
          from[7].k = ez;
          from[7].loc = BACK_DOWN_RIGHT;
          from[7].c = 0;
          from[8].i = er;
          from[8].j = ephi;
          from[8].k = ez;
          from[8].loc = FRONT_UP_RIGHT;
          from[8].c = 1;
          from[9].i = er;
          from[9].j = ephi;
          from[9].k = ez;
          from[9].loc = FRONT_UP_LEFT;
          from[9].c = 1;
          from[10].i = er;
          from[10].j = ephi;
          from[10].k = ez;
          from[10].loc = FRONT_DOWN_LEFT;
          from[10].c = 1;
          from[11].i = er;
          from[11].j = ephi;
          from[11].k = ez;
          from[11].loc = FRONT_DOWN_RIGHT;
          from[11].c = 1;
          from[12].i = er;
          from[12].j = ephi;
          from[12].k = ez;
          from[12].loc = BACK_UP_RIGHT;
          from[12].c = 1;
          from[13].i = er;
          from[13].j = ephi;
          from[13].k = ez;
          from[13].loc = BACK_UP_LEFT;
          from[13].c = 1;
          from[14].i = er;
          from[14].j = ephi;
          from[14].k = ez;
          from[14].loc = BACK_DOWN_LEFT;
          from[14].c = 1;
          from[15].i = er;
          from[15].j = ephi;
          from[15].k = ez;
          from[15].loc = BACK_DOWN_RIGHT;
          from[15].c = 1;
          from[16].i = er;
          from[16].j = ephi;
          from[16].k = ez;
          from[16].loc = FRONT_UP_RIGHT;
          from[16].c = 2;
          from[17].i = er;
          from[17].j = ephi;
          from[17].k = ez;
          from[17].loc = FRONT_UP_LEFT;
          from[17].c = 2;
          from[18].i = er;
          from[18].j = ephi;
          from[18].k = ez;
          from[18].loc = FRONT_DOWN_LEFT;
          from[18].c = 2;
          from[19].i = er;
          from[19].j = ephi;
          from[19].k = ez;
          from[19].loc = FRONT_DOWN_RIGHT;
          from[19].c = 2;
          from[20].i = er;
          from[20].j = ephi;
          from[20].k = ez;
          from[20].loc = BACK_UP_RIGHT;
          from[20].c = 2;
          from[21].i = er;
          from[21].j = ephi;
          from[21].k = ez;
          from[21].loc = BACK_UP_LEFT;
          from[21].c = 2;
          from[22].i = er;
          from[22].j = ephi;
          from[22].k = ez;
          from[22].loc = BACK_DOWN_LEFT;
          from[22].c = 2;
          from[23].i = er;
          from[23].j = ephi;
          from[23].k = ez;
          from[23].loc = BACK_DOWN_RIGHT;
          from[23].c = 2;
          DMStagVecGetValuesStencil(da, VLocal, 24, from, valFrom);
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = ELEMENT;
          to[0].c = 0;
          valTo[0] = 0.125 * (valFrom[0] + valFrom[1] + valFrom[2] + valFrom[3] + valFrom[4] + valFrom[5] + valFrom[6] + valFrom[7]);
          to[1].i = er;
          to[1].j = ephi;
          to[1].k = ez;
          to[1].loc = ELEMENT;
          to[1].c = 1;
          valTo[1] = 0.125 * (valFrom[8] + valFrom[9] + valFrom[10] + valFrom[11] + valFrom[12] + valFrom[13] + valFrom[14] + valFrom[15]);
          to[2].i = er;
          to[2].j = ephi;
          to[2].k = ez;
          to[2].loc = ELEMENT;
          to[2].c = 2;
          valTo[2] = 0.125 * (valFrom[16] + valFrom[17] + valFrom[18] + valFrom[19] + valFrom[20] + valFrom[21] + valFrom[22] + valFrom[23]);
          DMStagVecSetValuesStencil(da, C, 3, to, valTo, INSERT_VALUES);

        }
      }
    }
    VecAssemblyBegin(C);
    VecAssemblyEnd(C);
    DMRestoreLocalVector(da,&VLocal);

    PetscLogEventEnd(USER_EVENT,0,0,0,0);

    return(0);
}

PetscErrorCode VertexToEdgeReconstructionMat(TS ts, Mat M, void *ptr)
{
    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;

    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    MatZeroEntries(M);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil row, col[2];
          PetscScalar valJ[2];
          PetscInt nEntries = 2;

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_UP_LEFT;
          col[1].c = 0;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_LEFT;
          col[0].c = 1;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_LEFT;
          col[1].c = 1;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_RIGHT;
          col[1].c = 0;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_RIGHT;
          col[0].c = 1;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_RIGHT;
          col[1].c = 1;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP_RIGHT;
          col[1].c = 0;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_LEFT;
          col[0].c = 1;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_DOWN_LEFT;
          col[1].c = 1;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_DOWN_RIGHT;
          col[1].c = 0;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN_RIGHT;
          col[0].c = 1;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP_RIGHT;
          col[1].c = 1;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_RIGHT;
          col[0].c = 2;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_UP_RIGHT;
          col[1].c = 2;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_LEFT;
          col[0].c = 2;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_UP_LEFT;
          col[1].c = 2;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN_RIGHT;
          col[0].c = 2;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_RIGHT;
          col[1].c = 2;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN_LEFT;
          col[0].c = 2;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_LEFT;
          col[1].c = 2;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

        }
      }
    }
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    return(0);
}

PetscErrorCode VertexToFaceReconstructionMat(TS ts, Mat M, void *ptr)
{
    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;

    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    MatZeroEntries(M);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil row, col[4];
          PetscScalar valJ[4];
          PetscInt nEntries = 4;

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_RIGHT;
          col[0].c = 1;
          valJ[0] = 0.25;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP_LEFT;
          col[1].c = 1;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez;
          col[2].loc = FRONT_UP_RIGHT;
          col[2].c = 1;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = FRONT_UP_LEFT;
          col[3].c = 1;
          valJ[3] = 0.25;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN_LEFT;
          col[0].c = 0;
          valJ[0] = 0.25;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP_LEFT;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez;
          col[2].loc = FRONT_UP_LEFT;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = FRONT_DOWN_LEFT;
          col[3].c = 0;
          valJ[3] = 0.25;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN_LEFT;
          col[0].c = 1;
          valJ[0] = 0.25;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_DOWN_RIGHT;
          col[1].c = 1;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez;
          col[2].loc = FRONT_DOWN_RIGHT;
          col[2].c = 1;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = FRONT_DOWN_LEFT;
          col[3].c = 1;
          valJ[3] = 0.25;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.25;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_DOWN_RIGHT;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez;
          col[2].loc = FRONT_DOWN_RIGHT;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = FRONT_UP_RIGHT;
          col[3].c = 0;
          valJ[3] = 0.25;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_RIGHT;
          col[0].c = 2;
          valJ[0] = 0.25;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_DOWN_RIGHT;
          col[1].c = 2;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez;
          col[2].loc = BACK_DOWN_LEFT;
          col[2].c = 2;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = BACK_UP_LEFT;
          col[3].c = 2;
          valJ[3] = 0.25;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_RIGHT;
          col[0].c = 2;
          valJ[0] = 0.25;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_RIGHT;
          col[1].c = 2;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez;
          col[2].loc = FRONT_DOWN_LEFT;
          col[2].c = 2;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = FRONT_UP_LEFT;
          col[3].c = 2;
          valJ[3] = 0.25;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

        }
      }
    }
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    return(0);
}

PetscErrorCode EdgeToCellReconstructionMat(TS ts, Mat M, void *ptr)
{
    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;

    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    MatZeroEntries(M);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil row, col[4];
          PetscScalar valJ[4];
          PetscInt nEntries = 4;

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = ELEMENT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP;
          col[0].c = 0;
          valJ[0] = 0.25;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez;
          col[2].loc = BACK_UP;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = BACK_DOWN;
          col[3].c = 0;
          valJ[3] = 0.25;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = ELEMENT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.25;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_LEFT;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez;
          col[2].loc = BACK_RIGHT;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = BACK_LEFT;
          col[3].c = 0;
          valJ[3] = 0.25;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = ELEMENT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.25;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = DOWN_RIGHT;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez;
          col[2].loc = UP_LEFT;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = DOWN_LEFT;
          col[3].c = 0;
          valJ[3] = 0.25;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

        }
      }
    }
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    return(0);
}

PetscErrorCode FaceToCellReconstructionMat(TS ts, Mat M, void *ptr)
{
    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;

    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    MatZeroEntries(M);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil row, col[2];
          PetscScalar valJ[2];
          PetscInt nEntries = 2;

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = ELEMENT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP;
          col[0].c = 0;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = DOWN;
          col[1].c = 0;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = ELEMENT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT;
          col[0].c = 0;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK;
          col[1].c = 0;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = ELEMENT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;

          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = LEFT;
          col[1].c = 0;
          valJ[1] = 0.5;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

        }
      }
    }
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    return(0);
}

PetscErrorCode FaceToVertexProjectionMat(TS ts, Mat M, void *ptr)
{
    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;

    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    MatZeroEntries(M);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil row, col[4];
          PetscScalar valJ[4];
          PetscInt nEntries;

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_LEFT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(er!=0 && ez!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = DOWN;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez-1;
          col[2].loc = DOWN;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er-1;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = DOWN;
          col[3].c = 0;
          valJ[3] = 0.25;

          nEntries = 4;
          }
          else if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = DOWN;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else if(er!=0) {
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = DOWN;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_LEFT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK;
          col[0].c = 0;
          valJ[0] = 0.25;

          if ((user -> phibtype && !(er == 0)) || ((!(er == 0 || ephi == 0)) && !(user -> phibtype))) {
          col[1].i = er-1;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = BACK;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi-1;
          col[2].k = ez;
          col[2].loc = BACK;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er-1;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = BACK;
          col[3].c = 0;
          valJ[3] = 0.25;

          nEntries = 4;
          }
          else if(er!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = BACK;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else {
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = LEFT;
          col[0].c = 0;
          valJ[0] = 0.25;

          if ((user -> phibtype && !(ez == 0)) || ((!(ez == 0 || ephi == 0)) && !(user -> phibtype))) {
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez-1;
          col[1].loc = LEFT;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi-1;
          col[2].k = ez;
          col[2].loc = LEFT;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez-1;
          col[3].loc = LEFT;
          col[3].c = 0;
          valJ[3] = 0.25;

          nEntries = 4;
          }
          else if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = LEFT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = LEFT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          if(er==N[0]-1){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_RIGHT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = DOWN;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_RIGHT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK;
          col[0].c = 0;
          valJ[0] = 0.25;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = BACK;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else {
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = RIGHT;
          col[0].c = 0;
          valJ[0] = 0.25;

          if ((user -> phibtype && !(ez == 0)) || ((!(ez == 0 || ephi == 0)) && !(user -> phibtype))) {
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez-1;
          col[1].loc = RIGHT;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi-1;
          col[2].k = ez;
          col[2].loc = RIGHT;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er;
          col[3].j = ephi;
          col[3].k = ez-1;
          col[3].loc = RIGHT;
          col[3].c = 0;
          valJ[3] = 0.25;

          nEntries = 4;
          }
          else if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = RIGHT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = RIGHT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(ez==N[2]-1){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_LEFT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(er!=0) {
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = DOWN;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = LEFT;
          col[0].c = 0;
          valJ[0] = 0.25;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = LEFT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_LEFT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT;
          col[0].c = 0;
          valJ[0] = 0.25;

          if ((user -> phibtype && !(er == 0)) || ((!(er == 0 || ephi == 0)) && !(user -> phibtype))) {
          col[1].i = er-1;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = FRONT;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi-1;
          col[2].k = ez;
          col[2].loc = FRONT;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er-1;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = FRONT;
          col[3].c = 0;
          valJ[3] = 0.25;

          nEntries = 4;
          }
          else if(er!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = FRONT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else {
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(ephi==N[1]-1 && !user->phibtype){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_LEFT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(er!=0 && ez!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = UP;
          col[1].c = 0;
          valJ[1] = 0.25;

          col[2].i = er;
          col[2].j = ephi;
          col[2].k = ez-1;
          col[2].loc = UP;
          col[2].c = 0;
          valJ[2] = 0.25;

          col[3].i = er-1;
          col[3].j = ephi;
          col[3].k = ez;
          col[3].loc = UP;
          col[3].c = 0;
          valJ[3] = 0.25;

          nEntries = 4;
          }
          else if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = UP;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else if(er!=0) {
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = UP;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = LEFT;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = LEFT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_LEFT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(er!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(er==N[0]-1 && ephi==N[1]-1 && !user->phibtype){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_RIGHT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = UP;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = RIGHT;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = RIGHT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_RIGHT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK;
          col[0].c = 0;
          valJ[0] = 0.25;
          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_LEFT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(er!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = UP;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_LEFT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT;
          col[0].c = 0;
          valJ[0] = 0.25;

          if(er!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = LEFT;
          col[0].c = 0;
          valJ[0] = 0.25;
          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(er==N[0]-1 && ez==N[2]-1){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_RIGHT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN;
          col[0].c = 0;
          valJ[0] = 0.25;
          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = RIGHT;
          col[0].c = 0;
          valJ[0] = 0.25;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = RIGHT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else {
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_RIGHT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT;
          col[0].c = 0;
          valJ[0] = 0.25;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = FRONT;
          col[1].c = 0;
          valJ[1] = 0.25;

          nEntries = 2;
          }
          else {
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(er==N[0]-1 && ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_RIGHT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT;
          col[0].c = 0;
          valJ[0] = 0.25;
          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = RIGHT;
          col[0].c = 0;
          valJ[0] = 0.25;
          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_RIGHT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP;
          col[0].c = 0;
          valJ[0] = 0.25;
          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

        }
      }
    }
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    return(0);
}

PetscErrorCode EdgeToVertexProjection_Original(TS ts, Vec E, Vec V, void *ptr)
{
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("EdgeToVertexProjection_Original",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;
    Vec            ELocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er = 0, ephi = 0, ez = 0, n = 1;

    PetscErrorCode ierr = 0;

    VecZeroEntries(V);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
    DMGetLocalVector(da,&ELocal);
    DMGlobalToLocal(da,E,INSERT_VALUES,ELocal);

    DMStagStencil from[24], to[36];
    PetscScalar valFrom[24], valTo[36];
    PetscInt nEntries, nValues = 0;
    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      if(ez==N[2]-1){continue;}
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {

          if(er==N[0]-1){continue;}

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = BACK_LEFT;
          from[0].c = 0;
          nEntries = 1;
          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
            from[1].i = er;
            from[1].j = ephi-1;
            from[1].k = ez;
            from[1].loc = BACK_LEFT;
            from[1].c = 0;
          }
          else{
            from[1].i = er;
            from[1].j = ephi;
            from[1].k = ez;
            from[1].loc = BACK_LEFT;
            from[1].c = 0;
          }
          nEntries += 1;

          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 1;
          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          else{
            valTo[0] = 0.5 * valFrom[0];
          }
          nValues += 1;


          from[2].i = er;
          from[2].j = ephi;
          from[2].k = ez;
          from[2].loc = DOWN_LEFT;
          from[2].c = 0;
          nEntries += 1;
          if(ez!=0){
            from[3].i = er;
            from[3].j = ephi;
            from[3].k = ez-1;
            from[3].loc = DOWN_LEFT;
            from[3].c = 0;
          }
          else{
            from[3].i = er;
            from[3].j = ephi;
            from[3].k = ez;
            from[3].loc = DOWN_LEFT;
            from[3].c = 0;
          }
          nEntries += 1;

          to[1].i = er;
          to[1].j = ephi;
          to[1].k = ez;
          to[1].loc = BACK_DOWN_LEFT;
          to[1].c = 2;
          if(ez!=0){
            valTo[1] = 0.5 * (valFrom[2] + valFrom[3]);
          }
          else{
            valTo[0] = 0.5 * valFrom[2];
          }
          nValues += 1;


          from[4].i = er;
          from[4].j = ephi;
          from[4].k = ez;
          from[4].loc = BACK_DOWN;
          from[4].c = 0;
          nEntries += 1;
          if(er!=0){
            from[5].i = er-1;
            from[5].j = ephi;
            from[5].k = ez;
            from[5].loc = BACK_DOWN;
            from[5].c = 0;
          }
          else{
            from[5].i = er;
            from[5].j = ephi;
            from[5].k = ez;
            from[5].loc = BACK_DOWN;
            from[5].c = 0;
          }
          nEntries += 1;

          to[2].i = er;
          to[2].j = ephi;
          to[2].k = ez;
          to[2].loc = BACK_DOWN_LEFT;
          to[2].c = 0;
          if(er!=0){
            valTo[2] = 0.5 * (valFrom[5] + valFrom[4]);
          }
          else{
            valTo[2] = 0.5 * valFrom[4];
          }
          nValues += 1;

        }
      }
    }

    if(er==N[0]-1){
      //skip ez == N[2]-1
      //loop over ephi and ez
      from[6].i = er;
      from[6].j = ephi;
      from[6].k = ez;
      from[6].loc = BACK_RIGHT;
      from[6].c = 0;
      nEntries += 1;
      if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
        from[7].i = er;
        from[7].j = ephi-1;
        from[7].k = ez;
        from[7].loc = BACK_RIGHT;
        from[7].c = 0;
      }
      else{
        from[7].i = er;
        from[7].j = ephi;
        from[7].k = ez;
        from[7].loc = BACK_RIGHT;
        from[7].c = 0;
      }
      nEntries += 1;

      to[3].i = er;
      to[3].j = ephi;
      to[3].k = ez;
      to[3].loc = BACK_DOWN_RIGHT;
      to[3].c = 1;
      if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
        valTo[3] = 0.5 * (valFrom[7] + valFrom[6]);
      }
      else{
        valTo[3] = 0.5 * valFrom[6];
      }
      nValues += 1;


      from[8].i = er;
      from[8].j = ephi;
      from[8].k = ez;
      from[8].loc = DOWN_RIGHT;
      from[8].c = 0;
      nEntries += 1;
      if(ez!=0){
        from[9].i = er;
        from[9].j = ephi;
        from[9].k = ez-1;
        from[9].loc = DOWN_RIGHT;
        from[9].c = 0;
      }
      else{
        from[9].i = er;
        from[9].j = ephi;
        from[9].k = ez;
        from[9].loc = DOWN_RIGHT;
        from[9].c = 0;
      }
      nEntries += 1;

      to[4].i = er;
      to[4].j = ephi;
      to[4].k = ez;
      to[4].loc = BACK_DOWN_RIGHT;
      to[4].c = 2;
      if(ez!=0){
        valTo[4] = 0.5 * (valFrom[9] + valFrom[8]);
      }
      else{
        valTo[4] = 0.5 * valFrom[8];
      }
      nValues += 1;


      from[10].i = er;
      from[10].j = ephi;
      from[10].k = ez;
      from[10].loc = BACK_DOWN;
      from[10].c = 0;
      nEntries += 1;

      to[5].i = er;
      to[5].j = ephi;
      to[5].k = ez;
      to[5].loc = BACK_DOWN_RIGHT;
      to[5].c = 0;
      valTo[5] = 0.5 * valFrom[0];
      nValues += 1;
    }



    if(ez==N[2]-1){
      //skip er = N[0]-1
      //loop over er and ephi
      from[0].i = er;
      from[0].j = ephi;
      from[0].k = ez;
      from[0].loc = FRONT_LEFT;
      from[0].c = 0;
      nEntries = 1;
      if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
        from[1].i = er;
        from[1].j = ephi-1;
        from[1].k = ez;
        from[1].loc = FRONT_LEFT;
        from[1].c = 0;
        nEntries = 2;
      }
      DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom);
      to[0].i = er;
      to[0].j = ephi;
      to[0].k = ez;
      to[0].loc = FRONT_DOWN_LEFT;
      to[0].c = 1;
      if(nEntries==1){
        valTo[0] = 0.5 * valFrom[0];
      }
      else{
        valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
      }
      DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES);

      from[0].i = er;
      from[0].j = ephi;
      from[0].k = ez;
      from[0].loc = FRONT_DOWN;
      from[0].c = 0;
      nEntries = 1;
      if(er!=0){
        from[1].i = er-1;
        from[1].j = ephi;
        from[1].k = ez;
        from[1].loc = FRONT_DOWN;
        from[1].c = 0;
        nEntries = 2;
      }
      DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom);
      to[0].i = er;
      to[0].j = ephi;
      to[0].k = ez;
      to[0].loc = FRONT_DOWN_LEFT;
      to[0].c = 0;
      if(nEntries==1){
        valTo[0] = 0.5 * valFrom[0];
      }
      else{
        valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
      }
      DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES);

      from[0].i = er;
      from[0].j = ephi;
      from[0].k = ez;
      from[0].loc = DOWN_LEFT;
      from[0].c = 0;
      nEntries = 1;
      DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom);
      to[0].i = er;
      to[0].j = ephi;
      to[0].k = ez;
      to[0].loc = FRONT_DOWN_LEFT;
      to[0].c = 2;
      if(nEntries==1){
        valTo[0] = 0.5 * valFrom[0];
      }
      else{
        valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
      }
      DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES);
    }


    if(er==N[0]-1 && ez==N[2]-1){
      //loop over ephi
      from[0].i = er;
      from[0].j = ephi;
      from[0].k = ez;
      from[0].loc = FRONT_RIGHT;
      from[0].c = 0;
      nEntries = 1;
      if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
        from[1].i = er;
        from[1].j = ephi-1;
        from[1].k = ez;
        from[1].loc = FRONT_RIGHT;
        from[1].c = 0;
        nEntries = 2;
      }
      DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom);
      to[0].i = er;
      to[0].j = ephi;
      to[0].k = ez;
      to[0].loc = FRONT_DOWN_RIGHT;
      to[0].c = 1;
      if(nEntries==1){
        valTo[0] = 0.5 * valFrom[0];
      }
      else{
        valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
      }
      DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES);

      from[0].i = er;
      from[0].j = ephi;
      from[0].k = ez;
      from[0].loc = FRONT_DOWN;
      from[0].c = 0;
      nEntries = 1;
      DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom);
      to[0].i = er;
      to[0].j = ephi;
      to[0].k = ez;
      to[0].loc = FRONT_DOWN_RIGHT;
      to[0].c = 0;
      if(nEntries==1){
        valTo[0] = 0.5 * valFrom[0];
      }
      else{
        valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
      }
      DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES);

      from[0].i = er;
      from[0].j = ephi;
      from[0].k = ez;
      from[0].loc = DOWN_RIGHT;
      from[0].c = 0;
      nEntries = 1;
      DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom);
      to[0].i = er;
      to[0].j = ephi;
      to[0].k = ez;
      to[0].loc = FRONT_DOWN_RIGHT;
      to[0].c = 2;
      if(nEntries==1){
        valTo[0] = 0.5 * valFrom[0];
      }
      else{
        valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
      }
      DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES);
    }

    VecAssemblyBegin(V);
    VecAssemblyEnd(V);
    DMRestoreLocalVector(da,&ELocal);

  PetscLogEventEnd(USER_EVENT,0,0,0,0);
    return(0);
}

PetscErrorCode EdgeToVertexProjectionMat(TS ts, Mat M, void *ptr)
{
    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;

    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    MatZeroEntries(M);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil row, col[2];
          PetscScalar valJ[2];
          PetscInt nEntries;

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_LEFT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
            col[1].i = er;
            col[1].j = ephi-1;
            col[1].k = ez;
            col[1].loc = BACK_LEFT;
            col[1].c = 0;
            valJ[1] = 0.5;

            nEntries = 2;
          }
          else {
            nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_LEFT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(ez!=0){
            col[1].i = er;
            col[1].j = ephi;
            col[1].k = ez-1;
            col[1].loc = DOWN_LEFT;
            col[1].c = 0;
            valJ[1] = 0.5;

            nEntries = 2;
          }
          else {
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(er!=0){
            col[1].i = er-1;
            col[1].j = ephi;
            col[1].k = ez;
            col[1].loc = BACK_DOWN;
            col[1].c = 0;
            valJ[1] = 0.5;

            nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          if(er==N[0]-1){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_RIGHT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = BACK_RIGHT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_RIGHT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = DOWN_RIGHT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else {
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN;
          col[0].c = 0;
          valJ[0] = 0.5;

          nEntries = 1;


          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(ez==N[2]-1){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_LEFT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype) {
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = FRONT_LEFT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(er!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_LEFT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(ephi==N[1]-1 && !user->phibtype){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_LEFT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(er!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_LEFT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = UP_LEFT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(er==N[0]-1 && ephi==N[1]-1 && !user->phibtype){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_RIGHT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;

          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP;
          col[0].c = 0;
          valJ[0] = 0.5;

          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_RIGHT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(ez!=0){
            col[1].i = er;
            col[1].j = ephi;
            col[1].k = ez-1;
            col[1].loc = UP_RIGHT;
            col[1].c = 0;
            valJ[1] = 0.5;

            nEntries = 2;
          }
          else{
            nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_LEFT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_LEFT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP_LEFT;
          col[0].c = 0;
          valJ[0] = 0.5;

          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(er!=0){
            col[1].i = er-1;
            col[1].j = ephi;
            col[1].k = ez;
            col[1].loc = FRONT_UP;
            col[1].c = 0;
            valJ[1] = 0.5;

            nEntries = 2;
          }
          else{
            nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(er==N[0]-1 && ez==N[2]-1){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_RIGHT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype) {
            col[1].i = er;
            col[1].j = ephi-1;
            col[1].k = ez;
            col[1].loc = FRONT_RIGHT;
            col[1].c = 0;
            valJ[1] = 0.5;

            nEntries = 2;
          }
          else{
            nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN;
          col[0].c = 0;
          valJ[0] = 0.5;

          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_RIGHT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;

          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

          if(er==N[0]-1 && ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_RIGHT;
          row.c = 2;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;
          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP;
          col[0].c = 0;
          valJ[0] = 0.5;
          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_RIGHT;
          row.c = 1;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_RIGHT;
          col[0].c = 0;
          valJ[0] = 0.5;
          nEntries = 1;

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);
          }

        }
      }
    }
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    return(0);
}

PetscErrorCode CellToFaceProjectionMat(TS ts, Mat M, void *ptr)
{
    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;

    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    MatZeroEntries(M);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil row, col[2];
          PetscScalar valJ[2];
          PetscInt nEntries = 2;

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = ELEMENT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if((ephi != N[1]-1 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi+1;
          col[1].k = ez;
          col[1].loc = ELEMENT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = ELEMENT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(er!=0){
          col[1].i = er-1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = ELEMENT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = ELEMENT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          col[1].i = er;
          col[1].j = ephi-1;
          col[1].k = ez;
          col[1].loc = ELEMENT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = ELEMENT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(er != N[0]-1){
          col[1].i = er+1;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = ELEMENT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = ELEMENT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(ez!=0){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez-1;
          col[1].loc = ELEMENT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = ELEMENT;
          col[0].c = 0;
          valJ[0] = 0.5;

          if(ez != N[2]-1){
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez+1;
          col[1].loc = ELEMENT;
          col[1].c = 0;
          valJ[1] = 0.5;

          nEntries = 2;
          }
          else{
          nEntries = 1;
          }

          DMStagMatSetValuesStencil(da, M, 1, & row, nEntries, col, valJ, INSERT_VALUES);

        }
      }
    }
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
    return(0);
}

PetscErrorCode FromPetscVecToArray(TS ts, Vec X, PetscScalar *gf_BR, PetscScalar *gf_BP, PetscScalar *gf_BZ, PetscScalar *g_R, PetscScalar *g_P, PetscScalar *g_Z, void *ptr)
{
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FromPetscVecToArray",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User           *user = (User*)ptr;
  DM             da, coordDA = user->coorda;
  DM             dmFr, dmFphi,dmFz, daFr, daFphi,daFz;
  DM             dmBr, dmBphi,dmBz, daBr, daBphi,daBz;
  DM             dmCoord,dmCoorda;
  PetscInt       startr,startphi,startz,nr,nphi,nz;
  Vec            vecFr, vecFphi, vecFz, F_r, F_r2, F_phi, F_phi2, F_z, F_z2, vecBr, vecBphi, vecBz, B_r, B_r2, B_phi, B_phi2, B_z, B_z2;
  Vec            B_rLocal, B_phiLocal, B_zLocal, F_rLocal, F_phiLocal, F_zLocal, XLocal;
  Vec            coordLocal, coordaLocal;
  PetscInt       N[3],er,ephi,ez,d;
  PetscInt       icBrp[3],icBphip[3],icBzp[3],icBrm[3],icBphim[3],icBzm[3];
  PetscInt       ivBrp,ivBphip,ivBzp,ivBrm,ivBphim,ivBzm;
  PetscScalar       ****arrCoord,****arrCoorda,****arrX,****arrCr,****arrCphi,****arrCz,****arrFr,****arrFphi,****arrFz,****arrVr,****arrVphi,****arrVz;
  const PetscScalar *array;
  int            len;

    TSGetDM(ts,&da);

    DMStagCreateCompatibleDMStag(da,0,0,3,0,&dmFr); /* 3 dofs per face */
    DMSetUp(dmFr);
    DMStagSetUniformCoordinatesExplicit(dmFr,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmFr,&F_r2);

    DMStagCreateCompatibleDMStag(da,0,0,3,0,&dmFphi); /* 3 dofs per face */
    DMSetUp(dmFphi);
    DMStagSetUniformCoordinatesExplicit(dmFphi,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmFphi,&F_phi2);

    DMStagCreateCompatibleDMStag(da,0,0,3,0,&dmFz); /* 3 dofs per face */
    DMSetUp(dmFz);
    DMStagSetUniformCoordinatesExplicit(dmFz,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmFz,&F_z2);

    DMStagCreateCompatibleDMStag(da,0,0,1,0,&dmBr); /* 1 dof per face */
    DMSetUp(dmBr);
    DMStagSetUniformCoordinatesExplicit(dmBr,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmBr,&B_r2);

    DMStagCreateCompatibleDMStag(da,0,0,1,0,&dmBphi); /* 1 dof per face */
    DMSetUp(dmBphi);
    DMStagSetUniformCoordinatesExplicit(dmBphi,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmBphi,&B_phi2);

    DMStagCreateCompatibleDMStag(da,0,0,1,0,&dmBz); /* 1 dof per face */
    DMSetUp(dmBz);
    DMStagSetUniformCoordinatesExplicit(dmBz,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmBz,&B_z2);

    DMGetLocalVector(da, & XLocal);
    DMGlobalToLocalBegin(da, X, INSERT_VALUES, XLocal);
    DMGlobalToLocalEnd(da, X, INSERT_VALUES, XLocal);


    //PetscPrintf(PETSC_COMM_WORLD,"Before copying B_r values\n");
    DMStagGetCorners(dmBr,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmBr,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[1];
                PetscScalar   valFrom[1];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = LEFT;    from[0].c = 0;
                DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom);
                DMStagVecSetValuesStencil(dmBr,B_r2,1,from,valFrom,INSERT_VALUES);
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = RIGHT;    from[0].c = 0;
                    DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom);
                    DMStagVecSetValuesStencil(dmBr,B_r2,1,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(B_r2);
    VecAssemblyEnd(B_r2);

    DMStagVecSplitToDMDA(dmBr,B_r2,LEFT,-1,&daBr,&vecBr); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecBr,"rFace_center_values");


    //PetscPrintf(PETSC_COMM_WORLD,"Before copying B_phi values\n");
    DMStagGetCorners(dmBphi,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmBphi,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[1];
                PetscScalar   valFrom[1];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN;    from[0].c = 0;
                DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom);
                DMStagVecSetValuesStencil(dmBphi,B_phi2,1,from,valFrom,INSERT_VALUES);
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP; from[0].c = 0;
                    DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom);
                    DMStagVecSetValuesStencil(dmBphi,B_phi2,1,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(B_phi2);
    VecAssemblyEnd(B_phi2);

    DMStagVecSplitToDMDA(dmBphi,B_phi2,DOWN,-1,&daBphi,&vecBphi); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecBphi,"phiFace_center_values");


    //PetscPrintf(PETSC_COMM_WORLD,"Before copying B_z values\n");
    DMStagGetCorners(dmBz,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmBz,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK;    from[0].c = 0;
                DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom);
                DMStagVecSetValuesStencil(dmBz,B_z2,1,from,valFrom,INSERT_VALUES);
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT; from[0].c = 0;
                    DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom);
                    DMStagVecSetValuesStencil(dmBz,B_z2,1,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(B_z2);
    VecAssemblyEnd(B_z2);

    DMStagVecSplitToDMDA(dmBz,B_z2,BACK,-1,&daBz,&vecBz); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecBz,"zFace_center_values");


    //PetscPrintf(PETSC_COMM_WORLD,"Before copying F_r coordinates\n");
    DMGetCoordinatesLocal(dmFr, &F_r);
    DMStagGetCorners(dmFr,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmFr,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = LEFT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = LEFT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = LEFT;    from[2].c = 2;
                DMStagVecGetValuesStencil(dmFr,F_r,3,from,valFrom);
                DMStagVecSetValuesStencil(dmFr,F_r2,3,from,valFrom,INSERT_VALUES);
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = RIGHT;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = RIGHT;    from[2].c = 2;
                    DMStagVecGetValuesStencil(dmFr,F_r,3,from,valFrom);
                    DMStagVecSetValuesStencil(dmFr,F_r2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(F_r2);
    VecAssemblyEnd(F_r2);

    DMStagVecSplitToDMDA(dmFr,F_r2,LEFT,-3,&daFr,&vecFr); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecFr,"rFace_center_coordinates");


    //PetscPrintf(PETSC_COMM_WORLD,"Before copying F_phi coordinates\n");
    DMGetCoordinatesLocal(dmFphi, &F_phi);
    DMStagGetCorners(dmFphi,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmFphi,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = DOWN;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = DOWN;    from[2].c = 2;
                DMStagVecGetValuesStencil(dmFphi,F_phi,3,from,valFrom);
                DMStagVecSetValuesStencil(dmFphi,F_phi2,3,from,valFrom,INSERT_VALUES);
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = UP;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = UP;    from[2].c = 2;
                    DMStagVecGetValuesStencil(dmFphi,F_phi,3,from,valFrom);
                    DMStagVecSetValuesStencil(dmFphi,F_phi2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(F_phi2);
    VecAssemblyEnd(F_phi2);

    DMStagVecSplitToDMDA(dmFphi,F_phi2,DOWN,-3,&daFphi,&vecFphi); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecFphi,"phiFace_center_coordinates");


    //PetscPrintf(PETSC_COMM_WORLD,"Before copying F_z coordinates\n");
    DMGetCoordinatesLocal(dmFz, &F_z);
    DMStagGetCorners(dmFz,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmFz,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK;    from[2].c = 2;
                DMStagVecGetValuesStencil(dmFz,F_z,3,from,valFrom);
                DMStagVecSetValuesStencil(dmFz,F_z2,3,from,valFrom,INSERT_VALUES);
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT;    from[2].c = 2;
                    DMStagVecGetValuesStencil(dmFz,F_z,3,from,valFrom);
                    DMStagVecSetValuesStencil(dmFz,F_z2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(F_z2);
    VecAssemblyEnd(F_z2);

    DMStagVecSplitToDMDA(dmFz,F_z2,BACK,-3,&daFz,&vecFz); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecFz,"zFace_center_coordinates");


    PetscMPIInt rank;
    MPI_Comm    comm;
    VecScatter  scat;
    Vec         Xseq, naturalX;

    DMDACreateNaturalVector(daFr,&naturalX);
    DMDAGlobalToNaturalBegin(daFr, vecFr, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daFr, vecFr, INSERT_VALUES, naturalX);

    /* create scater to zero */
    //VecScatterCreateToZero(naturalX, &scat, &Xseq);
    VecScatterCreateToAll(naturalX, &scat, &Xseq);
    VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* Only rank == 0 has the entries of the patch, so run code only at that rank */
    if (rank == 0 || 1) {
      PetscInt sizeX;
      VecGetSize(Xseq, &sizeX);
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*(user->Nr+1)*user->Nz);
      VecGetArrayRead(Xseq, &array);
      len = sizeof(array)/sizeof(*array); // finding size of array
      //PetscPrintf(PETSC_COMM_SELF, "len = %d\n", len);
      memcpy(g_R, array, 3*(user->Nr+1)*user->Nphi*user->Nz*(sizeof(PetscScalar)));
      VecRestoreArrayRead(Xseq, &array);
    }

    VecDestroy(&Xseq);
    VecScatterDestroy(&scat);
    VecDestroy(&naturalX);


    DMDACreateNaturalVector(daFphi,&naturalX);
    DMDAGlobalToNaturalBegin(daFphi, vecFphi, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daFphi, vecFphi, INSERT_VALUES, naturalX);

    /* create scater to zero */
    //VecScatterCreateToZero(naturalX, &scat, &Xseq);
    VecScatterCreateToAll(naturalX, &scat, &Xseq);
    VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* Only rank == 0 has the entries of the patch, so run code only at that rank */
    if (rank == 0 || 1) {
      PetscInt sizeX;
      VecGetSize(Xseq, &sizeX);
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,(user->Nphi+1)*user->Nr*user->Nz);
      VecGetArrayRead(Xseq, &array);
      memcpy(g_P, array, 3*(user->Nphi+1)*user->Nr*user->Nz*(sizeof(PetscScalar)));
      VecRestoreArrayRead(Xseq, &array);
    }

    VecDestroy(&Xseq);
    VecScatterDestroy(&scat);
    VecDestroy(&naturalX);


    DMDACreateNaturalVector(daFz,&naturalX);
    DMDAGlobalToNaturalBegin(daFz, vecFz, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daFz, vecFz, INSERT_VALUES, naturalX);

    /* create scater to zero */
    //VecScatterCreateToZero(naturalX, &scat, &Xseq);
    VecScatterCreateToAll(naturalX, &scat, &Xseq);
    VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* Only rank == 0 has the entries of the patch, so run code only at that rank */
    if (rank == 0 || 1) {
      PetscInt sizeX;
      VecGetSize(Xseq, &sizeX);
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*user->Nr*(user->Nz+1));
      VecGetArrayRead(Xseq, &array);
      memcpy(g_Z, array, 3*(user->Nz+1)*user->Nr*user->Nphi*(sizeof(PetscScalar)));
      VecRestoreArrayRead(Xseq, &array);
    }

    VecDestroy(&Xseq);
    VecScatterDestroy(&scat);
    VecDestroy(&naturalX);


    DMDACreateNaturalVector(daBr,&naturalX);
    DMDAGlobalToNaturalBegin(daBr, vecBr, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daBr, vecBr, INSERT_VALUES, naturalX);

    /* create scater to zero */
    //VecScatterCreateToZero(naturalX, &scat, &Xseq);
    VecScatterCreateToAll(naturalX, &scat, &Xseq);
    VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* Only rank == 0 has the entries of the patch, so run code only at that rank */
    if (rank == 0 || 1) {
      PetscInt sizeX;
      VecGetSize(Xseq, &sizeX);
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*(user->Nr+1)*user->Nz);
      VecGetArrayRead(Xseq, &array);
      memcpy(gf_BR, array, (user->Nr+1)*user->Nz*user->Nphi*(sizeof(PetscScalar)));
      VecRestoreArrayRead(Xseq, &array);
    }

    VecDestroy(&Xseq);
    VecScatterDestroy(&scat);
    VecDestroy(&naturalX);


    DMDACreateNaturalVector(daBphi,&naturalX);
    DMDAGlobalToNaturalBegin(daBphi, vecBphi, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daBphi, vecBphi, INSERT_VALUES, naturalX);

    /* create scater to zero */
    //VecScatterCreateToZero(naturalX, &scat, &Xseq);
    VecScatterCreateToAll(naturalX, &scat, &Xseq);
    VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* Only rank == 0 has the entries of the patch, so run code only at that rank */
    if (rank == 0 || 1) {
      PetscInt sizeX;
      VecGetSize(Xseq, &sizeX);
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,(user->Nphi+1)*user->Nr*user->Nz);
      VecGetArrayRead(Xseq, &array);
      memcpy(gf_BP, array, (user->Nphi+1)*user->Nz*user->Nr*(sizeof(PetscScalar)));
      VecRestoreArrayRead(Xseq, &array);
    }

    VecDestroy(&Xseq);
    VecScatterDestroy(&scat);
    VecDestroy(&naturalX);


    DMDACreateNaturalVector(daBz,&naturalX);
    DMDAGlobalToNaturalBegin(daBz, vecBz, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daBz, vecBz, INSERT_VALUES, naturalX);

    /* create scater to zero */
    //VecScatterCreateToZero(naturalX, &scat, &Xseq);
    VecScatterCreateToAll(naturalX, &scat, &Xseq);
    VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* Only rank == 0 has the entries of the patch, so run code only at that rank */
    if (rank == 0 || 1) {
      PetscInt sizeX;
      VecGetSize(Xseq, &sizeX);
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*user->Nr*(user->Nz+1));
      VecGetArrayRead(Xseq, &array);
      memcpy(gf_BZ, array, (user->Nz+1)*user->Nphi*user->Nr*(sizeof(PetscScalar)));
      VecRestoreArrayRead(Xseq, &array);
    }

    VecDestroy(&Xseq);
    VecScatterDestroy(&scat);
    VecDestroy(&naturalX);


    DMDestroy(&dmFr);
    DMDestroy(&dmFphi);
    DMDestroy(&dmFz);

    DMDestroy(&daFr);
    DMDestroy(&daFphi);
    DMDestroy(&daFz);

    DMDestroy(&dmBr);
    DMDestroy(&dmBphi);
    DMDestroy(&dmBz);

    DMDestroy(&daBr);
    DMDestroy(&daBphi);
    DMDestroy(&daBz);

    VecDestroy(&vecFr);
    VecDestroy(&vecFphi);
    VecDestroy(&vecFz);
    VecDestroy(&vecBr);
    VecDestroy(&vecBphi);
    VecDestroy(&vecBz);

    VecDestroy(&F_r2);
    VecDestroy(&F_phi2);
    VecDestroy(&F_z2);
    VecDestroy(&B_r2);
    VecDestroy(&B_phi2);
    VecDestroy(&B_z2);

  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return(0);
}

PetscErrorCode CellCoordArrays(TS ts, PetscScalar *vecCR, PetscScalar* vecCZ, void *ptr)
{
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("CellCoordArrays",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;
    DM             daCr, daCz, dmC;

    PetscInt       startr,startphi,startz,nr,nphi,nz;

    Vec            vecCr, Cr2, vecCz, C, Cz2;
    Vec            coordLocal, coordaLocal;
    PetscInt       N[3],er,ephi,ez,d;
    PetscInt       icp[3];
    DM             dmCoord,dmCoorda;


    TSGetDM(ts,&da);

    DMStagCreateCompatibleDMStag(da,0,0,0,3,&dmC); /* 1 dof per cell */
    DMSetUp(dmC);
    DMStagSetUniformCoordinatesExplicit(dmC,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmC,&Cr2);
    DMCreateGlobalVector(dmC,&Cz2);

    //PetscPrintf(PETSC_COMM_WORLD,"Before copying cell coordinates\n");
    DMGetCoordinatesLocal(dmC, &C);
    DMStagGetCorners(dmC,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[1],to[1];
                PetscScalar   valFrom[1];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = ELEMENT;    from[0].c = 0;
                DMStagVecGetValuesStencil(dmC,C,1,from,valFrom);
                DMStagVecSetValuesStencil(dmC,Cr2,1,from,valFrom,INSERT_VALUES);
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = ELEMENT;    from[0].c = 2;
                to[0].i = er; to[0].j = ephi; to[0].k = ez; to[0].loc = ELEMENT;    to[0].c = 0;
                DMStagVecGetValuesStencil(dmC,C,1,from,valFrom);
                DMStagVecSetValuesStencil(dmC,Cz2,1,to,valFrom,INSERT_VALUES);
            }
        }
    }
    VecAssemblyBegin(Cr2);
    VecAssemblyEnd(Cr2);
    VecAssemblyBegin(Cz2);
    VecAssemblyEnd(Cz2);

    DMStagVecSplitToDMDA(dmC,Cr2,ELEMENT,-1,&daCr,&vecCr); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecCr,"Cell_center_r_coordinates");
    DMStagVecSplitToDMDA(dmC,Cz2,ELEMENT,-1,&daCz,&vecCz); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecCr,"Cell_center_z_coordinates");

    PetscMPIInt rank;
    MPI_Comm    comm;
    VecScatter  scat;
    Vec         Xseq, naturalX;
    const PetscScalar *array;

    DMDACreateNaturalVector(daCr,&naturalX);
    DMDAGlobalToNaturalBegin(daCr, vecCr, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daCr, vecCr, INSERT_VALUES, naturalX);

    /* create scater to zero */
    //VecScatterCreateToZero(naturalX, &scat, &Xseq);
    VecScatterCreateToAll(naturalX, &scat, &Xseq);
    VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* Only rank == 0 has the entries of the patch, so run code only at that rank */
    if (rank == 0 || 1) {
      PetscInt sizeX;
      VecGetSize(Xseq, &sizeX);
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*user->Nr*user->Nz);
      VecGetArrayRead(Xseq, &array);
      memcpy(vecCR, array, user->Nr*user->Nz*user->Nphi*(sizeof(PetscScalar)));
      VecRestoreArrayRead(Xseq, &array);
    }

    VecDestroy(&Xseq);
    VecScatterDestroy(&scat);
    VecDestroy(&naturalX);

    DMDACreateNaturalVector(daCz,&naturalX);
    DMDAGlobalToNaturalBegin(daCz, vecCz, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daCz, vecCz, INSERT_VALUES, naturalX);

    /* create scater to zero */
    //VecScatterCreateToZero(naturalX, &scat, &Xseq);
    VecScatterCreateToAll(naturalX, &scat, &Xseq);
    VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* Only rank == 0 has the entries of the patch, so run code only at that rank */
    if (rank == 0 || 1) {
      PetscInt sizeX;
      VecGetSize(Xseq, &sizeX);
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*user->Nr*user->Nz);
      VecGetArrayRead(Xseq, &array);
      memcpy(vecCZ, array, user->Nr*user->Nz*user->Nphi*(sizeof(PetscScalar)));
      VecRestoreArrayRead(Xseq, &array);
    }

    VecDestroy(&Xseq);
    VecScatterDestroy(&scat);
    VecDestroy(&naturalX);


    DMDestroy(&dmC);

    DMDestroy(&daCr);
    DMDestroy(&daCz);

    VecDestroy(&vecCr);
    VecDestroy(&vecCz);

    VecDestroy(&Cr2);
    VecDestroy(&Cz2);

  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return(0);
}

PetscErrorCode ScatterTest(TS ts, void *ptr)
{
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("ScatterTest",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;
    DM             daC, dmC;
    DM             dmFr, dmFphi,dmFz, daFr, daFphi,daFz;
    DM             dmEr, dmEphi,dmEz, daEr, daEphi,daEz;
    DM             dmV, daV;

    PetscInt       startr,startphi,startz,nr,nphi,nz;

    Vec            vecC, C, C2, vecFr, vecFphi, vecFz, F_r, F_r2, F_phi, F_phi2, F_z, F_z2, vecEr, vecEphi, vecEz, E_r, E_r2, E_phi, E_phi2, E_z, E_z2, vecV, V, V2;

    Vec            C_rLocal, C_phiLocal, C_zLocal, F_rLocal, F_phiLocal, F_zLocal, V_rLocal, V_phiLocal, V_zLocal, xLocal;
    Vec            coordLocal, coordaLocal;
    PetscInt       N[3],er,ephi,ez,d;

    PetscInt       icBrp[3],icBphip[3],icBzp[3],icBrm[3],icBphim[3],icBzm[3];
    PetscInt       icErmzm[3],icErmzp[3],icErpzm[3],icErpzp[3];
    PetscInt       icEphimzm[3],icEphipzm[3],icEphimzp[3],icEphipzp[3];
    PetscInt       icErmphim[3],icErpphim[3],icErmphip[3],icErpphip[3];
    PetscInt       icrmphimzm[3],icrpphimzm[3],icrmphipzm[3],icrpphipzm[3];
    PetscInt       icrmphimzp[3],icrpphimzp[3],icrmphipzp[3],icrpphipzp[3];
    PetscInt       icp[3];


    PetscInt        ivBrp,ivBphip,ivBzp,ivBrm,ivBphim,ivBzm;

    PetscInt        ivErmzm,ivErmzp,ivErpzm,ivErpzp;
    PetscInt        ivEphimzm,ivEphipzm,ivEphimzp,ivEphipzp;
    PetscInt        ivErmphim,ivErpphim,ivErmphip,ivErpphip;
    DM              dmCoord,dmCoorda;
    PetscScalar       ****arrCoord,****arrCoorda,****arrX,****arrCr,****arrCphi,****arrCz,****arrFr,****arrFphi,****arrFz,****arrVr,****arrVphi,****arrVz;

    TSGetDM(ts,&da);

    DMStagCreateCompatibleDMStag(da,0,0,3,0,&dmFr); /* 3 dofs per face */
    DMSetUp(dmFr);
    DMStagSetUniformCoordinatesExplicit(dmFr,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmFr,&F_r2);

    DMStagCreateCompatibleDMStag(da,0,0,3,0,&dmFphi); /* 3 dofs per face */
    DMSetUp(dmFphi);
    DMStagSetUniformCoordinatesExplicit(dmFphi,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmFphi,&F_phi2);

    DMStagCreateCompatibleDMStag(da,0,0,3,0,&dmFz); /* 3 dofs per face */
    DMSetUp(dmFz);
    DMStagSetUniformCoordinatesExplicit(dmFz,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmFz,&F_z2);

    DMStagCreateCompatibleDMStag(da,0,0,0,3,&dmC); /* 3 dofs per cell */
    DMSetUp(dmC);
    DMStagSetUniformCoordinatesExplicit(dmC,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax);
    DMCreateGlobalVector(dmC,&C2);


    PetscPrintf(PETSC_COMM_WORLD,"Before copying F_r coordinates\n");
    DMGetCoordinatesLocal(dmFr, &F_r);
    DMStagGetCorners(dmFr,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmFr,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = LEFT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = LEFT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = LEFT;    from[2].c = 2;
                DMStagVecGetValuesStencil(dmFr,F_r,3,from,valFrom);
                DMStagVecSetValuesStencil(dmFr,F_r2,3,from,valFrom,INSERT_VALUES);
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = RIGHT;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = RIGHT;    from[2].c = 2;
                    DMStagVecGetValuesStencil(dmFr,F_r,3,from,valFrom);
                    DMStagVecSetValuesStencil(dmFr,F_r2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(F_r2);
    VecAssemblyEnd(F_r2);

    DMStagVecSplitToDMDA(dmFr,F_r2,LEFT,-3,&daFr,&vecFr); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecFr,"rFace_center_coordinates");


    PetscPrintf(PETSC_COMM_WORLD,"Before copying F_phi coordinates\n");
    DMGetCoordinatesLocal(dmFphi, &F_phi);
    DMStagGetCorners(dmFphi,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmFphi,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = DOWN;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = DOWN;    from[2].c = 2;
                DMStagVecGetValuesStencil(dmFphi,F_phi,3,from,valFrom);
                DMStagVecSetValuesStencil(dmFphi,F_phi2,3,from,valFrom,INSERT_VALUES);
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = UP;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = UP;    from[2].c = 2;
                    DMStagVecGetValuesStencil(dmFphi,F_phi,3,from,valFrom);
                    DMStagVecSetValuesStencil(dmFphi,F_phi2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(F_phi2);
    VecAssemblyEnd(F_phi2);

    DMStagVecSplitToDMDA(dmFphi,F_phi2,DOWN,-3,&daFphi,&vecFphi); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecFphi,"phiFace_center_coordinates");


    PetscPrintf(PETSC_COMM_WORLD,"Before copying F_z coordinates\n");
    DMGetCoordinatesLocal(dmFz, &F_z);
    DMStagGetCorners(dmFz,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    DMStagGetGlobalSizes(dmFz,&N[0],&N[1],&N[2]);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK;    from[2].c = 2;
                DMStagVecGetValuesStencil(dmFz,F_z,3,from,valFrom);
                DMStagVecSetValuesStencil(dmFz,F_z2,3,from,valFrom,INSERT_VALUES);
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT;    from[2].c = 2;
                    DMStagVecGetValuesStencil(dmFz,F_z,3,from,valFrom);
                    DMStagVecSetValuesStencil(dmFz,F_z2,3,from,valFrom,INSERT_VALUES);
                }
            }
        }
    }
    VecAssemblyBegin(F_z2);
    VecAssemblyEnd(F_z2);

    DMStagVecSplitToDMDA(dmFz,F_z2,BACK,-3,&daFz,&vecFz); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecFz,"zFace_center_coordinates");


    PetscPrintf(PETSC_COMM_WORLD,"Before copying cell coordinates\n");
    DMGetCoordinatesLocal(dmC, &C);
    DMStagGetCorners(dmC,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL);
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = ELEMENT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = ELEMENT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = ELEMENT;    from[2].c = 2;
                DMStagVecGetValuesStencil(dmC,C,3,from,valFrom);
                DMStagVecSetValuesStencil(dmC,C2,3,from,valFrom,INSERT_VALUES);
            }
        }
    }
    VecAssemblyBegin(C2);
    VecAssemblyEnd(C2);

    DMStagVecSplitToDMDA(dmC,C2,ELEMENT,-3,&daC,&vecC); /* note -3 : pad with zero */
    PetscObjectSetName((PetscObject)vecC,"Cell_center_coordinates");

    PetscMPIInt rank;
    MPI_Comm    comm;
    VecScatter  scat;
    Vec         Xseq, naturalX;
    const PetscScalar *array;

    DMDACreateNaturalVector(daC,&naturalX);
    DMDAGlobalToNaturalBegin(daC, vecC, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daC, vecC, INSERT_VALUES, naturalX);

    /* create scater to zero */
    VecScatterCreateToZero(naturalX, &scat, &Xseq);
    VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    /* Only rank == 0 has the entries of the patch, so run code only at that rank */
    if (rank == 0) {
      PetscInt sizeX;
      VecGetSize(Xseq, &sizeX);
      PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*user->Nr*user->Nz);
      VecGetArrayRead(Xseq, &array);
      PetscInt      i, j, k;

      /* Loop over the patch of the entire domain */
      for (k = 0; k < user->Nz; k++) {
        for (j = 0; j < user->Nphi; j++) {
          for (i = 0; i < user->Nr; i++) {
            PetscPrintf(PETSC_COMM_SELF, "The (r,phi,z) cell center coordinates at global index %d + %d*Nr + %d*Nr*Nphi are (%g,%g,%g)\n", i, j, k, (double)array[3*(k*user->Nr*user->Nphi+j*user->Nr+i)], (double)array[1+3*(k*user->Nr*user->Nphi+j*user->Nr+i)], (double)array[2+3*(k*user->Nr*user->Nphi+j*user->Nr+i)]);
          }
        }
      }

      VecRestoreArrayRead(Xseq, &array);
    }

    VecDestroy(&Xseq);
    VecDestroy(&naturalX);

    DMDestroy(&dmFr);
    DMDestroy(&dmFphi);
    DMDestroy(&dmFz);
    DMDestroy(&dmC);

    DMDestroy(&daFr);
    DMDestroy(&daFphi);
    DMDestroy(&daFz);
    DMDestroy(&daC);

    VecDestroy(&vecC);
    VecDestroy(&vecFr);
    VecDestroy(&vecFphi);
    VecDestroy(&vecFz);

    VecDestroy(&C2);
    VecDestroy(&F_r2);
    VecDestroy(&F_phi2);
    VecDestroy(&F_z2);

  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return(0);
}

PetscErrorCode ComputeIsEBoundary(TS ts, IS * isE_boundary, void * ptr) {

  User * user = (User * ) ptr;
  DM da;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  PetscInt N[3], er, ephi, ez, n = 1;

  PetscReal t;
  Vec dummyX;
  Mat Jpre;
  PC dummypc;
  KSP dummyksp, * subksp;
  PetscErrorCode ierr = 0;

  TSGetDM(ts, & da);
  DMCreateGlobalVector(da, & dummyX);
  FormInitialSolution(ts, dummyX, user);
  TSGetTime(ts, & t);
  DMCreateMatrix(da, & Jpre);
  MatZeroEntries(Jpre);

  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col[1];
        PetscScalar valJ[1];
        PetscInt nEntries;

        /* The edges are oriented in the directions of unit vectors: e_r, e_phi and e_z */

        /* E field part */

        if (er == 0 || ez == 0) {
          /* Equation on the back or left boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if ((user -> phibtype && (er == 0)) || (((er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on the left or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if ((user -> phibtype && (ez == 0)) || (((ephi == 0 || ez == 0)) && !(user -> phibtype))) {
          /* Equation on the back or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (er == N[0] - 1) {
          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_LEFT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);

          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP_RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && ez == N[2] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Equation on the front and right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_RIGHT;
          col[0].c = 0;
          valJ[0] = 1.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if ((user -> phibtype && !(er == 0)) || ((!(er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal down left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = DOWN_LEFT;
          col[0].c = 0;
          valJ[0] = 0.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if ((user -> phibtype && !(ez == 0)) || ((!(ez == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal back down edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN;
          col[0].c = 0;
          valJ[0] = 0.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (!(er == 0 || ez == 0)) {
          /* Equation on internal back left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_LEFT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_LEFT;
          col[0].c = 0;
          valJ[0] = 0.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        /* B field part */

        /*if (!(er == 0)){*/
        /* Equation on left face */
        nEntries = 1;
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = LEFT;
        row.c = 0;

        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = LEFT;
        col[0].c = 0;
        valJ[0] = 0.0;

        DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        /*}*/

        /*if(user->phibtype || (ephi != 0 && !(user->phibtype)) ){*/
        /* Equation on down face */
        nEntries = 1;
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN;
        row.c = 0;

        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = DOWN;
        col[0].c = 0;
        valJ[0] = 0.0;

        DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        /*}*/

        /*if (!(ez == 0)){*/
        /* Equation on back face */
        nEntries = 1;
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK;
        row.c = 0;

        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = BACK;
        col[0].c = 0;
        valJ[0] = 0.0;

        DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        /*}*/

        if (er == N[0] - 1) {
          /* Equation on right boundary face */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = RIGHT;
          col[0].c = 0;
          valJ[0] = 0.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on up boundary face */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = UP;
          col[0].c = 0;
          valJ[0] = 0.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on front boundary face */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;

          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT;
          col[0].c = 0;
          valJ[0] = 0.0;

          DMStagMatSetValuesStencil(da, Jpre, 1, & row, nEntries, col, valJ, INSERT_VALUES);
        }

      }
    }
  }

  MatAssemblyBegin(Jpre, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Jpre, MAT_FINAL_ASSEMBLY);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Jpre:\n");
    MatView(Jpre, PETSC_VIEWER_STDOUT_WORLD);
  }

  KSPCreate(PETSC_COMM_WORLD, & dummyksp);
  KSPSetFromOptions(dummyksp);
  KSPSetOperators(dummyksp, Jpre, Jpre);
  PetscBarrier((PetscObject)Jpre);
  MatDestroy( & Jpre);
  KSPGetPC(dummyksp, & dummypc);
  PCSetType(dummypc, PCFIELDSPLIT);
  PCFieldSplitSetDetectSaddlePoint(dummypc, PETSC_TRUE);
  PCSetUp(dummypc);
  PCFieldSplitSetSchurFactType(dummypc, PC_FIELDSPLIT_SCHUR_FACT_FULL);
  PCFieldSplitSetSchurPre(dummypc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
  KSPSetUp(dummyksp);
  KSPSetTolerances(dummyksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 outer iteration for the dummy solve
  PCFieldSplitGetSubKSP(dummypc, & n, & subksp);
  KSPSetTolerances(subksp[1], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 inner iteration maximum for the dummy solve
  KSPSolve(dummyksp, dummyX, dummyX);
  PetscBarrier((PetscObject)dummyX);
  VecDestroy( & dummyX);
  IS isdup;
  PCFieldSplitGetISByIndex(dummypc, 0, & isdup);
  ISDuplicate(isdup, isE_boundary);
  PetscBarrier((PetscObject)dummyksp);
  KSPDestroy( & dummyksp);

  return (0);
}

PetscErrorCode VertexToEdgeReconstruction_scalar(TS ts, Vec V, Vec E, void *ptr)
{
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("VertexToEdgeReconstruction_scalar",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

    User           *user = (User*)ptr;
    DM             da, coordDA = user->coorda;
    Vec            VLocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez, n = 1;

    PetscErrorCode ierr = 0;

    VecZeroEntries(E);
    TSGetDM(ts, & da);
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
    DMGetLocalVector(da,&VLocal);
    DMGlobalToLocal(da,V,INSERT_VALUES,VLocal);

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[8], to[12];
          PetscScalar valFrom[8], valTo[12];

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = FRONT_UP_RIGHT;
          from[0].c = 3;
          from[1].i = er;
          from[1].j = ephi;
          from[1].k = ez;
          from[1].loc = FRONT_UP_LEFT;
          from[1].c = 3;
          from[2].i = er;
          from[2].j = ephi;
          from[2].k = ez;
          from[2].loc = FRONT_DOWN_LEFT;
          from[2].c = 3;
          from[3].i = er;
          from[3].j = ephi;
          from[3].k = ez;
          from[3].loc = FRONT_DOWN_RIGHT;
          from[3].c = 3;
          from[4].i = er;
          from[4].j = ephi;
          from[4].k = ez;
          from[4].loc = BACK_UP_RIGHT;
          from[4].c = 3;
          from[5].i = er;
          from[5].j = ephi;
          from[5].k = ez;
          from[5].loc = BACK_UP_LEFT;
          from[5].c = 3;
          from[6].i = er;
          from[6].j = ephi;
          from[6].k = ez;
          from[6].loc = BACK_DOWN_LEFT;
          from[6].c = 3;
          from[7].i = er;
          from[7].j = ephi;
          from[7].k = ez;
          from[7].loc = BACK_DOWN_RIGHT;
          from[7].c = 3;
          DMStagVecGetValuesStencil(da, VLocal, 8, from, valFrom);
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = FRONT_UP;
          to[0].c = 0;
          valTo[0] = 0.5 * (valFrom[0] + valFrom[1]);
          to[1].i = er;
          to[1].j = ephi;
          to[1].k = ez;
          to[1].loc = FRONT_LEFT;
          to[1].c = 0;
          valTo[1] = 0.5 * (valFrom[1] + valFrom[2]);
          to[2].i = er;
          to[2].j = ephi;
          to[2].k = ez;
          to[2].loc = FRONT_DOWN;
          to[2].c = 0;
          valTo[2] = 0.5 * (valFrom[2] + valFrom[3]);
          to[3].i = er;
          to[3].j = ephi;
          to[3].k = ez;
          to[3].loc = FRONT_RIGHT;
          to[3].c = 0;
          valTo[3] = 0.5 * (valFrom[0] + valFrom[3]);
          to[4].i = er;
          to[4].j = ephi;
          to[4].k = ez;
          to[4].loc = BACK_UP;
          to[4].c = 0;
          valTo[4] = 0.5 * (valFrom[4] + valFrom[5]);
          to[5].i = er;
          to[5].j = ephi;
          to[5].k = ez;
          to[5].loc = BACK_LEFT;
          to[5].c = 0;
          valTo[5] = 0.5 * (valFrom[5] + valFrom[6]);
          to[6].i = er;
          to[6].j = ephi;
          to[6].k = ez;
          to[6].loc = BACK_DOWN;
          to[6].c = 0;
          valTo[6] = 0.5 * (valFrom[6] + valFrom[7]);
          to[7].i = er;
          to[7].j = ephi;
          to[7].k = ez;
          to[7].loc = BACK_RIGHT;
          to[7].c = 0;
          valTo[7] = 0.5 * (valFrom[4] + valFrom[7]);
          to[8].i = er;
          to[8].j = ephi;
          to[8].k = ez;
          to[8].loc = UP_RIGHT;
          to[8].c = 0;
          valTo[8] = 0.5 * (valFrom[0] + valFrom[4]);
          to[9].i = er;
          to[9].j = ephi;
          to[9].k = ez;
          to[9].loc = UP_LEFT;
          to[9].c = 0;
          valTo[9] = 0.5 * (valFrom[1] + valFrom[5]);
          to[10].i = er;
          to[10].j = ephi;
          to[10].k = ez;
          to[10].loc = DOWN_RIGHT;
          to[10].c = 0;
          valTo[10] = 0.5 * (valFrom[3] + valFrom[7]);
          to[11].i = er;
          to[11].j = ephi;
          to[11].k = ez;
          to[11].loc = DOWN_LEFT;
          to[11].c = 0;
          valTo[11] = 0.5 * (valFrom[2] + valFrom[6]);
          DMStagVecSetValuesStencil(da, E, 12, to, valTo, INSERT_VALUES);

        }
      }
    }
    VecAssemblyBegin(E);
    VecAssemblyEnd(E);
    DMRestoreLocalVector(da,&VLocal);

    PetscLogEventEnd(USER_EVENT,0,0,0,0);
    return(0);
}

PetscErrorCode getJArray(TS ts, Vec X, PetscScalar *gf_V, void *ptr) {
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("getVArray",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  User           *user = (User*)ptr;
  DM             da, dmV, daV;
  PetscInt       startr,startphi,startz,nr,nphi,nz;
  Vec            vecV, V, X_local;
  PetscInt       er,ephi,ez,d;
  const PetscScalar *array;
  int            len;

  PetscCall(TSGetDM(ts,& da));

  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmV); /* 3 dofs per element */
  PetscCall(DMSetUp(dmV));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmV, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));
  PetscCall(DMCreateGlobalVector(dmV, & V));
  PetscCall(DMGetLocalVector(da, & X_local));
  PetscCall(DMGlobalToLocal(da, X, INSERT_VALUES, X_local));

  PetscCall(DMStagGetCorners(dmV, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));

  for (ez = startz; ez < startz + nz; ++ez)
  {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi)
    {
      for (er = startr; er < startr + nr; ++er)
      {
        DMStagStencil from[24], to[3];
        PetscScalar valFrom[24], valTo[3];
        for (PetscInt comp = 0; comp < 3; ++comp)
        {
          for (PetscInt index = 0; index < 8; ++index)
          {
            from[index + comp*8].i = er;
            from[index + comp*8].j = ephi;
            from[index + comp*8].k = ez;
            from[index + comp*8].c = comp;
          }
          from[0 + comp*8].loc = BACK_DOWN_LEFT;
          from[1 + comp*8].loc = BACK_DOWN_RIGHT;
          from[2 + comp*8].loc = BACK_UP_LEFT;
          from[3 + comp*8].loc = BACK_UP_RIGHT;
          from[4 + comp*8].loc = FRONT_DOWN_LEFT;
          from[5 + comp*8].loc = FRONT_DOWN_RIGHT;
          from[6 + comp*8].loc = FRONT_UP_LEFT;
          from[7 + comp*8].loc = FRONT_UP_RIGHT;
        }

        PetscCall(DMStagVecGetValuesStencil(da, X_local, 24, from, valFrom));

        for (PetscInt comp = 0; comp < 3; ++comp)
        {
          to[comp].i = er;
          to[comp].j = ephi;
          to[comp].k = ez;
          to[comp].loc = ELEMENT;
          to[comp].c = comp;
          valTo[comp] = 0.0;
          for (PetscInt index = 0; index < 8; ++index)
            valTo[comp] += valFrom[index + comp*8];
          valTo[comp] /= 8.0;
        }

        PetscCall(DMStagVecSetValuesStencil(dmV, V, 3, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(V));
  PetscCall(VecAssemblyEnd(V));

  DMStagVecSplitToDMDA(dmV, V, ELEMENT, -3, & daV, & vecV); /* note -3 : pad with zero in 2D case */
  PetscCall(PetscObjectSetName((PetscObject) vecV, "Velocity"));

  PetscMPIInt rank;
  MPI_Comm    comm;
  VecScatter  scat;
  Vec         Xseq, naturalX;


  DMDACreateNaturalVector(daV,&naturalX);
  DMDAGlobalToNaturalBegin(daV, vecV, INSERT_VALUES, naturalX);
  DMDAGlobalToNaturalEnd(daV, vecV, INSERT_VALUES, naturalX);

  /* create scater to zero */
  //VecScatterCreateToZero(naturalX, &scat, &Xseq);
  VecScatterCreateToAll(naturalX, &scat, &Xseq);
  VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD);

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  /* Only rank == 0 has the entries of the patch, so run code only at that rank */
  if (rank == 0 || 1) {
    PetscInt sizeX;
    PetscCall(VecGetSize(Xseq, &sizeX));
    //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*(user->Nr+1)*user->Nz);
    PetscCall(VecGetArrayRead(Xseq, &array));
    memcpy(gf_V, array, 3*user->Nr*user->Nz*user->Nphi*(sizeof(PetscScalar)));
    PetscCall(VecRestoreArrayRead(Xseq, &array));
  }

  PetscCall(VecDestroy(&naturalX));


  /* Destroy DMDAs and Vecs */
  PetscCall(VecDestroy( & vecV));
  PetscCall(DMDestroy( & daV));
  PetscCall(VecDestroy( & V));
  PetscCall(DMDestroy( & dmV));


  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

int isInDomain(const double * R, const double * Z, void * ptr){
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("isInDomain",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  int value = -3;
  int loopbr = 0;
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, d, N[3], er, ephi, ez;
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  DM dmCoorda, coordDA = user -> coorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoord;

  PetscCall(DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
  PetscCall(DMGetCoordinateDM(coordDA, & dmCoorda));
  PetscCall(DMGetCoordinatesLocal(coordDA, & coordaLocal));
  PetscCall(DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoord));
  for (d = 0; d < 3; ++d) {
    /* Face coordinates */
    PetscCall(DMStagGetLocationSlot(dmCoorda, LEFT, d, & icBrm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, DOWN, d, & icBphim[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, BACK, d, & icBzm[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, RIGHT, d, & icBrp[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, UP, d, & icBphip[d]));
    PetscCall(DMStagGetLocationSlot(dmCoorda, FRONT, d, & icBzp[d]));
  }
  PetscCall(DMStagGetGlobalSizes(user -> coorda, & N[0], & N[1], & N[2]));

  for (ez = startz; ez < startz + nz; ++ez) {
    if(loopbr){
      break;
    }
    ephi = 0;
    for (er = startr; er < startr + nr; ++er) {
      if((*Z >= arrCoord[ez][ephi][er][icBzm[2]]) && (*Z <= arrCoord[ez][ephi][er][icBzp[2]]) && (*R >= arrCoord[ez][ephi][er][icBrm[0]]) && (*R <= arrCoord[ez][ephi][er][icBrp[0]])) {
        value = (int) (user->dataC[er + ez * N[1] * N[0]]);
        loopbr = 1;
        break;
      }
    }
  }

  if(value == -3){
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "cell indices not found for (R,Z) = (%f,%f)\n", *R,*Z));
  }
  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  return value;
}
