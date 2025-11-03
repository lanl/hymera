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

#include <mfd_config.h>

#include <ts_functions.h>

#include <monitor_functions.h>

#include <geometry.h>

#include <mass_matrix_coefficients.h>

#include <mimetic_operators.h>

PetscErrorCode FormMaterialPropertiesMatrix(TS ts, Vec MPdiag, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec tempMPdiag;
  Mat MP;
  PetscInt N[3], er, ephi, ez, d;

  PetscPrintf(PETSC_COMM_WORLD, "Before Creating MP\n");
  DMCreateMatrix(coordDA, & MP);
  MatZeroEntries(MP);
  DMCreateGlobalVector(coordDA, & tempMPdiag);
  VecZeroEntries(tempMPdiag);

  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col;
        PetscScalar val;
        PetscInt nEntries = 1;

        /* The edges are oriented in the directions of unit vectors: e_r, e_phi and e_z */

        /* Cell material properties matrix: MP */
        /* Equation on element */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = ELEMENT;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = ELEMENT;
        col.c = 0;

        if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {val = user -> etaplasma ;}
        else if((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12)) {val = user -> etawall ;}
        else if((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12)) {val = user -> etawall ;}
        else{val = user -> etaout ;}
        val = val / user->mu0;
        DMStagMatSetValuesStencil(coordDA, MP, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);

      }
    }
  }

  MatAssemblyBegin(MP, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(MP, MAT_FINAL_ASSEMBLY);
  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "MP matrix: \n");
    MatView(MP, PETSC_VIEWER_STDOUT_WORLD);
  }

  MatGetDiagonal(MP, tempMPdiag);
  VecCopy(tempMPdiag, MPdiag);

  PetscBarrier((PetscObject)MP);
  PetscBarrier((PetscObject)tempMPdiag);

  MatDestroy( & MP);
  VecDestroy( & tempMPdiag);
  return (0);
}


PetscErrorCode FormFaceMassMatrix(TS ts, Vec Mfdiag, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec tempMfdiag, tempInteriorIndicatordiag;
  Mat Mf, tempMf, InteriorIndicator, tempInteriorIndicator;
  PetscInt N[3], er, ephi, ez, d;

  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  DMCreateMatrix(da, & Mf);
  MatZeroEntries(Mf);

  DMCreateMatrix(da, & InteriorIndicator);
  MatZeroEntries(InteriorIndicator);

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col;
        PetscScalar val;
        PetscInt nEntries = 1;

        /* The edges are oriented in the directions of unit vectors: e_r, e_phi and e_z */

        /* Face mass matrix: Mf */
        /* Equation on left face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = LEFT;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = LEFT;
        col.c = 0;
        val = betaf_wmp(er, ephi, ez, LEFT, user);
        DMStagMatSetValuesStencil(da, Mf, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
          val = 1.0;
          DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        if (er == N[0] - 1) {
          /* Equation on right face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = RIGHT;
          col.c = 0;
          val = betaf_wmp(er, ephi, ez, RIGHT, user);
          DMStagMatSetValuesStencil(da, Mf, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
        }
        /* Equation on back face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK;
          col.c = 0;
          val = betaf_wmp(er, ephi, ez, BACK, user);
          DMStagMatSetValuesStencil(da, Mf, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        if (ez == N[2] - 1) {
          /* Equation on front face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT;
          col.c = 0;
          val = betaf_wmp(er, ephi, ez, FRONT, user);
          DMStagMatSetValuesStencil(da, Mf, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
        }
        /* Equation on down face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = DOWN;
        col.c = 0;
        val = betaf_wmp(er, ephi, ez, DOWN, user);
        DMStagMatSetValuesStencil(da, Mf, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
          val = 1.0;
          DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        if (ephi == N[1] - 1) {
          /* Equation on up face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP;
          col.c = 0;
          val = betaf_wmp(er, ephi, ez, UP, user);
          DMStagMatSetValuesStencil(da, Mf, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
        }
      }
    }
  }

  MatAssemblyBegin(Mf, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Mf, MAT_FINAL_ASSEMBLY);
  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Mf matrix: \n");
    MatView(Mf, PETSC_VIEWER_STDOUT_WORLD);
  }

  MatCreateSubMatrix(Mf, user -> isB, user -> isB, MAT_INITIAL_MATRIX, & tempMf);
  MatCreateVecs(tempMf, & tempMfdiag, NULL);

  PetscErrorCode ierr = 0;
  PetscPrintf(PETSC_COMM_SELF, "Before MatGetDiagonal matrix: \n");
  ierr = MatGetDiagonal(tempMf, tempMfdiag);
  CHKERRQ(ierr);

  VecDuplicate(tempMfdiag, & Mfdiag);
  VecCopy(tempMfdiag, Mfdiag);
  PetscBarrier((PetscObject)tempMf);
  MatDestroy( & tempMf);
  MatDestroy( & Mf);
  VecDestroy( & tempMfdiag);

  if(0){
    Vec dummysol, dummyB1, dummyX;
    DMCreateGlobalVector(da, & dummyX);
    VecZeroEntries(dummyX);
    VecDuplicate(dummyX, & dummysol);
    VecZeroEntries(dummysol);
    VecGetSubVector(dummysol, user->isB, & dummyB1);
    VecCopy(Mfdiag, dummyB1);
    VecRestoreSubVector(dummysol, user->isB, & dummyB1);

    DumpError(ts, 7, dummysol, user);
    PetscBarrier((PetscObject)dummysol);
    VecDestroy( & dummysol);
    return(0);
  }

  MatAssemblyBegin(InteriorIndicator, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(InteriorIndicator, MAT_FINAL_ASSEMBLY);
  MatCreateSubMatrix(InteriorIndicator, user -> isB, user -> isB, MAT_INITIAL_MATRIX, & tempInteriorIndicator);
  MatCreateVecs(tempInteriorIndicator, & tempInteriorIndicatordiag, NULL);
  ierr = MatGetDiagonal(tempInteriorIndicator, tempInteriorIndicatordiag);
  CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "======Mfdiag Vector======\n");
  VecView(Mfdiag, PETSC_VIEWER_STDOUT_WORLD);

  PetscViewer Matviewer;
  if (user -> etawall == 0.000001) {
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Mf_eta_1_6.m", & Matviewer);
  } else if (user -> etawall == 0.001) {
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Mf_eta_1_3.m", & Matviewer);
  } else {
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Mf_eta_1_inf.m", & Matviewer);
  }
  PetscViewerPushFormat(Matviewer, PETSC_VIEWER_ASCII_MATLAB);
  VecView(Mfdiag, Matviewer);
  PetscViewerPopFormat(Matviewer);

  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "InteriorfIndicator.m", & Matviewer);
  PetscViewerPushFormat(Matviewer, PETSC_VIEWER_ASCII_MATLAB);
  VecView(tempInteriorIndicatordiag, Matviewer);
  PetscViewerPopFormat(Matviewer);

  PetscBarrier((PetscObject)tempInteriorIndicator);
  MatDestroy( & tempInteriorIndicator);
  MatDestroy( & InteriorIndicator);
  VecDestroy( & tempInteriorIndicatordiag);

  PetscViewerDestroy( & Matviewer);
  PetscBarrier((PetscObject)Mfdiag);
  VecDestroy( & Mfdiag);
  return (0);
}


PetscErrorCode FormEdgeMassMatrix(TS ts, Vec Mediag, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec tempMediag, tempInteriorIndicatordiag;
  Mat Me, tempMe, InteriorIndicator, tempInteriorIndicator;
  PetscInt N[3], er, ephi, ez, d;

  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  DMCreateMatrix(da, & Me);
  MatZeroEntries(Me);

  DMCreateMatrix(da, & InteriorIndicator);
  MatZeroEntries(InteriorIndicator);

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col;
        PetscScalar val;
        PetscInt nEntries = 1;

        /* The edges are oriented in the directions of unit vectors: e_r, e_phi and e_z */

        /* Edge mass matrix: Me */
        /* Equation on down left edge */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN_LEFT;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = DOWN_LEFT;
        col.c = 0;
        val = betae(er, ephi, ez, DOWN_LEFT, user);
        DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
          val = 1.0;
          DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        if (er == N[0] - 1) {
          /* Equation on down right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_RIGHT;
          col.c = 0;
          val = betae(er, ephi, ez, DOWN_RIGHT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
          /* Equation on back right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_RIGHT;
          col.c = 0;
          val = betae(er, ephi, ez, BACK_RIGHT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
        }
        /* Equation on back left edge */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_LEFT;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = BACK_LEFT;
        col.c = 0;
        val = betae(er, ephi, ez, BACK_LEFT, user);
        DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
          val = 1.0;
          DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        if (ez == N[2] - 1) {
          /* Equation on front left edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_LEFT;
          col.c = 0;
          val = betae(er, ephi, ez, FRONT_LEFT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
          /* Equation on front down edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_DOWN;
          col.c = 0;
          val = betae(er, ephi, ez, FRONT_DOWN, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
        }
        /* Equation on back down edge */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_DOWN;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = BACK_DOWN;
        col.c = 0;
        val = betae(er, ephi, ez, BACK_DOWN, user);
        DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
          val = 1.0;
          DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        if (ephi == N[1] - 1) {
          /* Equation on back up edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_UP;
          col.c = 0;
          val = betae(er, ephi, ez, BACK_UP, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
          /* Equation on up left edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_LEFT;
          col.c = 0;
          val = betae(er, ephi, ez, UP_LEFT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
        }
        if (er == N[0] - 1 && ephi == N[1] - 1) {
          /* Equation on up right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_RIGHT;
          col.c = 0;
          val = betae(er, ephi, ez, UP_RIGHT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
        }
        if (ephi == N[1] - 1 && ez == N[2] - 1) {
          /* Equation on front up edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_UP;
          col.c = 0;
          val = betae(er, ephi, ez, FRONT_UP, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Equation on front right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_RIGHT;
          col.c = 0;
          val = betae(er, ephi, ez, FRONT_RIGHT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if ((ephi > -1 && ephi < N[1] && fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0)< 1e-12) || (ephi == -1 && fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) || (ephi == N[1] && fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12)) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, InteriorIndicator, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          }
        }
      }
    }
  }

  MatAssemblyBegin(Me, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Me, MAT_FINAL_ASSEMBLY);
  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Me matrix: \n");
    MatView(Me, PETSC_VIEWER_STDOUT_WORLD);
  }

  MatCreateSubMatrix(Me, user -> istau, user -> istau, MAT_INITIAL_MATRIX, & tempMe);
  MatCreateVecs(tempMe, & tempMediag, NULL);

  PetscErrorCode ierr = 0;
  PetscPrintf(PETSC_COMM_SELF, "Before MatGetDiagonal matrix: \n");
  ierr = MatGetDiagonal(tempMe, tempMediag);
  CHKERRQ(ierr);

  VecDuplicate(tempMediag, & Mediag);
  VecCopy(tempMediag, Mediag);
  PetscBarrier((PetscObject)tempMe);
  MatDestroy( & tempMe);
  MatDestroy( & Me);
  VecDestroy( & tempMediag);

  MatAssemblyBegin(InteriorIndicator, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(InteriorIndicator, MAT_FINAL_ASSEMBLY);
  MatCreateSubMatrix(InteriorIndicator, user -> istau, user -> istau, MAT_INITIAL_MATRIX, & tempInteriorIndicator);
  MatCreateVecs(tempInteriorIndicator, & tempInteriorIndicatordiag, NULL);
  ierr = MatGetDiagonal(tempInteriorIndicator, tempInteriorIndicatordiag);
  CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "======Mediag Vector======\n");
  VecView(Mediag, PETSC_VIEWER_STDOUT_WORLD);

  PetscViewer Matviewer;
  if (user -> etawall == 0.000001) {
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Me_eta_1_6.m", & Matviewer);
  } else if (user -> etawall == 0.001) {
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Me_eta_1_3.m", & Matviewer);
  } else {
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, "Me_eta_1_inf.m", & Matviewer);
  }
  PetscViewerPushFormat(Matviewer, PETSC_VIEWER_ASCII_MATLAB);
  VecView(Mediag, Matviewer);
  PetscViewerPopFormat(Matviewer);

  PetscViewerASCIIOpen(PETSC_COMM_WORLD, "InteriorIndicator.m", & Matviewer);
  PetscViewerPushFormat(Matviewer, PETSC_VIEWER_ASCII_MATLAB);
  VecView(tempInteriorIndicatordiag, Matviewer);
  PetscViewerPopFormat(Matviewer);

  PetscBarrier((PetscObject)tempInteriorIndicator);
  MatDestroy( & tempInteriorIndicator);
  MatDestroy( & InteriorIndicator);
  VecDestroy( & tempInteriorIndicatordiag);

  PetscViewerDestroy( & Matviewer);
  PetscBarrier((PetscObject)Mediag);
  VecDestroy( & Mediag);
  return (0);
}

PetscErrorCode FormDiscreteDivergence(TS ts, DM newda, Mat D, Vec X, Vec F, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec coordLocal, fLocal, xLocal;
  PetscInt N[3], er, ephi, ez, d;
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm, ivC;
  DM dmCoord;
  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX;

  VecZeroEntries(F);
  MatZeroEntries(D);
  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  DMGetLocalVector(newda, & fLocal);
  DMStagVecGetArray(newda, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);

  /* Cell locations */
  DMStagGetLocationSlot(newda, ELEMENT, 0, & ivC);

  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);

  for (d = 0; d < 3; ++d) {
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
  }

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col[6];
        PetscScalar cellvolume, valD[6];
        PetscInt nEntries = 6;

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
        } else {
          cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
            PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
        }

        arrF[ez][ephi][er][ivC] = (-surface(er, ephi, ez, LEFT, user) / cellvolume) * arrX[ez][ephi][er][ivBrm] + (surface(er, ephi, ez, RIGHT, user) / cellvolume) * arrX[ez][ephi][er][ivBrp] + (-surface(er, ephi, ez, DOWN, user) / cellvolume) * arrX[ez][ephi][er][ivBphim] + (surface(er, ephi, ez, UP, user) / cellvolume) * arrX[ez][ephi][er][ivBphip] + (-surface(er, ephi, ez, BACK, user) / cellvolume) * arrX[ez][ephi][er][ivBzm] + (surface(er, ephi, ez, FRONT, user) / cellvolume) * arrX[ez][ephi][er][ivBzp];

        if (user -> debug) {
          PetscPrintf(PETSC_COMM_WORLD, "div(B,%d,%d,%d)= %g\n", er, ephi, ez, arrF[ez][ephi][er][ivC]);
          PetscPrintf(PETSC_COMM_WORLD, "1st term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double)(-surface(er, ephi, ez, LEFT, user) / cellvolume) * arrX[ez][ephi][er][ivBrm]);
          PetscPrintf(PETSC_COMM_WORLD, "B in 1st term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double) arrX[ez][ephi][er][ivBrm]);
          PetscPrintf(PETSC_COMM_WORLD, "div coefficient in 1st term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double)(-surface(er, ephi, ez, LEFT, user) / cellvolume));
          PetscPrintf(PETSC_COMM_WORLD, "2nd term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double)(surface(er, ephi, ez, RIGHT, user) / cellvolume) * arrX[ez][ephi][er][ivBrp]);
          PetscPrintf(PETSC_COMM_WORLD, "B in 2nd term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double) arrX[ez][ephi][er][ivBrp]);
          PetscPrintf(PETSC_COMM_WORLD, "div coefficient in 2nd term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double)(surface(er, ephi, ez, RIGHT, user) / cellvolume));

          PetscPrintf(PETSC_COMM_WORLD, "3rd term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double)(-surface(er, ephi, ez, DOWN, user) / cellvolume) * arrX[ez][ephi][er][ivBphim]);
          PetscPrintf(PETSC_COMM_WORLD, "4th term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double)(surface(er, ephi, ez, UP, user) / cellvolume) * arrX[ez][ephi][er][ivBphip]);
          PetscPrintf(PETSC_COMM_WORLD, "5th term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double)(-surface(er, ephi, ez, BACK, user) / cellvolume) * arrX[ez][ephi][er][ivBzm]);
          PetscPrintf(PETSC_COMM_WORLD, "6th term in div(B,%d,%d,%d)= %g\n", er, ephi, ez, (double)(surface(er, ephi, ez, FRONT, user) / cellvolume) * arrX[ez][ephi][er][ivBzp]);
        }

        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = ELEMENT;
        row.c = 0;
        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = FRONT;
        col[0].c = 0;
        valD[0] = surface(er, ephi, ez, FRONT, user) / cellvolume;
        col[1].i = er;
        col[1].j = ephi;
        col[1].k = ez;
        col[1].loc = BACK;
        col[1].c = 0;
        valD[1] = -surface(er, ephi, ez, BACK, user) / cellvolume;
        col[2].i = er;
        col[2].j = ephi;
        col[2].k = ez;
        col[2].loc = UP;
        col[2].c = 0;
        valD[2] = surface(er, ephi, ez, UP, user) / cellvolume;
        col[3].i = er;
        col[3].j = ephi;
        col[3].k = ez;
        col[3].loc = DOWN;
        col[3].c = 0;
        valD[3] = -surface(er, ephi, ez, DOWN, user) / cellvolume;
        col[4].i = er;
        col[4].j = ephi;
        col[4].k = ez;
        col[4].loc = RIGHT;
        col[4].c = 0;
        valD[4] = surface(er, ephi, ez, RIGHT, user) / cellvolume;
        col[5].i = er;
        col[5].j = ephi;
        col[5].k = ez;
        col[5].loc = LEFT;
        col[5].c = 0;
        valD[5] = -surface(er, ephi, ez, LEFT, user) / cellvolume;
        DMStagMatSetValuesStencil(coordDA, D, 1, & row, nEntries, col, valD, INSERT_VALUES);

      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);
  DMStagVecRestoreArray(newda, fLocal, & arrF);
  DMLocalToGlobal(newda, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(newda, & fLocal);

  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);
  MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);
  return (0);
}

PetscErrorCode FormGradDerivedDivergence(TS ts, Mat GD, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Mat tempGD, extGD, D, G, dummymat;
  Vec dummyX;
  PetscInt N[3], er, ephi, ez;
  KSP dummyksp;
  PC dummypc;
  IS isBE; /* indexing PETSc object for edge and face dofs */
  IS isOthers; /* indexing PETSc object for other dofs */

  TSGetDM(ts, & da);
  //MatCreate(PETSC_COMM_WORLD,&G);
  DMCreateMatrix(coordDA, & G);

  PetscPrintf(PETSC_COMM_WORLD, "Before FormDiscreteGradient\n");
  FormDiscreteGradient(ts, G, user);
  PetscPrintf(PETSC_COMM_WORLD, "After FormDiscreteGradient\n");

  //MatCreate(PETSC_COMM_WORLD,&D);
  DMCreateMatrix(coordDA, & D);

  PetscPrintf(PETSC_COMM_WORLD, "Before FormDerivedDivergence\n");
  FormDerivedDivergence(ts, D, user);
  PetscPrintf(PETSC_COMM_WORLD, "After FormDerivedDivergence\n");

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Derived_divergence matrix D:\n");
    MatView(D, PETSC_VIEWER_STDOUT_WORLD);
  }

  //MatCreate(PETSC_COMM_WORLD,&extGD);
  DMCreateMatrix(coordDA, & extGD);
  MatMatMult(G, D, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & extGD);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======extGD = G*D on coordDA======\n");
    MatView(extGD, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscBarrier((PetscObject)G);
  PetscBarrier((PetscObject)D);
  MatDestroy( & G);
  MatDestroy( & D);

  DMCreateGlobalVector(coordDA, & dummyX);
  VecSet(dummyX, 0.0);
  KSPCreate(PETSC_COMM_WORLD, & dummyksp);

  DMCreateMatrix(coordDA, & dummymat);
  //MatAssemblyBegin(dummymat,MAT_FINAL_ASSEMBLY);
  //MatAssemblyEnd(dummymat,MAT_FINAL_ASSEMBLY);
  //MatConvert(extGD, MATSAME, MAT_INITIAL_MATRIX, &dummymat);
  //MatCopy(extGD, dummymat, DIFFERENT_NONZERO_PATTERN);
  //MatDuplicate(extGD,MAT_COPY_VALUES,&dummymat);

  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col;
        PetscScalar val;
        PetscInt nEntries = 1;

        /* B field part */

        /* Equation on left face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = LEFT;
        row.c = 0;

        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = LEFT;
        col.c = 0;
        val = 1.0;

        DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

        /* Equation on down face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN;
        row.c = 0;

        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = DOWN;
        col.c = 0;
        val = 1.0;

        DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

        /* Equation on back face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK;
        row.c = 0;

        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = BACK;
        col.c = 0;
        val = 1.0;

        DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

        if (er == N[0] - 1) {
          /* Equation on right boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on up boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on front boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
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

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && (er == 0)) || (((er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on the left or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && (ez == 0)) || (((ephi == 0 || ez == 0)) && !(user -> phibtype))) {
          /* Equation on the left or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_DOWN;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (er == N[0] - 1) {
          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_UP;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_DOWN;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && ez == N[2] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_UP;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Equation on the front and right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && !(er == 0)) || ((!(er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal down left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && !(ez == 0)) || ((!(ez == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal back down edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_DOWN;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (!(er == 0 || ez == 0)) {
          /* Equation on internal back left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

      }
    }
  }

  MatAssemblyBegin(dummymat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dummymat, MAT_FINAL_ASSEMBLY);

  KSPSetOptionsPrefix(dummyksp, "dummy_");
  KSPSetOperators(dummyksp, dummymat, dummymat);
  KSPGetPC(dummyksp, & dummypc);
  PCSetType(dummypc, PCFIELDSPLIT);
  //PCFieldSplitSchurPrecondition(dummypc,PC_FIELDSPLIT_SCHUR_PRE_SELF,NULL);
  PCFieldSplitSetDetectSaddlePoint(dummypc, PETSC_TRUE);
  PCFieldSplitSetSchurPre(dummypc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
  KSPSetFromOptions(dummyksp);

  KSPSetUp(dummyksp);
  KSPSetTolerances(dummyksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 iteration for the dummy solve

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======dummymat Matrix======\n");
    MatView(dummymat, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "======dummyX Vector======\n");
    VecView(dummyX, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscPrintf(PETSC_COMM_WORLD, "======Before KSPSolve dummyksp ======\n");
  KSPSolve(dummyksp, dummyX, dummyX);
  PetscPrintf(PETSC_COMM_WORLD, "======After KSPSolve dummyksp ======\n");
  PetscBarrier((PetscObject)dummyX);
  PetscBarrier((PetscObject)dummymat);
  VecDestroy( & dummyX);
  MatDestroy( & dummymat);
  PCFieldSplitGetISByIndex(dummypc, 0, & isBE);
  PCFieldSplitGetISByIndex(dummypc, 1, & isOthers);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======edge & face dofs======\n");
    ISView(isBE, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "======other dofs======\n");
    ISView(isOthers, PETSC_VIEWER_STDOUT_WORLD);
  }

  if (user -> etawall == 1 && user -> etaplasma == 1) {
    PetscPrintf(PETSC_COMM_WORLD, "======Before MatGetSubVector DiagMe1======\n");
    PetscErrorCode ierr = 0;
    ierr = VecGetSubVector(user -> DiagMe1, isBE, & (user -> MeVec1));
    CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "======After MatGetSubVector DiagMe1======\n");
    VecAssemblyBegin(user -> MeVec1);
    VecAssemblyEnd(user -> MeVec1);
  }
  PetscPrintf(PETSC_COMM_WORLD, "======Before MatGetSubVector======\n");
  PetscErrorCode ierr = 0;
  ierr = VecGetSubVector(user -> DiagMe, isBE, & (user -> MeVec));
  CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "======After MatGetSubVector======\n");
  VecAssemblyBegin(user -> MeVec);
  VecAssemblyEnd(user -> MeVec);

  DMCreateMatrix(da, & tempGD);

  PetscPrintf(PETSC_COMM_WORLD, "======Before MatCreateSubMatrix======\n");
  MatCreateSubMatrix(extGD, isBE, isBE, MAT_INITIAL_MATRIX, & tempGD);
  PetscPrintf(PETSC_COMM_WORLD, "======After MatCreateSubMatrix======\n");
  MatAssemblyBegin(tempGD, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(tempGD, MAT_FINAL_ASSEMBLY);

  MatCopy(tempGD, GD, DIFFERENT_NONZERO_PATTERN);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======GD = MatCreateSubMatrix(extGD)======\n");
    MatView(GD, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscBarrier((PetscObject)extGD);
  PetscBarrier((PetscObject)tempGD);
  PetscBarrier((PetscObject)dummyksp);
  MatDestroy( & extGD);
  MatDestroy( & tempGD);
  KSPDestroy( & dummyksp);
  return (0);
}

PetscErrorCode FormDerivedGradDivergence(TS ts, Mat GD, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda, newda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Mat tempGD, extGD, D, G, dummymat;
  Vec dummyX, div, X;
  PetscInt N[3], er, ephi, ez, procr, procphi, procz;
  KSP dummyksp;
  PC dummypc;
  IS isBE; /* indexing PETSc object for edge and face dofs */
  IS isOthers; /* indexing PETSc object for other dofs */
  PetscErrorCode ierr = 0;

  TSGetDM(ts, & da);
  DMStagGetNumRanks(da, & procr, & procphi, & procz);
  if (user -> phibtype) {
    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE, user -> Nr, user -> Nphi, user -> Nz, procr, procphi, procz, 0, 0, 0, 1, DMSTAG_STENCIL_BOX, 1, NULL, NULL, NULL, & newda);
  } else {
    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, user -> Nr, user -> Nphi, user -> Nz, procr, procphi, procz, 0, 0, 0, 1, DMSTAG_STENCIL_BOX, 1, NULL, NULL, NULL, & newda);
  }
  DMSetFromOptions(newda);
  DMSetUp(newda);
  DMStagSetUniformCoordinatesExplicit(newda, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);

  ierr = DMCreateGlobalVector(newda, & div);
  CHKERRQ(ierr);
  DMCreateMatrix(user -> coorda, & D);
  DMCreateGlobalVector(da, & X);
  VecSet(X, 0.0);
  PetscPrintf(PETSC_COMM_WORLD, "Before FormDiscreteDivergence\n");
  FormDiscreteDivergence(ts, newda, D, X, div, user);
  PetscPrintf(PETSC_COMM_WORLD, "After FormDiscreteDivergence\n");
  PetscBarrier((PetscObject)newda);
  DMDestroy( & newda);
  PetscBarrier((PetscObject)div);
  VecDestroy( & div);
  PetscBarrier((PetscObject)X);
  VecDestroy( & X);
  DMCreateMatrix(user -> coorda, & G);
  MatZeroEntries(G);
  MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);
  PetscPrintf(PETSC_COMM_WORLD, "Before FormDerivedGradient\n");
  FormDerivedGradient(ts, D, G, user);
  PetscPrintf(PETSC_COMM_WORLD, "After FormDerivedGradient\n");

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Derived Gradient matrix G:\n");
    MatView(G, PETSC_VIEWER_STDOUT_WORLD);
  }

  //MatCreate(PETSC_COMM_WORLD,&extGD);
  DMCreateMatrix(user -> coorda, & extGD);

  if(0){/* Cancelling the image of the primary divergence on boundary cells */
   IS isC_boundary;
   PetscErrorCode ierr = 0;
   ComputeIsEBoundary(ts, & isC_boundary, user);
   ierr = MatZeroRowsIS(D, isC_boundary, 0.0, NULL, NULL); CHKERRQ(ierr);
  }

  MatMatMult(G, D, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & extGD);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======extGD = G*D on coordDA======\n");
    MatView(extGD, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscBarrier((PetscObject)G);
  PetscBarrier((PetscObject)D);
  MatDestroy( & G);
  MatDestroy( & D);

  DMCreateGlobalVector(user -> coorda, & dummyX);
  VecSet(dummyX, 0.0);
  KSPCreate(PETSC_COMM_WORLD, & dummyksp);

  DMCreateMatrix(user -> coorda, & dummymat);

  DMStagGetGlobalSizes(user -> coorda, & N[0], & N[1], & N[2]);
  DMStagGetCorners(user -> coorda, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col;
        PetscScalar val;
        PetscInt nEntries = 1;

        /* B field part */

        /* Equation on left face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = LEFT;
        row.c = 0;

        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = LEFT;
        col.c = 0;
        val = 1.0;

        DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

        /* Equation on down face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN;
        row.c = 0;

        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = DOWN;
        col.c = 0;
        val = 1.0;

        DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

        /* Equation on back face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK;
        row.c = 0;

        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = BACK;
        col.c = 0;
        val = 1.0;

        DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

        if (er == N[0] - 1) {
          /* Equation on right boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on up boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on front boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
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

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && (er == 0)) || (((er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on the left or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && (ez == 0)) || (((ephi == 0 || ez == 0)) && !(user -> phibtype))) {
          /* Equation on the left or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_DOWN;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (er == N[0] - 1) {
          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_UP;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_DOWN;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && ez == N[2] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_UP;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Equation on the front and right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && !(er == 0)) || ((!(er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal down left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && !(ez == 0)) || ((!(ez == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal back down edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_DOWN;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (!(er == 0 || ez == 0)) {
          /* Equation on internal back left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

      }
    }
  }

  MatAssemblyBegin(dummymat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dummymat, MAT_FINAL_ASSEMBLY);

  KSPSetOptionsPrefix(dummyksp, "dummy_");
  KSPSetOperators(dummyksp, dummymat, dummymat);
  KSPGetPC(dummyksp, & dummypc);
  PCSetType(dummypc, PCFIELDSPLIT);
  //PCFieldSplitSchurPrecondition(dummypc,PC_FIELDSPLIT_SCHUR_PRE_SELF,NULL);
  PCFieldSplitSetDetectSaddlePoint(dummypc, PETSC_TRUE);
  PCFieldSplitSetSchurPre(dummypc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
  KSPSetFromOptions(dummyksp);

  KSPSetUp(dummyksp);
  KSPSetTolerances(dummyksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 iteration for the dummy solve

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======dummymat Matrix======\n");
    MatView(dummymat, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "======dummyX Vector======\n");
    VecView(dummyX, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscPrintf(PETSC_COMM_WORLD, "======Before KSPSolve dummyksp ======\n");
  KSPSolve(dummyksp, dummyX, dummyX);
  PetscPrintf(PETSC_COMM_WORLD, "======After KSPSolve dummyksp ======\n");
  PetscBarrier((PetscObject)dummyX);
  PetscBarrier((PetscObject)dummymat);
  VecDestroy( & dummyX);
  MatDestroy( & dummymat);
  PCFieldSplitGetISByIndex(dummypc, 0, & isBE);
  PCFieldSplitGetISByIndex(dummypc, 1, & isOthers);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======edge & face dofs======\n");
    ISView(isBE, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "======other dofs======\n");
    ISView(isOthers, PETSC_VIEWER_STDOUT_WORLD);
  }

  DMCreateMatrix(da, & tempGD);

  PetscPrintf(PETSC_COMM_WORLD, "======Before MatCreateSubMatrix======\n");
  MatCreateSubMatrix(extGD, isBE, isBE, MAT_INITIAL_MATRIX, & tempGD);
  PetscPrintf(PETSC_COMM_WORLD, "======After MatCreateSubMatrix======\n");
  MatAssemblyBegin(tempGD, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(tempGD, MAT_FINAL_ASSEMBLY);

  MatCopy(tempGD, GD, DIFFERENT_NONZERO_PATTERN);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======GD = MatCreateSubMatrix(extGD)======\n");
    MatView(GD, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscBarrier((PetscObject)extGD);
  PetscBarrier((PetscObject)tempGD);
  PetscBarrier((PetscObject)dummyksp);
  MatDestroy( & extGD);
  MatDestroy( & tempGD);
  KSPDestroy( & dummyksp);
  return (0);
}

PetscErrorCode FormDerivedGradEtaDivergence(TS ts, Mat GD, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda, newda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Mat tempGD, extGD, D, G, dummymat;
  Vec dummyX, div, X, MPdiag;
  PetscInt N[3], er, ephi, ez, procr, procphi, procz;
  KSP dummyksp;
  PC dummypc;
  IS isBE; /* indexing PETSc object for edge and face dofs */
  IS isOthers; /* indexing PETSc object for other dofs */
  PetscErrorCode ierr = 0;

  TSGetDM(ts, & da);
  DMStagGetNumRanks(da, & procr, & procphi, & procz);
  if (user -> phibtype) {
    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE, user -> Nr, user -> Nphi, user -> Nz, procr, procphi, procz, 0, 0, 0, 1, DMSTAG_STENCIL_BOX, 1, NULL, NULL, NULL, & newda);
  } else {
    DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, user -> Nr, user -> Nphi, user -> Nz, procr, procphi, procz, 0, 0, 0, 1, DMSTAG_STENCIL_BOX, 1, NULL, NULL, NULL, & newda);
  }
  DMSetFromOptions(newda);
  DMSetUp(newda);
  DMStagSetUniformCoordinatesExplicit(newda, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);

  ierr = DMCreateGlobalVector(newda, & div);
  CHKERRQ(ierr);
  DMCreateMatrix(user -> coorda, & D);
  DMCreateGlobalVector(da, & X);
  VecSet(X, 0.0);
  PetscPrintf(PETSC_COMM_WORLD, "Before FormDiscreteDivergence\n");
  FormDiscreteDivergence(ts, newda, D, X, div, user);
  PetscPrintf(PETSC_COMM_WORLD, "After FormDiscreteDivergence\n");
  PetscBarrier((PetscObject)newda);
  DMDestroy( & newda);
  PetscBarrier((PetscObject)div);
  VecDestroy( & div);
  PetscBarrier((PetscObject)X);
  VecDestroy( & X);
  DMCreateMatrix(user -> coorda, & G);
  MatZeroEntries(G);
  MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);

  // Store previous resistivities
  PetscReal plasma_Resist, out_Resist, wall_Resist, permeability;
  permeability = user -> mu0;
  plasma_Resist = user -> etaplasma;
  wall_Resist = user -> etawall;
  out_Resist = user -> etaout;
  // Change resistivities in user context to one
  user -> mu0 = 1.0;
  user -> etaplasma = 1.0;
  user -> etawall = 1.0;
  user -> etaout = 1.0;

  PetscPrintf(PETSC_COMM_WORLD, "Before FormDerivedGradient\n");
  FormDerivedGradient(ts, D, G, user);
  PetscPrintf(PETSC_COMM_WORLD, "After FormDerivedGradient\n");

  // Switch back to original resistivities
  user -> mu0 = permeability;
  user -> etaplasma = plasma_Resist;
  user -> etawall = wall_Resist;
  user -> etaout = out_Resist;

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Derived Gradient matrix G:\n");
    MatView(G, PETSC_VIEWER_STDOUT_WORLD);
  }

  DMCreateGlobalVector(coordDA, & MPdiag);
  VecZeroEntries(MPdiag);
  FormMaterialPropertiesMatrix(ts, MPdiag, user);

  DMCreateMatrix(user -> coorda, & extGD);

  if(0){/* Cancelling the image of the primary divergence on boundary cells */
   IS isC_boundary;
   PetscErrorCode ierr = 0;
   ComputeIsEBoundary(ts, & isC_boundary, user);
   ierr = MatZeroRowsIS(D, isC_boundary, 0.0, NULL, NULL); CHKERRQ(ierr);
  }

  MatDiagonalScale( D, MPdiag, NULL); // D := MP * D
  MatMatMult(G, D, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & extGD); // extGD := G * MP * D

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======extGD = G*D on coordDA======\n");
    MatView(extGD, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscBarrier((PetscObject)G);
  PetscBarrier((PetscObject)D);
  MatDestroy( & G);
  MatDestroy( & D);
  VecDestroy( & MPdiag);

  DMCreateGlobalVector(user -> coorda, & dummyX);
  VecSet(dummyX, 0.0);
  KSPCreate(PETSC_COMM_WORLD, & dummyksp);

  DMCreateMatrix(user -> coorda, & dummymat);

  DMStagGetGlobalSizes(user -> coorda, & N[0], & N[1], & N[2]);
  DMStagGetCorners(user -> coorda, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col;
        PetscScalar val;
        PetscInt nEntries = 1;

        /* B field part */

        /* Equation on left face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = LEFT;
        row.c = 0;

        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = LEFT;
        col.c = 0;
        val = 1.0;

        DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

        /* Equation on down face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN;
        row.c = 0;

        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = DOWN;
        col.c = 0;
        val = 1.0;

        DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

        /* Equation on back face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK;
        row.c = 0;

        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = BACK;
        col.c = 0;
        val = 1.0;

        DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

        if (er == N[0] - 1) {
          /* Equation on right boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on up boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on front boundary face */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
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

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && (er == 0)) || (((er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on the left or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && (ez == 0)) || (((ephi == 0 || ez == 0)) && !(user -> phibtype))) {
          /* Equation on the left or down boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_DOWN;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (er == N[0] - 1) {
          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

          /* Equation on the right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_UP;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

          /* Equation on the up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);

          /* Equation on the front boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_DOWN;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && ez == N[2] - 1 && !(user -> phibtype)) {
          /* Equation on the right and up boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_UP;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Equation on the front and right boundary */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_RIGHT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && !(er == 0)) || ((!(er == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal down left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if ((user -> phibtype && !(ez == 0)) || ((!(ez == 0 || ephi == 0)) && !(user -> phibtype))) {
          /* Equation on internal back down edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_DOWN;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        if (!(er == 0 || ez == 0)) {
          /* Equation on internal back left edge */
          nEntries = 1;
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_LEFT;
          row.c = 0;

          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_LEFT;
          col.c = 0;
          val = 1.0;

          DMStagMatSetValuesStencil(coordDA, dummymat, 1, & row, nEntries, & col, & val, INSERT_VALUES);
        }

      }
    }
  }

  MatAssemblyBegin(dummymat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(dummymat, MAT_FINAL_ASSEMBLY);

  KSPSetOptionsPrefix(dummyksp, "dummy_");
  KSPSetOperators(dummyksp, dummymat, dummymat);
  KSPGetPC(dummyksp, & dummypc);
  PCSetType(dummypc, PCFIELDSPLIT);
  //PCFieldSplitSchurPrecondition(dummypc,PC_FIELDSPLIT_SCHUR_PRE_SELF,NULL);
  PCFieldSplitSetDetectSaddlePoint(dummypc, PETSC_TRUE);
  PCFieldSplitSetSchurPre(dummypc, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
  KSPSetFromOptions(dummyksp);

  KSPSetUp(dummyksp);
  KSPSetTolerances(dummyksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 iteration for the dummy solve

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======dummymat Matrix======\n");
    MatView(dummymat, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "======dummyX Vector======\n");
    VecView(dummyX, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscPrintf(PETSC_COMM_WORLD, "======Before KSPSolve dummyksp ======\n");
  KSPSolve(dummyksp, dummyX, dummyX);
  PetscPrintf(PETSC_COMM_WORLD, "======After KSPSolve dummyksp ======\n");
  PetscBarrier((PetscObject)dummyX);
  PetscBarrier((PetscObject)dummymat);
  VecDestroy( & dummyX);
  MatDestroy( & dummymat);
  PCFieldSplitGetISByIndex(dummypc, 0, & isBE);
  PCFieldSplitGetISByIndex(dummypc, 1, & isOthers);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======edge & face dofs======\n");
    ISView(isBE, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "======other dofs======\n");
    ISView(isOthers, PETSC_VIEWER_STDOUT_WORLD);
  }

  DMCreateMatrix(da, & tempGD);

  PetscPrintf(PETSC_COMM_WORLD, "======Before MatCreateSubMatrix======\n");
  MatCreateSubMatrix(extGD, isBE, isBE, MAT_INITIAL_MATRIX, & tempGD);
  PetscPrintf(PETSC_COMM_WORLD, "======After MatCreateSubMatrix======\n");
  MatAssemblyBegin(tempGD, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(tempGD, MAT_FINAL_ASSEMBLY);

  MatCopy(tempGD, GD, DIFFERENT_NONZERO_PATTERN);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "======GD = MatCreateSubMatrix(extGD)======\n");
    MatView(GD, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscBarrier((PetscObject)extGD);
  PetscBarrier((PetscObject)tempGD);
  PetscBarrier((PetscObject)dummyksp);
  MatDestroy( & extGD);
  MatDestroy( & tempGD);
  KSPDestroy( & dummyksp);
  return (0);
}


PetscErrorCode FormDerivedDivergence(TS ts, Mat D, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec coordLocal, coordaLocal, Mediag, MnInvdiag, BoundaryCancelDiag, dummyX, dummyF;
  Mat G, Me, MnInv, BoundaryCancel;
  PetscInt N[3], er, ephi, ez, d;
  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];
  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrpphimzm[3], icrmphipzm[3], icrpphipzm[3];
  PetscInt icrmphimzp[3], icrpphimzp[3], icrmphipzp[3], icrpphipzp[3];
  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;
  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  DM dmCoord, dmCoorda;
  PetscScalar ** ** arrCoorda, ** ** arrCoord;

  PetscPrintf(PETSC_COMM_WORLD, "Before Creating MnInv and Me\n");
  DMCreateMatrix(da, & MnInv);
  DMCreateMatrix(da, & Me);
  DMCreateMatrix(da, & BoundaryCancel);
  MatZeroEntries(MnInv);
  MatZeroEntries(Me);
  MatZeroEntries(BoundaryCancel);

  PetscPrintf(PETSC_COMM_WORLD, "Before Creating MnInvdiag and Mediag\n");
  DMCreateGlobalVector(da, & MnInvdiag);
  DMCreateGlobalVector(da, & Mediag);
  DMCreateGlobalVector(da, & BoundaryCancelDiag);
  //VecCreate(PETSC_COMM_WORLD,&MnInvdiag);
  //VecCreate(PETSC_COMM_WORLD,&Mediag);

  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);

  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);

  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }
  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col;
        PetscScalar rmzmedgelength, rmzpedgelength, rmphimedgelength, rmphipedgelength, rpzmedgelength, rpzpedgelength, rpphimedgelength, rpphipedgelength, phipzmedgelength, phimzmedgelength, phimzpedgelength, phipzpedgelength, val;
        PetscInt nEntries = 1;

        /* The edges are oriented in the directions of unit vectors: e_r, e_phi and e_z */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]); /* front right = rpzp */
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* Edge mass matrix: Me */
        /* Equation on down left edge */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN_LEFT;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = DOWN_LEFT;
        col.c = 0;
        val = betaenomp(er, ephi, ez, DOWN_LEFT, user);
        DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        if (er == N[0] - 1) {
          /* Equation on down right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = DOWN_RIGHT;
          col.c = 0;
          val = betaenomp(er, ephi, ez, DOWN_RIGHT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          /* Equation on back right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_RIGHT;
          col.c = 0;
          val = betaenomp(er, ephi, ez, BACK_RIGHT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        /* Equation on back left edge */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_LEFT;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = BACK_LEFT;
        col.c = 0;
        val = betaenomp(er, ephi, ez, BACK_LEFT, user);
        DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        if (ez == N[2] - 1) {
          /* Equation on front left edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_LEFT;
          col.c = 0;
          val = betaenomp(er, ephi, ez, FRONT_LEFT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          /* Equation on front down edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_DOWN;
          col.c = 0;
          val = betaenomp(er, ephi, ez, FRONT_DOWN, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        /* Equation on back down edge */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_DOWN;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = BACK_DOWN;
        col.c = 0;
        val = betaenomp(er, ephi, ez, BACK_DOWN, user);
        DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        if (ephi == N[1] - 1) {
          /* Equation on back up edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_UP;
          col.c = 0;
          val = betaenomp(er, ephi, ez, BACK_UP, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          /* Equation on up left edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_LEFT;
          col.c = 0;
          val = betaenomp(er, ephi, ez, UP_LEFT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        if (er == N[0] - 1 && ephi == N[1] - 1) {
          /* Equation on up right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP_RIGHT;
          col.c = 0;
          val = betaenomp(er, ephi, ez, UP_RIGHT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        if (ephi == N[1] - 1 && ez == N[2] - 1) {
          /* Equation on front up edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_UP;
          col.c = 0;
          val = betaenomp(er, ephi, ez, FRONT_UP, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Equation on front right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_RIGHT;
          col.c = 0;
          val = betaenomp(er, ephi, ez, FRONT_RIGHT, user);
          DMStagMatSetValuesStencil(da, Me, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        /* The inverse of vertex mass matrix: MnInv = M_v^{-1} */
        /* Equation on back down left vertex */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_DOWN_LEFT;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = BACK_DOWN_LEFT;
        col.c = 0;
        val = 1.0 / betavnomp(er, ephi, ez, BACK_DOWN_LEFT, user);
        DMStagMatSetValuesStencil(da, MnInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        if (er != 0 && ez != 0 && (user -> phibtype || ephi != 0)) {
          val = 1.0;
          DMStagMatSetValuesStencil(da, BoundaryCancel, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if (user -> debug) printf("INTERIOR VERTEX DETECTED: BACK_DOWN_LEFT VERTEX OF CELL(%d,%d,%d)\n", er, ephi, ez);
        }

        /* Equation on back down right vertex */
        if (er == N[0] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_DOWN_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_DOWN_RIGHT;
          col.c = 0;
          val = 1.0 / betavnomp(er, ephi, ez, BACK_DOWN_RIGHT, user);
          DMStagMatSetValuesStencil(da, MnInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        /* Equation on back up left vertex */
        if (ephi == N[1] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_LEFT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_UP_LEFT;
          col.c = 0;
          val = 1.0 / betavnomp(er, ephi, ez, BACK_UP_LEFT, user);
          DMStagMatSetValuesStencil(da, MnInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
          if (er != 0 && ez != 0 && user -> phibtype) {
            val = 1.0;
            DMStagMatSetValuesStencil(da, BoundaryCancel, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
            if (user -> debug) printf("INTERIOR VERTEX DETECTED: BACK_UP_LEFT VERTEX OF CELL(%d,%d,%d)\n", er, ephi, ez);
          }
        }

        /* Equation on front down left vertex */
        if (ez == N[2] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_LEFT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_DOWN_LEFT;
          col.c = 0;
          val = 1.0 / betavnomp(er, ephi, ez, FRONT_DOWN_LEFT, user);
          DMStagMatSetValuesStencil(da, MnInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        /* Equation on back up right vertex */
        if (er == N[0] - 1 && ephi == N[1] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = BACK_UP_RIGHT;
          col.c = 0;
          val = 1.0 / betavnomp(er, ephi, ez, BACK_UP_RIGHT, user);
          DMStagMatSetValuesStencil(da, MnInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        /* Equation on front down right vertex */
        if (er == N[0] - 1 && ez == N[2] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_DOWN_RIGHT;
          col.c = 0;
          val = 1.0 / betavnomp(er, ephi, ez, FRONT_DOWN_RIGHT, user);
          DMStagMatSetValuesStencil(da, MnInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        /* Equation on front up left vertex */
        if (ephi == N[1] - 1 && ez == N[2] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_LEFT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_UP_LEFT;
          col.c = 0;
          val = 1.0 / betavnomp(er, ephi, ez, FRONT_UP_LEFT, user);
          DMStagMatSetValuesStencil(da, MnInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        /* Equation on front up right vertex */
        if (er == N[0] - 1 && ephi == N[1] - 1 && ez == N[2] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_UP_RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT_UP_RIGHT;
          col.c = 0;
          val = 1.0 / betavnomp(er, ephi, ez, FRONT_UP_RIGHT, user);
          DMStagMatSetValuesStencil(da, MnInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);
  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscPrintf(PETSC_COMM_WORLD, "Before assembling Me\n");
  MatAssemblyBegin(Me, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Me, MAT_FINAL_ASSEMBLY);
  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Me matrix: \n");
    MatView(Me, PETSC_VIEWER_STDOUT_WORLD);
  }
  PetscPrintf(PETSC_COMM_WORLD, "Before assembling MnInv\n");
  MatAssemblyBegin(MnInv, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(MnInv, MAT_FINAL_ASSEMBLY);
  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "MnInv matrix: \n");
    MatView(MnInv, PETSC_VIEWER_STDOUT_WORLD);
  }
  MatAssemblyBegin(BoundaryCancel, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(BoundaryCancel, MAT_FINAL_ASSEMBLY);
  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "BoundaryCancel matrix: \n");
    MatView(BoundaryCancel, PETSC_VIEWER_STDOUT_WORLD);
  }
  //MatCreate(PETSC_COMM_WORLD,&G);
  DMCreateMatrix(da, & G);

  DMCreateGlobalVector(da, & dummyX);
  DMCreateGlobalVector(da, & dummyF);
  VecSet(dummyX, 0.0);
  VecSet(dummyF, 0.0);
  FormDiscreteGradientEP(ts, G, dummyX, dummyF, user);
  VecDestroy( & dummyX);
  VecDestroy( & dummyF);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Primary Mimetic Gradient matrix: \n");
    MatView(G, PETSC_VIEWER_STDOUT_WORLD);
  }
  MatTranspose(G, MAT_INPLACE_MATRIX, & G); // G := G^T
  /*
  MatMatMatMult(MnInv, G, Me, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &D); // D := MnInv * G * Me
  MatZeroEntries(Me); // Me := 0
  MatAYPX(D, -1.0, Me, DIFFERENT_NONZERO_PATTERN); // D := -D + Me
  */
  MatGetDiagonal(Me, Mediag);

  PetscPrintf(PETSC_COMM_WORLD, "Before VecCopy user->DiagMe = Mediag \n");
  VecCopy(Mediag, user -> DiagMe); //user->DiagMe = Mediag;
  PetscPrintf(PETSC_COMM_WORLD, "After VecCopy user->DiagMe = Mediag \n");

  MatGetDiagonal(MnInv, MnInvdiag);
  MatDiagonalScale(G, MnInvdiag, Mediag);
  MatScale(G, -1.0);
  //MatConvert(G, MATSAME, MAT_INITIAL_MATRIX, &D);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Before MatCopy G, D \n");
    MatView(G, PETSC_VIEWER_STDOUT_WORLD);
  }
  //MatCopy(G, D, DIFFERENT_NONZERO_PATTERN);
  MatZeroEntries(D);
  MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);

  MatGetDiagonal(BoundaryCancel, BoundaryCancelDiag);
  MatDiagonalScale(G, BoundaryCancelDiag, NULL);

  MatAXPY(D, 1.0, G, DIFFERENT_NONZERO_PATTERN);
  //MatMatMult(BoundaryCancel,G,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&D); // D := BoundaryCancel * G this cancels the derived_divergence on boundary vertices

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "D := BoundaryCancel * G matrix\n");
    MatView(D, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscBarrier((PetscObject)G);
  PetscBarrier((PetscObject)MnInv);
  PetscBarrier((PetscObject)Me);
  PetscBarrier((PetscObject)BoundaryCancel);
  PetscBarrier((PetscObject)MnInvdiag);
  PetscBarrier((PetscObject)Mediag);
  PetscBarrier((PetscObject)BoundaryCancelDiag);

  MatDestroy( & G);
  MatDestroy( & MnInv);
  MatDestroy( & Me);
  MatDestroy( & BoundaryCancel);
  VecDestroy( & MnInvdiag);
  VecDestroy( & Mediag);
  VecDestroy( & BoundaryCancelDiag);
  return (0);
}

PetscErrorCode ApplyDerivedDivergence(TS ts, Vec X, Vec F, void * ptr) {
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("ApplyDerivedDivergence",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec fLocal, xLocal;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];

  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;
  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  PetscInt ivEPrmphimzm, ivEPrmphimzp, ivEPrmphipzm;
  PetscInt ivEPrpphimzm, ivEPrpphimzp, ivEPrpphipzm;
  PetscInt ivEPrmphipzp, ivEPrpphipzp;

  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;

  TSGetDM(ts, & da);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Vertex locations */
  DMStagGetLocationSlot(da, BACK_DOWN_LEFT, 3, & ivEPrmphimzm);
  DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, 3, & ivEPrpphimzm);
  DMStagGetLocationSlot(da, BACK_UP_LEFT, 3, & ivEPrmphipzm);
  DMStagGetLocationSlot(da, BACK_UP_RIGHT, 3, & ivEPrpphipzm);
  DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, 3, & ivEPrmphimzp);
  DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, 3, & ivEPrpphimzp);
  DMStagGetLocationSlot(da, FRONT_UP_LEFT, 3, & ivEPrmphipzp);
  DMStagGetLocationSlot(da, FRONT_UP_RIGHT, 3, & ivEPrpphipzp);
  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);
  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Compute function over the locally owned part of the grid */
  /* f(E) = derived_mimetic_divergence(E) */
  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }
        /*PetscPrintf(PETSC_COMM_WORLD,"dphi = %g\n",(double)user->dphi);
        PetscPrintf(PETSC_COMM_WORLD,"Phip - Phim at CELL(%d,%d,%d) = %g\n",(double)(arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]),er,ephi,ez);*/

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* f(E) = - primary_mimetic_gradient^T (beta_e_nomp E)/beta_v_nomp ≡ - M_v^{-1} Grad^T M_e E  */
        /* Inner Back Down Left Vertex */
        if (er != 0 && ez != 0 && (ephi != 0 || user -> phibtype)) {
          arrF[ez][ephi][er][ivEPrmphimzm] = - (- arrX[ez][ephi][er][ivErmphim] * betaenomp(er, ephi, ez, DOWN_LEFT, user) / rmphimedgelength + arrX[ez-1][ephi][er][ivErmphim] * betaenomp(er, ephi, ez-1, DOWN_LEFT, user) / rmphimedgelength - arrX[ez][ephi][er][ivErmzm] * betaenomp(er, ephi, ez, BACK_LEFT, user) / rmzmedgelength + arrX[ez][ephi-1][er][ivErmzm] * betaenomp(er, ephi-1, ez, BACK_LEFT, user) / rmzmedgelength - arrX[ez][ephi][er][ivEphimzm] * betaenomp(er, ephi, ez, BACK_DOWN, user) / phimzmedgelength + arrX[ez][ephi][er-1][ivEphimzm] * betaenomp(er-1, ephi, ez, BACK_DOWN, user) / phimzmedgelength) / betavnomp(er,ephi,ez,BACK_DOWN_LEFT,user);
        }
      }
    }
  }

  /* Restore vectors */
  /* DMStagVecRestoreArray(da,F,&arrF); */

  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode ApplyDeltastar(TS ts, Vec X, Vec F, void * ptr) {
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("ApplyDeltastar",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  Vec Fcopy;

  /* f(EP) = -(1/r) Delta*(EP) = - der_Div ((1/r) prim_Grad(EP)) = (1/r^2) \partial EP / \partial r - (1/r) \partial^2 EP / \partial^2 r - (1/r) \partial^2 EP / \partial^2 z*/

  /* Compute 1/r \tilde{\nabla} (EP) */
  VecDuplicate(X, & Fcopy);
  VecCopy(X, Fcopy);
  FormDiscreteGradientEP_tilde(ts, X, Fcopy, user);

  ApplyDerivedDivergence(ts, Fcopy, F, user);

  /* Restore vectors */
  VecDestroy( & Fcopy);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode ApplyDeltastar2(TS ts, Vec X, Vec F, void * ptr) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("ApplyDeltastar2",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec coordLocal, fLocal, xLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;
  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX;

  //VecZeroEntries(F);
  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  for (d = 0; d < 4; ++d) {
    /* Vertex locations */
    DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]);
    DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]);
  }

  for (d = 0; d < 3; ++d) {
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        if(er>0 && ez>0 && ez < N[2]-1 && er < N[0]-1){
          arrF[ez][ephi][er][ivVrmphimzm[3]] = (arrX[ez][ephi][er+1][ivVrmphimzm[3]] - arrX[ez][ephi][er-1][ivVrmphimzm[3]]) / (2 * (user -> dr) * PetscSqr(arrCoorda[ez][ephi][er][icrmphimzm[0]])) - (arrX[ez][ephi][er+1][ivVrmphimzm[3]] - 2 * arrX[ez][ephi][er][ivVrmphimzm[3]] + arrX[ez][ephi][er-1][ivVrmphimzm[3]]) / (PetscSqr(user -> dr) * arrCoorda[ez][ephi][er][icrmphimzm[0]]) - (arrX[ez+1][ephi][er][ivVrmphimzm[3]] - 2 * arrX[ez][ephi][er][ivVrmphimzm[3]] + arrX[ez-1][ephi][er][ivVrmphimzm[3]]) / (PetscSqr(user -> dz) * arrCoorda[ez][ephi][er][icrmphimzm[0]]);
        }
        if(er>0 && ez == N[2]-1 && er < N[0]-1){
          arrF[ez][ephi][er][ivVrmphimzm[3]] = (arrX[ez][ephi][er+1][ivVrmphimzm[3]] - arrX[ez][ephi][er-1][ivVrmphimzm[3]]) / (2 * (user -> dr) * PetscSqr(arrCoorda[ez][ephi][er][icrmphimzm[0]])) - (arrX[ez][ephi][er+1][ivVrmphimzm[3]] - 2 * arrX[ez][ephi][er][ivVrmphimzm[3]] + arrX[ez][ephi][er-1][ivVrmphimzm[3]]) / (PetscSqr(user -> dr) * arrCoorda[ez][ephi][er][icrmphimzm[0]]) - (arrX[ez][ephi][er][ivVrmphimzp[3]] - 2 * arrX[ez][ephi][er][ivVrmphimzm[3]] + arrX[ez-1][ephi][er][ivVrmphimzm[3]]) / (PetscSqr(user -> dz) * arrCoorda[ez][ephi][er][icrmphimzm[0]]);
        }
        if(ez>0 && er == N[0]-1 && ez < N[2]-1){
          arrF[ez][ephi][er][ivVrmphimzm[3]] = (arrX[ez][ephi][er][ivVrpphimzm[3]] - arrX[ez][ephi][er-1][ivVrmphimzm[3]]) / (2 * (user -> dr) * PetscSqr(arrCoorda[ez][ephi][er][icrmphimzm[0]])) - (arrX[ez][ephi][er][ivVrpphimzm[3]] - 2 * arrX[ez][ephi][er][ivVrmphimzm[3]] + arrX[ez][ephi][er-1][ivVrmphimzm[3]]) / (PetscSqr(user -> dr) * arrCoorda[ez][ephi][er][icrmphimzm[0]]) - (arrX[ez+1][ephi][er][ivVrmphimzm[3]] - 2 * arrX[ez][ephi][er][ivVrmphimzm[3]] + arrX[ez-1][ephi][er][ivVrmphimzm[3]]) / (PetscSqr(user -> dz) * arrCoorda[ez][ephi][er][icrmphimzm[0]]);
        }
        if(er == N[0]-1 && ez == N[2]-1){
          arrF[ez][ephi][er][ivVrmphimzm[3]] = (arrX[ez][ephi][er][ivVrpphimzm[3]] - arrX[ez][ephi][er-1][ivVrmphimzm[3]]) / (2 * (user -> dr) * PetscSqr(arrCoorda[ez][ephi][er][icrmphimzm[0]])) - (arrX[ez][ephi][er][ivVrpphimzm[3]] - 2 * arrX[ez][ephi][er][ivVrmphimzm[3]] + arrX[ez][ephi][er-1][ivVrmphimzm[3]]) / (PetscSqr(user -> dr) * arrCoorda[ez][ephi][er][icrmphimzm[0]]) - (arrX[ez][ephi][er][ivVrmphimzp[3]] - 2 * arrX[ez][ephi][er][ivVrmphimzm[3]] + arrX[ez-1][ephi][er][ivVrmphimzm[3]]) / (PetscSqr(user -> dz) * arrCoorda[ez][ephi][er][icrmphimzm[0]]);
        }
      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode ApplyVectorLaplacian(TS ts, Vec X, Vec F, void * ptr) {
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("ApplyVectorLaplacian",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec fLocal, xLocal, F1, F2, F3, F1Local, F2Local, F3Local;
  Vec coordLocal, GradV1, GradV2, GradV1Local, GradV2Local;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];

  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;
  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  PetscInt ivVrmphimzm[3], ivVrmphimzp[3], ivVrmphipzm[3];
  PetscInt ivVrpphimzm[3], ivVrpphimzp[3], ivVrpphipzm[3];
  PetscInt ivVrmphipzp[3], ivVrpphipzp[3];

  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, ** ** arrF3, ** ** arrF1, ** ** arrF2, ** ** arrGradV1, ** ** arrGradV2, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;

  {
    /* Compute the gradient of V */
    VecDuplicate(X, & F1);
    VecCopy(X, F1);
    VecDuplicate(X, & F2);
    VecCopy(X, F2);
    VecDuplicate(X, & F3);
    VecCopy(X, F3);
    FormDiscreteGradientVectorField(ts, X, F1, F2, F3, user);
  }
  TSGetDM(ts, & da);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Vertex locations */
  for (d = 0; d < 3; ++d) {
    /* Vertex locations */
    DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]);
    DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]);
  }

  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);
  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Compute function over the locally owned part of the grid */
  /* f(V) = Vector_Laplacian(V) */
  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  DMGetLocalVector(da, & F1Local);
  DMGlobalToLocalBegin(da, F1, INSERT_VALUES, F1Local);
  DMGlobalToLocalEnd(da, F1, INSERT_VALUES, F1Local);
  DMStagVecGetArray(da, F1Local, & arrF1);

  DMGetLocalVector(da, & F2Local);
  DMGlobalToLocalBegin(da, F2, INSERT_VALUES, F2Local);
  DMGlobalToLocalEnd(da, F2, INSERT_VALUES, F2Local);
  DMStagVecGetArray(da, F2Local, & arrF2);

  DMGetLocalVector(da, & F3Local);
  DMGlobalToLocalBegin(da, F3, INSERT_VALUES, F3Local);
  DMGlobalToLocalEnd(da, F3, INSERT_VALUES, F3Local);
  DMStagVecGetArray(da, F3Local, & arrF3);

  /* P_{e->v}(prim_grad(V)) */
  DMCreateGlobalVector(da, & GradV1);
  EdgeToVertexProjection(ts, F1, GradV1, user);
  DMGetLocalVector(da, & GradV1Local);
  DMGlobalToLocalBegin(da, GradV1, INSERT_VALUES, GradV1Local);
  DMGlobalToLocalEnd(da, GradV1, INSERT_VALUES, GradV1Local);
  DMStagVecGetArrayRead(da, GradV1Local, & arrGradV1);

  DMCreateGlobalVector(da, & GradV2);
  EdgeToVertexProjection(ts, F2, GradV2, user);
  DMGetLocalVector(da, & GradV2Local);
  DMGlobalToLocalBegin(da, GradV2, INSERT_VALUES, GradV2Local);
  DMGlobalToLocalEnd(da, GradV2, INSERT_VALUES, GradV2Local);
  DMStagVecGetArrayRead(da, GradV2Local, & arrGradV2);

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* f_r(V) = - primary_mimetic_gradient^T (beta_e_nomp primary_mimetic_gradient(V_r))/beta_v_nomp - V_r/r^2 - (2/r) primary_mimetic_gradient_{phi,phi}(V) ≡ - M_v^{-1} Grad^T M_e primary_mimetic_gradient(V_r) - V_r/r^2 - (2/r) primary_mimetic_gradient_{phi,phi}(V)
           f_phi(V) = - primary_mimetic_gradient^T (beta_e_nomp primary_mimetic_gradient(V_phi))/beta_v_nomp - V_phi/r^2 - (2/r) primary_mimetic_gradient_{r,phi}(V) ≡ - M_v^{-1} Grad^T M_e primary_mimetic_gradient(V_phi) - V_phi/r^2 - (2/r) primary_mimetic_gradient_{r,phi}(V)
           f_z(V) = - primary_mimetic_gradient^T (beta_e_nomp primary_mimetic_gradient(V_z))/beta_v_nomp ≡ - M_v^{-1} Grad^T M_e primary_mimetic_gradient(V_z)
         */
        /* Inner Back Down Left Vertex */
        if (er != 0 && ez != 0 && (ephi != 0 || user -> phibtype)) {
          arrF[ez][ephi][er][ivVrmphimzm[0]] = - (- arrF1[ez][ephi][er][ivErmphim] * betaenomp(er, ephi, ez, DOWN_LEFT, user) / rmphimedgelength + arrF1[ez-1][ephi][er][ivErmphim] * betaenomp(er, ephi, ez-1, DOWN_LEFT, user) / rmphimedgelength - arrF1[ez][ephi][er][ivErmzm] * betaenomp(er, ephi, ez, BACK_LEFT, user) / rmzmedgelength + arrF1[ez][ephi-1][er][ivErmzm] * betaenomp(er, ephi-1, ez, BACK_LEFT, user) / rmzmedgelength - arrF1[ez][ephi][er][ivEphimzm] * betaenomp(er, ephi, ez, BACK_DOWN, user) / phimzmedgelength + arrF1[ez][ephi][er-1][ivEphimzm] * betaenomp(er-1, ephi, ez, BACK_DOWN, user) / phimzmedgelength) / betavnomp(er,ephi,ez,BACK_DOWN_LEFT,user) - arrX[ez][ephi][er][ivVrmphimzm[0]] / PetscSqr(arrCoorda[ez][ephi][er][icrmphimzm[0]]) - 2.0 * arrGradV2[ez][ephi][er][ivVrmphimzm[1]] / arrCoorda[ez][ephi][er][icrmphimzm[0]];
          arrF[ez][ephi][er][ivVrmphimzm[1]] = - (- arrF2[ez][ephi][er][ivErmphim] * betaenomp(er, ephi, ez, DOWN_LEFT, user) / rmphimedgelength + arrF2[ez-1][ephi][er][ivErmphim] * betaenomp(er, ephi, ez-1, DOWN_LEFT, user) / rmphimedgelength - arrF2[ez][ephi][er][ivErmzm] * betaenomp(er, ephi, ez, BACK_LEFT, user) / rmzmedgelength + arrF2[ez][ephi-1][er][ivErmzm] * betaenomp(er, ephi-1, ez, BACK_LEFT, user) / rmzmedgelength - arrF2[ez][ephi][er][ivEphimzm] * betaenomp(er, ephi, ez, BACK_DOWN, user) / phimzmedgelength + arrF2[ez][ephi][er-1][ivEphimzm] * betaenomp(er-1, ephi, ez, BACK_DOWN, user) / phimzmedgelength) / betavnomp(er,ephi,ez,BACK_DOWN_LEFT,user) - arrX[ez][ephi][er][ivVrmphimzm[1]] / PetscSqr(arrCoorda[ez][ephi][er][icrmphimzm[0]]) + 2.0 * arrGradV1[ez][ephi][er][ivVrmphimzm[1]] / arrCoorda[ez][ephi][er][icrmphimzm[0]] ;
          arrF[ez][ephi][er][ivVrmphimzm[2]] = - (- arrF3[ez][ephi][er][ivErmphim] * betaenomp(er, ephi, ez, DOWN_LEFT, user) / rmphimedgelength + arrF3[ez-1][ephi][er][ivErmphim] * betaenomp(er, ephi, ez-1, DOWN_LEFT, user) / rmphimedgelength - arrF3[ez][ephi][er][ivErmzm] * betaenomp(er, ephi, ez, BACK_LEFT, user) / rmzmedgelength + arrF3[ez][ephi-1][er][ivErmzm] * betaenomp(er, ephi-1, ez, BACK_LEFT, user) / rmzmedgelength - arrF3[ez][ephi][er][ivEphimzm] * betaenomp(er, ephi, ez, BACK_DOWN, user) / phimzmedgelength + arrF3[ez][ephi][er-1][ivEphimzm] * betaenomp(er-1, ephi, ez, BACK_DOWN, user) / phimzmedgelength) / betavnomp(er,ephi,ez,BACK_DOWN_LEFT,user);
        }
      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArray(da, F1Local, & arrF1);
  DMRestoreLocalVector(da, & F1Local);

  DMStagVecRestoreArray(da, F2Local, & arrF2);
  DMRestoreLocalVector(da, & F2Local);

  DMStagVecRestoreArray(da, F3Local, & arrF3);
  DMRestoreLocalVector(da, & F3Local);

  DMStagVecRestoreArrayRead(da, GradV1Local, & arrGradV1);
  DMRestoreLocalVector(da, & GradV1Local);

  DMStagVecRestoreArrayRead(da, GradV2Local, & arrGradV2);
  DMRestoreLocalVector(da, & GradV2Local);

  VecDestroy(& GradV1);
  VecDestroy(& GradV2);
  VecDestroy(& F1);
  VecDestroy(& F2);
  VecDestroy(& F3);

  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode FormDerivedGradient(TS ts, Mat D, Mat G, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec coordLocal, coordaLocal, Mcdiag, MfInvdiag;
  Mat Mc, MfInv;
  PetscInt N[3], er, ephi, ez, d;
  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];
  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrpphimzm[3], icrmphipzm[3], icrpphipzm[3];
  PetscInt icrmphimzp[3], icrpphimzp[3], icrmphipzp[3], icrpphipzp[3];
  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;
  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  DM dmCoord, dmCoorda;
  PetscScalar ** ** arrCoorda, ** ** arrCoord;

  PetscPrintf(PETSC_COMM_WORLD, "Before Creating MnInv and Me\n");
  DMCreateMatrix(coordDA, & MfInv);
  DMCreateMatrix(coordDA, & Mc);
  MatZeroEntries(Mc);
  MatZeroEntries(MfInv);

  PetscPrintf(PETSC_COMM_WORLD, "Before Creating MfInvdiag and Mcdiag\n");
  DMCreateGlobalVector(coordDA, & MfInvdiag);
  DMCreateGlobalVector(coordDA, & Mcdiag);

  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);

  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);

  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }
  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col;
        PetscScalar rmzmedgelength, rmzpedgelength, rmphimedgelength, rmphipedgelength, rpzmedgelength, rpzpedgelength, rpphimedgelength, rpphipedgelength, phipzmedgelength, phimzmedgelength, phimzpedgelength, phipzpedgelength, val;
        PetscInt nEntries = 1;

        /* The edges are oriented in the directions of unit vectors: e_r, e_phi and e_z */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]); /* front right = rpzp */
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* Cell mass matrix: Mc */
        /* Equation on element */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = ELEMENT;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = ELEMENT;
        col.c = 0;
        val = alphac(er, ephi, ez, user);
        DMStagMatSetValuesStencil(coordDA, Mc, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);

        /* The inverse of face mass matrix with material properties: MfInv = \widehat{M}_f^{-1} */
        /* Equation on left face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = LEFT;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = LEFT;
        col.c = 0;
        val = 1.0 / betaf_wmp(er, ephi, ez, LEFT, user);
        DMStagMatSetValuesStencil(coordDA, MfInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);

        /* Equation on down face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = DOWN;
        col.c = 0;
        val = 1.0 / betaf_wmp(er, ephi, ez, DOWN, user);
        DMStagMatSetValuesStencil(coordDA, MfInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);

        /* Equation on back face */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK;
        row.c = 0;
        col.i = er;
        col.j = ephi;
        col.k = ez;
        col.loc = BACK;
        col.c = 0;
        val = 1.0 / betaf_wmp(er, ephi, ez, BACK, user);
        DMStagMatSetValuesStencil(coordDA, MfInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);

        /* Equation on right face */
        if (er == N[0] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = RIGHT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = RIGHT;
          col.c = 0;
          val = 1.0 / betaf_wmp(er, ephi, ez, RIGHT, user);
          DMStagMatSetValuesStencil(coordDA, MfInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        /* Equation on up face */
        if (ephi == N[1] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = UP;
          col.c = 0;
          val = 1.0 / betaf_wmp(er, ephi, ez, UP, user);
          DMStagMatSetValuesStencil(coordDA, MfInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

        /* Equation on front face */
        if (ez == N[2] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT;
          row.c = 0;
          col.i = er;
          col.j = ephi;
          col.k = ez;
          col.loc = FRONT;
          col.c = 0;
          val = 1.0 / betaf_wmp(er, ephi, ez, FRONT, user);
          DMStagMatSetValuesStencil(coordDA, MfInv, nEntries, & row, nEntries, & col, & val, INSERT_VALUES);
        }

      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);
  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscPrintf(PETSC_COMM_WORLD, "Before assembling Mc\n");
  MatAssemblyBegin(Mc, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Mc, MAT_FINAL_ASSEMBLY);
  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Mc matrix: \n");
    MatView(Mc, PETSC_VIEWER_STDOUT_WORLD);
  }
  PetscPrintf(PETSC_COMM_WORLD, "Before assembling widehat{M}_f^{-1}\n");
  MatAssemblyBegin(MfInv, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(MfInv, MAT_FINAL_ASSEMBLY);
  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "widehat{M}_f^{-1} matrix: \n");
    MatView(MfInv, PETSC_VIEWER_STDOUT_WORLD);
  }

  MatCopy(D, G, DIFFERENT_NONZERO_PATTERN);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Primary Mimetic Divergence matrix: \n");
    MatView(G, PETSC_VIEWER_STDOUT_WORLD);
  }
  MatTranspose(G, MAT_INPLACE_MATRIX, & G); // G := G^T
  MatGetDiagonal(Mc, Mcdiag);

  MatGetDiagonal(MfInv, MfInvdiag);
  MatDiagonalScale(G, MfInvdiag, Mcdiag);
  MatScale(G, -1.0);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Derived Gradient matrix\n");
    MatView(G, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscBarrier((PetscObject)MfInv);
  PetscBarrier((PetscObject)Mc);
  PetscBarrier((PetscObject)MfInvdiag);
  PetscBarrier((PetscObject)Mcdiag);

  MatDestroy( & MfInv);
  MatDestroy( & Mc);
  VecDestroy( & MfInvdiag);
  VecDestroy( & Mcdiag);
  return (0);
}

PetscErrorCode FormDiscreteGradient(TS ts, Mat G, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec coordLocal, coordaLocal;
  PetscInt N[3], er, ephi, ez, d;
  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];
  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrpphimzm[3], icrmphipzm[3], icrpphipzm[3];
  PetscInt icrmphimzp[3], icrpphimzp[3], icrmphipzp[3], icrpphipzp[3];
  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;
  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  DM dmCoord, dmCoorda;
  PetscScalar ** ** arrCoorda, ** ** arrCoord;

  MatZeroEntries(G);
  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);

  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);

  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col[2];
        PetscScalar rmzmedgelength, rmzpedgelength, rmphimedgelength, rmphipedgelength,
        rpzmedgelength, rpzpedgelength, rpphimedgelength, rpphipedgelength,
        phipzmedgelength, phimzmedgelength, phimzpedgelength, phipzpedgelength,
        valG[2];
        PetscInt nEntries = 2;

        /* The edges are oriented in the directions of unit vectors: e_r, e_phi and e_z */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]); /* front right = rpzp */
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* Equation on down left edge */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN_LEFT;
        row.c = 0;
        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = FRONT_DOWN_LEFT;
        col[0].c = 0;
        valG[0] = 1.0 / rmphimedgelength;
        col[1].i = er;
        col[1].j = ephi;
        col[1].k = ez;
        col[1].loc = BACK_DOWN_LEFT;
        col[1].c = 0;
        valG[1] = -1.0 / rmphimedgelength;
        DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        if (er == N[0] - 1) {
          /* Equation on down right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN_RIGHT;
          col[0].c = 0;
          valG[0] = 1.0 / rpphimedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_DOWN_RIGHT;
          col[1].c = 0;
          valG[1] = -1.0 / rpphimedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
          /* Equation on back right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_RIGHT;
          col[0].c = 0;
          valG[0] = 1.0 / rpzmedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_DOWN_RIGHT;
          col[1].c = 0;
          valG[1] = -1.0 / rpzmedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }
        /* Equation on back left edge */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_LEFT;
        row.c = 0;
        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = BACK_UP_LEFT;
        col[0].c = 0;
        valG[0] = 1.0 / rmzmedgelength;
        col[1].i = er;
        col[1].j = ephi;
        col[1].k = ez;
        col[1].loc = BACK_DOWN_LEFT;
        col[1].c = 0;
        valG[1] = -1.0 / rmzmedgelength;
        DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        if (ez == N[2] - 1) {
          /* Equation on front left edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_LEFT;
          col[0].c = 0;
          valG[0] = 1.0 / rmzpedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_LEFT;
          col[1].c = 0;
          valG[1] = -1.0 / rmzpedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
          /* Equation on front down edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN_RIGHT;
          col[0].c = 0;
          valG[0] = 1.0 / phimzpedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_LEFT;
          col[1].c = 0;
          valG[1] = -1.0 / phimzpedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }
        /* Equation on back down edge */
        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_DOWN;
        row.c = 0;
        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = BACK_DOWN_RIGHT;
        col[0].c = 0;
        valG[0] = 1.0 / phimzmedgelength;
        col[1].i = er;
        col[1].j = ephi;
        col[1].k = ez;
        col[1].loc = BACK_DOWN_LEFT;
        col[1].c = 0;
        valG[1] = -1.0 / phimzmedgelength;
        DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        if (ephi == N[1] - 1) {
          /* Equation on back up edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_RIGHT;
          col[0].c = 0;
          valG[0] = 1.0 / phipzmedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP_LEFT;
          col[1].c = 0;
          valG[1] = -1.0 / phipzmedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
          /* Equation on up left edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_LEFT;
          col[0].c = 0;
          valG[0] = 1.0 / rmphipedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP_LEFT;
          col[1].c = 0;
          valG[1] = -1.0 / rmphipedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }
        if (er == N[0] - 1 && ephi == N[1] - 1) {
          /* Equation on up right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_RIGHT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_RIGHT;
          col[0].c = 0;
          valG[0] = 1.0 / rpphipedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP_RIGHT;
          col[1].c = 0;
          valG[1] = -1.0 / rpphipedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }
        if (ephi == N[1] - 1 && ez == N[2] - 1) {
          /* Equation on front up edge */
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
          valG[0] = 1.0 / phipzpedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_UP_LEFT;
          col[1].c = 0;
          valG[1] = -1.0 / phipzpedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Equation on front right edge */
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_UP_RIGHT;
          col[0].c = 0;
          valG[0] = 1.0 / rpzpedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_RIGHT;
          col[1].c = 0;
          valG[1] = -1.0 / rpzpedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }
      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);
  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);
  return (0);

}

PetscErrorCode FormPrimaryCurl(TS ts, Vec X, Vec F, void * ptr) {
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FormPrimaryCurl",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec fLocal, xLocal;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];

  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];
  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;

  TSGetDM(ts, & da);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);
  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);
  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Compute function over the locally owned part of the grid */
  /* f1(B,E) = primary_mimetic_curl(E) */
  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal );
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }
        /*PetscPrintf(PETSC_COMM_WORLD,"dphi = %g\n",(double)user->dphi);
        PetscPrintf(PETSC_COMM_WORLD,"Phip - Phim at CELL(%d,%d,%d) = %g\n",(double)(arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]),er,ephi,ez);*/

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* f1(B,E) = primary_mimetic_curl(E) */
        arrF[ez][ephi][er][ivBrm] = (rmzmedgelength * arrX[ez][ephi][er][ivErmzm] - rmzpedgelength * arrX[ez][ephi][er][ivErmzp] + rmphipedgelength * arrX[ez][ephi][er][ivErmphip] - rmphimedgelength * arrX[ez][ephi][er][ivErmphim]) / surface(er, ephi, ez, LEFT, user); /* Left face */

        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivBrp] = (rpzmedgelength * arrX[ez][ephi][er][ivErpzm] - rpzpedgelength * arrX[ez][ephi][er][ivErpzp] + rpphipedgelength * arrX[ez][ephi][er][ivErpphip] - rpphimedgelength * arrX[ez][ephi][er][ivErpphim]) / surface(er, ephi, ez, RIGHT, user); /* Right face */
        }

        arrF[ez][ephi][er][ivBphim] = (-phimzmedgelength * arrX[ez][ephi][er][ivEphimzm] + phimzpedgelength * arrX[ez][ephi][er][ivEphimzp] - rpphimedgelength * arrX[ez][ephi][er][ivErpphim] + rmphimedgelength * arrX[ez][ephi][er][ivErmphim]) / surface(er, ephi, ez, DOWN, user); /* Face down */

        if (ephi == N[1] - 1) {
          arrF[ez][ephi][er][ivBphip] = (-phipzmedgelength * arrX[ez][ephi][er][ivEphipzm] + phipzpedgelength * arrX[ez][ephi][er][ivEphipzp] - rpphipedgelength * arrX[ez][ephi][er][ivErpphip] + rmphipedgelength * arrX[ez][ephi][er][ivErmphip]) / surface(er, ephi, ez, UP, user); /* Face up */
        }

        arrF[ez][ephi][er][ivBzm] = (phimzmedgelength * arrX[ez][ephi][er][ivEphimzm] - phipzmedgelength * arrX[ez][ephi][er][ivEphipzm] + rpzmedgelength * arrX[ez][ephi][er][ivErpzm] - rmzmedgelength * arrX[ez][ephi][er][ivErmzm]) / surface(er, ephi, ez, BACK, user); /* Back face */

        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivBzp] = (phimzpedgelength * arrX[ez][ephi][er][ivEphimzp] - phipzpedgelength * arrX[ez][ephi][er][ivEphipzp] + rpzpedgelength * arrX[ez][ephi][er][ivErpzp] - rmzpedgelength * arrX[ez][ephi][er][ivErmzp]) / surface(er, ephi, ez, FRONT, user); /* Front face */
        }

      }
    }
  }
  /* End of triple for loop */

  /* Restore vectors */
  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);
  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);
  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "F = \n");
    VecView(F, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode FormDerivedCurl(TS ts, Vec X, Vec F, void * ptr) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FormDerivedCurl",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec fLocal, xLocal;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];

  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];
  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;

  TSGetDM(ts, & da);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);
  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);
  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Compute function over the locally owned part of the grid */
  /* f(B) = derived_mimetic_curl(B) */
  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }
        /*PetscPrintf(PETSC_COMM_WORLD,"dphi = %g\n",(double)user->dphi);
        PetscPrintf(PETSC_COMM_WORLD,"Phip - Phim at CELL(%d,%d,%d) = %g\n",(double)(arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]),er,ephi,ez);*/

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* f(B) = primary_mimetic_curl^T (beta_f B)/beta_e ≡ M_e^{-1} Curl^T M_f B  */
        /* Back Left edge */
        //if (er == 0 && ez == 0) {arrF[ez][ephi][er][ivErmzm] = ( arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * rmzmedgelength / betae(er, ephi, ez, BACK_LEFT, user);} /* Back Left boundary edge */
        //else if (er == 0) {arrF[ez][ephi][er][ivErmzm] = ( arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) - arrX[ez-1][ephi][er][ivBrm] * betaf(er, ephi, ez-1, LEFT, user)/surface(er, ephi, ez-1, LEFT, user)) * rmzmedgelength / betae(er, ephi, ez, BACK_LEFT, user);} /* Left boundary edge */
        //else if (ez == 0) {arrF[ez][ephi][er][ivErmzm] = ( arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) + arrX[ez][ephi][er-1][ivBzm] * betaf(er-1, ephi, ez, BACK, user)/surface(er-1, ephi, ez, BACK, user)) * rmzmedgelength / betae(er, ephi, ez, BACK_LEFT, user);} /* Back boundary edge */
        if (er != 0 && ez != 0) {
          arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
            arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
            arrX[ez - 1][ephi][er][ivBrm] * betaf(er, ephi, ez - 1, LEFT, user) / surface(er, ephi, ez - 1, LEFT, user) +
            arrX[ez][ephi][er - 1][ivBzm] * betaf(er - 1, ephi, ez, BACK, user) / surface(er - 1, ephi, ez, BACK, user)) * rmzmedgelength / betae(er, ephi, ez, BACK_LEFT, user);
        } /* Internal edge */

        if (!(user -> phibtype)) {
          /* Back Down edge */
          //if (ephi == 0 && ez == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);} /* Back Down boundary edge */
          //else if (ephi == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) + arrX[ez-1][ephi][er][ivBphim] * betaf(er, ephi, ez-1, DOWN, user)/surface(er, ephi, ez-1, DOWN, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);} /* Down boundary edge */
          //else if (ez == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) - arrX[ez][ephi-1][er][ivBzm] * betaf(er, ephi-1, ez, BACK, user)/surface(er, ephi-1, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);} /* Back boundary edge */
          if (ephi != 0 && ez != 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
              arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);
          } /* Internal edge */

          /* Down Left edge */
          //if (er == 0 && ephi == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);} /* Down Left boundary edge */
          //else if (er == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi-1][er][ivBrm] * betaf(er, ephi-1, ez, LEFT, user)/surface(er, ephi-1, ez, LEFT, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);} /* Left boundary edge */
          //else if (ephi == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er-1][ivBphim] * betaf(er-1, ephi, ez, DOWN, user)/surface(er-1, ephi, ez, DOWN, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);} /* Down boundary edge */
          if (ephi != 0 && ez != 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
              arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);
          } /* Internal edge */
        } else {
          //if (ez == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) - arrX[ez][ephi-1][er][ivBzm] * betaf(er, ephi-1, ez, BACK, user)/surface(er, ephi-1, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);} /* Back boundary edge */
          if (ez != 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
              arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);
          } /* Internal edge */

          //if (er == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi-1][er][ivBrm] * betaf(er, ephi-1, ez, LEFT, user)/surface(er, ephi-1, ez, LEFT, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);} /* Left boundary edge */
          if (er != 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
              arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);
          } /* Internal edge */
        }

        if (er == N[0] - 1) {
          /* Right boundary Back edge */
          //if (ez == 0) {arrF[ez][ephi][er][ivErpzm] = ( arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * rpzmedgelength / betae(er, ephi, ez, BACK_RIGHT, user);} /* Back Right boundary edge */
          //else {arrF[ez][ephi][er][ivErpzm] = ( arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) - arrX[ez-1][ephi][er][ivBrp] * betaf(er, ephi, ez-1, RIGHT, user)/surface(er, ephi, ez-1, RIGHT, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * rpzmedgelength / betae(er, ephi, ez, BACK_RIGHT, user);} /* Right boundary edge */

          /* Right boundary Down edge */
          //if (ephi == 0 && !(user->phibtype)) {arrF[ez][ephi][er][ivErpphim] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user)) * rpphimedgelength / betae(er, ephi, ez, DOWN_RIGHT, user);} /* Down Right boundary edge */
          //else {arrF[ez][ephi][er][ivErpphim] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi-1][er][ivBrp] * betaf(er, ephi-1, ez, RIGHT, user)/surface(er, ephi-1, ez, RIGHT, user)) * rpphimedgelength / betae(er, ephi, ez, DOWN_RIGHT, user);} /* Right boundary edge */
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Up boundary Back edge */
          //if (ez == 0) {arrF[ez][ephi][er][ivEphipzm] = ( -arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * phipzmedgelength / betae(er, ephi, ez, BACK_UP, user);} /* Back Up boundary edge */
          //else {arrF[ez][ephi][er][ivEphipzm] = ( -arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) + arrX[ez-1][ephi][er][ivBphip] * betaf(er, ephi, ez-1, UP, user)/surface(er, ephi, ez-1, UP, user)) * phipzmedgelength / betae(er, ephi, ez, BACK_UP, user);} /* Up boundary edge */

          /* Up boundary Left edge */
          //if (er == 0) {arrF[ez][ephi][er][ivErmphip] = ( arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user)) * rmphipedgelength / betae(er, ephi, ez, UP_LEFT, user);} /* Up Left boundary edge */
          //else {arrF[ez][ephi][er][ivErmphip] = ( arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er-1][ivBphip] * betaf(er-1, ephi, ez, UP, user)/surface(er-1, ephi, ez, UP, user)) * rmphipedgelength / betae(er, ephi, ez, UP_LEFT, user);} /* Up boundary edge */
        }

        if (ez == N[2] - 1) {
          /* Front boundary Left edge */
          //if (er == 0) {arrF[ez][ephi][er][ivErmzp] = ( -arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user)) * rmzpedgelength / betae(er, ephi, ez, FRONT_LEFT, user);} /* Front Left boundary edge */
          //else {arrF[ez][ephi][er][ivErmzp] = ( -arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er-1][ivBzp] * betaf(er-1, ephi, ez, FRONT, user)/surface(er-1, ephi, ez, FRONT, user)) * rmzpedgelength / betae(er, ephi, ez, FRONT_LEFT, user);} /* Front boundary edge */

          /* Front boundary Down edge */
          //if (ephi == 0 && !(user->phibtype)) {arrF[ez][ephi][er][ivEphimzp] = ( arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user)) * phimzpedgelength / betae(er, ephi, ez, FRONT_DOWN, user);} /* Front Down boundary edge */
          //else {arrF[ez][ephi][er][ivEphimzp] = ( arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi-1][er][ivBzp] * betaf(er, ephi-1, ez, FRONT, user)/surface(er, ephi-1, ez, FRONT, user)) * phimzpedgelength / betae(er, ephi, ez, FRONT_DOWN, user);} /* Front boundary edge */
        }

        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Up Right boundary edge */
          //arrF[ez][ephi][er][ivErpphip] = ( -arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user)) * rpphipedgelength / betae(er, ephi, ez, UP_RIGHT, user);
        }
        if (ephi == N[1] - 1 && ez == N[2] - 1 && !(user -> phibtype)) {
          /* Front Up boundary edge */
          //arrF[ez][ephi][er][ivEphipzp] = ( arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user)) * phipzpedgelength / betae(er, ephi, ez, FRONT_UP, user);
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Front Right boundary edge */
          //arrF[ez][ephi][er][ivErpzp] = ( -arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user)) * rpzpedgelength / betae(er, ephi, ez, FRONT_RIGHT, user);
        }
      }
    }
  }

  /* Restore vectors */
  /* DMStagVecRestoreArray(da,F,&arrF); */

  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode FormDerivedCurlnores(TS ts, Vec X, Vec F, void * ptr) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FormDerivedCurlnores",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec fLocal, xLocal;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];

  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];
  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;

  TSGetDM(ts, & da);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);
  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);
  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Compute function over the locally owned part of the grid */
  /* f(B) = derived_mimetic_curl(B) */
  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }
        /*PetscPrintf(PETSC_COMM_WORLD,"dphi = %g\n",(double)user->dphi);
        PetscPrintf(PETSC_COMM_WORLD,"Phip - Phim at CELL(%d,%d,%d) = %g\n",(double)(arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]),er,ephi,ez);*/

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* f(B) = primary_mimetic_curl^T (beta_f B)/beta_e ≡ M_e^{-1} Curl^T M_f B  */
        /* Back Left edge */
        //if (er == 0 && ez == 0) {arrF[ez][ephi][er][ivErmzm] = ( arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * rmzmedgelength / betaenores(er, ephi, ez, BACK_LEFT, user);} /* Back Left boundary edge */
        //else if (er == 0) {arrF[ez][ephi][er][ivErmzm] = ( arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) - arrX[ez-1][ephi][er][ivBrm] * betaf(er, ephi, ez-1, LEFT, user)/surface(er, ephi, ez-1, LEFT, user)) * rmzmedgelength / betaenores(er, ephi, ez, BACK_LEFT, user);} /* Left boundary edge */
        //else if (ez == 0) {arrF[ez][ephi][er][ivErmzm] = ( arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) + arrX[ez][ephi][er-1][ivBzm] * betaf(er-1, ephi, ez, BACK, user)/surface(er-1, ephi, ez, BACK, user)) * rmzmedgelength / betaenores(er, ephi, ez, BACK_LEFT, user);} /* Back boundary edge */
        if (er != 0 && ez != 0) {
          arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
            arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
            arrX[ez - 1][ephi][er][ivBrm] * betaf(er, ephi, ez - 1, LEFT, user) / surface(er, ephi, ez - 1, LEFT, user) +
            arrX[ez][ephi][er - 1][ivBzm] * betaf(er - 1, ephi, ez, BACK, user) / surface(er - 1, ephi, ez, BACK, user)) * rmzmedgelength / betaenores(er, ephi, ez, BACK_LEFT, user);
        } /* Internal edge */

        if (!(user -> phibtype)) {
          /* Back Down edge */
          //if (ephi == 0 && ez == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * phimzmedgelength / betaenores(er, ephi, ez, BACK_DOWN, user);} /* Back Down boundary edge */
          //else if (ephi == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) + arrX[ez-1][ephi][er][ivBphim] * betaf(er, ephi, ez-1, DOWN, user)/surface(er, ephi, ez-1, DOWN, user)) * phimzmedgelength / betaenores(er, ephi, ez, BACK_DOWN, user);} /* Down boundary edge */
          //else if (ez == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) - arrX[ez][ephi-1][er][ivBzm] * betaf(er, ephi-1, ez, BACK, user)/surface(er, ephi-1, ez, BACK, user)) * phimzmedgelength / betaenores(er, ephi, ez, BACK_DOWN, user);} /* Back boundary edge */
          if (ephi != 0 && ez != 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
              arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betaenores(er, ephi, ez, BACK_DOWN, user);
          } /* Internal edge */

          /* Down Left edge */
          //if (er == 0 && ephi == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user)) * rmphimedgelength / betaenores(er, ephi, ez, DOWN_LEFT, user);} /* Down Left boundary edge */
          //else if (er == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi-1][er][ivBrm] * betaf(er, ephi-1, ez, LEFT, user)/surface(er, ephi-1, ez, LEFT, user)) * rmphimedgelength / betaenores(er, ephi, ez, DOWN_LEFT, user);} /* Left boundary edge */
          //else if (ephi == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er-1][ivBphim] * betaf(er-1, ephi, ez, DOWN, user)/surface(er-1, ephi, ez, DOWN, user)) * rmphimedgelength / betaenores(er, ephi, ez, DOWN_LEFT, user);} /* Down boundary edge */
          if (ephi != 0 && ez != 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
              arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betaenores(er, ephi, ez, DOWN_LEFT, user);
          } /* Internal edge */
        } else {
          //if (ez == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) - arrX[ez][ephi-1][er][ivBzm] * betaf(er, ephi-1, ez, BACK, user)/surface(er, ephi-1, ez, BACK, user)) * phimzmedgelength / betaenores(er, ephi, ez, BACK_DOWN, user);} /* Back boundary edge */
          if (ez != 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
              arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betaenores(er, ephi, ez, BACK_DOWN, user);
          } /* Internal edge */

          //if (er == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi-1][er][ivBrm] * betaf(er, ephi-1, ez, LEFT, user)/surface(er, ephi-1, ez, LEFT, user)) * rmphimedgelength / betaenores(er, ephi, ez, DOWN_LEFT, user);} /* Left boundary edge */
          if (er != 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
              arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betaenores(er, ephi, ez, DOWN_LEFT, user);
          } /* Internal edge */
        }

        if (er == N[0] - 1) {
          /* Right boundary Back edge */
          //if (ez == 0) {arrF[ez][ephi][er][ivErpzm] = ( arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * rpzmedgelength / betaenores(er, ephi, ez, BACK_RIGHT, user);} /* Back Right boundary edge */
          //else {arrF[ez][ephi][er][ivErpzm] = ( arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) - arrX[ez-1][ephi][er][ivBrp] * betaf(er, ephi, ez-1, RIGHT, user)/surface(er, ephi, ez-1, RIGHT, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * rpzmedgelength / betaenores(er, ephi, ez, BACK_RIGHT, user);} /* Right boundary edge */

          /* Right boundary Down edge */
          //if (ephi == 0 && !(user->phibtype)) {arrF[ez][ephi][er][ivErpphim] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user)) * rpphimedgelength / betaenores(er, ephi, ez, DOWN_RIGHT, user);} /* Down Right boundary edge */
          //else {arrF[ez][ephi][er][ivErpphim] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi-1][er][ivBrp] * betaf(er, ephi-1, ez, RIGHT, user)/surface(er, ephi-1, ez, RIGHT, user)) * rpphimedgelength / betaenores(er, ephi, ez, DOWN_RIGHT, user);} /* Right boundary edge */
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Up boundary Back edge */
          //if (ez == 0) {arrF[ez][ephi][er][ivEphipzm] = ( -arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * phipzmedgelength / betaenores(er, ephi, ez, BACK_UP, user);} /* Back Up boundary edge */
          //else {arrF[ez][ephi][er][ivEphipzm] = ( -arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) + arrX[ez-1][ephi][er][ivBphip] * betaf(er, ephi, ez-1, UP, user)/surface(er, ephi, ez-1, UP, user)) * phipzmedgelength / betaenores(er, ephi, ez, BACK_UP, user);} /* Up boundary edge */

          /* Up boundary Left edge */
          //if (er == 0) {arrF[ez][ephi][er][ivErmphip] = ( arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user)) * rmphipedgelength / betaenores(er, ephi, ez, UP_LEFT, user);} /* Up Left boundary edge */
          //else {arrF[ez][ephi][er][ivErmphip] = ( arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er-1][ivBphip] * betaf(er-1, ephi, ez, UP, user)/surface(er-1, ephi, ez, UP, user)) * rmphipedgelength / betaenores(er, ephi, ez, UP_LEFT, user);} /* Up boundary edge */
        }

        if (ez == N[2] - 1) {
          /* Front boundary Left edge */
          //if (er == 0) {arrF[ez][ephi][er][ivErmzp] = ( -arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user)) * rmzpedgelength / betaenores(er, ephi, ez, FRONT_LEFT, user);} /* Front Left boundary edge */
          //else {arrF[ez][ephi][er][ivErmzp] = ( -arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er-1][ivBzp] * betaf(er-1, ephi, ez, FRONT, user)/surface(er-1, ephi, ez, FRONT, user)) * rmzpedgelength / betaenores(er, ephi, ez, FRONT_LEFT, user);} /* Front boundary edge */

          /* Front boundary Down edge */
          //if (ephi == 0 && !(user->phibtype)) {arrF[ez][ephi][er][ivEphimzp] = ( arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user)) * phimzpedgelength / betaenores(er, ephi, ez, FRONT_DOWN, user);} /* Front Down boundary edge */
          //else {arrF[ez][ephi][er][ivEphimzp] = ( arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi-1][er][ivBzp] * betaf(er, ephi-1, ez, FRONT, user)/surface(er, ephi-1, ez, FRONT, user)) * phimzpedgelength / betaenores(er, ephi, ez, FRONT_DOWN, user);} /* Front boundary edge */
        }

        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Up Right boundary edge */
          //arrF[ez][ephi][er][ivErpphip] = ( -arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user)) * rpphipedgelength / betaenores(er, ephi, ez, UP_RIGHT, user);
        }
        if (ephi == N[1] - 1 && ez == N[2] - 1 && !(user -> phibtype)) {
          /* Front Up boundary edge */
          //arrF[ez][ephi][er][ivEphipzp] = ( arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user)) * phipzpedgelength / betaenores(er, ephi, ez, FRONT_UP, user);
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Front Right boundary edge */
          //arrF[ez][ephi][er][ivErpzp] = ( -arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user)) * rpzpedgelength / betaenores(er, ephi, ez, FRONT_RIGHT, user);
        }
      }
    }
  }

  /* Restore vectors */
  /* DMStagVecRestoreArray(da,F,&arrF); */

  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode FormDerivedCurlnomp(TS ts, Vec X, Vec F, void * ptr) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FormDerivedCurlnomp",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec fLocal, xLocal;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];

  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];
  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;

  TSGetDM(ts, & da);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);
  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);
  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Compute function over the locally owned part of the grid */
  /* f(B) = derived_mimetic_curl(B) */
  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }
        /*PetscPrintf(PETSC_COMM_WORLD,"dphi = %g\n",(double)user->dphi);
        PetscPrintf(PETSC_COMM_WORLD,"Phip - Phim at CELL(%d,%d,%d) = %g\n",(double)(arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]),er,ephi,ez);*/

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* f(B) = primary_mimetic_curl^T (beta_f B)/beta_e ≡ M_e^{-1} Curl^T M_f B  */
        /* Back Left edge */
        //if (er == 0 && ez == 0) {arrF[ez][ephi][er][ivErmzm] = ( arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * rmzmedgelength / betaenomp(er, ephi, ez, BACK_LEFT, user);} /* Back Left boundary edge */
        //else if (er == 0) {arrF[ez][ephi][er][ivErmzm] = ( arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) - arrX[ez-1][ephi][er][ivBrm] * betaf(er, ephi, ez-1, LEFT, user)/surface(er, ephi, ez-1, LEFT, user)) * rmzmedgelength / betaenomp(er, ephi, ez, BACK_LEFT, user);} /* Left boundary edge */
        //else if (ez == 0) {arrF[ez][ephi][er][ivErmzm] = ( arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) + arrX[ez][ephi][er-1][ivBzm] * betaf(er-1, ephi, ez, BACK, user)/surface(er-1, ephi, ez, BACK, user)) * rmzmedgelength / betaenomp(er, ephi, ez, BACK_LEFT, user);} /* Back boundary edge */
        if (er != 0 && ez != 0) {
          arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
            arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
            arrX[ez - 1][ephi][er][ivBrm] * betaf(er, ephi, ez - 1, LEFT, user) / surface(er, ephi, ez - 1, LEFT, user) +
            arrX[ez][ephi][er - 1][ivBzm] * betaf(er - 1, ephi, ez, BACK, user) / surface(er - 1, ephi, ez, BACK, user)) * rmzmedgelength / betaenomp(er, ephi, ez, BACK_LEFT, user);
        } /* Internal edge */

        if (!(user -> phibtype)) {
          /* Back Down edge */
          //if (ephi == 0 && ez == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * phimzmedgelength / betaenomp(er, ephi, ez, BACK_DOWN, user);} /* Back Down boundary edge */
          //else if (ephi == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) + arrX[ez-1][ephi][er][ivBphim] * betaf(er, ephi, ez-1, DOWN, user)/surface(er, ephi, ez-1, DOWN, user)) * phimzmedgelength / betaenomp(er, ephi, ez, BACK_DOWN, user);} /* Down boundary edge */
          //else if (ez == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) - arrX[ez][ephi-1][er][ivBzm] * betaf(er, ephi-1, ez, BACK, user)/surface(er, ephi-1, ez, BACK, user)) * phimzmedgelength / betaenomp(er, ephi, ez, BACK_DOWN, user);} /* Back boundary edge */
          if (ephi != 0 && ez != 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
              arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betaenomp(er, ephi, ez, BACK_DOWN, user);
          } /* Internal edge */

          /* Down Left edge */
          //if (er == 0 && ephi == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user)) * rmphimedgelength / betaenomp(er, ephi, ez, DOWN_LEFT, user);} /* Down Left boundary edge */
          //else if (er == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi-1][er][ivBrm] * betaf(er, ephi-1, ez, LEFT, user)/surface(er, ephi-1, ez, LEFT, user)) * rmphimedgelength / betaenomp(er, ephi, ez, DOWN_LEFT, user);} /* Left boundary edge */
          //else if (ephi == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er-1][ivBphim] * betaf(er-1, ephi, ez, DOWN, user)/surface(er-1, ephi, ez, DOWN, user)) * rmphimedgelength / betaenomp(er, ephi, ez, DOWN_LEFT, user);} /* Down boundary edge */
          if (ephi != 0 && ez != 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
              arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betaenomp(er, ephi, ez, DOWN_LEFT, user);
          } /* Internal edge */
        } else {
          //if (ez == 0) {arrF[ez][ephi][er][ivEphimzm] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) - arrX[ez][ephi-1][er][ivBzm] * betaf(er, ephi-1, ez, BACK, user)/surface(er, ephi-1, ez, BACK, user)) * phimzmedgelength / betaenomp(er, ephi, ez, BACK_DOWN, user);} /* Back boundary edge */
          if (ez != 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
              arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betaenomp(er, ephi, ez, BACK_DOWN, user);
          } /* Internal edge */

          //if (er == 0) {arrF[ez][ephi][er][ivErmphim] = ( -arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi-1][er][ivBrm] * betaf(er, ephi-1, ez, LEFT, user)/surface(er, ephi-1, ez, LEFT, user)) * rmphimedgelength / betaenomp(er, ephi, ez, DOWN_LEFT, user);} /* Left boundary edge */
          if (er != 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
              arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betaenomp(er, ephi, ez, DOWN_LEFT, user);
          } /* Internal edge */
        }

        if (er == N[0] - 1) {
          /* Right boundary Back edge */
          //if (ez == 0) {arrF[ez][ephi][er][ivErpzm] = ( arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * rpzmedgelength / betaenomp(er, ephi, ez, BACK_RIGHT, user);} /* Back Right boundary edge */
          //else {arrF[ez][ephi][er][ivErpzm] = ( arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) - arrX[ez-1][ephi][er][ivBrp] * betaf(er, ephi, ez-1, RIGHT, user)/surface(er, ephi, ez-1, RIGHT, user) + arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * rpzmedgelength / betaenomp(er, ephi, ez, BACK_RIGHT, user);} /* Right boundary edge */

          /* Right boundary Down edge */
          //if (ephi == 0 && !(user->phibtype)) {arrF[ez][ephi][er][ivErpphim] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user)) * rpphimedgelength / betaenomp(er, ephi, ez, DOWN_RIGHT, user);} /* Down Right boundary edge */
          //else {arrF[ez][ephi][er][ivErpphim] = ( -arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi-1][er][ivBrp] * betaf(er, ephi-1, ez, RIGHT, user)/surface(er, ephi-1, ez, RIGHT, user)) * rpphimedgelength / betaenomp(er, ephi, ez, DOWN_RIGHT, user);} /* Right boundary edge */
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Up boundary Back edge */
          //if (ez == 0) {arrF[ez][ephi][er][ivEphipzm] = ( -arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user)) * phipzmedgelength / betaenomp(er, ephi, ez, BACK_UP, user);} /* Back Up boundary edge */
          //else {arrF[ez][ephi][er][ivEphipzm] = ( -arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user)/surface(er, ephi, ez, BACK, user) + arrX[ez-1][ephi][er][ivBphip] * betaf(er, ephi, ez-1, UP, user)/surface(er, ephi, ez-1, UP, user)) * phipzmedgelength / betaenomp(er, ephi, ez, BACK_UP, user);} /* Up boundary edge */

          /* Up boundary Left edge */
          //if (er == 0) {arrF[ez][ephi][er][ivErmphip] = ( arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user)) * rmphipedgelength / betaenomp(er, ephi, ez, UP_LEFT, user);} /* Up Left boundary edge */
          //else {arrF[ez][ephi][er][ivErmphip] = ( arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er-1][ivBphip] * betaf(er-1, ephi, ez, UP, user)/surface(er-1, ephi, ez, UP, user)) * rmphipedgelength / betaenomp(er, ephi, ez, UP_LEFT, user);} /* Up boundary edge */
        }

        if (ez == N[2] - 1) {
          /* Front boundary Left edge */
          //if (er == 0) {arrF[ez][ephi][er][ivErmzp] = ( -arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user)) * rmzpedgelength / betaenomp(er, ephi, ez, FRONT_LEFT, user);} /* Front Left boundary edge */
          //else {arrF[ez][ephi][er][ivErmzp] = ( -arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user)/surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er-1][ivBzp] * betaf(er-1, ephi, ez, FRONT, user)/surface(er-1, ephi, ez, FRONT, user)) * rmzpedgelength / betaenomp(er, ephi, ez, FRONT_LEFT, user);} /* Front boundary edge */

          /* Front boundary Down edge */
          //if (ephi == 0 && !(user->phibtype)) {arrF[ez][ephi][er][ivEphimzp] = ( arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user)) * phimzpedgelength / betaenomp(er, ephi, ez, FRONT_DOWN, user);} /* Front Down boundary edge */
          //else {arrF[ez][ephi][er][ivEphimzp] = ( arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user)/surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi-1][er][ivBzp] * betaf(er, ephi-1, ez, FRONT, user)/surface(er, ephi-1, ez, FRONT, user)) * phimzpedgelength / betaenomp(er, ephi, ez, FRONT_DOWN, user);} /* Front boundary edge */
        }

        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Up Right boundary edge */
          //arrF[ez][ephi][er][ivErpphip] = ( -arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user)) * rpphipedgelength / betaenomp(er, ephi, ez, UP_RIGHT, user);
        }
        if (ephi == N[1] - 1 && ez == N[2] - 1 && !(user -> phibtype)) {
          /* Front Up boundary edge */
          //arrF[ez][ephi][er][ivEphipzp] = ( arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user)/surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user)) * phipzpedgelength / betaenomp(er, ephi, ez, FRONT_UP, user);
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Front Right boundary edge */
          //arrF[ez][ephi][er][ivErpzp] = ( -arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user)/surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user)/surface(er, ephi, ez, FRONT, user)) * rpzpedgelength / betaenomp(er, ephi, ez, FRONT_RIGHT, user);
        }
      }
    }
  }

  /* Restore vectors */
  /* DMStagVecRestoreArray(da,F,&arrF); */

  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode FormDerivedCurlExt(Vec X, Vec F, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec fLocal, xLocal;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];

  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];
  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  DM dmCoord;
  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda;

  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;

  da = coordDA;
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);
  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);
  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Compute function over the locally owned part of the grid */
  /* f(B) = derived_mimetic_curl(B) */
  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }
        /*PetscPrintf(PETSC_COMM_WORLD,"dphi = %g\n",(double)user->dphi);
        PetscPrintf(PETSC_COMM_WORLD,"Phip - Phim at CELL(%d,%d,%d) = %g\n",(double)(arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]),er,ephi,ez);*/

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        /* f(B) = primary_mimetic_curl^T (beta_f B)/beta_e ≡ M_e^{-1} Curl^T M_f B  */
        /* Back Left edge */
        if (er == 0 && ez == 0) {
          arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
            arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user)) * rmzmedgelength / betae(er, ephi, ez, BACK_LEFT, user);
        } /* Back Left boundary edge */
        else if (er == 0) {
          arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
            arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
            arrX[ez - 1][ephi][er][ivBrm] * betaf(er, ephi, ez - 1, LEFT, user) / surface(er, ephi, ez - 1, LEFT, user)) * rmzmedgelength / betae(er, ephi, ez, BACK_LEFT, user);
        } /* Left boundary edge */
        else if (ez == 0) {
          arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
            arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
            arrX[ez][ephi][er - 1][ivBzm] * betaf(er - 1, ephi, ez, BACK, user) / surface(er - 1, ephi, ez, BACK, user)) * rmzmedgelength / betae(er, ephi, ez, BACK_LEFT, user);
        } /* Back boundary edge */
        else {
          arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) -
            arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
            arrX[ez - 1][ephi][er][ivBrm] * betaf(er, ephi, ez - 1, LEFT, user) / surface(er, ephi, ez - 1, LEFT, user) +
            arrX[ez][ephi][er - 1][ivBzm] * betaf(er - 1, ephi, ez, BACK, user) / surface(er - 1, ephi, ez, BACK, user)) * rmzmedgelength / betae(er, ephi, ez, BACK_LEFT, user);
        } /* Internal edge */

        if (!(user -> phibtype)) {
          /* Back Down edge */
          if (ephi == 0 && ez == 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);
          } /* Back Down boundary edge */
          else if (ephi == 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
              arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);
          } /* Down boundary edge */
          else if (ez == 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);
          } /* Back boundary edge */
          else {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
              arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);
          } /* Internal edge */

          /* Down Left edge */
          if (er == 0 && ephi == 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);
          } /* Down Left boundary edge */
          else if (er == 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);
          } /* Left boundary edge */
          else if (ephi == 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) -
              arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);
          } /* Down boundary edge */
          else {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
              arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);
          } /* Internal edge */
        } else {
          if (ez == 0) {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);
          } /* Back boundary edge */
          else {
            arrF[ez][ephi][er][ivEphimzm] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) +
              arrX[ez - 1][ephi][er][ivBphim] * betaf(er, ephi, ez - 1, DOWN, user) / surface(er, ephi, ez - 1, DOWN, user) -
              arrX[ez][ephi - 1][er][ivBzm] * betaf(er, ephi - 1, ez, BACK, user) / surface(er, ephi - 1, ez, BACK, user)) * phimzmedgelength / betae(er, ephi, ez, BACK_DOWN, user);
          } /* Internal edge */

          if (er == 0) {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);
          } /* Left boundary edge */
          else {
            arrF[ez][ephi][er][ivErmphim] = (-arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) +
              arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) +
              arrX[ez][ephi - 1][er][ivBrm] * betaf(er, ephi - 1, ez, LEFT, user) / surface(er, ephi - 1, ez, LEFT, user) -
              arrX[ez][ephi][er - 1][ivBphim] * betaf(er - 1, ephi, ez, DOWN, user) / surface(er - 1, ephi, ez, DOWN, user)) * rmphimedgelength / betae(er, ephi, ez, DOWN_LEFT, user);
          } /* Internal edge */
        }

        if (er == N[0] - 1) {
          /* Right boundary Back edge */
          if (ez == 0) {
            arrF[ez][ephi][er][ivErpzm] = (arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user) / surface(er, ephi, ez, RIGHT, user) +
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user)) * rpzmedgelength / betae(er, ephi, ez, BACK_RIGHT, user);
          } /* Back Right boundary edge */
          else {
            arrF[ez][ephi][er][ivErpzm] = (arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user) / surface(er, ephi, ez, RIGHT, user) -
              arrX[ez - 1][ephi][er][ivBrp] * betaf(er, ephi, ez - 1, RIGHT, user) / surface(er, ephi, ez - 1, RIGHT, user) -
              arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user)) * rpzmedgelength / betae(er, ephi, ez, BACK_RIGHT, user);
          } /* Right boundary edge */

          /* Right boundary Down edge */
          if (ephi == 0 && !(user -> phibtype)) {
            arrF[ez][ephi][er][ivErpphim] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) -
              arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user) / surface(er, ephi, ez, RIGHT, user)) * rpphimedgelength / betae(er, ephi, ez, DOWN_RIGHT, user);
          } /* Down Right boundary edge */
          else {
            arrF[ez][ephi][er][ivErpphim] = (-arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) - arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user) / surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi - 1][er][ivBrp] * betaf(er, ephi - 1, ez, RIGHT, user) / surface(er, ephi - 1, ez, RIGHT, user)) * rpphimedgelength / betae(er, ephi, ez, DOWN_RIGHT, user);
          } /* Right boundary edge */
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Up boundary Back edge */
          if (ez == 0) {
            arrF[ez][ephi][er][ivEphipzm] = (-arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user) / surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user)) * phipzmedgelength / betae(er, ephi, ez, BACK_UP, user);
          } /* Back Up boundary edge */
          else {
            arrF[ez][ephi][er][ivEphipzm] = (-arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user) / surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzm] * betaf(er, ephi, ez, BACK, user) / surface(er, ephi, ez, BACK, user) + arrX[ez - 1][ephi][er][ivBphip] * betaf(er, ephi, ez - 1, UP, user) / surface(er, ephi, ez - 1, UP, user)) * phipzmedgelength / betae(er, ephi, ez, BACK_UP, user);
          } /* Up boundary edge */

          /* Up boundary Left edge */
          if (er == 0) {
            arrF[ez][ephi][er][ivErmphip] = (arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user) / surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user)) * rmphipedgelength / betae(er, ephi, ez, UP_LEFT, user);
          } /* Up Left boundary edge */
          else {
            arrF[ez][ephi][er][ivErmphip] = (arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user) / surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) - arrX[ez][ephi][er - 1][ivBphip] * betaf(er - 1, ephi, ez, UP, user) / surface(er - 1, ephi, ez, UP, user)) * rmphipedgelength / betae(er, ephi, ez, UP_LEFT, user);
          } /* Up boundary edge */
        }

        if (ez == N[2] - 1) {
          /* Front boundary Left edge */
          if (er == 0) {
            arrF[ez][ephi][er][ivErmzp] = (-arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user) / surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user)) * rmzpedgelength / betae(er, ephi, ez, FRONT_LEFT, user);
          } /* Front Left boundary edge */
          else {
            arrF[ez][ephi][er][ivErmzp] = (-arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user) / surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi][er][ivBrm] * betaf(er, ephi, ez, LEFT, user) / surface(er, ephi, ez, LEFT, user) + arrX[ez][ephi][er - 1][ivBzp] * betaf(er - 1, ephi, ez, FRONT, user) / surface(er - 1, ephi, ez, FRONT, user)) * rmzpedgelength / betae(er, ephi, ez, FRONT_LEFT, user);
          } /* Front boundary edge */

          /* Front boundary Down edge */
          if (ephi == 0 && !(user -> phibtype)) {
            arrF[ez][ephi][er][ivEphimzp] = (arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user) / surface(er, ephi, ez, FRONT, user)) * phimzpedgelength / betae(er, ephi, ez, FRONT_DOWN, user);
          } /* Front Down boundary edge */
          else {
            arrF[ez][ephi][er][ivEphimzp] = (arrX[ez][ephi][er][ivBphim] * betaf(er, ephi, ez, DOWN, user) / surface(er, ephi, ez, DOWN, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user) / surface(er, ephi, ez, FRONT, user) - arrX[ez][ephi - 1][er][ivBzp] * betaf(er, ephi - 1, ez, FRONT, user) / surface(er, ephi - 1, ez, FRONT, user)) * phimzpedgelength / betae(er, ephi, ez, FRONT_DOWN, user);
          } /* Front boundary edge */
        }

        if (er == N[0] - 1 && ephi == N[1] - 1 && !(user -> phibtype)) {
          /* Up Right boundary edge */
          arrF[ez][ephi][er][ivErpphip] = (-arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user) / surface(er, ephi, ez, UP, user) + arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user) / surface(er, ephi, ez, RIGHT, user)) * rpphipedgelength / betae(er, ephi, ez, UP_RIGHT, user);
        }
        if (ephi == N[1] - 1 && ez == N[2] - 1 && !(user -> phibtype)) {
          /* Front Up boundary edge */
          arrF[ez][ephi][er][ivEphipzp] = (arrX[ez][ephi][er][ivBphip] * betaf(er, ephi, ez, UP, user) / surface(er, ephi, ez, UP, user) - arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user) / surface(er, ephi, ez, FRONT, user)) * phipzpedgelength / betae(er, ephi, ez, FRONT_UP, user);
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          /* Front Right boundary edge */
          arrF[ez][ephi][er][ivErpzp] = (-arrX[ez][ephi][er][ivBrp] * betaf(er, ephi, ez, RIGHT, user) / surface(er, ephi, ez, RIGHT, user) + arrX[ez][ephi][er][ivBzp] * betaf(er, ephi, ez, FRONT, user) / surface(er, ephi, ez, FRONT, user)) * rpzpedgelength / betae(er, ephi, ez, FRONT_RIGHT, user);
        }
      }
    }
  }

  /* Restore vectors */
  /* DMStagVecRestoreArray(da,F,&arrF); */

  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  return (0);
}

PetscErrorCode FormSourceTermPotential(TS ts, PetscReal time, Vec P, void * ptr) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FormSourceTermPotential",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da;
  PetscInt startr, startphi, startz, nr, nphi, nz;

  Vec pLocal;
  Vec coordLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];

  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];
  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];

  PetscInt ivBrp, ivBphip, ivBzp, ivBrm, ivBphim, ivBzm;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  DM dmCoord;
  PetscScalar ** ** arrCoord, ** ** arrP;
  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  /*DMCreateGlobalVector(da,&X);*/
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);
  /* Face locations */
  DMStagGetLocationSlot(da, LEFT, 0, & ivBrm);
  DMStagGetLocationSlot(da, DOWN, 0, & ivBphim);
  DMStagGetLocationSlot(da, BACK, 0, & ivBzm);
  DMStagGetLocationSlot(da, RIGHT, 0, & ivBrp);
  DMStagGetLocationSlot(da, UP, 0, & ivBphip);
  DMStagGetLocationSlot(da, FRONT, 0, & ivBzp);
  DMGetCoordinateDM(da, & dmCoord);
  DMGetCoordinatesLocal(da, & coordLocal);
  DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);
  for (d = 0; d < 3; ++d) {
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoord, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoord, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoord, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoord, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT, d, & icBzp[d]);
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoord, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoord, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoord, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoord, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoord, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoord, FRONT_UP, d, & icEphipzp[d]);
  }
  /* Compute function over the locally owned part of the grid */
  DMGetLocalVector(da, & pLocal);
  DMStagVecGetArray(da, pLocal, & arrP);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        if (user -> ictype == 1) {
          /* Nothing on faces */
          arrP[ez][ephi][er][ivBrm] = 0.0;
          arrP[ez][ephi][er][ivBphim] = 0.0;
          arrP[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t) = (cos(pi*t)*(2*z-1)/r) e_phi - (cos(pi*t)*cos(phi)/r^2) e_z
         A(r, phi, z, t) = pi*sin(pi*t)*(cos(phi)-phi*(z-0.5)^2) e_z
         P = A + E */
          arrP[ez][ephi][er][ivErmzm] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErmzm[2]] - 1.0) / arrCoord[ez][ephi][er][icErmzm[0]];
          arrP[ez][ephi][er][ivEphimzm] = 0.0;
          arrP[ez][ephi][er][ivErmphim] = -PetscCosReal(arrCoord[ez][ephi][er][icErmphim[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErmphim[0]]) +
            PETSC_PI * PetscSinReal(PETSC_PI * time) * (PetscCosReal(arrCoord[ez][ephi][er][icErmphim[1]]) - arrCoord[ez][ephi][er][icErmphim[1]] * PetscSqr(arrCoord[ez][ephi][er][icErmphim[2]] - 0.5));

          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivErpzm] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErpzm[2]] - 1.0) / arrCoord[ez][ephi][er][icErpzm[0]];
            arrP[ez][ephi][er][ivErpphim] = -PetscCosReal(arrCoord[ez][ephi][er][icErpphim[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErpphim[0]]) +
              PETSC_PI * PetscSinReal(PETSC_PI * time) * (PetscCosReal(arrCoord[ez][ephi][er][icErpphim[1]]) - arrCoord[ez][ephi][er][icErpphim[1]] * PetscSqr(arrCoord[ez][ephi][er][icErpphim[2]] - 0.5));
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivEphipzm] = 0.0;
            arrP[ez][ephi][er][ivErmphip] = -PetscCosReal(arrCoord[ez][ephi][er][icErmphip[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErmphip[0]]) +
              PETSC_PI * PetscSinReal(PETSC_PI * time) * (PetscCosReal(arrCoord[ez][ephi][er][icErmphip[1]]) - arrCoord[ez][ephi][er][icErmphip[1]] * PetscSqr(arrCoord[ez][ephi][er][icErmphip[2]] - 0.5));
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErmzp] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErmzp[2]] - 1.0) / arrCoord[ez][ephi][er][icErmzp[0]];
            arrP[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivErpphip] = -PetscCosReal(arrCoord[ez][ephi][er][icErpphip[1]]) * PetscCosReal(PETSC_PI * time) / PetscSqr(arrCoord[ez][ephi][er][icErpphip[0]]) +
              PETSC_PI * PetscSinReal(PETSC_PI * time) * (PetscCosReal(arrCoord[ez][ephi][er][icErpphip[1]]) - arrCoord[ez][ephi][er][icErpphip[1]] * PetscSqr(arrCoord[ez][ephi][er][icErpphip[2]] - 0.5));
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErpzp] = PetscCosReal(PETSC_PI * time) * (2.0 * arrCoord[ez][ephi][er][icErpzp[2]] - 1.0) / arrCoord[ez][ephi][er][icErpzp[0]];
          }
        } else if (user -> ictype == 2) {
          /* Nothing on faces */
          arrP[ez][ephi][er][ivBrm] = 0.0;
          arrP[ez][ephi][er][ivBphim] = 0.0;
          arrP[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t) = (pi*cos(pi*r) + sin(pi*r)/r) * exp(-t) e_z
       A(r, phi, z, t) = - (cos(pi*r)/pi) * exp(-t) e_z
       P = A + E */
          arrP[ez][ephi][er][ivErmzm] = 0.0;
          arrP[ez][ephi][er][ivEphimzm] = 0.0;
          arrP[ez][ephi][er][ivErmphim] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphim[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphim[0]]) / arrCoord[ez][ephi][er][icErmphim[0]]) * PetscExpReal(-time) -
            PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphim[0]]) * PetscExpReal(-time) / PETSC_PI;
          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivErpzm] = 0.0;
            arrP[ez][ephi][er][ivErpphim] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphim[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphim[0]]) / arrCoord[ez][ephi][er][icErpphim[0]]) * PetscExpReal(-time) -
              PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphim[0]]) * PetscExpReal(-time) / PETSC_PI;
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivEphipzm] = 0.0;
            arrP[ez][ephi][er][ivErmphip] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphip[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphip[0]]) / arrCoord[ez][ephi][er][icErmphip[0]]) * PetscExpReal(-time) -
              PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErmphip[0]]) * PetscExpReal(-time) / PETSC_PI;
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErmzp] = 0.0;
            arrP[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivErpphip] = (PETSC_PI * PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphip[0]]) + PetscSinReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphip[0]]) / arrCoord[ez][ephi][er][icErpphip[0]]) * PetscExpReal(-time) -
              PetscCosReal(PETSC_PI * arrCoord[ez][ephi][er][icErpphip[0]]) * PetscExpReal(-time) / PETSC_PI;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErpzp] = 0.0;
          }
        } else if (user -> ictype == 7) {
          /* Nothing on faces */
          arrP[ez][ephi][er][ivBrm] = 0.0;
          arrP[ez][ephi][er][ivBphim] = 0.0;
          arrP[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t=0) = (cos(phi)/r) e_r
             A(r, phi, z, t) = - (pi * r * sin(pi*t)/2) e_phi
             P = A + E */
          arrP[ez][ephi][er][ivErmzm] = -PETSC_PI * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSinReal(PETSC_PI * time) / 2.0;
          arrP[ez][ephi][er][ivEphimzm] = PetscCosReal(arrCoord[ez][ephi][er][icEphimzm[1]]) / arrCoord[ez][ephi][er][icEphimzm[0]];
          arrP[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivErpzm] = -PETSC_PI * arrCoord[ez][ephi][er][icErpzm[0]] * PetscSinReal(PETSC_PI * time) / 2.0;
            arrP[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivEphipzm] = PetscCosReal(arrCoord[ez][ephi][er][icEphipzm[1]]) / arrCoord[ez][ephi][er][icEphipzm[0]];
            arrP[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErmzp] = -PETSC_PI * arrCoord[ez][ephi][er][icErmzp[0]] * PetscSinReal(PETSC_PI * time) / 2.0;
            arrP[ez][ephi][er][ivEphimzp] = PetscCosReal(arrCoord[ez][ephi][er][icEphimzp[1]]) / arrCoord[ez][ephi][er][icEphimzp[0]];
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivEphipzp] = PetscCosReal(arrCoord[ez][ephi][er][icEphipzp[1]]) / arrCoord[ez][ephi][er][icEphipzp[0]];
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErpzp] = -PETSC_PI * arrCoord[ez][ephi][er][icErpzp[0]] * PetscSinReal(PETSC_PI * time) / 2.0;
          }
        } else if (user -> ictype == 3 || user -> ictype == 4 || user -> ictype == 5 || user -> ictype == 8 || user -> ictype == 9 || user -> ictype == 10 || user -> ictype == 12 || user -> ictype == 13) {
          /* Nothing on faces */
          arrP[ez][ephi][er][ivBrm] = 0.0;
          arrP[ez][ephi][er][ivBphim] = 0.0;
          arrP[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t) = 0
           A(r, phi, z, t) = 0
           P = A + E */
          arrP[ez][ephi][er][ivErmzm] = 0.0;
          arrP[ez][ephi][er][ivEphimzm] = 0.0;
          arrP[ez][ephi][er][ivErmphim] = 0.0;
          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivErpzm] = 0.0;
            arrP[ez][ephi][er][ivErpphim] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivEphipzm] = 0.0;
            arrP[ez][ephi][er][ivErmphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErmzp] = 0.0;
            arrP[ez][ephi][er][ivEphimzp] = 0.0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivErpphip] = 0.0;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivEphipzp] = 0.0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErpzp] = 0.0;
          }
        } else if (user -> ictype == 11) {
          /* Nothing on faces */
          arrP[ez][ephi][er][ivBrm] = 0.0;
          arrP[ez][ephi][er][ivBphim] = 0.0;
          arrP[ez][ephi][er][ivBzm] = 0.0;
          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivBrp] = 0.0;
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivBphip] = 0.0;
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivBzp] = 0.0;
          }
          /* E(r, phi, z, t) = (1/V_A)*(eta/mu0)*(cos(eta/mu0 * t)/r)*(e_z + [1 / sqrt(2 * ln(2*r_max/r))] e_phi)
           A(r, phi, z, t) = - (eta/mu0)*(sin(eta/mu0 * t) * (z - r*phi*sqrt(2 * ln(2*r_max/r))) e_r
           P = A + E */
          arrP[ez][ephi][er][ivErmzm] = condu(er, ephi, ez, BACK_LEFT, user) * PetscCosReal((condu(er, ephi, ez, BACK_LEFT, user) / user->mu0) * time) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzm[0]])));
          arrP[ez][ephi][er][ivEphimzm] = - condu(er, ephi, ez, BACK_DOWN, user) * PetscSinReal((condu(er, ephi, ez, BACK_DOWN, user) / user->mu0) * time) * (arrCoord[ez][ephi][er][icEphimzm[2]] - arrCoord[ez][ephi][er][icEphimzm[0]] * arrCoord[ez][ephi][er][icEphimzm[1]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icEphimzm[0]]))) / user->mu0;
          arrP[ez][ephi][er][ivErmphim] = condu(er, ephi, ez, DOWN_LEFT, user) * PetscCosReal((condu(er, ephi, ez, DOWN_LEFT, user) / user->mu0) * time) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphim[0]]);
          if (er == N[0] - 1) {
            arrP[ez][ephi][er][ivErpzm] = condu(er, ephi, ez, BACK_RIGHT, user) * PetscCosReal((condu(er, ephi, ez, BACK_RIGHT, user) / user->mu0) * time) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzm[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzm[0]])));
            arrP[ez][ephi][er][ivErpphim] = condu(er, ephi, ez, DOWN_RIGHT, user) * PetscCosReal((condu(er, ephi, ez, DOWN_RIGHT, user) / user->mu0) * time) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphim[0]]);
          }
          if (ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivEphipzm] = - condu(er, ephi, ez, BACK_UP, user) * PetscSinReal((condu(er, ephi, ez, BACK_UP, user) / user->mu0) * time) * (arrCoord[ez][ephi][er][icEphipzm[2]] - arrCoord[ez][ephi][er][icEphipzm[0]] * arrCoord[ez][ephi][er][icEphipzm[1]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icEphipzm[0]]))) / user->mu0;
            arrP[ez][ephi][er][ivErmphip] = condu(er, ephi, ez, UP_LEFT, user) * PetscCosReal((condu(er, ephi, ez, UP_LEFT, user) / user->mu0) * time) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmphip[0]]);
          }
          if (ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErmzp] = condu(er, ephi, ez, FRONT_LEFT, user) * PetscCosReal((condu(er, ephi, ez, FRONT_LEFT, user) / user->mu0) * time) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErmzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErmzp[0]])));
            arrP[ez][ephi][er][ivEphimzp] = - condu(er, ephi, ez, FRONT_DOWN, user) * PetscSinReal((condu(er, ephi, ez, FRONT_DOWN, user) / user->mu0) * time) * (arrCoord[ez][ephi][er][icEphimzp[2]] - arrCoord[ez][ephi][er][icEphimzp[0]] * arrCoord[ez][ephi][er][icEphimzp[1]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icEphimzp[0]]))) / user->mu0;
          }
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrP[ez][ephi][er][ivErpphip] = condu(er, ephi, ez, UP_RIGHT, user) * PetscCosReal((condu(er, ephi, ez, UP_RIGHT, user) / user->mu0) * time) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivEphipzp] = - condu(er, ephi, ez, FRONT_UP, user) * PetscSinReal((condu(er, ephi, ez, FRONT_UP, user) / user->mu0) * time) * (arrCoord[ez][ephi][er][icEphipzp[2]] - arrCoord[ez][ephi][er][icEphipzp[0]] * arrCoord[ez][ephi][er][icEphipzp[1]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icEphipzp[0]]))) / user->mu0;
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrP[ez][ephi][er][ivErpzp] = condu(er, ephi, ez, FRONT_RIGHT, user) * PetscCosReal((condu(er, ephi, ez, FRONT_RIGHT, user) / user->mu0) * time) / (11000000.0 * user->mu0 * arrCoord[ez][ephi][er][icErpzp[0]] * PetscSqrtScalar(2.0 * PetscLogReal(2.0 * user->rmax / arrCoord[ez][ephi][er][icErpzp[0]])));
          }
        } else {
          PetscPrintf(PETSC_COMM_WORLD, "Test case not set\n");
        }
      }
    }
  }
  /* Restore vectors */
  DMStagVecRestoreArray(da, pLocal, & arrP);
  DMLocalToGlobal(da, pLocal, INSERT_VALUES, P);
  DMRestoreLocalVector(da, & pLocal);
  DMStagVecRestoreArrayRead(dmCoord, coordLocal, & arrCoord);

  if (user -> debug) {
    /*This print is just for debugging*/
    PetscPrintf(PETSC_COMM_WORLD, "Source term potential vector\n");
    VecView(P, PETSC_VIEWER_STDOUT_WORLD);
  }

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode FormDiscreteGradientEP(TS ts, Mat G, Vec X, Vec F, void * ptr) {
  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec coordLocal, fLocal, xLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;
  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX;

  VecZeroEntries(F);
  MatZeroEntries(G);
  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  for (d = 0; d < 4; ++d) {
    /* Vertex locations */
    DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]);
    DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]);
  }
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);

  for (d = 0; d < 3; ++d) {
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil row, col[2];
        PetscScalar valG[2];
        PetscInt nEntries = 2;

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivVrmphipzm[3]] - arrX[ez][ephi][er][ivVrmphimzm[3]]) / rmzmedgelength;
        arrF[ez][ephi][er][ivErmphim] = (arrX[ez][ephi][er][ivVrmphimzp[3]] - arrX[ez][ephi][er][ivVrmphimzm[3]]) / rmphimedgelength;
        arrF[ez][ephi][er][ivEphimzm] = (arrX[ez][ephi][er][ivVrpphimzm[3]] - arrX[ez][ephi][er][ivVrmphimzm[3]]) / phimzmedgelength;
        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivErpzm] = (arrX[ez][ephi][er][ivVrpphipzm[3]] - arrX[ez][ephi][er][ivVrpphimzm[3]]) / rpzmedgelength;
          arrF[ez][ephi][er][ivErpphim] = (arrX[ez][ephi][er][ivVrpphimzp[3]] - arrX[ez][ephi][er][ivVrpphimzm[3]]) / rpphimedgelength;
        }
        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivErmzp] = (arrX[ez][ephi][er][ivVrmphipzp[3]] - arrX[ez][ephi][er][ivVrmphimzp[3]]) / rmzpedgelength;
          arrF[ez][ephi][er][ivEphimzp] = (arrX[ez][ephi][er][ivVrpphimzp[3]] - arrX[ez][ephi][er][ivVrmphimzp[3]]) / phimzpedgelength;
        }
        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivEphipzm] = (arrX[ez][ephi][er][ivVrpphipzm[3]] - arrX[ez][ephi][er][ivVrmphipzm[3]]) / phipzmedgelength;
          arrF[ez][ephi][er][ivErmphip] = (arrX[ez][ephi][er][ivVrmphipzp[3]] - arrX[ez][ephi][er][ivVrmphipzm[3]]) / rmphipedgelength;
        }
        if (!(user -> phibtype)) {
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrF[ez][ephi][er][ivErpphip] = (arrX[ez][ephi][er][ivVrpphipzp[3]] - arrX[ez][ephi][er][ivVrpphipzm[3]]) / rpphipedgelength;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrF[ez][ephi][er][ivEphipzp] = (arrX[ez][ephi][er][ivVrpphipzp[3]] - arrX[ez][ephi][er][ivVrmphipzp[3]]) / phipzpedgelength;
          }
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrF[ez][ephi][er][ivErpzp] = (arrX[ez][ephi][er][ivVrpphipzp[3]] - arrX[ez][ephi][er][ivVrpphimzp[3]]) / rpzpedgelength;
        }

        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_LEFT;
        row.c = 0;
        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = BACK_DOWN_LEFT;
        col[0].c = 3;
        valG[0] = - 1.0 / rmzmedgelength;
        col[1].i = er;
        col[1].j = ephi;
        col[1].k = ez;
        col[1].loc = BACK_UP_LEFT;
        col[1].c = 3;
        valG[1] = 1.0 / rmzmedgelength;
        DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);

        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = DOWN_LEFT;
        row.c = 0;
        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = BACK_DOWN_LEFT;
        col[0].c = 3;
        valG[0] = - 1.0 / rmphimedgelength;
        col[1].i = er;
        col[1].j = ephi;
        col[1].k = ez;
        col[1].loc = FRONT_DOWN_LEFT;
        col[1].c = 3;
        valG[1] = 1.0 / rmphimedgelength;
        DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);

        row.i = er;
        row.j = ephi;
        row.k = ez;
        row.loc = BACK_DOWN;
        row.c = 0;
        col[0].i = er;
        col[0].j = ephi;
        col[0].k = ez;
        col[0].loc = BACK_DOWN_LEFT;
        col[0].c = 3;
        valG[0] = - 1.0 / phimzmedgelength;
        col[1].i = er;
        col[1].j = ephi;
        col[1].k = ez;
        col[1].loc = BACK_DOWN_RIGHT;
        col[1].c = 3;
        valG[1] = 1.0 / phimzmedgelength;
        DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);

        if (er == N[0] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_RIGHT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN_RIGHT;
          col[0].c = 3;
          valG[0] = - 1.0 / rpzmedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP_RIGHT;
          col[1].c = 3;
          valG[1] = 1.0 / rpzmedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = DOWN_RIGHT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_DOWN_RIGHT;
          col[0].c = 3;
          valG[0] = - 1.0 / rpphimedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_RIGHT;
          col[1].c = 3;
          valG[1] = 1.0 / rpphimedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }

        if (ez == N[2] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_LEFT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN_LEFT;
          col[0].c = 3;
          valG[0] = - 1.0 / rmzpedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_UP_LEFT;
          col[1].c = 3;
          valG[1] = 1.0 / rmzpedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_DOWN;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN_LEFT;
          col[0].c = 3;
          valG[0] = - 1.0 / phimzpedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_DOWN_RIGHT;
          col[1].c = 3;
          valG[1] = 1.0 / phimzpedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }

        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = BACK_UP;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_LEFT;
          col[0].c = 3;
          valG[0] = - 1.0 / phipzmedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = BACK_UP_RIGHT;
          col[1].c = 3;
          valG[1] = 1.0 / phipzmedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);

          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = UP_LEFT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = BACK_UP_LEFT;
          col[0].c = 3;
          valG[0] = - 1.0 / rmphipedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_UP_LEFT;
          col[1].c = 3;
          valG[1] = 1.0 / rmphipedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }

        if (!(user -> phibtype)) {
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            row.i = er;
            row.j = ephi;
            row.k = ez;
            row.loc = UP_RIGHT;
            row.c = 0;
            col[0].i = er;
            col[0].j = ephi;
            col[0].k = ez;
            col[0].loc = BACK_UP_RIGHT;
            col[0].c = 3;
            valG[0] = - 1.0 / rpphipedgelength;
            col[1].i = er;
            col[1].j = ephi;
            col[1].k = ez;
            col[1].loc = FRONT_UP_RIGHT;
            col[1].c = 3;
            valG[1] = 1.0 / rpphipedgelength;
            DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            row.i = er;
            row.j = ephi;
            row.k = ez;
            row.loc = FRONT_UP;
            row.c = 0;
            col[0].i = er;
            col[0].j = ephi;
            col[0].k = ez;
            col[0].loc = FRONT_UP_LEFT;
            col[0].c = 3;
            valG[0] = - 1.0 / phipzpedgelength;
            col[1].i = er;
            col[1].j = ephi;
            col[1].k = ez;
            col[1].loc = FRONT_UP_RIGHT;
            col[1].c = 3;
            valG[1] = 1.0 / phipzpedgelength;
            DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
          }
        }

        if (er == N[0] - 1 && ez == N[2] - 1) {
          row.i = er;
          row.j = ephi;
          row.k = ez;
          row.loc = FRONT_RIGHT;
          row.c = 0;
          col[0].i = er;
          col[0].j = ephi;
          col[0].k = ez;
          col[0].loc = FRONT_DOWN_RIGHT;
          col[0].c = 3;
          valG[0] = - 1.0 / rpzpedgelength;
          col[1].i = er;
          col[1].j = ephi;
          col[1].k = ez;
          col[1].loc = FRONT_UP_RIGHT;
          col[1].c = 3;
          valG[1] = 1.0 / rpzpedgelength;
          DMStagMatSetValuesStencil(coordDA, G, 1, & row, nEntries, col, valG, INSERT_VALUES);
        }
      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);
  MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);
  return (0);
}

PetscErrorCode FormDiscreteGradientEP_noMat(TS ts, Vec X, Vec F, void * ptr) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FormDiscreteGradientEP_noMat",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec coordLocal, fLocal, xLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;
  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX;

  VecZeroEntries(F);
  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  for (d = 0; d < 4; ++d) {
    /* Vertex locations */
    DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]);
    DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]);
  }
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);

  for (d = 0; d < 3; ++d) {
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        arrF[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivVrmphipzm[3]] - arrX[ez][ephi][er][ivVrmphimzm[3]]) / rmzmedgelength;
        arrF[ez][ephi][er][ivErmphim] = (arrX[ez][ephi][er][ivVrmphimzp[3]] - arrX[ez][ephi][er][ivVrmphimzm[3]]) / rmphimedgelength;
        arrF[ez][ephi][er][ivEphimzm] = (arrX[ez][ephi][er][ivVrpphimzm[3]] - arrX[ez][ephi][er][ivVrmphimzm[3]]) / phimzmedgelength;
        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivErpzm] = (arrX[ez][ephi][er][ivVrpphipzm[3]] - arrX[ez][ephi][er][ivVrpphimzm[3]]) / rpzmedgelength;
          arrF[ez][ephi][er][ivErpphim] = (arrX[ez][ephi][er][ivVrpphimzp[3]] - arrX[ez][ephi][er][ivVrpphimzm[3]]) / rpphimedgelength;
        }
        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivErmzp] = (arrX[ez][ephi][er][ivVrmphipzp[3]] - arrX[ez][ephi][er][ivVrmphimzp[3]]) / rmzpedgelength;
          arrF[ez][ephi][er][ivEphimzp] = (arrX[ez][ephi][er][ivVrpphimzp[3]] - arrX[ez][ephi][er][ivVrmphimzp[3]]) / phimzpedgelength;
        }
        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivEphipzm] = (arrX[ez][ephi][er][ivVrpphipzm[3]] - arrX[ez][ephi][er][ivVrmphipzm[3]]) / phipzmedgelength;
          arrF[ez][ephi][er][ivErmphip] = (arrX[ez][ephi][er][ivVrmphipzp[3]] - arrX[ez][ephi][er][ivVrmphipzm[3]]) / rmphipedgelength;
        }
        if (!(user -> phibtype)) {
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrF[ez][ephi][er][ivErpphip] = (arrX[ez][ephi][er][ivVrpphipzp[3]] - arrX[ez][ephi][er][ivVrpphipzm[3]]) / rpphipedgelength;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrF[ez][ephi][er][ivEphipzp] = (arrX[ez][ephi][er][ivVrpphipzp[3]] - arrX[ez][ephi][er][ivVrmphipzp[3]]) / phipzpedgelength;
          }
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrF[ez][ephi][er][ivErpzp] = (arrX[ez][ephi][er][ivVrpphipzp[3]] - arrX[ez][ephi][er][ivVrpphimzp[3]]) / rpzpedgelength;
        }


      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode FormDiscreteGradientEP_tilde(TS ts, Vec X, Vec F, void * ptr) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FormDiscreteGradientEP_tilde",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec coordLocal, fLocal, xLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];
  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];


  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;
  PetscScalar ** ** arrCoord, ** ** arrF, ** ** arrX;

  VecZeroEntries(F);
  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  DMGetLocalVector(da, & fLocal);
  DMGlobalToLocalBegin(da, F, INSERT_VALUES, fLocal);
  DMGlobalToLocalEnd(da, F, INSERT_VALUES, fLocal);
  DMStagVecGetArray(da, fLocal, & arrF);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  for (d = 0; d < 4; ++d) {
    /* Vertex locations */
    DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]);
    DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]);
  }
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);

  for (d = 0; d < 3; ++d) {
    /* Edge coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_LEFT, d, & icErmzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN, d, & icEphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_RIGHT, d, & icErpzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP, d, & icEphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, DOWN_LEFT, d, & icErmphim[d]);
    DMStagGetLocationSlot(dmCoorda, DOWN_RIGHT, d, & icErpphim[d]);
    DMStagGetLocationSlot(dmCoorda, UP_LEFT, d, & icErmphip[d]);
    DMStagGetLocationSlot(dmCoorda, UP_RIGHT, d, & icErpphip[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN, d, & icEphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_LEFT, d, & icErmzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_RIGHT, d, & icErpzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP, d, & icEphipzp[d]);
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* f(EP) = - (1/r) \tilde{nabla}(EP) = - (1/r) \tilde{prim_Grad}(EP) = (1/r^2) \partial EP / \partial r - (1/r) \partial^2 EP / \partial^2 r - (1/r) \partial^2 EP / \partial^2 z*/

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        arrF[ez][ephi][er][ivErmphim] = - (arrX[ez][ephi][er][ivVrmphimzp[3]] - arrX[ez][ephi][er][ivVrmphimzm[3]]) / (rmphimedgelength * arrCoorda[ez][ephi][er][icErmphim[0]]);
        arrF[ez][ephi][er][ivEphimzm] = - (arrX[ez][ephi][er][ivVrpphimzm[3]] - arrX[ez][ephi][er][ivVrmphimzm[3]]) / (phimzmedgelength * arrCoorda[ez][ephi][er][icEphimzm[0]]);
        if (er == N[0] - 1) {
          arrF[ez][ephi][er][ivErpphim] = - (arrX[ez][ephi][er][ivVrpphimzp[3]] - arrX[ez][ephi][er][ivVrpphimzm[3]]) / (rpphimedgelength * arrCoorda[ez][ephi][er][icErpphim[0]]);
        }
        if (ez == N[2] - 1) {
          arrF[ez][ephi][er][ivEphimzp] = - (arrX[ez][ephi][er][ivVrpphimzp[3]] - arrX[ez][ephi][er][ivVrmphimzp[3]]) / (phimzpedgelength * arrCoorda[ez][ephi][er][icEphimzp[0]]);
        }
        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF[ez][ephi][er][ivEphipzm] = - (arrX[ez][ephi][er][ivVrpphipzm[3]] - arrX[ez][ephi][er][ivVrmphipzm[3]]) / (phipzmedgelength * arrCoorda[ez][ephi][er][icEphipzm[0]]);
          arrF[ez][ephi][er][ivErmphip] = - (arrX[ez][ephi][er][ivVrmphipzp[3]] - arrX[ez][ephi][er][ivVrmphipzm[3]]) / (rmphipedgelength * arrCoorda[ez][ephi][er][icErmphip[0]]);
        }
        if (!(user -> phibtype)) {
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrF[ez][ephi][er][ivErpphip] = - (arrX[ez][ephi][er][ivVrpphipzp[3]] - arrX[ez][ephi][er][ivVrpphipzm[3]]) / (rpphipedgelength * arrCoorda[ez][ephi][er][icErpphip[0]]);
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrF[ez][ephi][er][ivEphipzp] = - (arrX[ez][ephi][er][ivVrpphipzp[3]] - arrX[ez][ephi][er][ivVrmphipzp[3]]) / (phipzpedgelength * arrCoorda[ez][ephi][er][icEphipzp[0]]);
          }
        }


      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArray(da, fLocal, & arrF);
  DMLocalToGlobal(da, fLocal, INSERT_VALUES, F);
  DMRestoreLocalVector(da, & fLocal);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode FormDiscreteGradientVectorField(TS ts, Vec X, Vec F1, Vec F2, Vec F3, void * ptr) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FormDiscreteGradientVectorField",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;
  DM da, coordDA = user -> coorda;
  PetscInt startr, startphi, startz, nr, nphi, nz;
  Vec coordLocal, f1Local, f2Local, f3Local, xLocal;
  PetscInt N[3], er, ephi, ez, d;

  PetscInt ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;

  PetscInt ivVrmphimzm[4], ivVrmphimzp[4], ivVrmphipzm[4], ivVrmphipzp[4];
  PetscInt ivVrpphimzm[4], ivVrpphipzm[4], ivVrpphipzp[4], ivVrpphimzp[4];

  PetscInt icrmphimzm[3], icrmphimzp[3], icrmphipzm[3], icrmphipzp[3];
  PetscInt icrpphimzm[3], icrpphimzp[3], icrpphipzm[3], icrpphipzp[3];

  DM dmCoorda;
  Vec coordaLocal;
  PetscScalar ** ** arrCoorda, rmzmedgelength, rmphimedgelength, rmzpedgelength, rmphipedgelength, phimzmedgelength, phimzpedgelength, rpphimedgelength, rpzmedgelength, phipzmedgelength, rpzpedgelength, rpphipedgelength, phipzpedgelength;
  PetscScalar ** ** arrCoord, ** ** arrF1, ** ** arrF2, ** ** arrF3, ** ** arrX;

  VecZeroEntries(F1);
  VecZeroEntries(F2);
  VecZeroEntries(F3);
  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  DMGetLocalVector(da, & f1Local);
  DMGlobalToLocalBegin(da, F1, INSERT_VALUES, f1Local);
  DMGlobalToLocalEnd(da, F1, INSERT_VALUES, f1Local);
  DMStagVecGetArray(da, f1Local, & arrF1);

  DMGetLocalVector(da, & f2Local);
  DMGlobalToLocalBegin(da, F2, INSERT_VALUES, f2Local);
  DMGlobalToLocalEnd(da, F2, INSERT_VALUES, f2Local);
  DMStagVecGetArray(da, f2Local, & arrF2);

  DMGetLocalVector(da, & f3Local);
  DMGlobalToLocalBegin(da, F3, INSERT_VALUES, f3Local);
  DMGlobalToLocalEnd(da, F3, INSERT_VALUES, f3Local);
  DMStagVecGetArray(da, f3Local, & arrF3);

  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArrayRead(da, xLocal, & arrX);

  DMGetCoordinateDM(coordDA, & dmCoorda);
  DMGetCoordinatesLocal(coordDA, & coordaLocal);
  DMStagVecGetArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  for (d = 0; d < 4; ++d) {
    /* Vertex locations */
    DMStagGetLocationSlot(da, BACK_DOWN_LEFT, d, & ivVrmphimzm[d]);
    DMStagGetLocationSlot(da, BACK_DOWN_RIGHT, d, & ivVrpphimzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_LEFT, d, & ivVrmphipzm[d]);
    DMStagGetLocationSlot(da, BACK_UP_RIGHT, d, & ivVrpphipzm[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_LEFT, d, & ivVrmphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_DOWN_RIGHT, d, & ivVrpphimzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_LEFT, d, & ivVrmphipzp[d]);
    DMStagGetLocationSlot(da, FRONT_UP_RIGHT, d, & ivVrpphipzp[d]);
  }
  /* Edge locations */
  DMStagGetLocationSlot(da, BACK_LEFT, 0, & ivErmzm);
  DMStagGetLocationSlot(da, BACK_DOWN, 0, & ivEphimzm);
  DMStagGetLocationSlot(da, BACK_RIGHT, 0, & ivErpzm);
  DMStagGetLocationSlot(da, BACK_UP, 0, & ivEphipzm);
  DMStagGetLocationSlot(da, DOWN_LEFT, 0, & ivErmphim);
  DMStagGetLocationSlot(da, DOWN_RIGHT, 0, & ivErpphim);
  DMStagGetLocationSlot(da, UP_LEFT, 0, & ivErmphip);
  DMStagGetLocationSlot(da, UP_RIGHT, 0, & ivErpphip);
  DMStagGetLocationSlot(da, FRONT_DOWN, 0, & ivEphimzp);
  DMStagGetLocationSlot(da, FRONT_LEFT, 0, & ivErmzp);
  DMStagGetLocationSlot(da, FRONT_RIGHT, 0, & ivErpzp);
  DMStagGetLocationSlot(da, FRONT_UP, 0, & ivEphipzp);

  for (d = 0; d < 3; ++d) {
    /* Vertex coordinates */
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_LEFT, d, & icrmphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_DOWN_RIGHT, d, & icrpphimzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_LEFT, d, & icrmphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, BACK_UP_RIGHT, d, & icrpphipzm[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_LEFT, d, & icrmphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_DOWN_RIGHT, d, & icrpphimzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_LEFT, d, & icrmphipzp[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT_UP_RIGHT, d, & icrpphipzp[d]);
  }

  /* Loop over all local elements */
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzm[0]];
        } else {
          rmzmedgelength = arrCoorda[ez][ephi][er][icrmphimzm[0]] * (arrCoorda[ez][ephi][er][icrmphipzm[1]] - arrCoorda[ez][ephi][er][icrmphimzm[1]]); /* back left = rmzm */
        }

        rmphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]]); /* down left = rmphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rmzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrmphimzp[0]];
        } else {
          rmzpedgelength = arrCoorda[ez][ephi][er][icrmphimzp[0]] * (arrCoorda[ez][ephi][er][icrmphipzp[1]] - arrCoorda[ez][ephi][er][icrmphimzp[1]]); /* front left = rmzp */
        }

        rmphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]]); /* up left = rmphip */

        phimzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzm[0]], arrCoorda[ez][ephi][er][icrmphimzm[1]], arrCoorda[ez][ephi][er][icrmphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]]); /* back down = phimzm */

        phimzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphimzp[0]], arrCoorda[ez][ephi][er][icrmphimzp[1]], arrCoorda[ez][ephi][er][icrmphimzp[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* front down = phimzp */

        rpphimedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphimzm[0]], arrCoorda[ez][ephi][er][icrpphimzm[1]], arrCoorda[ez][ephi][er][icrpphimzm[2]], arrCoorda[ez][ephi][er][icrpphimzp[0]], arrCoorda[ez][ephi][er][icrpphimzp[1]], arrCoorda[ez][ephi][er][icrpphimzp[2]]); /* down right = rpphim */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzmedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzm[0]];
        } else {
          rpzmedgelength = arrCoorda[ez][ephi][er][icrpphimzm[0]] * (arrCoorda[ez][ephi][er][icrpphipzm[1]] - arrCoorda[ez][ephi][er][icrpphimzm[1]]); /* back right = rpzm */
        }

        phipzmedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzm[0]], arrCoorda[ez][ephi][er][icrmphipzm[1]], arrCoorda[ez][ephi][er][icrmphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]]); /* back up = phipzm */

        if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
          rpzpedgelength = user -> dphi * arrCoorda[ez][ephi][er][icrpphimzp[0]];
        } else {
          rpzpedgelength = arrCoorda[ez][ephi][er][icrpphimzp[0]] * (arrCoorda[ez][ephi][er][icrpphipzp[1]] - arrCoorda[ez][ephi][er][icrpphimzp[1]]);
        }

        rpphipedgelength = cyldistance(arrCoorda[ez][ephi][er][icrpphipzm[0]], arrCoorda[ez][ephi][er][icrpphipzm[1]], arrCoorda[ez][ephi][er][icrpphipzm[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* up right = rpphip */

        phipzpedgelength = cyldistance(arrCoorda[ez][ephi][er][icrmphipzp[0]], arrCoorda[ez][ephi][er][icrmphipzp[1]], arrCoorda[ez][ephi][er][icrmphipzp[2]], arrCoorda[ez][ephi][er][icrpphipzp[0]], arrCoorda[ez][ephi][er][icrpphipzp[1]], arrCoorda[ez][ephi][er][icrpphipzp[2]]); /* front up = phipzp */

        // F1: r component derivatives
        arrF1[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivVrmphipzm[0]] - arrX[ez][ephi][er][ivVrmphimzm[0]]) / rmzmedgelength;
        arrF1[ez][ephi][er][ivErmphim] = (arrX[ez][ephi][er][ivVrmphimzp[0]] - arrX[ez][ephi][er][ivVrmphimzm[0]]) / rmphimedgelength;
        arrF1[ez][ephi][er][ivEphimzm] = (arrX[ez][ephi][er][ivVrpphimzm[0]] - arrX[ez][ephi][er][ivVrmphimzm[0]]) / phimzmedgelength;
        if (er == N[0] - 1) {
          arrF1[ez][ephi][er][ivErpzm] = (arrX[ez][ephi][er][ivVrpphipzm[0]] - arrX[ez][ephi][er][ivVrpphimzm[0]]) / rpzmedgelength;
          arrF1[ez][ephi][er][ivErpphim] = (arrX[ez][ephi][er][ivVrpphimzp[0]] - arrX[ez][ephi][er][ivVrpphimzm[0]]) / rpphimedgelength;
        }
        if (ez == N[2] - 1) {
          arrF1[ez][ephi][er][ivErmzp] = (arrX[ez][ephi][er][ivVrmphipzp[0]] - arrX[ez][ephi][er][ivVrmphimzp[0]]) / rmzpedgelength;
          arrF1[ez][ephi][er][ivEphimzp] = (arrX[ez][ephi][er][ivVrpphimzp[0]] - arrX[ez][ephi][er][ivVrmphimzp[0]]) / phimzpedgelength;
        }
        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF1[ez][ephi][er][ivEphipzm] = (arrX[ez][ephi][er][ivVrpphipzm[0]] - arrX[ez][ephi][er][ivVrmphipzm[0]]) / phipzmedgelength;
          arrF1[ez][ephi][er][ivErmphip] = (arrX[ez][ephi][er][ivVrmphipzp[0]] - arrX[ez][ephi][er][ivVrmphipzm[0]]) / rmphipedgelength;
        }
        if (!(user -> phibtype)) {
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrF1[ez][ephi][er][ivErpphip] = (arrX[ez][ephi][er][ivVrpphipzp[0]] - arrX[ez][ephi][er][ivVrpphipzm[0]]) / rpphipedgelength;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrF1[ez][ephi][er][ivEphipzp] = (arrX[ez][ephi][er][ivVrpphipzp[0]] - arrX[ez][ephi][er][ivVrmphipzp[0]]) / phipzpedgelength;
          }
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrF1[ez][ephi][er][ivErpzp] = (arrX[ez][ephi][er][ivVrpphipzp[0]] - arrX[ez][ephi][er][ivVrpphimzp[0]]) / rpzpedgelength;
        }

        // F2: phi component derivatives
        arrF2[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivVrmphipzm[1]] - arrX[ez][ephi][er][ivVrmphimzm[1]]) / rmzmedgelength;
        arrF2[ez][ephi][er][ivErmphim] = (arrX[ez][ephi][er][ivVrmphimzp[1]] - arrX[ez][ephi][er][ivVrmphimzm[1]]) / rmphimedgelength;
        arrF2[ez][ephi][er][ivEphimzm] = (arrX[ez][ephi][er][ivVrpphimzm[1]] - arrX[ez][ephi][er][ivVrmphimzm[1]]) / phimzmedgelength;
        if (er == N[0] - 1) {
          arrF2[ez][ephi][er][ivErpzm] = (arrX[ez][ephi][er][ivVrpphipzm[1]] - arrX[ez][ephi][er][ivVrpphimzm[1]]) / rpzmedgelength;
          arrF2[ez][ephi][er][ivErpphim] = (arrX[ez][ephi][er][ivVrpphimzp[1]] - arrX[ez][ephi][er][ivVrpphimzm[1]]) / rpphimedgelength;
        }
        if (ez == N[2] - 1) {
          arrF2[ez][ephi][er][ivErmzp] = (arrX[ez][ephi][er][ivVrmphipzp[1]] - arrX[ez][ephi][er][ivVrmphimzp[1]]) / rmzpedgelength;
          arrF2[ez][ephi][er][ivEphimzp] = (arrX[ez][ephi][er][ivVrpphimzp[1]] - arrX[ez][ephi][er][ivVrmphimzp[1]]) / phimzpedgelength;
        }
        if (ephi == N[1] - 1 && !(user -> phibtype)) {
          arrF2[ez][ephi][er][ivEphipzm] = (arrX[ez][ephi][er][ivVrpphipzm[1]] - arrX[ez][ephi][er][ivVrmphipzm[1]]) / phipzmedgelength;
          arrF2[ez][ephi][er][ivErmphip] = (arrX[ez][ephi][er][ivVrmphipzp[1]] - arrX[ez][ephi][er][ivVrmphipzm[1]]) / rmphipedgelength;
        }
        if (!(user -> phibtype)) {
          if (er == N[0] - 1 && ephi == N[1] - 1) {
            arrF2[ez][ephi][er][ivErpphip] = (arrX[ez][ephi][er][ivVrpphipzp[1]] - arrX[ez][ephi][er][ivVrpphipzm[1]]) / rpphipedgelength;
          }
          if (ephi == N[1] - 1 && ez == N[2] - 1) {
            arrF2[ez][ephi][er][ivEphipzp] = (arrX[ez][ephi][er][ivVrpphipzp[1]] - arrX[ez][ephi][er][ivVrmphipzp[1]]) / phipzpedgelength;
          }
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrF2[ez][ephi][er][ivErpzp] = (arrX[ez][ephi][er][ivVrpphipzp[1]] - arrX[ez][ephi][er][ivVrpphimzp[1]]) / rpzpedgelength;
        }

        // F3: z component derivatives
          arrF3[ez][ephi][er][ivErmzm] = (arrX[ez][ephi][er][ivVrmphipzm[2]] - arrX[ez][ephi][er][ivVrmphimzm[2]]) / rmzmedgelength;
          arrF3[ez][ephi][er][ivErmphim] = (arrX[ez][ephi][er][ivVrmphimzp[2]] - arrX[ez][ephi][er][ivVrmphimzm[2]]) / rmphimedgelength;
          arrF3[ez][ephi][er][ivEphimzm] = (arrX[ez][ephi][er][ivVrpphimzm[2]] - arrX[ez][ephi][er][ivVrmphimzm[2]]) / phimzmedgelength;
          if (er == N[0] - 1) {
            arrF3[ez][ephi][er][ivErpzm] = (arrX[ez][ephi][er][ivVrpphipzm[2]] - arrX[ez][ephi][er][ivVrpphimzm[2]]) / rpzmedgelength;
            arrF3[ez][ephi][er][ivErpphim] = (arrX[ez][ephi][er][ivVrpphimzp[2]] - arrX[ez][ephi][er][ivVrpphimzm[2]]) / rpphimedgelength;
          }
          if (ez == N[2] - 1) {
            arrF3[ez][ephi][er][ivErmzp] = (arrX[ez][ephi][er][ivVrmphipzp[2]] - arrX[ez][ephi][er][ivVrmphimzp[2]]) / rmzpedgelength;
            arrF3[ez][ephi][er][ivEphimzp] = (arrX[ez][ephi][er][ivVrpphimzp[2]] - arrX[ez][ephi][er][ivVrmphimzp[2]]) / phimzpedgelength;
          }
          if (ephi == N[1] - 1 && !(user -> phibtype)) {
            arrF3[ez][ephi][er][ivEphipzm] = (arrX[ez][ephi][er][ivVrpphipzm[2]] - arrX[ez][ephi][er][ivVrmphipzm[2]]) / phipzmedgelength;
            arrF3[ez][ephi][er][ivErmphip] = (arrX[ez][ephi][er][ivVrmphipzp[2]] - arrX[ez][ephi][er][ivVrmphipzm[2]]) / rmphipedgelength;
          }
          if (!(user -> phibtype)) {
            if (er == N[0] - 1 && ephi == N[1] - 1) {
              arrF3[ez][ephi][er][ivErpphip] = (arrX[ez][ephi][er][ivVrpphipzp[2]] - arrX[ez][ephi][er][ivVrpphipzm[2]]) / rpphipedgelength;
            }
            if (ephi == N[1] - 1 && ez == N[2] - 1) {
              arrF3[ez][ephi][er][ivEphipzp] = (arrX[ez][ephi][er][ivVrpphipzp[2]] - arrX[ez][ephi][er][ivVrmphipzp[2]]) / phipzpedgelength;
            }
          }
          if (er == N[0] - 1 && ez == N[2] - 1) {
            arrF3[ez][ephi][er][ivErpzp] = (arrX[ez][ephi][er][ivVrpphipzp[2]] - arrX[ez][ephi][er][ivVrpphimzp[2]]) / rpzpedgelength;
          }

      }
    }
  }

  /* Restore vectors */
  DMStagVecRestoreArray(da, f1Local, & arrF1);
  DMLocalToGlobal(da, f1Local, INSERT_VALUES, F1);
  DMRestoreLocalVector(da, & f1Local);

  DMStagVecRestoreArray(da, f2Local, & arrF2);
  DMLocalToGlobal(da, f2Local, INSERT_VALUES, F2);
  DMRestoreLocalVector(da, & f2Local);

  DMStagVecRestoreArray(da, f3Local, & arrF3);
  DMLocalToGlobal(da, f3Local, INSERT_VALUES, F3);
  DMRestoreLocalVector(da, & f3Local);

  DMStagVecRestoreArrayRead(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, & arrCoorda);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}

PetscErrorCode FormElectricField(TS ts, Vec X, Vec F, void * ptr) {

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;
  PetscLogDouble user_event_flops;

  PetscClassIdRegister("class name",&classid);
  PetscLogEventRegister("FormElectricField",classid,&USER_EVENT);
  PetscLogEventBegin(USER_EVENT,0,0,0,0);

  User * user = (User * ) ptr;

  FormDiscreteGradientEP_noMat(ts, X, F, user);
  VecScale(F,-1.0);
  VecAXPY(F,1.0,X);

  PetscLogFlops(user_event_flops);
  PetscLogEventEnd(USER_EVENT,0,0,0,0);

  return (0);
}
