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
// objects (verified via nm) as of this commit. They were moved out of
// monitor_functions.c verbatim (byte-identical bodies) to shrink that file
// while keeping the code in-tree and greppable for physics review.
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

PetscErrorCode DumpSolution(TS ts, PetscInt step, Vec X, void * ptr) {
  User * user = (User * ) ptr;
  DM da, dmC, daC, dmBAvg, daBAvg, dmEAvg, daEAvg, dmJAvg, daJAvg, dmV, daV, dmEP, daEP;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec J, J_local, X_local, vecC, C, vecBAvg, BAvg, vecEAvg, EAvg, vecJAvg, JAvg, vecV, V, vecEP, EP;
  PetscReal time = 0.0;

  TSGetDM(ts, & da);
  TSGetTime(ts, & time);

  //PetscPrintf(PETSC_COMM_WORLD,"Current time: t = %f\n", time);

  DMStagCreateCompatibleDMStag(da, 1, 0, 0, 0, & dmEP); /* 1 dof per vertex */
  DMStagCreateCompatibleDMStag(da, 3, 0, 0, 0, & dmV); /* 3 dofs per vertex */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmBAvg); /* 3 dof per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmEAvg); /* 3 dof per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmJAvg); /* 3 dof per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 1, & dmC); /* 1 dof per element */

  DMSetUp(dmEP);
  DMSetUp(dmV);
  DMSetUp(dmBAvg);
  DMSetUp(dmEAvg);
  DMSetUp(dmJAvg);
  DMSetUp(dmC);

  DMStagSetUniformCoordinatesExplicit(dmEP, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
  DMStagSetUniformCoordinatesExplicit(dmV, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
  DMStagSetUniformCoordinatesExplicit(dmBAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
  DMStagSetUniformCoordinatesExplicit(dmEAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
  DMStagSetUniformCoordinatesExplicit(dmJAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
  DMStagSetUniformCoordinatesExplicit(dmC, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);

  DMCreateGlobalVector(dmEP, & EP);
  DMCreateGlobalVector(dmV, & V);
  DMCreateGlobalVector(dmBAvg, & BAvg);
  DMCreateGlobalVector(dmEAvg, & EAvg);
  DMCreateGlobalVector(dmJAvg, & JAvg);
  DMCreateGlobalVector(dmC, & C);

  DMGetLocalVector(da, & X_local);
  DMGlobalToLocal(da, X, INSERT_VALUES, X_local);

  // Compute \tilde{J}:= curl(\tilde{B})
  VecDuplicate(X,&J);
  VecCopy(X, J);
  FormDerivedCurlnomp(ts, X, J, user);
  DMGetLocalVector(da, & J_local);
  DMGlobalToLocal(da, J, INSERT_VALUES, J_local);

  DMStagGetCorners(dmEP, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[8], to[8];
        PetscScalar valFrom[8], valTo[8];

        from[0].i = er;
        from[0].j = ephi;
        from[0].k = ez;
        from[0].loc = BACK_DOWN_LEFT;
        from[0].c = 3;
        from[1].i = er;
        from[1].j = ephi;
        from[1].k = ez;
        from[1].loc = BACK_DOWN_RIGHT;
        from[1].c = 3;
        from[2].i = er;
        from[2].j = ephi;
        from[2].k = ez;
        from[2].loc = BACK_UP_LEFT;
        from[2].c = 3;
        from[3].i = er;
        from[3].j = ephi;
        from[3].k = ez;
        from[3].loc = BACK_UP_RIGHT;
        from[3].c = 3;
        from[4].i = er;
        from[4].j = ephi;
        from[4].k = ez;
        from[4].loc = FRONT_DOWN_LEFT;
        from[4].c = 3;
        from[5].i = er;
        from[5].j = ephi;
        from[5].k = ez;
        from[5].loc = FRONT_DOWN_RIGHT;
        from[5].c = 3;
        from[6].i = er;
        from[6].j = ephi;
        from[6].k = ez;
        from[6].loc = FRONT_UP_LEFT;
        from[6].c = 3;
        from[7].i = er;
        from[7].j = ephi;
        from[7].k = ez;
        from[7].loc = FRONT_UP_RIGHT;
        from[7].c = 3;
        DMStagVecGetValuesStencil(da, X_local, 8, from, valFrom);

        to[0].i = er;
        to[0].j = ephi;
        to[0].k = ez;
        to[0].loc = BACK_DOWN_LEFT;
        to[0].c = 0;
        valTo[0] = valFrom[0];
        to[1].i = er;
        to[1].j = ephi;
        to[1].k = ez;
        to[1].loc = BACK_DOWN_RIGHT;
        to[1].c = 0;
        valTo[1] = valFrom[1];
        to[2].i = er;
        to[2].j = ephi;
        to[2].k = ez;
        to[2].loc = BACK_UP_LEFT;
        to[2].c = 0;
        valTo[2] = valFrom[2];
        to[3].i = er;
        to[3].j = ephi;
        to[3].k = ez;
        to[3].loc = BACK_UP_RIGHT;
        to[3].c = 0;
        valTo[3] = valFrom[3];
        to[4].i = er;
        to[4].j = ephi;
        to[4].k = ez;
        to[4].loc = FRONT_DOWN_LEFT;
        to[4].c = 0;
        valTo[4] = valFrom[4];
        to[5].i = er;
        to[5].j = ephi;
        to[5].k = ez;
        to[5].loc = FRONT_DOWN_RIGHT;
        to[5].c = 0;
        valTo[5] = valFrom[5];
        to[6].i = er;
        to[6].j = ephi;
        to[6].k = ez;
        to[6].loc = FRONT_UP_LEFT;
        to[6].c = 0;
        valTo[6] = valFrom[6];
        to[7].i = er;
        to[7].j = ephi;
        to[7].k = ez;
        to[7].loc = FRONT_UP_RIGHT;
        to[7].c = 0;
        valTo[7] = valFrom[7];

        DMStagVecSetValuesStencil(dmEP, EP, 8, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(EP);
  VecAssemblyEnd(EP);

  DMStagGetCorners(dmV, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[24], to[24];
        PetscScalar valFrom[24], valTo[24];

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = BACK_DOWN_LEFT;
          from[0].c = 0;
          from[1].i = er;
          from[1].j = ephi;
          from[1].k = ez;
          from[1].loc = BACK_DOWN_RIGHT;
          from[1].c = 0;
          from[2].i = er;
          from[2].j = ephi;
          from[2].k = ez;
          from[2].loc = BACK_UP_LEFT;
          from[2].c = 0;
          from[3].i = er;
          from[3].j = ephi;
          from[3].k = ez;
          from[3].loc = BACK_UP_RIGHT;
          from[3].c = 0;
          from[4].i = er;
          from[4].j = ephi;
          from[4].k = ez;
          from[4].loc = FRONT_DOWN_LEFT;
          from[4].c = 0;
          from[5].i = er;
          from[5].j = ephi;
          from[5].k = ez;
          from[5].loc = FRONT_DOWN_RIGHT;
          from[5].c = 0;
          from[6].i = er;
          from[6].j = ephi;
          from[6].k = ez;
          from[6].loc = FRONT_UP_LEFT;
          from[6].c = 0;
          from[7].i = er;
          from[7].j = ephi;
          from[7].k = ez;
          from[7].loc = FRONT_UP_RIGHT;
          from[7].c = 0;
          from[8].i = er;
          from[8].j = ephi;
          from[8].k = ez;
          from[8].loc = BACK_DOWN_LEFT;
          from[8].c = 1;
          from[9].i = er;
          from[9].j = ephi;
          from[9].k = ez;
          from[9].loc = BACK_DOWN_RIGHT;
          from[9].c = 1;
          from[10].i = er;
          from[10].j = ephi;
          from[10].k = ez;
          from[10].loc = BACK_UP_LEFT;
          from[10].c = 1;
          from[11].i = er;
          from[11].j = ephi;
          from[11].k = ez;
          from[11].loc = BACK_UP_RIGHT;
          from[11].c = 1;
          from[12].i = er;
          from[12].j = ephi;
          from[12].k = ez;
          from[12].loc = FRONT_DOWN_LEFT;
          from[12].c = 1;
          from[13].i = er;
          from[13].j = ephi;
          from[13].k = ez;
          from[13].loc = FRONT_DOWN_RIGHT;
          from[13].c = 1;
          from[14].i = er;
          from[14].j = ephi;
          from[14].k = ez;
          from[14].loc = FRONT_UP_LEFT;
          from[14].c = 1;
          from[15].i = er;
          from[15].j = ephi;
          from[15].k = ez;
          from[15].loc = FRONT_UP_RIGHT;
          from[15].c = 1;
          from[16].i = er;
          from[16].j = ephi;
          from[16].k = ez;
          from[16].loc = BACK_DOWN_LEFT;
          from[16].c = 2;
          from[17].i = er;
          from[17].j = ephi;
          from[17].k = ez;
          from[17].loc = BACK_DOWN_RIGHT;
          from[17].c = 2;
          from[18].i = er;
          from[18].j = ephi;
          from[18].k = ez;
          from[18].loc = BACK_UP_LEFT;
          from[18].c = 2;
          from[19].i = er;
          from[19].j = ephi;
          from[19].k = ez;
          from[19].loc = BACK_UP_RIGHT;
          from[19].c = 2;
          from[20].i = er;
          from[20].j = ephi;
          from[20].k = ez;
          from[20].loc = FRONT_DOWN_LEFT;
          from[20].c = 2;
          from[21].i = er;
          from[21].j = ephi;
          from[21].k = ez;
          from[21].loc = FRONT_DOWN_RIGHT;
          from[21].c = 2;
          from[22].i = er;
          from[22].j = ephi;
          from[22].k = ez;
          from[22].loc = FRONT_UP_LEFT;
          from[22].c = 2;
          from[23].i = er;
          from[23].j = ephi;
          from[23].k = ez;
          from[23].loc = FRONT_UP_RIGHT;
          from[23].c = 2;
          DMStagVecGetValuesStencil(da, X_local, 24, from, valFrom);

        to[0].i = er;
        to[0].j = ephi;
        to[0].k = ez;
        to[0].loc = BACK_DOWN_LEFT;
        to[0].c = 0;
        valTo[0] = valFrom[0];
        to[1].i = er;
        to[1].j = ephi;
        to[1].k = ez;
        to[1].loc = BACK_DOWN_RIGHT;
        to[1].c = 0;
        valTo[1] = valFrom[1];
        to[2].i = er;
        to[2].j = ephi;
        to[2].k = ez;
        to[2].loc = BACK_UP_LEFT;
        to[2].c = 0;
        valTo[2] = valFrom[2];
        to[3].i = er;
        to[3].j = ephi;
        to[3].k = ez;
        to[3].loc = BACK_UP_RIGHT;
        to[3].c = 0;
        valTo[3] = valFrom[3];
        to[4].i = er;
        to[4].j = ephi;
        to[4].k = ez;
        to[4].loc = FRONT_DOWN_LEFT;
        to[4].c = 0;
        valTo[4] = valFrom[4];
        to[5].i = er;
        to[5].j = ephi;
        to[5].k = ez;
        to[5].loc = FRONT_DOWN_RIGHT;
        to[5].c = 0;
        valTo[5] = valFrom[5];
        to[6].i = er;
        to[6].j = ephi;
        to[6].k = ez;
        to[6].loc = FRONT_UP_LEFT;
        to[6].c = 0;
        valTo[6] = valFrom[6];
        to[7].i = er;
        to[7].j = ephi;
        to[7].k = ez;
        to[7].loc = FRONT_UP_RIGHT;
        to[7].c = 0;
        valTo[7] = valFrom[7];

        to[8].i = er;
        to[8].j = ephi;
        to[8].k = ez;
        to[8].loc = BACK_DOWN_LEFT;
        to[8].c = 1;
        valTo[8] = valFrom[8];
        to[9].i = er;
        to[9].j = ephi;
        to[9].k = ez;
        to[9].loc = BACK_DOWN_RIGHT;
        to[9].c = 1;
        valTo[9] = valFrom[9];
        to[10].i = er;
        to[10].j = ephi;
        to[10].k = ez;
        to[10].loc = BACK_UP_LEFT;
        to[10].c = 1;
        valTo[10] = valFrom[10];
        to[11].i = er;
        to[11].j = ephi;
        to[11].k = ez;
        to[11].loc = BACK_UP_RIGHT;
        to[11].c = 1;
        valTo[11] = valFrom[11];
        to[12].i = er;
        to[12].j = ephi;
        to[12].k = ez;
        to[12].loc = FRONT_DOWN_LEFT;
        to[12].c = 1;
        valTo[12] = valFrom[12];
        to[13].i = er;
        to[13].j = ephi;
        to[13].k = ez;
        to[13].loc = FRONT_DOWN_RIGHT;
        to[13].c = 1;
        valTo[13] = valFrom[13];
        to[14].i = er;
        to[14].j = ephi;
        to[14].k = ez;
        to[14].loc = FRONT_UP_LEFT;
        to[14].c = 1;
        valTo[14] = valFrom[14];
        to[15].i = er;
        to[15].j = ephi;
        to[15].k = ez;
        to[15].loc = FRONT_UP_RIGHT;
        to[15].c = 1;
        valTo[15] = valFrom[15];

        to[16].i = er;
        to[16].j = ephi;
        to[16].k = ez;
        to[16].loc = BACK_DOWN_LEFT;
        to[16].c = 2;
        valTo[16] = valFrom[16];
        to[17].i = er;
        to[17].j = ephi;
        to[17].k = ez;
        to[17].loc = BACK_DOWN_RIGHT;
        to[17].c = 2;
        valTo[17] = valFrom[17];
        to[18].i = er;
        to[18].j = ephi;
        to[18].k = ez;
        to[18].loc = BACK_UP_LEFT;
        to[18].c = 2;
        valTo[18] = valFrom[18];
        to[19].i = er;
        to[19].j = ephi;
        to[19].k = ez;
        to[19].loc = BACK_UP_RIGHT;
        to[19].c = 2;
        valTo[19] = valFrom[19];
        to[20].i = er;
        to[20].j = ephi;
        to[20].k = ez;
        to[20].loc = FRONT_DOWN_LEFT;
        to[20].c = 2;
        valTo[20] = valFrom[20];
        to[21].i = er;
        to[21].j = ephi;
        to[21].k = ez;
        to[21].loc = FRONT_DOWN_RIGHT;
        to[21].c = 2;
        valTo[21] = valFrom[21];
        to[22].i = er;
        to[22].j = ephi;
        to[22].k = ez;
        to[22].loc = FRONT_UP_LEFT;
        to[22].c = 2;
        valTo[22] = valFrom[22];
        to[23].i = er;
        to[23].j = ephi;
        to[23].k = ez;
        to[23].loc = FRONT_UP_RIGHT;
        to[23].c = 2;
        valTo[23] = valFrom[23];

        DMStagVecSetValuesStencil(dmV, V, 24, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);

  DMStagGetCorners(dmBAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[6], to[3];
        PetscScalar valFrom[6], valTo[3];

        from[0].i = er;
        from[0].j = ephi;
        from[0].k = ez;
        from[0].loc = UP;
        from[0].c = 0;
        from[1].i = er;
        from[1].j = ephi;
        from[1].k = ez;
        from[1].loc = DOWN;
        from[1].c = 0;
        from[2].i = er;
        from[2].j = ephi;
        from[2].k = ez;
        from[2].loc = LEFT;
        from[2].c = 0;
        from[3].i = er;
        from[3].j = ephi;
        from[3].k = ez;
        from[3].loc = RIGHT;
        from[3].c = 0;
        from[4].i = er;
        from[4].j = ephi;
        from[4].k = ez;
        from[4].loc = FRONT;
        from[4].c = 0;
        from[5].i = er;
        from[5].j = ephi;
        from[5].k = ez;
        from[5].loc = BACK;
        from[5].c = 0;

        DMStagVecGetValuesStencil(da, X_local, 6, from, valFrom);

        to[0].i = er;
        to[0].j = ephi;
        to[0].k = ez;
        to[0].loc = ELEMENT;
        to[0].c = 0;
        valTo[0] = 0.5 * (valFrom[2] + valFrom[3]);
        to[1].i = er;
        to[1].j = ephi;
        to[1].k = ez;
        to[1].loc = ELEMENT;
        to[1].c = 1;
        valTo[1] = 0.5 * (valFrom[0] + valFrom[1]);
        to[2].i = er;
        to[2].j = ephi;
        to[2].k = ez;
        to[2].loc = ELEMENT;
        to[2].c = 2;
        valTo[2] = 0.5 * (valFrom[4] + valFrom[5]);

        DMStagVecSetValuesStencil(dmBAvg, BAvg, 3, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(BAvg);
  VecAssemblyEnd(BAvg);

  DMStagGetCorners(dmEAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[12], to[3];
        PetscScalar valFrom[12], valTo[3];

        from[0].i = er;
        from[0].j = ephi;
        from[0].k = ez;
        from[0].loc = FRONT_UP;
        from[0].c = 0;
        from[1].i = er;
        from[1].j = ephi;
        from[1].k = ez;
        from[1].loc = BACK_UP;
        from[1].c = 0;
        from[2].i = er;
        from[2].j = ephi;
        from[2].k = ez;
        from[2].loc = FRONT_DOWN;
        from[2].c = 0;
        from[3].i = er;
        from[3].j = ephi;
        from[3].k = ez;
        from[3].loc = BACK_DOWN;
        from[3].c = 0;
        from[4].i = er;
        from[4].j = ephi;
        from[4].k = ez;
        from[4].loc = FRONT_RIGHT;
        from[4].c = 0;
        from[5].i = er;
        from[5].j = ephi;
        from[5].k = ez;
        from[5].loc = BACK_RIGHT;
        from[5].c = 0;
        from[6].i = er;
        from[6].j = ephi;
        from[6].k = ez;
        from[6].loc = FRONT_LEFT;
        from[6].c = 0;
        from[7].i = er;
        from[7].j = ephi;
        from[7].k = ez;
        from[7].loc = BACK_LEFT;
        from[7].c = 0;
        from[8].i = er;
        from[8].j = ephi;
        from[8].k = ez;
        from[8].loc = UP_RIGHT;
        from[8].c = 0;
        from[9].i = er;
        from[9].j = ephi;
        from[9].k = ez;
        from[9].loc = DOWN_RIGHT;
        from[9].c = 0;
        from[10].i = er;
        from[10].j = ephi;
        from[10].k = ez;
        from[10].loc = UP_LEFT;
        from[10].c = 0;
        from[11].i = er;
        from[11].j = ephi;
        from[11].k = ez;
        from[11].loc = DOWN_LEFT;
        from[11].c = 0;
        DMStagVecGetValuesStencil(da, X_local, 12, from, valFrom);
        to[0].i = er;
        to[0].j = ephi;
        to[0].k = ez;
        to[0].loc = ELEMENT;
        to[0].c = 0;
        valTo[0] = 0.25 * (valFrom[0] + valFrom[1] + valFrom[2] + valFrom[3]);
        to[1].i = er;
        to[1].j = ephi;
        to[1].k = ez;
        to[1].loc = ELEMENT;
        to[1].c = 1;
        valTo[1] = 0.25 * (valFrom[4] + valFrom[5] + valFrom[6] + valFrom[7]);
        to[2].i = er;
        to[2].j = ephi;
        to[2].k = ez;
        to[2].loc = ELEMENT;
        to[2].c = 2;
        valTo[2] = 0.25 * (valFrom[8] + valFrom[9] + valFrom[10] + valFrom[11]);
        DMStagVecSetValuesStencil(dmEAvg, EAvg, 3, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(EAvg);
  VecAssemblyEnd(EAvg);

  DMStagGetCorners(dmJAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[12], to[3];
        PetscScalar valFrom[12], valTo[3];

        from[0].i = er;
        from[0].j = ephi;
        from[0].k = ez;
        from[0].loc = FRONT_UP;
        from[0].c = 0;
        from[1].i = er;
        from[1].j = ephi;
        from[1].k = ez;
        from[1].loc = BACK_UP;
        from[1].c = 0;
        from[2].i = er;
        from[2].j = ephi;
        from[2].k = ez;
        from[2].loc = FRONT_DOWN;
        from[2].c = 0;
        from[3].i = er;
        from[3].j = ephi;
        from[3].k = ez;
        from[3].loc = BACK_DOWN;
        from[3].c = 0;
        from[4].i = er;
        from[4].j = ephi;
        from[4].k = ez;
        from[4].loc = FRONT_RIGHT;
        from[4].c = 0;
        from[5].i = er;
        from[5].j = ephi;
        from[5].k = ez;
        from[5].loc = BACK_RIGHT;
        from[5].c = 0;
        from[6].i = er;
        from[6].j = ephi;
        from[6].k = ez;
        from[6].loc = FRONT_LEFT;
        from[6].c = 0;
        from[7].i = er;
        from[7].j = ephi;
        from[7].k = ez;
        from[7].loc = BACK_LEFT;
        from[7].c = 0;
        from[8].i = er;
        from[8].j = ephi;
        from[8].k = ez;
        from[8].loc = UP_RIGHT;
        from[8].c = 0;
        from[9].i = er;
        from[9].j = ephi;
        from[9].k = ez;
        from[9].loc = DOWN_RIGHT;
        from[9].c = 0;
        from[10].i = er;
        from[10].j = ephi;
        from[10].k = ez;
        from[10].loc = UP_LEFT;
        from[10].c = 0;
        from[11].i = er;
        from[11].j = ephi;
        from[11].k = ez;
        from[11].loc = DOWN_LEFT;
        from[11].c = 0;
        DMStagVecGetValuesStencil(da, J_local, 12, from, valFrom);
        to[0].i = er;
        to[0].j = ephi;
        to[0].k = ez;
        to[0].loc = ELEMENT;
        to[0].c = 0;
        valTo[0] = 0.25 * (valFrom[0] + valFrom[1] + valFrom[2] + valFrom[3]);
        to[1].i = er;
        to[1].j = ephi;
        to[1].k = ez;
        to[1].loc = ELEMENT;
        to[1].c = 1;
        valTo[1] = 0.25 * (valFrom[4] + valFrom[5] + valFrom[6] + valFrom[7]);
        to[2].i = er;
        to[2].j = ephi;
        to[2].k = ez;
        to[2].loc = ELEMENT;
        to[2].c = 2;
        valTo[2] = 0.25 * (valFrom[8] + valFrom[9] + valFrom[10] + valFrom[11]);
        DMStagVecSetValuesStencil(dmJAvg, JAvg, 3, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(JAvg);
  VecAssemblyEnd(JAvg);

  DMStagGetCorners(dmC, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[1], to[1];
        PetscScalar valFrom[1], valTo[1];

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = ELEMENT;
          from[0].c = 0;
          DMStagVecGetValuesStencil(da, X_local, 1, from, valFrom);
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = ELEMENT;
          to[0].c = 0;
          valTo[0] = valFrom[0];
          DMStagVecSetValuesStencil(dmC, C, 1, to, valTo, INSERT_VALUES);
        }
      }
    }
    VecAssemblyBegin(C);
    VecAssemblyEnd(C);

  DMRestoreLocalVector(da, & X_local);
  DMRestoreLocalVector(da, & J_local);

  DMStagVecSplitToDMDA(dmEP, EP, BACK_DOWN_LEFT, -1, & daEP, & vecEP); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmV, V, BACK_DOWN_LEFT, -3, & daV, & vecV); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmBAvg, BAvg, ELEMENT, -3, & daBAvg, & vecBAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmEAvg, EAvg, ELEMENT, -3, & daEAvg, & vecEAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmJAvg, JAvg, ELEMENT, -3, & daJAvg, & vecJAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmC, C, ELEMENT, -1, & daC, & vecC); /* note -3 : pad with zero in 2D case */

  PetscObjectSetName((PetscObject) vecEP, "Electrostatic Potential");
  PetscObjectSetName((PetscObject) vecV, "Velocity");
  PetscObjectSetName((PetscObject) vecBAvg, "Magnetic Field (Averaged)");
  PetscObjectSetName((PetscObject) vecEAvg, "Divergence-free part of Electric Field (Averaged)");
  PetscObjectSetName((PetscObject) vecJAvg, "Curl of Magnetic Field (Averaged)");
  PetscObjectSetName((PetscObject) vecC, "Number Density of Ions");

  /* Dump element-based fields to a .vtr file and create a .pvd file */
  {
    PetscViewer viewerC, viewerB, viewerE, viewerJ, viewerV, viewerEP;
    char filename[PETSC_MAX_PATH_LEN];
    FILE * pvdfile;

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D.pvd", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz);
    if (time == 0.0) {
      PetscFOpen(PETSC_COMM_WORLD, filename, "a", & pvdfile);
      //pvdfile = fopen(filename, "a");
      //pvdfile.open (filename, ios::out | ios::app);
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "<\?xml version=\"1.0\"?>\n");
      //pvdfile << "<\?xml version=\"1.0\"?>\n";
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
      //pvdfile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "  <Collection>\n");
      //pvdfile << "  <Collection>\n";
      PetscFClose(PETSC_COMM_WORLD, pvdfile);
      //fclose(pvdfile);
      //pvdfile.close();
    }

    PetscFOpen(PETSC_COMM_WORLD, filename, "a", & pvdfile);
    PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"0\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_vavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"1\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_bavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"2\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_tavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"3\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_navg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"4\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_javg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"5\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_lvlset_ic%.1D_grid%.2Dx%.2Dx%.2D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> Nr, user -> Nphi, user -> Nz);
    PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"6\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_epavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscFClose(PETSC_COMM_WORLD, pvdfile);

    if (time >= user -> ftime) {
      PetscFOpen(PETSC_COMM_WORLD, filename, "a", & pvdfile);
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "  </Collection>\n");
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "</VTKFile>\n");
      PetscFClose(PETSC_COMM_WORLD, pvdfile);
    }

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_epavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daEP), filename, FILE_MODE_WRITE, & viewerEP);
    VecView(vecEP, viewerEP);
    PetscViewerDestroy( & viewerEP);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_vavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daV), filename, FILE_MODE_WRITE, & viewerV);
    VecView(vecV, viewerV);
    PetscViewerDestroy( & viewerV);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_bavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daBAvg), filename, FILE_MODE_WRITE, & viewerB);
    VecView(vecBAvg, viewerB);
    PetscViewerDestroy( & viewerB);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_tavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daEAvg), filename, FILE_MODE_WRITE, & viewerE);
    VecView(vecEAvg, viewerE);
    PetscViewerDestroy( & viewerE);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_javg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daJAvg), filename, FILE_MODE_WRITE, & viewerJ);
    VecView(vecJAvg, viewerJ);
    PetscViewerDestroy( & viewerJ);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_navg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daC), filename, FILE_MODE_WRITE, & viewerC);
    VecView(vecC, viewerC);
    PetscViewerDestroy( & viewerC);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);
  }

  /* Destroy DMDAs and Vecs */
  VecDestroy( & vecEP);
  DMDestroy( & daEP);
  VecDestroy( & EP);
  DMDestroy( & dmEP);

  VecDestroy( & vecV);
  DMDestroy( & daV);
  VecDestroy( & V);
  DMDestroy( & dmV);

  VecDestroy( & vecBAvg);
  DMDestroy( & daBAvg);
  VecDestroy( & BAvg);
  DMDestroy( & dmBAvg);

  VecDestroy( & vecEAvg);
  DMDestroy( & daEAvg);
  VecDestroy( & EAvg);
  DMDestroy( & dmEAvg);

  VecDestroy( & vecJAvg);
  DMDestroy( & daJAvg);
  VecDestroy( & JAvg);
  DMDestroy( & dmJAvg);

  VecDestroy( & vecC);
  DMDestroy( & daC);
  VecDestroy( & C);
  DMDestroy( & dmC);

  VecDestroy( & J);

  return (0);
}

PetscErrorCode Dump1stVertexField(TS ts, PetscInt step, Vec X, void * ptr) {
  User * user = (User * ) ptr;
  DM da, dmV, daV;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec X_local, vecV, V;
  PetscReal time = 0.0;

  TSGetDM(ts, & da);
  TSGetTime(ts, & time);

  DMStagCreateCompatibleDMStag(da, 3, 0, 0, 0, & dmV); /* 3 dofs per vertex */
  DMSetUp(dmV);
  DMStagSetUniformCoordinatesExplicit(dmV, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
  DMCreateGlobalVector(dmV, & V);
  DMGetLocalVector(da, & X_local);
  DMGlobalToLocal(da, X, INSERT_VALUES, X_local);

  DMStagGetCorners(dmV, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[24], to[24];
        PetscScalar valFrom[24], valTo[24];

        from[0].i = er;
        from[0].j = ephi;
        from[0].k = ez;
        from[0].loc = BACK_DOWN_LEFT;
        from[0].c = 0;
        from[1].i = er;
        from[1].j = ephi;
        from[1].k = ez;
        from[1].loc = BACK_DOWN_RIGHT;
        from[1].c = 0;
        from[2].i = er;
        from[2].j = ephi;
        from[2].k = ez;
        from[2].loc = BACK_UP_LEFT;
        from[2].c = 0;
        from[3].i = er;
        from[3].j = ephi;
        from[3].k = ez;
        from[3].loc = BACK_UP_RIGHT;
        from[3].c = 0;
        from[4].i = er;
        from[4].j = ephi;
        from[4].k = ez;
        from[4].loc = FRONT_DOWN_LEFT;
        from[4].c = 0;
        from[5].i = er;
        from[5].j = ephi;
        from[5].k = ez;
        from[5].loc = FRONT_DOWN_RIGHT;
        from[5].c = 0;
        from[6].i = er;
        from[6].j = ephi;
        from[6].k = ez;
        from[6].loc = FRONT_UP_LEFT;
        from[6].c = 0;
        from[7].i = er;
        from[7].j = ephi;
        from[7].k = ez;
        from[7].loc = FRONT_UP_RIGHT;
        from[7].c = 0;
        from[8].i = er;
        from[8].j = ephi;
        from[8].k = ez;
        from[8].loc = BACK_DOWN_LEFT;
        from[8].c = 1;
        from[9].i = er;
        from[9].j = ephi;
        from[9].k = ez;
        from[9].loc = BACK_DOWN_RIGHT;
        from[9].c = 1;
        from[10].i = er;
        from[10].j = ephi;
        from[10].k = ez;
        from[10].loc = BACK_UP_LEFT;
        from[10].c = 1;
        from[11].i = er;
        from[11].j = ephi;
        from[11].k = ez;
        from[11].loc = BACK_UP_RIGHT;
        from[11].c = 1;
        from[12].i = er;
        from[12].j = ephi;
        from[12].k = ez;
        from[12].loc = FRONT_DOWN_LEFT;
        from[12].c = 1;
        from[13].i = er;
        from[13].j = ephi;
        from[13].k = ez;
        from[13].loc = FRONT_DOWN_RIGHT;
        from[13].c = 1;
        from[14].i = er;
        from[14].j = ephi;
        from[14].k = ez;
        from[14].loc = FRONT_UP_LEFT;
        from[14].c = 1;
        from[15].i = er;
        from[15].j = ephi;
        from[15].k = ez;
        from[15].loc = FRONT_UP_RIGHT;
        from[15].c = 1;
        from[16].i = er;
        from[16].j = ephi;
        from[16].k = ez;
        from[16].loc = BACK_DOWN_LEFT;
        from[16].c = 2;
        from[17].i = er;
        from[17].j = ephi;
        from[17].k = ez;
        from[17].loc = BACK_DOWN_RIGHT;
        from[17].c = 2;
        from[18].i = er;
        from[18].j = ephi;
        from[18].k = ez;
        from[18].loc = BACK_UP_LEFT;
        from[18].c = 2;
        from[19].i = er;
        from[19].j = ephi;
        from[19].k = ez;
        from[19].loc = BACK_UP_RIGHT;
        from[19].c = 2;
        from[20].i = er;
        from[20].j = ephi;
        from[20].k = ez;
        from[20].loc = FRONT_DOWN_LEFT;
        from[20].c = 2;
        from[21].i = er;
        from[21].j = ephi;
        from[21].k = ez;
        from[21].loc = FRONT_DOWN_RIGHT;
        from[21].c = 2;
        from[22].i = er;
        from[22].j = ephi;
        from[22].k = ez;
        from[22].loc = FRONT_UP_LEFT;
        from[22].c = 2;
        from[23].i = er;
        from[23].j = ephi;
        from[23].k = ez;
        from[23].loc = FRONT_UP_RIGHT;
        from[23].c = 2;
        DMStagVecGetValuesStencil(da, X_local, 24, from, valFrom);

        to[0].i = er;
        to[0].j = ephi;
        to[0].k = ez;
        to[0].loc = BACK_DOWN_LEFT;
        to[0].c = 0;
        valTo[0] = valFrom[0];
        to[1].i = er;
        to[1].j = ephi;
        to[1].k = ez;
        to[1].loc = BACK_DOWN_RIGHT;
        to[1].c = 0;
        valTo[1] = valFrom[1];
        to[2].i = er;
        to[2].j = ephi;
        to[2].k = ez;
        to[2].loc = BACK_UP_LEFT;
        to[2].c = 0;
        valTo[2] = valFrom[2];
        to[3].i = er;
        to[3].j = ephi;
        to[3].k = ez;
        to[3].loc = BACK_UP_RIGHT;
        to[3].c = 0;
        valTo[3] = valFrom[3];
        to[4].i = er;
        to[4].j = ephi;
        to[4].k = ez;
        to[4].loc = FRONT_DOWN_LEFT;
        to[4].c = 0;
        valTo[4] = valFrom[4];
        to[5].i = er;
        to[5].j = ephi;
        to[5].k = ez;
        to[5].loc = FRONT_DOWN_RIGHT;
        to[5].c = 0;
        valTo[5] = valFrom[5];
        to[6].i = er;
        to[6].j = ephi;
        to[6].k = ez;
        to[6].loc = FRONT_UP_LEFT;
        to[6].c = 0;
        valTo[6] = valFrom[6];
        to[7].i = er;
        to[7].j = ephi;
        to[7].k = ez;
        to[7].loc = FRONT_UP_RIGHT;
        to[7].c = 0;
        valTo[7] = valFrom[7];

        to[8].i = er;
        to[8].j = ephi;
        to[8].k = ez;
        to[8].loc = BACK_DOWN_LEFT;
        to[8].c = 1;
        valTo[8] = valFrom[8];
        to[9].i = er;
        to[9].j = ephi;
        to[9].k = ez;
        to[9].loc = BACK_DOWN_RIGHT;
        to[9].c = 1;
        valTo[9] = valFrom[9];
        to[10].i = er;
        to[10].j = ephi;
        to[10].k = ez;
        to[10].loc = BACK_UP_LEFT;
        to[10].c = 1;
        valTo[10] = valFrom[10];
        to[11].i = er;
        to[11].j = ephi;
        to[11].k = ez;
        to[11].loc = BACK_UP_RIGHT;
        to[11].c = 1;
        valTo[11] = valFrom[11];
        to[12].i = er;
        to[12].j = ephi;
        to[12].k = ez;
        to[12].loc = FRONT_DOWN_LEFT;
        to[12].c = 1;
        valTo[12] = valFrom[12];
        to[13].i = er;
        to[13].j = ephi;
        to[13].k = ez;
        to[13].loc = FRONT_DOWN_RIGHT;
        to[13].c = 1;
        valTo[13] = valFrom[13];
        to[14].i = er;
        to[14].j = ephi;
        to[14].k = ez;
        to[14].loc = FRONT_UP_LEFT;
        to[14].c = 1;
        valTo[14] = valFrom[14];
        to[15].i = er;
        to[15].j = ephi;
        to[15].k = ez;
        to[15].loc = FRONT_UP_RIGHT;
        to[15].c = 1;
        valTo[15] = valFrom[15];

        to[16].i = er;
        to[16].j = ephi;
        to[16].k = ez;
        to[16].loc = BACK_DOWN_LEFT;
        to[16].c = 2;
        valTo[16] = valFrom[16];
        to[17].i = er;
        to[17].j = ephi;
        to[17].k = ez;
        to[17].loc = BACK_DOWN_RIGHT;
        to[17].c = 2;
        valTo[17] = valFrom[17];
        to[18].i = er;
        to[18].j = ephi;
        to[18].k = ez;
        to[18].loc = BACK_UP_LEFT;
        to[18].c = 2;
        valTo[18] = valFrom[18];
        to[19].i = er;
        to[19].j = ephi;
        to[19].k = ez;
        to[19].loc = BACK_UP_RIGHT;
        to[19].c = 2;
        valTo[19] = valFrom[19];
        to[20].i = er;
        to[20].j = ephi;
        to[20].k = ez;
        to[20].loc = FRONT_DOWN_LEFT;
        to[20].c = 2;
        valTo[20] = valFrom[20];
        to[21].i = er;
        to[21].j = ephi;
        to[21].k = ez;
        to[21].loc = FRONT_DOWN_RIGHT;
        to[21].c = 2;
        valTo[21] = valFrom[21];
        to[22].i = er;
        to[22].j = ephi;
        to[22].k = ez;
        to[22].loc = FRONT_UP_LEFT;
        to[22].c = 2;
        valTo[22] = valFrom[22];
        to[23].i = er;
        to[23].j = ephi;
        to[23].k = ez;
        to[23].loc = FRONT_UP_RIGHT;
        to[23].c = 2;
        valTo[23] = valFrom[23];

        DMStagVecSetValuesStencil(dmV, V, 24, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);

  DMRestoreLocalVector(da, & X_local);

  DMStagVecSplitToDMDA(dmV, V, BACK_DOWN_LEFT, -3, & daV, & vecV); /* note -3 : pad with zero in 2D case */

  PetscObjectSetName((PetscObject) vecV, "Velocity");

  /* Dump element-based fields to a .vtr file and create a .pvd file */
  {
    PetscViewer viewerV;
    char filename[PETSC_MAX_PATH_LEN];
    FILE * pvdfile;

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/ictype9_dVdt_curlBxB/mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D.pvd", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz);
    if (time == 0.0) {
      PetscFOpen(PETSC_COMM_WORLD, filename, "a", & pvdfile);
      //pvdfile = fopen(filename, "a");
      //pvdfile.open (filename, ios::out | ios::app);
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "<\?xml version=\"1.0\"?>\n");
      //pvdfile << "<\?xml version=\"1.0\"?>\n";
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
      //pvdfile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "  <Collection>\n");
      //pvdfile << "  <Collection>\n";
      PetscFClose(PETSC_COMM_WORLD, pvdfile);
      //fclose(pvdfile);
      //pvdfile.close();
    }

    PetscFOpen(PETSC_COMM_WORLD, filename, "a", & pvdfile);
    PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"0\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_vavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscFClose(PETSC_COMM_WORLD, pvdfile);

    if (time >= user -> ftime) {
      PetscFOpen(PETSC_COMM_WORLD, filename, "a", & pvdfile);
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "  </Collection>\n");
      PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "</VTKFile>\n");
      PetscFClose(PETSC_COMM_WORLD, pvdfile);
    }

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/ictype9_dVdt_curlBxB/mfd_vavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daV), filename, FILE_MODE_WRITE, & viewerV);
    VecView(vecV, viewerV);
    PetscViewerDestroy( & viewerV);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);
  }

  /* Destroy DMDAs and Vecs */
  VecDestroy( & vecV);
  DMDestroy( & daV);
  VecDestroy( & V);
  DMDestroy( & dmV);

  return (0);
}

PetscErrorCode DumpPsi_Cell(TS ts, PetscInt step, PetscScalar *psi, void * ptr) {
  User * user = (User * ) ptr;
  DM da, dmC, daC;
  PetscInt er, ephi, ez, startr, startphi, startz, nr = 0, nphi = 0, nz = 0;
  Vec X_local, vecC, C;
  PetscReal time = 0.0;
  PetscMPIInt rank;
  PetscInt      i, j, k;
  const PetscInt dof0 = 0,
      dof1 = 0,
      dof2 = 0,
      dof3 = 1; /* 0 dof on each vertex, edge and face center and 1 dof on each cell center (1 for the poloidal magnetic flux function) */
  const PetscInt stencilWidth = 0;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  /* Only rank == 0 has the entries of the patch, so run code only at that rank */
  //if (rank == 0) {
    TSGetDM(ts, & da);
    TSGetTime(ts, & time);

    DMStagCreateCompatibleDMStag(da, 0, 0, 0, 1, & dmC); /* 3 dofs per element */

    //if (user->phibtype) {
      //DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE, user->Nr, user->Nphi, user->Nz, 1, 1, 1, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_NONE, stencilWidth, NULL, NULL, NULL, & dmC);
    //} else {
      //DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, user->Nr, user->Nphi, user->Nz, 1, 1, 1, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_NONE, stencilWidth, NULL, NULL, NULL, & dmC);
    //}

    //DMStagGetNumRanks(dmC, & nr, & nphi, & nz);

    DMSetUp(dmC);
    DMStagSetUniformCoordinatesExplicit(dmC, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
    DMCreateGlobalVector(dmC, & C);
    VecZeroEntries(C);


    DMStagGetCorners(dmC, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

    //PetscPrintf(PETSC_COMM_SELF,"(startr,startphi,startz,nr,nphi,nz) = (%d, %d, %d, %d, %d, %d)\n", startr, startphi, startz, nr, nphi, nz);
    //return (0);

    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil to[1];
          PetscScalar valTo[1];

          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = ELEMENT;
          to[0].c = 0;
          valTo[0] = psi[ez*user->Nr*user->Nphi+ephi*user->Nr+er] ;
          DMStagVecSetValuesStencil(dmC, C, 1, to, valTo, INSERT_VALUES);
        }
      }
    }
    VecAssemblyBegin(C);
    VecAssemblyEnd(C);

    DMStagVecSplitToDMDA(dmC, C, ELEMENT, -1, & daC, & vecC); /* note -3 : pad with zero in 2D case */

    PetscObjectSetName((PetscObject) vecC, "Poloidal magnetic flux function");

  /* Dump element-based fields to a .vtr file and create a .pvd file */

    PetscViewer viewerC;
    char filename[PETSC_MAX_PATH_LEN];
    FILE * pvdfile;

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_psavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daC), filename, FILE_MODE_WRITE, & viewerC);
    VecView(vecC, viewerC);
    PetscViewerDestroy( & viewerC);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);

    /* Destroy DMDAs and Vecs */
    VecDestroy( & vecC);
    DMDestroy( & daC);
    VecDestroy( & C);
    DMDestroy( & dmC);
  //}
  return (0);
}

PetscErrorCode getHermiteDataFD(const PetscScalar *gf, const int nR, const int nPhi, const int nZ, PetscScalar *data, const int i, const int j, const int nStencilWidth, const int bIntersect)
{
  PetscScalar E4[5];
  PetscScalar D4[5];

  if (nStencilWidth == 5)
  {
    E4[0] = 1.0/12.0, E4[1] = -8.0/12.0, E4[2] =   0.0     , E4[3] =  8.0/12.0, E4[4] = -1.0/12.0;
    D4[0] =-1.0/12.0, D4[1] = 16.0/12.0, D4[2] = -30.0/12.0, D4[3] = 16.0/12.0, D4[4] = -1.0/12.0;
    for (int ii = 0; ii < nStencilWidth; ++ii)
    {
      if(bIntersect == 1)
      {
        E4[ii] *= 4.0;
        D4[ii] *= 8.0;
      }
      else
      {
        E4[ii] *= 5.0;
        D4[ii] *= 25.0/2.0;
      }
    }
    for(int ii = 0; ii < 9; ++ii)
    {
      data[ii] = 0.0;
    }
    data[0+0] = gf[j*(nPhi*nR) + i];
    for (int jj = -nStencilWidth/2; jj <= nStencilWidth/2; ++jj)
    {
      data[1+3*0]     += gf[j*(nPhi*nR) + i+jj] * E4[jj+nStencilWidth/2];
      data[2+3*0]     += gf[j*(nPhi*nR) + i+jj] * D4[jj+nStencilWidth/2];
      data[0+3*1] += gf[(j+jj)*(nPhi*nR) + i] * E4[jj+nStencilWidth/2];
      data[0+3*2] += gf[(j+jj)*(nPhi*nR) + i] * D4[jj+nStencilWidth/2];
      for (int ii = -nStencilWidth/2; ii <= nStencilWidth/2; ++ii)
      {
        data[1+3*1] += gf[(j+jj)*(nPhi*nR) + i+ii] * E4[ii+nStencilWidth/2]*E4[jj+nStencilWidth/2];
        data[2+3*1] += gf[(j+jj)*(nPhi*nR) + i+ii] * D4[ii+nStencilWidth/2]*E4[jj+nStencilWidth/2];
        data[1+3*2] += gf[(j+jj)*(nPhi*nR) + i+ii] * E4[ii+nStencilWidth/2]*D4[jj+nStencilWidth/2];
        data[2+3*2] += gf[(j+jj)*(nPhi*nR) + i+ii] * D4[ii+nStencilWidth/2]*D4[jj+nStencilWidth/2];
      }
    }
  }
  else if(nStencilWidth == 3)
  {
    E4[0] = -0.5, E4[1] =  0.0, E4[2] = 0.5;
    D4[0] =  1.0, D4[1] = -2.0, D4[2] = 1.0;
    for (int ii = 0; ii < nStencilWidth; ++ii)
    {
      if(bIntersect == 1)
      {
        E4[ii] *= 2.0;
        D4[ii] *= 2.0;
      }
      else
      {
        E4[ii] *= 3.0;
        D4[ii] *= 9.0/2.0;
      }
    }

    for(int ii = 0; ii < 9; ++ii)
    {
      data[ii] = 0.0;
    }
    data[0+0] = gf[j*(nPhi*nR) + i];
    for (int jj = -nStencilWidth/2; jj <= nStencilWidth/2; ++jj)
    {
      data[1+3*0] += gf[j*(nPhi*nR) + i+jj] * E4[jj+nStencilWidth/2];
      data[2+3*0] += gf[j*(nPhi*nR) + i+jj] * D4[jj+nStencilWidth/2];
      data[0+3*1] += gf[(j+jj)*(nPhi*nR) + i] * E4[jj+nStencilWidth/2];
      data[0+3*2] += gf[(j+jj)*(nPhi*nR) + i] * D4[jj+nStencilWidth/2];
      for (int ii = -nStencilWidth/2; ii <= nStencilWidth/2; ++ii)
      {
        data[1+3*1] += gf[(j+jj)*(nPhi*nR) + i+ii] * E4[ii+nStencilWidth/2]*E4[jj+nStencilWidth/2];
        data[2+3*1] += gf[(j+jj)*(nPhi*nR) + i+ii] * D4[ii+nStencilWidth/2]*E4[jj+nStencilWidth/2];
        data[1+3*2] += gf[(j+jj)*(nPhi*nR) + i+ii] * E4[ii+nStencilWidth/2]*D4[jj+nStencilWidth/2];
        data[2+3*2] += gf[(j+jj)*(nPhi*nR) + i+ii] * D4[ii+nStencilWidth/2]*D4[jj+nStencilWidth/2];
      }
    }
  }
  return (0);
}

PetscScalar*   createHermiteFD(int *hdata_nR, int *hdata_nZ, const PetscScalar *gf, const PetscScalar *g, const int nStencilWidth, const int bIntersect, const int nR, const int nP, const int nZ){
  *hdata_nR = 0;
  *hdata_nZ = 0;
  int stp = nStencilWidth;

  if (bIntersect == 1)
    stp = nStencilWidth - 1;
  for (int j = nStencilWidth/2; j < nZ-nStencilWidth/2; j+= stp)
  {
    (*hdata_nZ)++;
  }
  for (int i = nStencilWidth/2; i < nR-nStencilWidth/2; i+= stp)
  {
    (*hdata_nR)++;
  }

  PetscScalar * primal = (PetscScalar*) malloc(sizeof(PetscScalar)*12* *hdata_nR * *hdata_nZ);
  int jj = 0;
  for (int j = nStencilWidth/2; j < nZ-nStencilWidth/2; j+= stp)
  {
    int ii = 0;
    for (int i = nStencilWidth/2; i < nR-nStencilWidth/2; i+= stp)
    {
      getHermiteDataFD(gf,  nR, nP, nZ, primal + jj* *hdata_nR*12 + ii*12, i, j, nStencilWidth, bIntersect);
      primal[jj* *hdata_nR*12 + ii*12 +  9] = g[3*i + 3*j*nR*nP    ];
      primal[jj* *hdata_nR*12 + ii*12 + 10] = g[3*i + 3*j*nR*nP + 1];
      primal[jj* *hdata_nR*12 + ii*12 + 11] = g[3*i + 3*j*nR*nP + 2];
      ii++;
    }
    jj++;
  }
  return primal;
}

PetscErrorCode multiplybyR(PetscScalar *gf, const PetscScalar *g, const int N){
  for (int i = 0; i < N; ++i){
    gf[i] *= g[3*i];
  }
  return (0);
}
