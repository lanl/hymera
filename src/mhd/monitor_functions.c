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

#include "RunawayKineticSolverWrapper.h"


PetscErrorCode DumpVelocity_Cell(TS ts, PetscInt step, Vec X, char* prefix, void * ptr) {
  User * user = (User * ) ptr;
  DM da, dmC, daC, dmBAvg, daBAvg, dmEAvg, daEAvg, dmJAvg, daJAvg, dmV, daV, dmEP, daEP;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec J, J_local, X_local, vecC, C, vecBAvg, BAvg, vecEAvg, EAvg, vecJAvg, JAvg, vecV, V, vecEP, EP;
  PetscReal time = 0.0;

  TSGetDM(ts, & da);
  TSGetTime(ts, & time);

  //PetscPrintf(PETSC_COMM_WORLD,"Current time: t = %f\n", time);

  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmV); /* 3 dofs per element */

  DMSetUp(dmV);

  DMStagSetUniformCoordinatesExplicit(dmV, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);

  DMCreateGlobalVector(dmV, & V);

  DMGetLocalVector(da, & X_local);
  DMGlobalToLocal(da, X, INSERT_VALUES, X_local);

  DMStagGetCorners(dmV, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[24], to[3];
        PetscScalar valFrom[24], valTo[3];

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
        to[0].loc = ELEMENT;
        to[0].c = 0;
        valTo[0] = (valFrom[0]+valFrom[1]+valFrom[2]+valFrom[3]+valFrom[4]+valFrom[5]+valFrom[6]+valFrom[7]) / 8.0;
        to[1].i = er;
        to[1].j = ephi;
        to[1].k = ez;
        to[1].loc = ELEMENT;
        to[1].c = 1;
        valTo[1] = (valFrom[8]+valFrom[9]+valFrom[10]+valFrom[11]+valFrom[12]+valFrom[13]+valFrom[14]+valFrom[15]) / 8.0;
        to[2].i = er;
        to[2].j = ephi;
        to[2].k = ez;
        to[2].loc = ELEMENT;
        to[2].c = 2;
        valTo[2] = (valFrom[16]+valFrom[17]+valFrom[18]+valFrom[19]+valFrom[20]+valFrom[21]+valFrom[22]+valFrom[23]) / 8.0;

        DMStagVecSetValuesStencil(dmV, V, 3, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);

  DMRestoreLocalVector(da, & X_local);

  DMStagVecSplitToDMDA(dmV, V, ELEMENT, -3, & daV, & vecV); /* note -3 : pad with zero in 2D case */

  if(prefix[0] == 'd'){
    PetscObjectSetName((PetscObject) vecV, "Velocity time derivative");
  }
  else if(prefix[0] == 'l'){
    PetscObjectSetName((PetscObject) vecV, "Velocity laplacian");
  }
  else{
    PetscObjectSetName((PetscObject) vecV, "Velocity");
  }

  /* Dump element-based fields to a .vtr file and create a .pvd file */
  {
    PetscViewer viewerC, viewerB, viewerE, viewerJ, viewerV, viewerEP;
    char filename[PETSC_MAX_PATH_LEN];
    FILE * pvdfile;


    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_%savg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", prefix, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
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

PetscErrorCode DumpSolution_Cell(TS ts, PetscInt step, Vec X, void * ptr) {
  User * user = (User * ) ptr;
  DM da, dmC, daC, dmBAvg, daBAvg, dmEAvg, daEAvg, dmJAvg, daJAvg, dmV, daV, dmEP, daEP;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec J, J_local, X_local, vecC, C, vecBAvg, BAvg, vecEAvg, EAvg, vecJAvg, JAvg, vecV, V, vecEP, EP;
  PetscReal time = 0.0;

  TSGetDM(ts, & da);
  TSGetTime(ts, & time);

  PetscPrintf(PETSC_COMM_WORLD,"Current time: t = %f\n", time);

  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 1, & dmEP); /* 1 dof per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmV); /* 3 dofs per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmBAvg); /* 3 dofs per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmEAvg); /* 3 dofs per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmJAvg); /* 3 dofs per element */
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
  VecScale(J, 1.0/user->L0);
  DMGetLocalVector(da, & J_local);
  DMGlobalToLocal(da, J, INSERT_VALUES, J_local);

  DMStagGetCorners(dmEP, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[8], to[1];
        PetscScalar valFrom[8], valTo[1];

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
        to[0].loc = ELEMENT;
        to[0].c = 0;
        valTo[0] = (valFrom[0]+valFrom[1]+valFrom[2]+valFrom[3]+valFrom[4]+valFrom[5]+valFrom[6]+valFrom[7]) / 8.0;

        DMStagVecSetValuesStencil(dmEP, EP, 1, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(EP);
  VecAssemblyEnd(EP);
  // PetscPrintf(PETSC_COMM_WORLD,"line 340");

  DMStagGetCorners(dmV, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        DMStagStencil from[24], to[3];
        PetscScalar valFrom[24], valTo[3];

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
        to[0].loc = ELEMENT;
        to[0].c = 0;
        valTo[0] = (valFrom[0]+valFrom[1]+valFrom[2]+valFrom[3]+valFrom[4]+valFrom[5]+valFrom[6]+valFrom[7]) / 8.0;
        to[1].i = er;
        to[1].j = ephi;
        to[1].k = ez;
        to[1].loc = ELEMENT;
        to[1].c = 1;
        valTo[1] = (valFrom[8]+valFrom[9]+valFrom[10]+valFrom[11]+valFrom[12]+valFrom[13]+valFrom[14]+valFrom[15]) / 8.0;
        to[2].i = er;
        to[2].j = ephi;
        to[2].k = ez;
        to[2].loc = ELEMENT;
        to[2].c = 2;
        valTo[2] = (valFrom[16]+valFrom[17]+valFrom[18]+valFrom[19]+valFrom[20]+valFrom[21]+valFrom[22]+valFrom[23]) / 8.0;

        DMStagVecSetValuesStencil(dmV, V, 3, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(V);
  VecAssemblyEnd(V);
    // PetscPrintf(PETSC_COMM_WORLD,"line 496");
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
  // PetscPrintf(PETSC_COMM_WORLD,"line 562");
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
  // PetscPrintf(PETSC_COMM_WORLD,"line 777");
  DMStagVecSplitToDMDA(dmEP, EP, ELEMENT, -1, & daEP, & vecEP); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmV, V, ELEMENT, -3, & daV, & vecV); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmBAvg, BAvg, ELEMENT, -3, & daBAvg, & vecBAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmEAvg, EAvg, ELEMENT, -3, & daEAvg, & vecEAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmJAvg, JAvg, ELEMENT, -3, & daJAvg, & vecJAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmC, C, ELEMENT, -1, & daC, & vecC); /* note -3 : pad with zero in 2D case */

  PetscObjectSetName((PetscObject) vecEP, "Electrostatic Potential");
  PetscObjectSetName((PetscObject) vecV, "Velocity");
  PetscObjectSetName((PetscObject) vecBAvg, "Magnetic Field");
  PetscObjectSetName((PetscObject) vecEAvg, "Divergence-free part of Electric Field");
  PetscObjectSetName((PetscObject) vecJAvg, "Curl of Magnetic Field");
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

PetscErrorCode DumpError(TS ts, PetscInt step, Vec X, void * ptr) {
  User * user = (User * ) ptr;
  DM da, dmC, daC, dmBAvg, daBAvg, dmEAvg, daEAvg, dmV, daV;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec X_local, vecC, C, vecBAvg, BAvg, vecEAvg, EAvg, vecV, V;

  TSGetDM(ts, & da);
  DMStagCreateCompatibleDMStag(da, 3, 0, 0, 0, & dmV); /* 3 dofs per vertex */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmEAvg); /* 3 dof per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmBAvg); /* 3 dof per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 1, & dmC); /* 1 dof per element */

  DMSetUp(dmC);
  DMSetUp(dmBAvg);
  DMSetUp(dmEAvg);
  DMSetUp(dmV);

  DMStagSetUniformCoordinatesExplicit(dmC, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
  DMStagSetUniformCoordinatesExplicit(dmBAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
  DMStagSetUniformCoordinatesExplicit(dmEAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);
  DMStagSetUniformCoordinatesExplicit(dmV, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);

  DMCreateGlobalVector(dmC, & C);
  DMCreateGlobalVector(dmBAvg, & BAvg);
  DMCreateGlobalVector(dmEAvg, & EAvg);
  DMCreateGlobalVector(dmV, & V);

  DMGetLocalVector(da, & X_local);
  DMGlobalToLocal(da, X, INSERT_VALUES, X_local);

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

  DMStagVecSplitToDMDA(dmC, C, ELEMENT, -1, & daC, & vecC); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmBAvg, BAvg, ELEMENT, -3, & daBAvg, & vecBAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmEAvg, EAvg, ELEMENT, -3, & daEAvg, & vecEAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmV, V, BACK_DOWN_LEFT, -3, & daV, & vecV); /* note -3 : pad with zero in 2D case */

  PetscObjectSetName((PetscObject) vecV, "Velocity Error");
  PetscObjectSetName((PetscObject) vecBAvg, "Magnetic Field Error (Averaged)");
  PetscObjectSetName((PetscObject) vecEAvg, "Electric Field Error (Averaged)");
  PetscObjectSetName((PetscObject) vecC, "Ions' Number Density Error");

  /* Dump element-based fields to a .vtr file */
  {
    PetscViewer viewerC, viewerB, viewerE, viewerV;
    char filename[PETSC_MAX_PATH_LEN];

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_nerravg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daC), filename, FILE_MODE_WRITE, & viewerC);
    VecView(vecC, viewerC);
    PetscViewerDestroy( & viewerC);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_berravg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daBAvg), filename, FILE_MODE_WRITE, & viewerB);
    VecView(vecBAvg, viewerB);
    PetscViewerDestroy( & viewerB);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_eerravg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daEAvg), filename, FILE_MODE_WRITE, & viewerE);
    VecView(vecEAvg, viewerE);
    PetscViewerDestroy( & viewerE);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_verravg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
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

  VecDestroy( & vecBAvg);
  DMDestroy( & daBAvg);
  VecDestroy( & BAvg);
  DMDestroy( & dmBAvg);

  VecDestroy( & vecEAvg);
  DMDestroy( & daEAvg);
  VecDestroy( & EAvg);
  DMDestroy( & dmEAvg);

  VecDestroy( & vecC);
  DMDestroy( & daC);
  VecDestroy( & C);
  DMDestroy( & dmC);

  return (0);
}

PetscErrorCode DumpDivergence(TS ts, DM newda, PetscInt step, Vec X, void * ptr) {
  User * user = (User * ) ptr;
  DM dmBAvg, daBAvg;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec X_local, vecBAvg, BAvg;

  DMStagCreateCompatibleDMStag(newda, 0, 0, 0, 1, & dmBAvg); /* 3 dof per element */

  DMSetUp(dmBAvg);

  DMStagSetUniformCoordinatesExplicit(dmBAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);

  DMCreateGlobalVector(dmBAvg, & BAvg);

  DMGetLocalVector(newda, & X_local);
  DMGlobalToLocal(newda, X, INSERT_VALUES, X_local);

  DMStagGetCorners(dmBAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
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
        DMStagVecGetValuesStencil(newda, X_local, 1, from, valFrom);
        to[0].i = er;
        to[0].j = ephi;
        to[0].k = ez;
        to[0].loc = ELEMENT;
        to[0].c = 0;
        valTo[0] = valFrom[0];
        DMStagVecSetValuesStencil(dmBAvg, BAvg, 1, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(BAvg);
  VecAssemblyEnd(BAvg);

  DMRestoreLocalVector(newda, & X_local);

  DMStagVecSplitToDMDA(dmBAvg, BAvg, ELEMENT, -1, & daBAvg, & vecBAvg); /* note -3 : pad with zero in 2D case */

  //DMStagVecSplitToDMDA(newda,X,ELEMENT,-1,&daBAvg,&vecBAvg); /* note -3 : pad with zero in 2D case */

  PetscObjectSetName((PetscObject) vecBAvg, "Divergence of Magnetic Field");

  /* Dump element-based fields to a .vtr file */
  {
    PetscViewer viewerB;
    char filename[PETSC_MAX_PATH_LEN];

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_divb_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daBAvg), filename, FILE_MODE_WRITE, & viewerB);
    VecView(vecBAvg, viewerB);
    PetscViewerDestroy( & viewerB);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);
  }

  /* Destroy DMDAs and Vecs */
  VecDestroy( & vecBAvg);
  DMDestroy( & daBAvg);
  VecDestroy( & BAvg);
  DMDestroy( & dmBAvg);

  return (0);
}

PetscErrorCode DumpLevelSet(TS ts, void * ptr) {
  User * user = (User * ) ptr;
  DM da, dmBAvg, daBAvg;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz, N[3];
  Vec vecBAvg, BAvg;

  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);

  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 1, & dmBAvg); /* 1 dof per element */

  DMSetUp(dmBAvg);

  DMStagSetUniformCoordinatesExplicit(dmBAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);

  DMCreateGlobalVector(dmBAvg, & BAvg);

  DMStagGetCorners(dmBAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
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
        valTo[0] = user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]];
        DMStagVecSetValuesStencil(dmBAvg, BAvg, 1, to, valTo, INSERT_VALUES);
      }
    }
  }
  VecAssemblyBegin(BAvg);
  VecAssemblyEnd(BAvg);

  DMStagVecSplitToDMDA(dmBAvg, BAvg, ELEMENT, -1, & daBAvg, & vecBAvg); /* note -3 : pad with zero in 2D case */

  PetscObjectSetName((PetscObject) vecBAvg, "Levelset function");

  /* Dump element-based fields to a .vtr file */
  {
    PetscViewer viewerB;
    char filename[PETSC_MAX_PATH_LEN];

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_lvlset_ic%.1D_grid%.2Dx%.2Dx%.2D.vtr", user -> ictype, user -> Nr, user -> Nphi, user -> Nz);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daBAvg), filename, FILE_MODE_WRITE, & viewerB);
    VecView(vecBAvg, viewerB);
    PetscViewerDestroy( & viewerB);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);
  }

  /* Destroy DMDAs and Vecs */
  VecDestroy( & vecBAvg);
  DMDestroy( & daBAvg);
  VecDestroy( & BAvg);
  DMDestroy( & dmBAvg);

  return (0);
}

PetscErrorCode SaveIntermediateSolution(TS ts, PetscInt step, PetscReal time, Vec X, void * ptr) {
  User * user = (User * ) ptr;
  char filename[PETSC_MAX_PATH_LEN];
  PetscViewer viewerX;

  /* Write X in binary for later use */
  PetscSNPrintf(filename, sizeof(filename), "%s/X_ic%.2D_grid%.2Dx%.2Dx%.2D_step%.3D_time%5.7f.dat", user->input_folder, user -> ictype, user -> Nr, user -> Nphi, user -> Nz, (int) step, (double) time);
  PetscPrintf(PETSC_COMM_WORLD, "Writing X vector into file %s ...\n", filename);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, & viewerX);
  VecView(X, viewerX);
  PetscViewerDestroy( & viewerX);
  PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);
  return(0);
}

PetscErrorCode ComputeCurrent(TS ts, Vec X, void * ptr) {
  User * user = (User * ) ptr;
  DM da;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  PetscInt N[3], ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  Vec xLocal, JLocal, J, GradEP, GradEPLocal;
  PetscScalar ** ** arrGradEP, ** ** arrJ, ** ** arrX, Javg, I1, I2, I3;

  TSGetDM(ts, & da);
  DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
  VecDuplicate(X, & J);
  VecCopy(X, J);

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

  // Compute J:= (1/mu0) curl(B)
  FormDerivedCurlnores(ts, X, J, user);
  // Multiply J by B_0/L_0
  VecScale(J, user->B0/user->L0);

  DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  if(0){
  DMGetLocalVector(da, & xLocal);
  DMGlobalToLocalBegin(da, X, INSERT_VALUES, xLocal);
  DMGlobalToLocalEnd(da, X, INSERT_VALUES, xLocal);
  DMStagVecGetArray(da, xLocal, & arrX);

  DMGetLocalVector(da, & JLocal);
  DMStagVecGetArray(da, JLocal, & arrJ);
  {
  // Compute Grad(EP)
    Mat G;
    /* Compute the gradient of EP */
    VecDuplicate(X, & GradEP);
    VecCopy(X, GradEP);
    DMCreateMatrix(da, & G);
    MatZeroEntries(G);
    FormDiscreteGradientEP(ts, G, X, GradEP, user);
    MatDestroy( & G);
    DMGetLocalVector(da, & GradEPLocal);
    DMGlobalToLocalBegin(da, GradEP, INSERT_VALUES, GradEPLocal);
    DMGlobalToLocalEnd(da, GradEP, INSERT_VALUES, GradEPLocal);
    DMStagVecGetArrayRead(da, GradEPLocal, & arrGradEP);
  }

  // Compute J:= (1/resistivity) * (tau - Grad(EP))
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {
        arrJ[ez][ephi][er][ivErmzm] = rese(er, ephi, ez, BACK_LEFT, user) * (arrX[ez][ephi][er][ivErmzm] - arrGradEP[ez][ephi][er][ivErmzm]);
        arrJ[ez][ephi][er][ivEphimzm] = rese(er, ephi, ez, BACK_DOWN, user) * (arrX[ez][ephi][er][ivEphimzm] - arrGradEP[ez][ephi][er][ivEphimzm]);
        arrJ[ez][ephi][er][ivErmphim] = rese(er, ephi, ez, DOWN_LEFT, user) * (arrX[ez][ephi][er][ivErmphim] - arrGradEP[ez][ephi][er][ivErmphim]);

        if (er == N[0] - 1) {
          arrJ[ez][ephi][er][ivErpzm] = rese(er, ephi, ez, BACK_RIGHT, user) * (arrX[ez][ephi][er][ivErpzm] - arrGradEP[ez][ephi][er][ivErpzm]);
          arrJ[ez][ephi][er][ivErpphim] = rese(er, ephi, ez, DOWN_RIGHT, user) * (arrX[ez][ephi][er][ivErpphim] - arrGradEP[ez][ephi][er][ivErpphim]);
        }
        if (ephi == N[1] - 1) {
          arrJ[ez][ephi][er][ivEphipzm] = rese(er, ephi, ez, BACK_UP, user) * (arrX[ez][ephi][er][ivEphipzm] - arrGradEP[ez][ephi][er][ivEphipzm]);
          arrJ[ez][ephi][er][ivErmphip] = rese(er, ephi, ez, UP_LEFT, user) * (arrX[ez][ephi][er][ivErmphip] - arrGradEP[ez][ephi][er][ivErmphip]);
        }
        if (ez == N[2] - 1) {
          arrJ[ez][ephi][er][ivErmzp] = rese(er, ephi, ez, FRONT_LEFT, user) * (arrX[ez][ephi][er][ivErmzp] - arrGradEP[ez][ephi][er][ivErmzp]);
          arrJ[ez][ephi][er][ivEphimzp] = rese(er, ephi, ez, FRONT_DOWN, user) * (arrX[ez][ephi][er][ivEphimzp] - arrGradEP[ez][ephi][er][ivEphimzp]);
        }
        if (er == N[0] - 1 && ephi == N[1] - 1) {
          arrJ[ez][ephi][er][ivErpphip] = rese(er, ephi, ez, UP_RIGHT, user) * (arrX[ez][ephi][er][ivErpphip] - arrGradEP[ez][ephi][er][ivErpphip]);
        }
        if (ephi == N[1] - 1 && ez == N[2] - 1) {
          arrJ[ez][ephi][er][ivEphipzp] = rese(er, ephi, ez, FRONT_UP, user) * (arrX[ez][ephi][er][ivEphipzp] - arrGradEP[ez][ephi][er][ivEphipzp]);
        }
        if (er == N[0] - 1 && ez == N[2] - 1) {
          arrJ[ez][ephi][er][ivErpzp] = rese(er, ephi, ez, FRONT_RIGHT, user) * (arrX[ez][ephi][er][ivErpzp] - arrGradEP[ez][ephi][er][ivErpzp]);
        }
      }
    }
  }

  DMStagVecRestoreArray(da, JLocal, & arrJ);
  DMLocalToGlobal(da, JLocal, INSERT_VALUES, J);

  DMStagVecRestoreArray(da, xLocal, & arrX);
  DMRestoreLocalVector(da, & xLocal);

  DMStagVecRestoreArrayRead(da, GradEPLocal, & arrGradEP);
  DMRestoreLocalVector(da, & GradEPLocal);

  VecDestroy( & GradEP);
  }

  DMGetLocalVector(da, & JLocal);
  DMGlobalToLocalBegin(da, J, INSERT_VALUES, JLocal);
  DMGlobalToLocalEnd(da, J, INSERT_VALUES, JLocal);
  DMStagVecGetArray(da, JLocal, & arrJ);

  user -> Iphi1 = 0.0;
  user -> Iphi2 = 0.0;
  user -> Iphi3 = 0.0;
  I1 = 0.0;
  I2 = 0.0;
  I3 = 0.0;

  for (ez = startz; ez < startz + nz; ++ez) {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
      for (er = startr; er < startr + nr; ++er) {

        if (ephi == 1) {
          if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
            Javg = 0.25 * (arrJ[ez][ephi][er][ivErmzm] + arrJ[ez][ephi][er][ivErpzm] + arrJ[ez][ephi][er][ivErmzp] + arrJ[ez][ephi][er][ivErpzp]);
            I2 += Javg * user -> dr * user -> dz;
          }
          if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.5) < 0.7) {
            Javg = 0.25 * (arrJ[ez][ephi][er][ivErmzm] + arrJ[ez][ephi][er][ivErpzm] + arrJ[ez][ephi][er][ivErmzp] + arrJ[ez][ephi][er][ivErpzp]);
            I1 += Javg * user -> dr * user -> dz;
          }
          if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 0.1) {
            Javg = 0.25 * (arrJ[ez][ephi][er][ivErmzm] + arrJ[ez][ephi][er][ivErpzm] + arrJ[ez][ephi][er][ivErmzp] + arrJ[ez][ephi][er][ivErpzp]);
            I3 += Javg * user -> dr * user -> dz;
          }
        }

      }
    }
  }

  I1 *= user->L0 * user->L0;
  I2 *= user->L0 * user->L0;
  I3 *= user->L0 * user->L0;

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_SELF, "Local Current intensity inside plasma = %g\n", I1);
    PetscPrintf(PETSC_COMM_SELF, "Local Current intensity outside plasma = %g\n", I2);
  }

  MPI_Reduce( & I1, & (user -> Iphi1), 1, MPI_DOUBLE, MPI_SUM, 0,
    PETSC_COMM_WORLD);

  MPI_Reduce( & I2, & (user -> Iphi2), 1, MPI_DOUBLE, MPI_SUM, 0,
    PETSC_COMM_WORLD);

  MPI_Reduce( & I3, & (user -> Iphi3), 1, MPI_DOUBLE, MPI_SUM, 0,
    PETSC_COMM_WORLD);

  if (user -> debug) {
    PetscPrintf(PETSC_COMM_WORLD, "Current intensity inside plasma = %g\n", user -> Iphi1);
    PetscPrintf(PETSC_COMM_WORLD, "Current intensity outside plasma = %g\n", user -> Iphi2);
  }

  DMStagVecRestoreArray(da, JLocal, & arrJ);
  DMRestoreLocalVector(da, & JLocal);
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

PetscErrorCode DumpEdgeField(TS ts, PetscInt step, Vec X, void * ptr) {
  User * user = (User * ) ptr;
  DM da, dmC, daC, dmBAvg, daBAvg, dmEAvg, daEAvg, dmJAvg, daJAvg, dmV, daV, dmEP, daEP;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec X_local, vecC, C, vecBAvg, BAvg, vecEAvg, EAvg, vecJAvg, JAvg, vecV, V, vecEP, EP;
  PetscReal time = 0.0;

  TSGetDM(ts, & da);
  TSGetTime(ts, & time);

  // PetscPrintf(PETSC_COMM_WORLD,"Current time: t = %f\n", time);

  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmEAvg); /* 3 dof per element */

  DMSetUp(dmEAvg);

  DMStagSetUniformCoordinatesExplicit(dmEAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax);

  DMCreateGlobalVector(dmEAvg, & EAvg);

  DMGetLocalVector(da, & X_local);
  DMGlobalToLocal(da, X, INSERT_VALUES, X_local);

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

  DMRestoreLocalVector(da, & X_local);

  DMStagVecSplitToDMDA(dmEAvg, EAvg, ELEMENT, -3, & daEAvg, & vecEAvg); /* note -3 : pad with zero in 2D case */

  PetscObjectSetName((PetscObject) vecEAvg, "Edge Field (Averaged)");

  /* Dump element-based fields to a .vtr file and create a .pvd file */
  {
    PetscViewer viewerC, viewerB, viewerE, viewerJ, viewerV, viewerEP;
    char filename[PETSC_MAX_PATH_LEN];
    FILE * pvdfile;

    PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_edgeavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step);
    PetscViewerVTKOpen(PetscObjectComm((PetscObject) daEAvg), filename, FILE_MODE_WRITE, & viewerE);
    VecView(vecEAvg, viewerE);
    PetscViewerDestroy( & viewerE);
    PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);
  }

  /* Destroy DMDAs and Vecs */
  VecDestroy( & vecEAvg);
  DMDestroy( & daEAvg);
  VecDestroy( & EAvg);
  DMDestroy( & dmEAvg);

  return (0);
}

PetscErrorCode multiplybyR(PetscScalar *gf, const PetscScalar *g, const int N){
  for (int i = 0; i < N; ++i){
    gf[i] *= g[3*i];
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

PetscErrorCode interpolateHermite(PetscScalar *primal, PetscScalar *dual, int hdata_nR, int hdata_nZ, void * ptr) {
  User * user = (User * ) ptr;
  for (int j = 0; j < hdata_nZ-1; ++j)
  {
    for (int i = 0; i < hdata_nR-1; ++i)
    {
      interpolate2D(primal + j * 12 * hdata_nR + i*12, dual + j * 36*(hdata_nR-1) + i*36, hdata_nR, user);
    }
  }
  return (0);
}

PetscErrorCode interpolate2D(PetscScalar *primal, PetscScalar *dual, int nR, void * ptr) {
  User * user = (User * ) ptr;
  int idx = 0;
  int idy = 0;
  int m = 2;

  PetscScalar * blu = primal;
  PetscScalar * bru = primal + 12;
  PetscScalar * tlu = primal + 12*(nR);
  PetscScalar * tru = primal + 12*(nR)+12;

  PetscScalar cofsyt[(2+1)*2*(2+1)];
  PetscScalar cofsyb[(2+1)*2*(2+1)];

    for (idy = 0; idy < m+1; ++idy)
    {
        interpolatex(idy, m, blu, bru, cofsyb, user);
        interpolatex(idy, m, tlu, tru, cofsyt, user);
    }

    for (idx = 0; idx < 2*(m+1); ++idx)
    {
        interpolatey(idx, m, cofsyb, cofsyt, dual, user);
    }

  return (0);
}

PetscErrorCode interpolatey(const int idx, const int m, PetscScalar *bdata, PetscScalar *tdata, PetscScalar *ac, void * ptr) {
  User * user = (User * ) ptr;
  int i,idy, k,j;

  for (i = 0; i < m+1; ++i){
    for (idy = 0; idy < m+1 - i; ++idy){
      user -> NT[i][idy] = bdata[i*2*(m+1) + idx];
      user -> NT[i][idy + m+1] = tdata[i*2*(m+1) + idx];
    }
  }

  //Fill in missing values between known data.
  for (i = 1; i < m+1; ++i){
    for (idy = m+1 - i; idy < m+1; ++idy){
      user -> NT[i][idy] = user -> NT[i-1][idy+1] - user -> NT[i-1][idy];
    }
  }

  //Fill in final part of a table
  for (i = m+1; i < 2*(m+1); ++i){
    for (idy = 0; idy < 2*(m+1) - i; ++idy){
      user -> NT[i][idy] = user -> NT[i-1][idy+1] - user -> NT[i-1][idy];
    }
  }

    //Get coefficients
  for(int i = 0; i < 2*(m+1); ++i){
    ac[i*2*(m+1) + idx]= user -> NT[i][0];
  }

  //Change basis
  for (k = 2*m; k >= 0; --k){
    for (j = k; j < 2*m + 1; ++j){
      if (k < m+1)
        ac[j*2*(m+1) + idx] = ac[j*2*(m+1) + idx] + 0.5 * ac[(j+1)*2*(m+1) + idx];
      else
        ac[j*2*(m+1) + idx] = ac[j*2*(m+1) + idx] - 0.5 * ac[(j+1)*2*(m+1) + idx];
    }
  }
  return(0);
}

PetscErrorCode interpolatex(const int idy, const int m, PetscScalar *ldata, PetscScalar *rdata, PetscScalar *ac, void * ptr) {
  User * user = (User * ) ptr;
  int i,idx, k,j;

  for (i = 0; i < m+1; ++i){
    for (idx = 0; idx < m+1 - i; ++idx){
      user -> NT[i][idx] = ldata[idy*(m+1) + i];
      user -> NT[i][idx + m+1] = rdata[idy*(m+1) + i];
    }
  }

  //Fill in missing values between known data.
  for (i = 1; i < m+1; ++i){
    for (idx = m+1 - i; idx < m+1; ++idx){
      user -> NT[i][idx] = user -> NT[i-1][idx+1] - user -> NT[i-1][idx];
    }
  }

  //Fill in final part of a table
  for (i = m+1; i < 2*(m+1); ++i){
    for (idx = 0; idx < 2*(m+1) - i; ++idx){
      user -> NT[i][idx] = user -> NT[i-1][idx+1] - user -> NT[i-1][idx];
    }
  }
  //Get coefficients
  for(i = 0; i < 2*(m+1); ++i){
    ac[idy*2*(m+1) + i]= user -> NT[i][0];
  }

  // Change basis
  for (k = 2*m; k >= 0; --k){
    for (j = k; j < 2*m + 1; ++j){
      if (k < m+1)
        ac[idy*2*(m+1) + j] = ac[idy*2*(m+1) + j] + 0.5 * ac[idy*2*(m+1) + j+1];
      else
        ac[idy*2*(m+1) + j] = ac[idy*2*(m+1) + j] - 0.5 * ac[idy*2*(m+1) + j+1];
    }
  }
  return(0);
}

PetscErrorCode integrateHermite_1D_r(PetscScalar *hd, PetscScalar *hi, int nR, int nZ, const PetscScalar hR, const PetscScalar hZ){
  for (int i = 0; i < 7*(nR-1); ++i) hi[i] = 0.0;
  for (int i = 0; i < nR-1; ++i)
  {
    PetscScalar II = 0.0;
    PetscScalar zloc = 1.0;
    for (int idy = 0; idy < 2*(2+1); ++idy)
    {
      PetscScalar scl = 1.0;
      for (int idx = 0; idx < 2*(2+1); idx+=2)
      {
        II += hd[idx + idy * 6 + i * 36] * hR / (idx+1) * scl * zloc;
        scl /= 4.0;
      }
      zloc *= -0.375;
    }
    for (int ii = i+1; ii < nR-1; ++ii)
    {
      hi[ii*7] += II;
    }
    zloc = 1.0;
    for (int idy = 0; idy < 2*(2+1); ++idy)
    {
      PetscScalar rloc = 0.5;
      for (int idx = 0; idx < 2*(2+1); ++idx)
      {
        hi[i*7] += hd[idx + idy*6 + i*36] * rloc * hR / (idx+1) * zloc;
        hi[(idx+1) + i*7] += hd[idx + idy*6 + i*36] * hR / (idx+1) * zloc;
        rloc *= -0.5;
      }
      zloc *= -0.375;
    }
  }
  return (0);
}

PetscErrorCode integrateHermite_2D_z(PetscScalar *hd, PetscScalar *hi, int nR, int nZ, const PetscScalar hR, const PetscScalar hZ){
  for (int i = 0; i < 42*(nR-1)*(nZ-1); ++i) hi[i] = 0.0;
  PetscScalar II[2*(2+1)];
  PetscScalar scl = 0.0;

  for (int j = 0; j < nZ-1; ++j)
  {
    for (int i = 0; i < nR-1; ++i)
    {
      for (int idx = 0; idx < 2*(2+1); ++idx) II[idx] = 0.0;
      scl = 1.0;

      for (int idy = 0; idy < 2*(2+1); idy+=2)
      {
        for (int idx = 0; idx < 2*(2+1); ++idx)
        {
          II[idx] += hd[idx + idy*6 + i * 36 + j * 36 * (nR-1)]* hZ / (PetscScalar)(idy+1) * scl;
        }
        scl /= 4.0;
      }
      for (int jj = j+1; jj < nZ-1; ++jj)
      {
        for (int idx = 0; idx < 2*(2+1); ++idx)
        {
          hi[idx + i * 42  + jj * 42 * (nR-1)] += II[idx];
        }
      }
      scl = 0.5;
      for (int idy = 0; idy < 2*(2+1); ++idy)
      {
        for (int idx = 0; idx < 2*(2+1); ++idx)
        {
          hi[idx + i * 42 + j * 42 * (nR-1)] += hd[idx + idy * 6 + i * 36 + j * 36 * (nR-1)] * scl * hZ / (PetscScalar)(idy+1);
          hi[idx + (idy+1) * 6 + i * 42 + j * 42 * (nR-1)] += hd[idx + idy * 6 + i * 36 + j * 36 * (nR-1)] * hZ / (PetscScalar)(idy+1);
        }
        scl *= -0.5;
      }
    }
  }
  return (0);
}

PetscScalar evalHermite2D(const PetscScalar *gR, const PetscScalar *gZ, const PetscScalar R, const PetscScalar Z, const PetscScalar hR, const PetscScalar hZ, const PetscScalar *hd, const int dr, const int dz, const int mr, const int mz, const int nR, const int nZ){
  tLocalCoordinate loc = getLocalCoordinate(gR, gZ, R, Z, hR, hZ);
  if (0 == 1)
  {
    printf("FATAL: local index out of bounds!\n");
    printf("R = %le, Rmin = %le, Rmax = %le\n", R, gR[0], gR[nR-1]);
    printf("Z = %le, Zmin = %le, Zmax = %le\n", Z, gZ[0], gZ[nR-1]);
    exit(0);
  }

  if (loc.iR < 0){
    loc.iR = 0; loc.R = -0.5;
  }
  else if ( loc.iR > nR-2){
    loc.iR = nR-2; loc.R = 0.5;
  }
  if (loc.iZ < 0){
    loc.iZ = 0; loc.Z = -0.5;
  }
  else if ( loc.iZ > nZ-2){
    loc.iZ = nZ-2; loc.Z = 0.5;
  }
  PetscScalar val = 0.0;
  PetscScalar sclz = 1.0;
  PetscScalar sclr = 1.0;
  for (int idy = 1; idy <= dz; ++idy) sclz /= hZ;
  for (int idy = dz; idy < mz; ++idy){
    int cz = 1;
    for (int i = 0; i < dz; ++i) cz*=(idy-i);
    sclr = 1.0;
    for (int idx = 1; idx <= dr; ++idx) sclr /= hR;
    for (int idx = dr; idx < mr; ++idx){
      int cr = 1;
      for (int i = 0; i < dr; ++i) cr*=(idx-i);

      val += sclr * sclz * hd[idx + idy * mr + loc.iR * mr * mz + loc.iZ * mr * mz * (nR-1)] * cr * cz;
      sclr *= loc.R;
    }
    sclz *= loc.Z;
  }
  return val;
}

PetscScalar evalHermite1D(const PetscScalar *gR, const PetscScalar *gZ, const PetscScalar R, const PetscScalar hR, const PetscScalar *hd, const int dr, const int mr, const int nR){
  tLocalCoordinate loc = getLocalCoordinate(gR, gZ, R, 0.0, hR, 0.3);

  //if (loc.iR < 0 || loc.iR > nR-1){
    //printf("FATAL: local index out of bounds!\n");
    //printf("R = %le", R);
    //exit(0);
  //}

  if (loc.iR < 0){
    loc.iR = 0; loc.R = -0.5;
  }
  else if ( loc.iR > nR-2){
    loc.iR = nR-2; loc.R = 0.5;
  }

  PetscScalar val = 0.0;
  PetscScalar sclr = 1.0;
  for (int idx = 1; idx <= dr; ++idx) sclr /= hR;
  for (int idx = dr; idx < mr; ++idx){
    int cr = 1;
    for (int i = 0; i < dr; ++i) cr*=(idx-i);
    val += sclr * hd[idx + loc.iR * mr] * cr;
    sclr *= loc.R;
  }
  return val;
}

tLocalCoordinate getLocalCoordinate(const PetscScalar *gR, const PetscScalar *gZ, const  PetscScalar R, const  PetscScalar Z, const PetscScalar hR, const PetscScalar hZ){
  tLocalCoordinate loc;
  loc.iR = (int)ceil((R - gR[0])/hR)-1;
  loc.iZ = (int)ceil((Z - gZ[0])/hZ)-1;
  loc.R = (R - gR[loc.iR]-0.5*hR)/hR;
  loc.Z = (Z - gZ[loc.iZ]-0.5*hZ)/hZ;

  return loc;
}

PetscErrorCode evalPsi(PetscScalar *psi, const PetscScalar *R, const PetscScalar *Z, int N, tHermiteDivFreeFields f){
  PetscScalar hR = f.m_hp_RR[1] - f.m_hp_RR[0];
  PetscScalar hZ = f.m_hp_ZR[1] - f.m_hp_ZR[0];

  //PetscPrintf(PETSC_COMM_WORLD, "Inside call to evalPsi: Before for loop.\n");
  //tLocalCoordinate loc = getLocalCoordinate(f.m_hp_RZ, f.m_hp_ZZ, R[0], 0.0, hR, 0.3);
  //printf("loc.iR = %d\n", loc.iR);
  //printf("f.m_hp_RZ = %g\n", f.m_hp_RZ[0]);
  //printf("R[0] = %g\n", R[0]);
  //printf("(R[0] - f.m_hp_RZ[0])/hR = %g\n", (R[0] - f.m_hp_RZ[0])/hR);

  for (int i = 0; i < N; ++i){
    psi[i] = evalHermite1D(f.m_hp_RZ, f.m_hp_ZZ, R[i], hR, f.m_hd_psi_ir, 0, 7, f.m_nRZ)
           - evalHermite2D(f.m_hp_RR, f.m_hp_ZR, R[i], Z[i], hR, hZ, f.m_hd_psi_iz, 0, 0, 6, 7, f.m_nRR, f.m_nZR);
    //PetscPrintf(PETSC_COMM_WORLD, "Inside call to evalPsi: Inside for loop --> i = %d.\n", i);
  }
  //PetscPrintf(PETSC_COMM_WORLD, "Inside call to evalPsi: After for loop.\n");
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

PetscErrorCode TSAdaptChoose_user(TSAdapt adapt, TS ts, PetscReal h, PetscInt *next_sc, PetscReal *next_h, PetscBool *accept, PetscReal *wlte, PetscReal *wltea, PetscReal *wlter) {
   static PetscReal dtTarget = -1.0;
   PetscReal        dtInitial;
   DM               dm;
   User          *ctx;
   PetscInt         step, n_record, iter;

   PetscFunctionBeginUser;
   PetscCall(TSGetDM(ts, &dm));
   PetscCall(DMGetApplicationContext(dm, &ctx));
   PetscCall(TSGetStepNumber(ts, &step));

   PetscReal currentstepsize;
   TSGetTimeStep(ts,&currentstepsize);
   PetscPrintf(PETSC_COMM_WORLD, "Current step size = %g\n", currentstepsize );
   PetscPrintf(PETSC_COMM_WORLD, "Current h = %g\n", h);

  /* double RunawayCurrent = RuKS_getIre(ctx->runaway_solver); */
  SNES snes;
  TSGetSNES(ts, & snes);
  SNESConvergedReason reason = SNES_CONVERGED_ITERATING;
  PetscCall(SNESGetConvergedReason(snes, &reason));
  SNESGetIterationNumber(snes, &iter);

  PetscPrintf(PETSC_COMM_WORLD, "Converged reason : %d\n", reason);

  n_record = ctx->n_record;
  if (step>0 && reason>0) {n_record++;}
  else {n_record=0;}
  // ADD CRITERION FOR INCREASE OF STEP SIZE: NUMBER OF ITERATIONS, GENTLER INCREASE OF STEP SIZE
  if (step>0 && PetscAbsScalar(1-ctx->prev_current/ctx->present_current) <=0.10 && reason>0 && h< 10 * ctx->V_A*ctx->L0*ctx->mu0/(100.0*ctx->etaplasma) && iter<6 && ctx->CorrectorIdentifier == 0) {
    *accept = PETSC_TRUE;
    *next_h = h*2.0;
    n_record = 0;
    ctx->dt = h*2.0;
    ctx->n_record_Steady_jRE = step;
  }
  else if(step > ctx->n_record_Steady_jRE && reason > 0 && PetscAbsScalar(1-ctx->prev_current/ctx->present_current) > 0.10 && ctx->CorrectorIdentifier == 1){
    PetscPrintf(PETSC_COMM_WORLD, "DID NOT CONVERGE. Converged reason : %d\n", reason);
    *accept = PETSC_FALSE;
    // *next_h = PetscMin(h*0.5, 0.001 * ctx->V_A * ctx->L0 * ctx->mu0 / ctx->etaplasma);
    *next_h = h*0.5; //PetscMin(h*0.5, 0.001 * ctx->V_A * ctx->L0 * ctx->mu0 / ctx->etaplasma);
    ctx->dt= h*0.5;
  }

  // else if (reason<0) {
  //   *accept = PETSC_FALSE;
  //   *next_h = PetscMin(h*0.5, 0.001 * ctx->V_A * ctx->L0 * ctx->mu0 / ctx->etaplasma);
  // }
  else{
    *accept = PETSC_TRUE;
    *next_h = h;
  }

  ctx->n_record = n_record;

   PetscPrintf(PETSC_COMM_WORLD, "next_h = %g\n", *next_h);
   //PetscPrintf(PETSC_COMM_WORLD, "ctx->n_record = %d\n", ctx->n_record);
   if(step>0){
        PetscPrintf(PETSC_COMM_WORLD, "RunawayCurrent=%g, Previous Runaway Current = %g\n", ctx->present_current, ctx->prev_current);
        PetscPrintf(PETSC_COMM_WORLD, "|RunawayCurrent - PreviousRunawayCurrent| / RunawayCurrent = %g\n", PetscAbsScalar(1-ctx->prev_current/ctx->present_current));
   }

   *next_sc = 0;  /* Reuse the same order scheme */
   *wlte    = -1; /* Weighted local truncation error was not evaluated */
   *wltea   = -1; /* Weighted absolute local truncation error was not evaluated */
   *wlter   = -1; /* Weighted relative local truncation error was not evaluated */
   PetscFunctionReturn(0);
}
