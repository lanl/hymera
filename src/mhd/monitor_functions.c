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


PetscErrorCode DumpVelocity_Cell(TS ts, PetscInt step, Vec X, char* prefix, void * ptr) {
  PetscFunctionBeginUser;
  User * user = (User * ) ptr;
  DM da, dmC, daC, dmBAvg, daBAvg, dmEAvg, daEAvg, dmJAvg, daJAvg, dmV, daV, dmEP, daEP;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec J, J_local, X_local, vecC, C, vecBAvg, BAvg, vecEAvg, EAvg, vecJAvg, JAvg, vecV, V, vecEP, EP;
  PetscReal time = 0.0;

  PetscCall(TSGetDM(ts, & da));
  PetscCall(TSGetTime(ts, & time));

  //PetscPrintf(PETSC_COMM_WORLD,"Current time: t = %f\n", time);

  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmV); /* 3 dofs per element */

  PetscCall(DMSetUp(dmV));

  PetscCall(DMStagSetUniformCoordinatesExplicit(dmV, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));

  PetscCall(DMCreateGlobalVector(dmV, & V));

  PetscCall(DMGetLocalVector(da, & X_local));
  PetscCall(DMGlobalToLocal(da, X, INSERT_VALUES, X_local));

  PetscCall(DMStagGetCorners(dmV, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
          PetscCall(DMStagVecGetValuesStencil(da, X_local, 24, from, valFrom));

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

        PetscCall(DMStagVecSetValuesStencil(dmV, V, 3, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(V));
  PetscCall(VecAssemblyEnd(V));

  PetscCall(DMRestoreLocalVector(da, & X_local));

  DMStagVecSplitToDMDA(dmV, V, ELEMENT, -3, & daV, & vecV); /* note -3 : pad with zero in 2D case */

  if(prefix[0] == 'd'){
    PetscCall(PetscObjectSetName((PetscObject) vecV, "Velocity time derivative"));
  }
  else if(prefix[0] == 'l'){
    PetscCall(PetscObjectSetName((PetscObject) vecV, "Velocity laplacian"));
  }
  else{
    PetscCall(PetscObjectSetName((PetscObject) vecV, "Velocity"));
  }

  /* Dump element-based fields to a .vtr file and create a .pvd file */
  {
    PetscViewer viewerC, viewerB, viewerE, viewerJ, viewerV, viewerEP;
    char filename[PETSC_MAX_PATH_LEN];
    FILE * pvdfile;


    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_%savg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", prefix, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daV), filename, FILE_MODE_WRITE, & viewerV));
    PetscCall(VecView(vecV, viewerV));
    PetscCall(PetscViewerDestroy( & viewerV));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));
  }

  /* Destroy DMDAs and Vecs */
  PetscCall(VecDestroy( & vecV));
  PetscCall(DMDestroy( & daV));
  PetscCall(VecDestroy( & V));
  PetscCall(DMDestroy( & dmV));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DumpSolution_Cell(TS ts, PetscInt step, Vec X, void * ptr) {
  PetscFunctionBeginUser;
  User * user = (User * ) ptr;
  DM da, dmC, daC, dmBAvg, daBAvg, dmEAvg, daEAvg, dmJAvg, daJAvg, dmV, daV, dmEP, daEP;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec J, J_local, X_local, vecC, C, vecBAvg, BAvg, vecEAvg, EAvg, vecJAvg, JAvg, vecV, V, vecEP, EP;
  PetscReal time = 0.0;

  PetscCall(TSGetDM(ts, & da));
  PetscCall(TSGetTime(ts, & time));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Current time: t = %f\n", time));

  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 1, & dmEP); /* 1 dof per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmV); /* 3 dofs per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmBAvg); /* 3 dofs per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmEAvg); /* 3 dofs per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmJAvg); /* 3 dofs per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 1, & dmC); /* 1 dof per element */

  PetscCall(DMSetUp(dmEP));
  PetscCall(DMSetUp(dmV));
  PetscCall(DMSetUp(dmBAvg));
  PetscCall(DMSetUp(dmEAvg));
  PetscCall(DMSetUp(dmJAvg));
  PetscCall(DMSetUp(dmC));

  PetscCall(DMStagSetUniformCoordinatesExplicit(dmEP, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmV, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmBAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmEAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmJAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmC, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));

  PetscCall(DMCreateGlobalVector(dmEP, & EP));
  PetscCall(DMCreateGlobalVector(dmV, & V));
  PetscCall(DMCreateGlobalVector(dmBAvg, & BAvg));
  PetscCall(DMCreateGlobalVector(dmEAvg, & EAvg));
  PetscCall(DMCreateGlobalVector(dmJAvg, & JAvg));
  PetscCall(DMCreateGlobalVector(dmC, & C));

  PetscCall(DMGetLocalVector(da, & X_local));
  PetscCall(DMGlobalToLocal(da, X, INSERT_VALUES, X_local));

  // Compute \tilde{J}:= curl(\tilde{B})
  PetscCall(VecDuplicate(X,&J));
  PetscCall(VecCopy(X, J));
  FormDerivedCurlnomp(ts, X, J, user);
  PetscCall(VecScale(J, 1.0/user->L0));
  PetscCall(DMGetLocalVector(da, & J_local));
  PetscCall(DMGlobalToLocal(da, J, INSERT_VALUES, J_local));

  PetscCall(DMStagGetCorners(dmEP, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
        PetscCall(DMStagVecGetValuesStencil(da, X_local, 8, from, valFrom));

        to[0].i = er;
        to[0].j = ephi;
        to[0].k = ez;
        to[0].loc = ELEMENT;
        to[0].c = 0;
        valTo[0] = (valFrom[0]+valFrom[1]+valFrom[2]+valFrom[3]+valFrom[4]+valFrom[5]+valFrom[6]+valFrom[7]) / 8.0;

        PetscCall(DMStagVecSetValuesStencil(dmEP, EP, 1, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(EP));
  PetscCall(VecAssemblyEnd(EP));
  // PetscPrintf(PETSC_COMM_WORLD,"line 340");

  PetscCall(DMStagGetCorners(dmV, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
          PetscCall(DMStagVecGetValuesStencil(da, X_local, 24, from, valFrom));

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

        PetscCall(DMStagVecSetValuesStencil(dmV, V, 3, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(V));
  PetscCall(VecAssemblyEnd(V));
    // PetscPrintf(PETSC_COMM_WORLD,"line 496");
  PetscCall(DMStagGetCorners(dmBAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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

        PetscCall(DMStagVecGetValuesStencil(da, X_local, 6, from, valFrom));

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

        PetscCall(DMStagVecSetValuesStencil(dmBAvg, BAvg, 3, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(BAvg));
  PetscCall(VecAssemblyEnd(BAvg));
  // PetscPrintf(PETSC_COMM_WORLD,"line 562");
  PetscCall(DMStagGetCorners(dmEAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
        PetscCall(DMStagVecGetValuesStencil(da, X_local, 12, from, valFrom));
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
        PetscCall(DMStagVecSetValuesStencil(dmEAvg, EAvg, 3, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(EAvg));
  PetscCall(VecAssemblyEnd(EAvg));

  PetscCall(DMStagGetCorners(dmJAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
        PetscCall(DMStagVecGetValuesStencil(da, J_local, 12, from, valFrom));
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
        PetscCall(DMStagVecSetValuesStencil(dmJAvg, JAvg, 3, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(JAvg));
  PetscCall(VecAssemblyEnd(JAvg));

  PetscCall(DMStagGetCorners(dmC, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
          PetscCall(DMStagVecGetValuesStencil(da, X_local, 1, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = ELEMENT;
          to[0].c = 0;
          valTo[0] = valFrom[0];
          PetscCall(DMStagVecSetValuesStencil(dmC, C, 1, to, valTo, INSERT_VALUES));
        }
      }
    }
    PetscCall(VecAssemblyBegin(C));
    PetscCall(VecAssemblyEnd(C));

  PetscCall(DMRestoreLocalVector(da, & X_local));
  PetscCall(DMRestoreLocalVector(da, & J_local));
  // PetscPrintf(PETSC_COMM_WORLD,"line 777");
  DMStagVecSplitToDMDA(dmEP, EP, ELEMENT, -1, & daEP, & vecEP); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmV, V, ELEMENT, -3, & daV, & vecV); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmBAvg, BAvg, ELEMENT, -3, & daBAvg, & vecBAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmEAvg, EAvg, ELEMENT, -3, & daEAvg, & vecEAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmJAvg, JAvg, ELEMENT, -3, & daJAvg, & vecJAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmC, C, ELEMENT, -1, & daC, & vecC); /* note -3 : pad with zero in 2D case */

  PetscCall(PetscObjectSetName((PetscObject) vecEP, "Electrostatic Potential"));
  PetscCall(PetscObjectSetName((PetscObject) vecV, "Velocity"));
  PetscCall(PetscObjectSetName((PetscObject) vecBAvg, "Magnetic Field"));
  PetscCall(PetscObjectSetName((PetscObject) vecEAvg, "Divergence-free part of Electric Field"));
  PetscCall(PetscObjectSetName((PetscObject) vecJAvg, "Curl of Magnetic Field"));
  PetscCall(PetscObjectSetName((PetscObject) vecC, "Number Density of Ions"));

  /* Dump element-based fields to a .vtr file and create a .pvd file */
  {
    PetscViewer viewerC, viewerB, viewerE, viewerJ, viewerV, viewerEP;
    char filename[PETSC_MAX_PATH_LEN];
    FILE * pvdfile;

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D.pvd", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz));
    if (time == 0.0) {
      PetscCall(PetscFOpen(PETSC_COMM_WORLD, filename, "a", & pvdfile));
      //pvdfile = fopen(filename, "a");
      //pvdfile.open (filename, ios::out | ios::app);
      PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "<\?xml version=\"1.0\"?>\n"));
      //pvdfile << "<\?xml version=\"1.0\"?>\n";
      PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n"));
      //pvdfile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
      PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "  <Collection>\n"));
      //pvdfile << "  <Collection>\n";
      PetscCall(PetscFClose(PETSC_COMM_WORLD, pvdfile));
      //fclose(pvdfile);
      //pvdfile.close();
    }

    PetscCall(PetscFOpen(PETSC_COMM_WORLD, filename, "a", & pvdfile));
    PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"0\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_vavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"1\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_bavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"2\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_tavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"3\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_navg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"4\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_javg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"5\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_lvlset_ic%.1D_grid%.2Dx%.2Dx%.2D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> Nr, user -> Nphi, user -> Nz));
    PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "    <DataSet timestep=\"%f\" group=\"\" part=\"6\" file=\"mfd_data_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D/mfd_epavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr\"/>\n", time, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscFClose(PETSC_COMM_WORLD, pvdfile));

    if (time >= user -> ftime) {
      PetscCall(PetscFOpen(PETSC_COMM_WORLD, filename, "a", & pvdfile));
      PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "  </Collection>\n"));
      PetscCall(PetscFPrintf(PETSC_COMM_WORLD, pvdfile, "</VTKFile>\n"));
      PetscCall(PetscFClose(PETSC_COMM_WORLD, pvdfile));
    }

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_epavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daEP), filename, FILE_MODE_WRITE, & viewerEP));
    PetscCall(VecView(vecEP, viewerEP));
    PetscCall(PetscViewerDestroy( & viewerEP));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_vavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daV), filename, FILE_MODE_WRITE, & viewerV));
    PetscCall(VecView(vecV, viewerV));
    PetscCall(PetscViewerDestroy( & viewerV));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_bavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daBAvg), filename, FILE_MODE_WRITE, & viewerB));
    PetscCall(VecView(vecBAvg, viewerB));
    PetscCall(PetscViewerDestroy( & viewerB));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_tavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daEAvg), filename, FILE_MODE_WRITE, & viewerE));
    PetscCall(VecView(vecEAvg, viewerE));
    PetscCall(PetscViewerDestroy( & viewerE));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_javg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daJAvg), filename, FILE_MODE_WRITE, & viewerJ));
    PetscCall(VecView(vecJAvg, viewerJ));
    PetscCall(PetscViewerDestroy( & viewerJ));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_navg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daC), filename, FILE_MODE_WRITE, & viewerC));
    PetscCall(VecView(vecC, viewerC));
    PetscCall(PetscViewerDestroy( & viewerC));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));
  }

  /* Destroy DMDAs and Vecs */
  PetscCall(VecDestroy( & vecEP));
  PetscCall(DMDestroy( & daEP));
  PetscCall(VecDestroy( & EP));
  PetscCall(DMDestroy( & dmEP));

  PetscCall(VecDestroy( & vecV));
  PetscCall(DMDestroy( & daV));
  PetscCall(VecDestroy( & V));
  PetscCall(DMDestroy( & dmV));

  PetscCall(VecDestroy( & vecBAvg));
  PetscCall(DMDestroy( & daBAvg));
  PetscCall(VecDestroy( & BAvg));
  PetscCall(DMDestroy( & dmBAvg));

  PetscCall(VecDestroy( & vecEAvg));
  PetscCall(DMDestroy( & daEAvg));
  PetscCall(VecDestroy( & EAvg));
  PetscCall(DMDestroy( & dmEAvg));

  PetscCall(VecDestroy( & vecJAvg));
  PetscCall(DMDestroy( & daJAvg));
  PetscCall(VecDestroy( & JAvg));
  PetscCall(DMDestroy( & dmJAvg));

  PetscCall(VecDestroy( & vecC));
  PetscCall(DMDestroy( & daC));
  PetscCall(VecDestroy( & C));
  PetscCall(DMDestroy( & dmC));

  PetscCall(VecDestroy( & J));

  PetscFunctionReturn(PETSC_SUCCESS);
}

#line 1762

PetscErrorCode DumpError(TS ts, PetscInt step, Vec X, void * ptr) {
  PetscFunctionBeginUser;
  User * user = (User * ) ptr;
  DM da, dmC, daC, dmBAvg, daBAvg, dmEAvg, daEAvg, dmV, daV;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec X_local, vecC, C, vecBAvg, BAvg, vecEAvg, EAvg, vecV, V;

  PetscCall(TSGetDM(ts, & da));
  DMStagCreateCompatibleDMStag(da, 3, 0, 0, 0, & dmV); /* 3 dofs per vertex */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmEAvg); /* 3 dof per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmBAvg); /* 3 dof per element */
  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 1, & dmC); /* 1 dof per element */

  PetscCall(DMSetUp(dmC));
  PetscCall(DMSetUp(dmBAvg));
  PetscCall(DMSetUp(dmEAvg));
  PetscCall(DMSetUp(dmV));

  PetscCall(DMStagSetUniformCoordinatesExplicit(dmC, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmBAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmEAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmV, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));

  PetscCall(DMCreateGlobalVector(dmC, & C));
  PetscCall(DMCreateGlobalVector(dmBAvg, & BAvg));
  PetscCall(DMCreateGlobalVector(dmEAvg, & EAvg));
  PetscCall(DMCreateGlobalVector(dmV, & V));

  PetscCall(DMGetLocalVector(da, & X_local));
  PetscCall(DMGlobalToLocal(da, X, INSERT_VALUES, X_local));

    PetscCall(DMStagGetCorners(dmC, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
            PetscCall(DMStagVecGetValuesStencil(da, X_local, 1, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = ELEMENT;
            to[0].c = 0;
            valTo[0] = valFrom[0];
            PetscCall(DMStagVecSetValuesStencil(dmC, C, 1, to, valTo, INSERT_VALUES));
          }
        }
      }
      PetscCall(VecAssemblyBegin(C));
      PetscCall(VecAssemblyEnd(C));

  PetscCall(DMStagGetCorners(dmBAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
        PetscCall(DMStagVecGetValuesStencil(da, X_local, 6, from, valFrom));
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
        PetscCall(DMStagVecSetValuesStencil(dmBAvg, BAvg, 3, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(BAvg));
  PetscCall(VecAssemblyEnd(BAvg));

  PetscCall(DMStagGetCorners(dmEAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
        PetscCall(DMStagVecGetValuesStencil(da, X_local, 12, from, valFrom));
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
        PetscCall(DMStagVecSetValuesStencil(dmEAvg, EAvg, 3, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(EAvg));
  PetscCall(VecAssemblyEnd(EAvg));

    PetscCall(DMStagGetCorners(dmV, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
            PetscCall(DMStagVecGetValuesStencil(da, X_local, 24, from, valFrom));

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

          PetscCall(DMStagVecSetValuesStencil(dmV, V, 24, to, valTo, INSERT_VALUES));
        }
      }
    }
    PetscCall(VecAssemblyBegin(V));
    PetscCall(VecAssemblyEnd(V));

  PetscCall(DMRestoreLocalVector(da, & X_local));

  DMStagVecSplitToDMDA(dmC, C, ELEMENT, -1, & daC, & vecC); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmBAvg, BAvg, ELEMENT, -3, & daBAvg, & vecBAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmEAvg, EAvg, ELEMENT, -3, & daEAvg, & vecEAvg); /* note -3 : pad with zero in 2D case */
  DMStagVecSplitToDMDA(dmV, V, BACK_DOWN_LEFT, -3, & daV, & vecV); /* note -3 : pad with zero in 2D case */

  PetscCall(PetscObjectSetName((PetscObject) vecV, "Velocity Error"));
  PetscCall(PetscObjectSetName((PetscObject) vecBAvg, "Magnetic Field Error (Averaged)"));
  PetscCall(PetscObjectSetName((PetscObject) vecEAvg, "Electric Field Error (Averaged)"));
  PetscCall(PetscObjectSetName((PetscObject) vecC, "Ions' Number Density Error"));

  /* Dump element-based fields to a .vtr file */
  {
    PetscViewer viewerC, viewerB, viewerE, viewerV;
    char filename[PETSC_MAX_PATH_LEN];

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_nerravg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daC), filename, FILE_MODE_WRITE, & viewerC));
    PetscCall(VecView(vecC, viewerC));
    PetscCall(PetscViewerDestroy( & viewerC));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_berravg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daBAvg), filename, FILE_MODE_WRITE, & viewerB));
    PetscCall(VecView(vecBAvg, viewerB));
    PetscCall(PetscViewerDestroy( & viewerB));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_eerravg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daEAvg), filename, FILE_MODE_WRITE, & viewerE));
    PetscCall(VecView(vecEAvg, viewerE));
    PetscCall(PetscViewerDestroy( & viewerE));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_verravg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daV), filename, FILE_MODE_WRITE, & viewerV));
    PetscCall(VecView(vecV, viewerV));
    PetscCall(PetscViewerDestroy( & viewerV));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));
  }

  /* Destroy DMDAs and Vecs */
  PetscCall(VecDestroy( & vecV));
  PetscCall(DMDestroy( & daV));
  PetscCall(VecDestroy( & V));
  PetscCall(DMDestroy( & dmV));

  PetscCall(VecDestroy( & vecBAvg));
  PetscCall(DMDestroy( & daBAvg));
  PetscCall(VecDestroy( & BAvg));
  PetscCall(DMDestroy( & dmBAvg));

  PetscCall(VecDestroy( & vecEAvg));
  PetscCall(DMDestroy( & daEAvg));
  PetscCall(VecDestroy( & EAvg));
  PetscCall(DMDestroy( & dmEAvg));

  PetscCall(VecDestroy( & vecC));
  PetscCall(DMDestroy( & daC));
  PetscCall(VecDestroy( & C));
  PetscCall(DMDestroy( & dmC));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DumpDivergence(TS ts, DM newda, PetscInt step, Vec X, void * ptr) {
  PetscFunctionBeginUser;
  User * user = (User * ) ptr;
  DM dmBAvg, daBAvg;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec X_local, vecBAvg, BAvg;

  DMStagCreateCompatibleDMStag(newda, 0, 0, 0, 1, & dmBAvg); /* 3 dof per element */

  PetscCall(DMSetUp(dmBAvg));

  PetscCall(DMStagSetUniformCoordinatesExplicit(dmBAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));

  PetscCall(DMCreateGlobalVector(dmBAvg, & BAvg));

  PetscCall(DMGetLocalVector(newda, & X_local));
  PetscCall(DMGlobalToLocal(newda, X, INSERT_VALUES, X_local));

  PetscCall(DMStagGetCorners(dmBAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
        PetscCall(DMStagVecGetValuesStencil(newda, X_local, 1, from, valFrom));
        to[0].i = er;
        to[0].j = ephi;
        to[0].k = ez;
        to[0].loc = ELEMENT;
        to[0].c = 0;
        valTo[0] = valFrom[0];
        PetscCall(DMStagVecSetValuesStencil(dmBAvg, BAvg, 1, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(BAvg));
  PetscCall(VecAssemblyEnd(BAvg));

  PetscCall(DMRestoreLocalVector(newda, & X_local));

  DMStagVecSplitToDMDA(dmBAvg, BAvg, ELEMENT, -1, & daBAvg, & vecBAvg); /* note -3 : pad with zero in 2D case */

  //DMStagVecSplitToDMDA(newda,X,ELEMENT,-1,&daBAvg,&vecBAvg); /* note -3 : pad with zero in 2D case */

  PetscCall(PetscObjectSetName((PetscObject) vecBAvg, "Divergence of Magnetic Field"));

  /* Dump element-based fields to a .vtr file */
  {
    PetscViewer viewerB;
    char filename[PETSC_MAX_PATH_LEN];

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_divb_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daBAvg), filename, FILE_MODE_WRITE, & viewerB));
    PetscCall(VecView(vecBAvg, viewerB));
    PetscCall(PetscViewerDestroy( & viewerB));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));
  }

  /* Destroy DMDAs and Vecs */
  PetscCall(VecDestroy( & vecBAvg));
  PetscCall(DMDestroy( & daBAvg));
  PetscCall(VecDestroy( & BAvg));
  PetscCall(DMDestroy( & dmBAvg));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DumpLevelSet(TS ts, void * ptr) {
  PetscFunctionBeginUser;
  User * user = (User * ) ptr;
  DM da, dmBAvg, daBAvg;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz, N[3];
  Vec vecBAvg, BAvg;

  PetscCall(TSGetDM(ts, & da));
  PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));

  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 1, & dmBAvg); /* 1 dof per element */

  PetscCall(DMSetUp(dmBAvg));

  PetscCall(DMStagSetUniformCoordinatesExplicit(dmBAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));

  PetscCall(DMCreateGlobalVector(dmBAvg, & BAvg));

  PetscCall(DMStagGetCorners(dmBAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
        PetscCall(DMStagVecSetValuesStencil(dmBAvg, BAvg, 1, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(BAvg));
  PetscCall(VecAssemblyEnd(BAvg));

  DMStagVecSplitToDMDA(dmBAvg, BAvg, ELEMENT, -1, & daBAvg, & vecBAvg); /* note -3 : pad with zero in 2D case */

  PetscCall(PetscObjectSetName((PetscObject) vecBAvg, "Levelset function"));

  /* Dump element-based fields to a .vtr file */
  {
    PetscViewer viewerB;
    char filename[PETSC_MAX_PATH_LEN];

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_lvlset_ic%.1D_grid%.2Dx%.2Dx%.2D.vtr", user -> ictype, user -> Nr, user -> Nphi, user -> Nz));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daBAvg), filename, FILE_MODE_WRITE, & viewerB));
    PetscCall(VecView(vecBAvg, viewerB));
    PetscCall(PetscViewerDestroy( & viewerB));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));
  }

  /* Destroy DMDAs and Vecs */
  PetscCall(VecDestroy( & vecBAvg));
  PetscCall(DMDestroy( & daBAvg));
  PetscCall(VecDestroy( & BAvg));
  PetscCall(DMDestroy( & dmBAvg));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode SaveIntermediateSolution(TS ts, PetscInt step, PetscReal time, Vec X, void * ptr) {
  PetscFunctionBeginUser;
  User * user = (User * ) ptr;
  char filename[PETSC_MAX_PATH_LEN];
  PetscViewer viewerX;

  /* Write X in binary for later use */
  PetscCall(PetscSNPrintf(filename, sizeof(filename), "%s/X_ic%.2D_grid%.2Dx%.2Dx%.2D_step%.3D_time%5.7f.dat", user->input_folder, user -> ictype, user -> Nr, user -> Nphi, user -> Nz, (int) step, (double) time));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Writing X vector into file %s ...\n", filename));
  PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, & viewerX));
  PetscCall(VecView(X, viewerX));
  PetscCall(PetscViewerDestroy( & viewerX));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode ComputeCurrent(TS ts, Vec X, void * ptr) {
  PetscFunctionBeginUser;
  User * user = (User * ) ptr;
  DM da;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  PetscInt N[3], ivErmzm, ivErmzp, ivErpzm, ivErpzp;
  PetscInt ivEphimzm, ivEphipzm, ivEphimzp, ivEphipzp;
  PetscInt ivErmphim, ivErpphim, ivErmphip, ivErpphip;
  Vec xLocal, JLocal, J, GradEP, GradEPLocal;
  PetscScalar ** ** arrGradEP, ** ** arrJ, ** ** arrX, Javg, I1, I2, I3;

  PetscCall(TSGetDM(ts, & da));
  PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
  PetscCall(VecDuplicate(X, & J));
  PetscCall(VecCopy(X, J));

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

  // Compute J:= (1/mu0) curl(B)
  FormDerivedCurlnores(ts, X, J, user);
  // Multiply J by B_0/L_0
  PetscCall(VecScale(J, user->B0/user->L0));

  PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));

#line 2582

  PetscCall(DMGetLocalVector(da, & JLocal));
  PetscCall(DMGlobalToLocalBegin(da, J, INSERT_VALUES, JLocal));
  PetscCall(DMGlobalToLocalEnd(da, J, INSERT_VALUES, JLocal));
  PetscCall(DMStagVecGetArray(da, JLocal, & arrJ));

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
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "Local Current intensity inside plasma = %g\n", I1));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "Local Current intensity outside plasma = %g\n", I2));
  }

  MPI_Reduce( & I1, & (user -> Iphi1), 1, MPI_DOUBLE, MPI_SUM, 0,
    PETSC_COMM_WORLD);

  MPI_Reduce( & I2, & (user -> Iphi2), 1, MPI_DOUBLE, MPI_SUM, 0,
    PETSC_COMM_WORLD);

  MPI_Reduce( & I3, & (user -> Iphi3), 1, MPI_DOUBLE, MPI_SUM, 0,
    PETSC_COMM_WORLD);

  if (user -> debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Current intensity inside plasma = %g\n", user -> Iphi1));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Current intensity outside plasma = %g\n", user -> Iphi2));
  }

  PetscCall(DMStagVecRestoreArray(da, JLocal, & arrJ));
  PetscCall(DMRestoreLocalVector(da, & JLocal));
  PetscCall(VecDestroy( & J));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#line 2996

PetscErrorCode DumpEdgeField(TS ts, PetscInt step, Vec X, void * ptr) {
  PetscFunctionBeginUser;
  User * user = (User * ) ptr;
  DM da, dmC, daC, dmBAvg, daBAvg, dmEAvg, daEAvg, dmJAvg, daJAvg, dmV, daV, dmEP, daEP;
  PetscInt er, ephi, ez, startr, startphi, startz, nr, nphi, nz;
  Vec X_local, vecC, C, vecBAvg, BAvg, vecEAvg, EAvg, vecJAvg, JAvg, vecV, V, vecEP, EP;
  PetscReal time = 0.0;

  PetscCall(TSGetDM(ts, & da));
  PetscCall(TSGetTime(ts, & time));

  // PetscPrintf(PETSC_COMM_WORLD,"Current time: t = %f\n", time);

  DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmEAvg); /* 3 dof per element */

  PetscCall(DMSetUp(dmEAvg));

  PetscCall(DMStagSetUniformCoordinatesExplicit(dmEAvg, user -> rmin, user -> rmax, user -> phimin, user -> phimax, user -> zmin, user -> zmax));

  PetscCall(DMCreateGlobalVector(dmEAvg, & EAvg));

  PetscCall(DMGetLocalVector(da, & X_local));
  PetscCall(DMGlobalToLocal(da, X, INSERT_VALUES, X_local));

  PetscCall(DMStagGetCorners(dmEAvg, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
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
        PetscCall(DMStagVecGetValuesStencil(da, X_local, 12, from, valFrom));
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
        PetscCall(DMStagVecSetValuesStencil(dmEAvg, EAvg, 3, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(EAvg));
  PetscCall(VecAssemblyEnd(EAvg));

  PetscCall(DMRestoreLocalVector(da, & X_local));

  DMStagVecSplitToDMDA(dmEAvg, EAvg, ELEMENT, -3, & daEAvg, & vecEAvg); /* note -3 : pad with zero in 2D case */

  PetscCall(PetscObjectSetName((PetscObject) vecEAvg, "Edge Field (Averaged)"));

  /* Dump element-based fields to a .vtr file and create a .pvd file */
  {
    PetscViewer viewerC, viewerB, viewerE, viewerJ, viewerV, viewerEP;
    char filename[PETSC_MAX_PATH_LEN];
    FILE * pvdfile;

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "vtrfiles/mfd_edgeavg_ic%.1D_ts%.1D_grid%.2Dx%.2Dx%.2D_step%.3D.vtr", user -> ictype, user -> tstype, user -> Nr, user -> Nphi, user -> Nz, step));
    PetscCall(PetscViewerVTKOpen(PetscObjectComm((PetscObject) daEAvg), filename, FILE_MODE_WRITE, & viewerE));
    PetscCall(VecView(vecEAvg, viewerE));
    PetscCall(PetscViewerDestroy( & viewerE));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename));
  }

  /* Destroy DMDAs and Vecs */
  PetscCall(VecDestroy( & vecEAvg));
  PetscCall(DMDestroy( & daEAvg));
  PetscCall(VecDestroy( & EAvg));
  PetscCall(DMDestroy( & dmEAvg));

  PetscFunctionReturn(PETSC_SUCCESS);
}

#line 3147

#line 3181

#line 3264


#line 3350

