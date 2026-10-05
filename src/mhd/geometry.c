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
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmpatch.h>
#include <petscsf.h>

#include "mhd.h"

#define PETSC_NULL_VEC PETSC_NULLPTR

PetscScalar cyldistance(PetscScalar r1, PetscScalar phi1, PetscScalar z1, PetscScalar r2, PetscScalar phi2, PetscScalar z2) {
  PetscScalar distance;
  distance = PetscSqrtScalar(PetscSqr(r1) + PetscSqr(r2) - 2.0 * r1 * r2 * PetscCosScalar(phi1 - phi2) + PetscSqr(z1 - z2));
  return distance;
}

PetscScalar surface(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, d, N[3];
  PetscInt icp[3];
  PetscInt icBrp[3], icBphip[3], icBzp[3], icBrm[3], icBphim[3], icBzm[3];
  PetscInt icErmzm[3], icErmzp[3], icErpzm[3], icErpzp[3];
  PetscInt icEphimzm[3], icEphipzm[3], icEphimzp[3], icEphipzp[3];
  PetscInt icErmphim[3], icErpphim[3], icErmphip[3], icErpphip[3];
  PetscInt icrmphimzm[3], icrpphimzm[3], icrmphipzm[3], icrpphipzm[3];
  PetscInt icrmphimzp[3], icrpphimzp[3], icrmphipzp[3], icrpphipzp[3];
  DM dmCoorda, coordDA = user -> coorda;
  PetscScalar ** ** arrCoord = user->arrCoord;
  PetscScalar surf;

  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  DMGetCoordinateDM(coordDA, & dmCoorda);
  for (d = 0; d < 3; ++d) {
    /* Element coordinates */
    DMStagGetLocationSlot(dmCoorda, ELEMENT, d, & icp[d]);
    /* Face coordinates */
    DMStagGetLocationSlot(dmCoorda, LEFT, d, & icBrm[d]);
    DMStagGetLocationSlot(dmCoorda, DOWN, d, & icBphim[d]);
    DMStagGetLocationSlot(dmCoorda, BACK, d, & icBzm[d]);
    DMStagGetLocationSlot(dmCoorda, RIGHT, d, & icBrp[d]);
    DMStagGetLocationSlot(dmCoorda, UP, d, & icBphip[d]);
    DMStagGetLocationSlot(dmCoorda, FRONT, d, & icBzp[d]);
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

  DMStagGetGlobalSizes(user -> coorda, & N[0], & N[1], & N[2]);

  /* Faces perpendicular to z direction */
  if (loc == BACK) {
    if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
      surf = user -> dphi * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icErpzm[0]]) - PetscSqr(arrCoord[ez][ephi][er][icErmzm[0]])) / 2.0; /* INT(r dphi dr) */
    } else {
      surf = PetscAbsReal(arrCoord[ez][ephi][er][icEphipzm[1]] - arrCoord[ez][ephi][er][icEphimzm[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icErpzm[0]]) - PetscSqr(arrCoord[ez][ephi][er][icErmzm[0]])) / 2.0; /* INT(r dphi dr) */
    }
  } else if (loc == FRONT) {
    if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
      surf = user -> dphi * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icErpzp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icErmzp[0]])) / 2.0; /* INT(r dphi dr) */
    } else {
      surf = PetscAbsReal(arrCoord[ez][ephi][er][icEphipzp[1]] - arrCoord[ez][ephi][er][icEphimzp[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icErpzp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icErmzp[0]])) / 2.0; /* INT(r dphi dr) */
    }
  }
  /* Faces perpendicular to phi direction */
  else if (loc == DOWN) {
    surf = PetscAbsReal(arrCoord[ez][ephi][er][icEphimzp[2]] - arrCoord[ez][ephi][er][icEphimzm[2]]) * PetscAbsReal(arrCoord[ez][ephi][er][icErpphim[0]] - arrCoord[ez][ephi][er][icErmphim[0]]); /* INT(dr dz) */

  } else if (loc == UP) {
    surf = PetscAbsReal(arrCoord[ez][ephi][er][icEphipzp[2]] - arrCoord[ez][ephi][er][icEphipzm[2]]) * PetscAbsReal(arrCoord[ez][ephi][er][icErpphip[0]] - arrCoord[ez][ephi][er][icErmphip[0]]); /* INT(dr dz) */
  }
  /* Faces perpendicular to r direction */
  else if (loc == LEFT) {

    if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
      surf = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icErmzp[2]] - arrCoord[ez][ephi][er][icErmzm[2]]) * arrCoord[ez][ephi][er][icBrm[0]]; /* INT(r dphi dz) */
    } else {
      surf = PetscAbsReal(arrCoord[ez][ephi][er][icErmzp[2]] - arrCoord[ez][ephi][er][icErmzm[2]]) * PetscAbsReal(arrCoord[ez][ephi][er][icErmphip[1]] - arrCoord[ez][ephi][er][icErmphim[1]]) * arrCoord[ez][ephi][er][icBrm[0]]; /* INT(r dphi dz) */
    }
    /* DEBUG PRINT */
    /*PetscPrintf(PETSC_COMM_WORLD,"Phip(%d,%d,%d) = %E\n",er,ephi,ez,(double)arrCoord[ez][ephi][er][icErmphip[1]]);
    PetscPrintf(PETSC_COMM_WORLD,"Phim(%d,%d,%d) = %E\n",er,ephi,ez,(double)arrCoord[ez][ephi][er][icErmphim[1]]);*/
  } else if (loc == RIGHT) {
    if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
      surf = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icErpzp[2]] - arrCoord[ez][ephi][er][icErpzm[2]]) * arrCoord[ez][ephi][er][icBrp[0]]; /* INT(r dphi dz) */
    } else {
      surf = PetscAbsReal(arrCoord[ez][ephi][er][icErpzp[2]] - arrCoord[ez][ephi][er][icErpzm[2]]) * PetscAbsReal(arrCoord[ez][ephi][er][icErpphip[1]] - arrCoord[ez][ephi][er][icErpphim[1]]) * arrCoord[ez][ephi][er][icBrp[0]]; /* INT(r dphi dz) */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in surface function");
  }
  return surf;
}

#line 555

#line 605

#line 1214

#line 1306

PetscErrorCode SaveSolution(TS ts, Vec X, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("SaveSolution",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

    User           *user = (User*)ptr;
    DM             da;
    DM             dmFr, dmFphi,dmFz, daFr, daFphi,daFz;

    PetscInt       startr,startphi,startz,nr,nphi,nz;

    Vec            vecFr, vecFphi, vecFz, F_r2, F_phi2, F_z2;

    Vec            XLocal;
    PetscInt       N[3],er,ephi,ez;

#line 1349
    PetscCall(TSGetDM(ts,&da));

    DMStagCreateCompatibleDMStag(da,0,0,1,0,&dmFr); /* 1 dof per face */
    PetscCall(DMSetUp(dmFr));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmFr,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmFr,&F_r2));

    DMStagCreateCompatibleDMStag(da,0,0,1,0,&dmFphi); /* 1 dof per face */
    PetscCall(DMSetUp(dmFphi));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmFphi,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmFphi,&F_phi2));

    DMStagCreateCompatibleDMStag(da,0,0,1,0,&dmFz); /* 1 dof per face */
    PetscCall(DMSetUp(dmFz));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmFz,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmFz,&F_z2));




    PetscCall(DMGetLocalVector(da, & XLocal));
    PetscCall(DMGlobalToLocalBegin(da, X, INSERT_VALUES, XLocal));
    PetscCall(DMGlobalToLocalEnd(da, X, INSERT_VALUES, XLocal));



    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying F_r values\n"));
    PetscCall(DMStagGetCorners(dmFr,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmFr,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[1];
                PetscScalar   valFrom[1];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = LEFT;    from[0].c = 0;
                PetscCall(DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmFr,F_r2,1,from,valFrom,INSERT_VALUES));
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = RIGHT;    from[0].c = 0;
                    PetscCall(DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmFr,F_r2,1,from,valFrom,INSERT_VALUES));
                }
            }
        }
    }
    PetscCall(VecAssemblyBegin(F_r2));
    PetscCall(VecAssemblyEnd(F_r2));

    DMStagVecSplitToDMDA(dmFr,F_r2,LEFT,-1,&daFr,&vecFr); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecFr,"rFace_center_values"));


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying F_phi values\n"));
    PetscCall(DMStagGetCorners(dmFphi,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmFphi,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[1];
                PetscScalar   valFrom[1];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN;    from[0].c = 0;
                PetscCall(DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmFphi,F_phi2,1,from,valFrom,INSERT_VALUES));
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP; from[0].c = 0;
                    PetscCall(DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmFphi,F_phi2,1,from,valFrom,INSERT_VALUES));
                }
            }
        }
    }
    PetscCall(VecAssemblyBegin(F_phi2));
    PetscCall(VecAssemblyEnd(F_phi2));

    DMStagVecSplitToDMDA(dmFphi,F_phi2,DOWN,-1,&daFphi,&vecFphi); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecFphi,"phiFace_center_values"));


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying F_z values\n"));
    PetscCall(DMStagGetCorners(dmFz,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmFz,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK;    from[0].c = 0;
                PetscCall(DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmFz,F_z2,1,from,valFrom,INSERT_VALUES));
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT; from[0].c = 0;
                    PetscCall(DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmFz,F_z2,1,from,valFrom,INSERT_VALUES));
                }
            }
        }
    }
    PetscCall(VecAssemblyBegin(F_z2));
    PetscCall(VecAssemblyEnd(F_z2));

    DMStagVecSplitToDMDA(dmFz,F_z2,BACK,-1,&daFz,&vecFz); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecFz,"zFace_center_values"));

    PetscViewer viewerD;
    char filename[PETSC_MAX_PATH_LEN];
    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecBr.m", user->input_folder));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before opening %s file\n", filename));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewerD));
    //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecBr.m",FILE_MODE_WRITE,&viewerD);
    //PetscViewerPushFormat(viewerD,PETSC_VIEWER_BINARY_MATLAB);
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecFr,viewerD));
    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecBphi.m", user->input_folder));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before opening %s file\n", filename));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewerD));
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecFphi,viewerD));
    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecBz.m", user->input_folder));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before opening %s file\n", filename));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewerD));
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecFz,viewerD));
    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscViewerDestroy(&viewerD));

    PetscCall(DMDestroy(&dmFr));
    PetscCall(DMDestroy(&dmFphi));
    PetscCall(DMDestroy(&dmFz));

    PetscCall(DMDestroy(&daFr));
    PetscCall(DMDestroy(&daFphi));
    PetscCall(DMDestroy(&daFz));

    PetscCall(VecDestroy(&vecFr));
    PetscCall(VecDestroy(&vecFphi));
    PetscCall(VecDestroy(&vecFz));

    PetscCall(VecDestroy(&F_r2));
    PetscCall(VecDestroy(&F_phi2));
    PetscCall(VecDestroy(&F_z2));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode SaveCoordinates(TS ts, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("SaveCoordinates",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

    User           *user = (User*)ptr;
    DM             da;
    DM             daC, dmC;
    DM             dmFr, dmFphi,dmFz, daFr, daFphi,daFz;
    DM             dmEr, dmEphi,dmEz, daEr, daEphi,daEz;

    PetscInt       startr,startphi,startz,nr,nphi,nz;

    Vec            vecC, C, C2, vecFr, vecFphi, vecFz, F_r, F_r2, F_phi, F_phi2, F_z, F_z2, vecEr, vecEphi, vecEz, E_r, E_r2, E_phi, E_phi2, E_z, E_z2;

    PetscInt       N[3],er,ephi,ez;

#line 1542
    PetscCall(TSGetDM(ts,&da));

    DMStagCreateCompatibleDMStag(da,0,0,3,0,&dmFr); /* 3 dofs per face */
    PetscCall(DMSetUp(dmFr));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmFr,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmFr,&F_r2));

    DMStagCreateCompatibleDMStag(da,0,0,3,0,&dmFphi); /* 3 dofs per face */
    PetscCall(DMSetUp(dmFphi));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmFphi,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmFphi,&F_phi2));

    DMStagCreateCompatibleDMStag(da,0,0,3,0,&dmFz); /* 3 dofs per face */
    PetscCall(DMSetUp(dmFz));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmFz,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmFz,&F_z2));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before creating compatible DMStag for E_r\n"));

    DMStagCreateCompatibleDMStag(da,0,3,0,0,&dmEr); /* 3 dofs per edge */
    PetscCall(DMSetUp(dmEr));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmEr,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmEr,&E_r2));

    DMStagCreateCompatibleDMStag(da,0,3,0,0,&dmEphi); /* 3 dofs per edge */
    PetscCall(DMSetUp(dmEphi));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmEphi,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmEphi,&E_phi2));

    DMStagCreateCompatibleDMStag(da,0,3,0,0,&dmEz); /* 3 dofs per edge */
    PetscCall(DMSetUp(dmEz));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmEz,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmEz,&E_z2));

    DMStagCreateCompatibleDMStag(da,0,0,0,3,&dmC); /* 3 dofs per cell */
    PetscCall(DMSetUp(dmC));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmC,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmC,&C2));
    //STOPPED HERE


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying E_r coordinates \n"));
    PetscCall(DMGetCoordinatesLocal(dmEr, &E_r));
    PetscCall(DMStagGetCorners(dmEr,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmEr,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK_DOWN;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK_DOWN;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK_DOWN;    from[2].c = 2;
                PetscCall(DMStagVecGetValuesStencil(dmEr,E_r,3,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmEr,E_r2,3,from,valFrom,INSERT_VALUES));
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK_UP;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK_UP;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK_UP;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmEr,E_r,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmEr,E_r2,1,from,valFrom,INSERT_VALUES));
                }
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT_DOWN;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT_DOWN;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT_DOWN;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmEr,E_r,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmEr,E_r2,3,from,valFrom,INSERT_VALUES));
                }
                if(ephi == N[1]-1 && ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT_UP;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT_UP;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT_UP;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmEr,E_r,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmEr,E_r2,3,from,valFrom,INSERT_VALUES));
                }
            }
        }
    }
    PetscCall(VecAssemblyBegin(E_r2));
    PetscCall(VecAssemblyEnd(E_r2));

    DMStagVecSplitToDMDA(dmEr,E_r2,BACK_DOWN,-3,&daEr,&vecEr); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecEr,"rEdge_center_coordinates"));


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying E_phi coordinates\n"));
    PetscCall(DMGetCoordinatesLocal(dmEphi, &E_phi));
    PetscCall(DMStagGetCorners(dmEphi,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmEphi,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK_LEFT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK_LEFT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK_LEFT;    from[2].c = 2;
                PetscCall(DMStagVecGetValuesStencil(dmEphi,E_phi,3,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmEphi,E_phi2,3,from,valFrom,INSERT_VALUES));
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK_RIGHT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK_RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK_RIGHT;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmEphi,E_phi,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmEphi,E_phi2,3,from,valFrom,INSERT_VALUES));
                }
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT_LEFT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT_LEFT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT_LEFT;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmEphi,E_phi,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmEphi,E_phi2,3,from,valFrom,INSERT_VALUES));
                }
                if(er == N[0]-1 && ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT_RIGHT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT_RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT_RIGHT;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmEphi,E_phi,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmEphi,E_phi2,3,from,valFrom,INSERT_VALUES));
                }
            }
        }
    }
    PetscCall(VecAssemblyBegin(E_phi2));
    PetscCall(VecAssemblyEnd(E_phi2));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before calling DMStagVecSplitToDMDA for E_phi coordinates\n"));
    DMStagVecSplitToDMDA(dmEphi,E_phi2,BACK_LEFT,-3,&daEphi,&vecEphi); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecEphi,"phiEdge_center_coordinates"));


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying E_z coordinates\n"));
    PetscCall(DMGetCoordinatesLocal(dmEz, &E_z));
    PetscCall(DMStagGetCorners(dmEz,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmEz,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN_LEFT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = DOWN_LEFT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = DOWN_LEFT;    from[2].c = 2;
                PetscCall(DMStagVecGetValuesStencil(dmEz,E_z,3,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmEz,E_z2,3,from,valFrom,INSERT_VALUES));
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN_RIGHT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = DOWN_RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = DOWN_RIGHT;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmEz,E_z,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmEz,E_z2,3,from,valFrom,INSERT_VALUES));
                }
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP_LEFT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = UP_LEFT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = UP_LEFT;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmEz,E_z,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmEz,E_z2,3,from,valFrom,INSERT_VALUES));
                }
                if(er == N[0]-1 && ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP_RIGHT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = UP_RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = UP_RIGHT;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmEz,E_z,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmEz,E_z2,3,from,valFrom,INSERT_VALUES));
                }
            }
        }
    }
    PetscCall(VecAssemblyBegin(E_z2));
    PetscCall(VecAssemblyEnd(E_z2));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before calling DMStagVecSplitToDMDA for E_z coordinates\n"));
    DMStagVecSplitToDMDA(dmEz,E_z2,DOWN_LEFT,-3,&daEz,&vecEz); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecEz,"zEdge_center_coordinates"));


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying F_r coordinates\n"));
    PetscCall(DMGetCoordinatesLocal(dmFr, &F_r));
    PetscCall(DMStagGetCorners(dmFr,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmFr,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = LEFT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = LEFT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = LEFT;    from[2].c = 2;
                PetscCall(DMStagVecGetValuesStencil(dmFr,F_r,3,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmFr,F_r2,3,from,valFrom,INSERT_VALUES));
                if(er == N[0]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = RIGHT;    from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = RIGHT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = RIGHT;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmFr,F_r,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmFr,F_r2,3,from,valFrom,INSERT_VALUES));
                }
            }
        }
    }
    PetscCall(VecAssemblyBegin(F_r2));
    PetscCall(VecAssemblyEnd(F_r2));

    DMStagVecSplitToDMDA(dmFr,F_r2,LEFT,-3,&daFr,&vecFr); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecFr,"rFace_center_coordinates"));


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying F_phi coordinates\n"));
    PetscCall(DMGetCoordinatesLocal(dmFphi, &F_phi));
    PetscCall(DMStagGetCorners(dmFphi,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmFphi,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = DOWN;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = DOWN;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = DOWN;    from[2].c = 2;
                PetscCall(DMStagVecGetValuesStencil(dmFphi,F_phi,3,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmFphi,F_phi2,3,from,valFrom,INSERT_VALUES));
                if(ephi == N[1]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = UP; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = UP;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = UP;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmFphi,F_phi,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmFphi,F_phi2,3,from,valFrom,INSERT_VALUES));
                }
            }
        }
    }
    PetscCall(VecAssemblyBegin(F_phi2));
    PetscCall(VecAssemblyEnd(F_phi2));

    DMStagVecSplitToDMDA(dmFphi,F_phi2,DOWN,-3,&daFphi,&vecFphi); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecFphi,"phiFace_center_coordinates"));


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying F_z coordinates\n"));
    PetscCall(DMGetCoordinatesLocal(dmFz, &F_z));
    PetscCall(DMStagGetCorners(dmFz,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmFz,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = BACK;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = BACK;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = BACK;    from[2].c = 2;
                PetscCall(DMStagVecGetValuesStencil(dmFz,F_z,3,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmFz,F_z2,3,from,valFrom,INSERT_VALUES));
                if(ez == N[2]-1){
                    from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = FRONT; from[0].c = 0;
                    from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = FRONT;    from[1].c = 1;
                    from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = FRONT;    from[2].c = 2;
                    PetscCall(DMStagVecGetValuesStencil(dmFz,F_z,3,from,valFrom));
                    PetscCall(DMStagVecSetValuesStencil(dmFz,F_z2,3,from,valFrom,INSERT_VALUES));
                }
            }
        }
    }
    PetscCall(VecAssemblyBegin(F_z2));
    PetscCall(VecAssemblyEnd(F_z2));

    DMStagVecSplitToDMDA(dmFz,F_z2,BACK,-3,&daFz,&vecFz); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecFz,"zFace_center_coordinates"));


    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before copying cell coordinates\n"));
    PetscCall(DMGetCoordinatesLocal(dmC, &C));
    PetscCall(DMStagGetCorners(dmC,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[3];
                PetscScalar   valFrom[3];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = ELEMENT;    from[0].c = 0;
                from[1].i = er; from[1].j = ephi; from[1].k = ez; from[1].loc = ELEMENT;    from[1].c = 1;
                from[2].i = er; from[2].j = ephi; from[2].k = ez; from[2].loc = ELEMENT;    from[2].c = 2;
                PetscCall(DMStagVecGetValuesStencil(dmC,C,3,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmC,C2,3,from,valFrom,INSERT_VALUES));
            }
        }
    }
    PetscCall(VecAssemblyBegin(C2));
    PetscCall(VecAssemblyEnd(C2));

    DMStagVecSplitToDMDA(dmC,C2,ELEMENT,-3,&daC,&vecC); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecC,"Cell_center_coordinates"));

    char filename[PETSC_MAX_PATH_LEN];
    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecC.m", user->input_folder));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Before opening %s file\n", filename));
    PetscViewer viewerD;
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewerD));
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecC,viewerD));

    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecFr.m", user->input_folder));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewerD));
    //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecFr.m",FILE_MODE_WRITE,&viewerD);
    //PetscViewerPushFormat(viewerD,PETSC_VIEWER_BINARY_MATLAB);
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecFr,viewerD));
    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecFphi.m", user->input_folder));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewerD));
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecFphi,viewerD));
    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecFz.m", user->input_folder));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewerD));
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecFz,viewerD));
    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecEr.m", user->input_folder));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewerD));
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecEr,viewerD));
    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecEphi.m", user->input_folder));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewerD));
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecEphi,viewerD));
    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscSNPrintf(filename, PETSC_MAX_PATH_LEN, "%s/vecEz.m", user->input_folder));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewerD));
    PetscCall(PetscViewerPushFormat(viewerD,PETSC_VIEWER_ASCII_MATLAB));
    PetscCall(VecView(vecEz,viewerD));
    PetscCall(PetscViewerPopFormat(viewerD));

    PetscCall(PetscBarrier((PetscObject)viewerD));
    PetscCall(PetscViewerDestroy(&viewerD));

    PetscCall(DMDestroy(&dmFr));
    PetscCall(DMDestroy(&dmFphi));
    PetscCall(DMDestroy(&dmFz));
    PetscCall(DMDestroy(&dmEr));
    PetscCall(DMDestroy(&dmEphi));
    PetscCall(DMDestroy(&dmEz));
    PetscCall(DMDestroy(&dmC));

    PetscCall(DMDestroy(&daFr));
    PetscCall(DMDestroy(&daFphi));
    PetscCall(DMDestroy(&daFz));
    PetscCall(DMDestroy(&daEr));
    PetscCall(DMDestroy(&daEphi));
    PetscCall(DMDestroy(&daEz));
    PetscCall(DMDestroy(&daC));

    PetscCall(VecDestroy(&vecC));
    PetscCall(VecDestroy(&vecEr));
    PetscCall(VecDestroy(&vecEphi));
    PetscCall(VecDestroy(&vecEz));
    PetscCall(VecDestroy(&vecFr));
    PetscCall(VecDestroy(&vecFphi));
    PetscCall(VecDestroy(&vecFz));

    PetscCall(VecDestroy(&C2));
    PetscCall(VecDestroy(&E_r2));
    PetscCall(VecDestroy(&E_phi2));
    PetscCall(VecDestroy(&E_z2));
    PetscCall(VecDestroy(&F_r2));
    PetscCall(VecDestroy(&F_phi2));
    PetscCall(VecDestroy(&F_z2));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

#line 2264

PetscErrorCode CellToVertexProjectionScalar(TS ts, Vec C, Vec V, void *ptr)
{
  PetscFunctionBeginUser;

  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("CellToVertexProjectionScalar",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

    User           *user = (User*)ptr;
    DM             da;
    Vec            CLocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;



    PetscCall(VecZeroEntries(V));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&CLocal));
    PetscCall(DMGlobalToLocal(da,C,INSERT_VALUES,CLocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[8], to[1];
          PetscScalar valFrom[8], valTo[1];

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = ELEMENT;
          from[0].c = 0;
          if(er>0 && ez>0 && (ephi>0 || user->phibtype)){
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
          PetscCall(DMStagVecGetValuesStencil(da, CLocal, 8, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 0;
          valTo[0] = 0.125 * (valFrom[0] + valFrom[1] + valFrom[2] + valFrom[3] + valFrom[4] + valFrom[5] + valFrom[6] + valFrom[7]);
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }
        }
      }
    }
    PetscCall(VecAssemblyBegin(V));
    PetscCall(VecAssemblyEnd(V));
    PetscCall(DMRestoreLocalVector(da,&CLocal));

    PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

    PetscFunctionReturn(PETSC_SUCCESS);
}

#line 2543

#line 2727

#line 1480

PetscErrorCode VertexToEdgeReconstruction(TS ts, Vec V, Vec E, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("VertexToEdgeReconstruction",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));


    DM             da;
    Vec            VLocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;
#line 1499
    PetscCall(VecZeroEntries(E));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&VLocal));
    PetscCall(DMGlobalToLocal(da,V,INSERT_VALUES,VLocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[24], to[12];
          PetscScalar valFrom[24], valTo[12];

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
          PetscCall(DMStagVecGetValuesStencil(da, VLocal, 24, from, valFrom));
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
          valTo[1] = 0.5 * (valFrom[9] + valFrom[10]);
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
          valTo[3] = 0.5 * (valFrom[8] + valFrom[11]);
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
          valTo[5] = 0.5 * (valFrom[13] + valFrom[14]);
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
          valTo[7] = 0.5 * (valFrom[12] + valFrom[15]);
          to[8].i = er;
          to[8].j = ephi;
          to[8].k = ez;
          to[8].loc = UP_RIGHT;
          to[8].c = 0;
          valTo[8] = 0.5 * (valFrom[16] + valFrom[20]);
          to[9].i = er;
          to[9].j = ephi;
          to[9].k = ez;
          to[9].loc = UP_LEFT;
          to[9].c = 0;
          valTo[9] = 0.5 * (valFrom[17] + valFrom[21]);
          to[10].i = er;
          to[10].j = ephi;
          to[10].k = ez;
          to[10].loc = DOWN_RIGHT;
          to[10].c = 0;
          valTo[10] = 0.5 * (valFrom[19] + valFrom[23]);
          to[11].i = er;
          to[11].j = ephi;
          to[11].k = ez;
          to[11].loc = DOWN_LEFT;
          to[11].c = 0;
          valTo[11] = 0.5 * (valFrom[18] + valFrom[22]);
          PetscCall(DMStagVecSetValuesStencil(da, E, 12, to, valTo, INSERT_VALUES));

        }
      }
    }
    PetscCall(VecAssemblyBegin(E));
    PetscCall(VecAssemblyEnd(E));
    PetscCall(DMRestoreLocalVector(da,&VLocal));

    PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));
    PetscFunctionReturn(PETSC_SUCCESS);
}

#line 3416

PetscErrorCode VertexToFaceReconstruction(TS ts, Vec V, Vec F, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("VertexToFaceReconstruction",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));


    DM             da;
    Vec            VLocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;



    PetscCall(VecZeroEntries(F));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&VLocal));
    PetscCall(DMGlobalToLocal(da,V,INSERT_VALUES,VLocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[24], to[6];
          PetscScalar valFrom[24], valTo[6];

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
          PetscCall(DMStagVecGetValuesStencil(da, VLocal, 24, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = UP;
          to[0].c = 0;
          valTo[0] = 0.25 * (valFrom[8] + valFrom[9] + valFrom[12] + valFrom[13]);
          to[1].i = er;
          to[1].j = ephi;
          to[1].k = ez;
          to[1].loc = LEFT;
          to[1].c = 0;
          valTo[1] = 0.25 * (valFrom[1] + valFrom[2] + valFrom[5] + valFrom[6]);
          to[2].i = er;
          to[2].j = ephi;
          to[2].k = ez;
          to[2].loc = DOWN;
          to[2].c = 0;
          valTo[2] = 0.25 * (valFrom[10] + valFrom[11] + valFrom[14] + valFrom[15]);
          to[3].i = er;
          to[3].j = ephi;
          to[3].k = ez;
          to[3].loc = RIGHT;
          to[3].c = 0;
          valTo[3] = 0.25 * (valFrom[0] + valFrom[3] + valFrom[4] + valFrom[7]);
          to[4].i = er;
          to[4].j = ephi;
          to[4].k = ez;
          to[4].loc = BACK;
          to[4].c = 0;
          valTo[4] = 0.25 * (valFrom[20] + valFrom[21] + valFrom[22] + valFrom[23]);
          to[5].i = er;
          to[5].j = ephi;
          to[5].k = ez;
          to[5].loc = FRONT;
          to[5].c = 0;
          valTo[5] = 0.25 * (valFrom[16] + valFrom[17] + valFrom[18] + valFrom[19]);
          PetscCall(DMStagVecSetValuesStencil(da, F, 6, to, valTo, INSERT_VALUES));

        }
      }
    }
    PetscCall(VecAssemblyBegin(F));
    PetscCall(VecAssemblyEnd(F));
    PetscCall(DMRestoreLocalVector(da,&VLocal));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));
    PetscFunctionReturn(PETSC_SUCCESS);
}

#line 3864

PetscErrorCode EdgeToCellReconstruction_r(TS ts, Vec E, Vec C, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("EdgeToCellReconstruction_r",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));


    DM             da;
    Vec            ELocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;



    PetscCall(VecZeroEntries(C));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&ELocal));
    PetscCall(DMGlobalToLocal(da,E,INSERT_VALUES,ELocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[4], to[1];
          PetscScalar valFrom[4], valTo[1];

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
          PetscCall(DMStagVecGetValuesStencil(da, ELocal, 4, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = ELEMENT;
          to[0].c = 0;
          valTo[0] = 0.25 * (valFrom[0] + valFrom[1] + valFrom[2] + valFrom[3]);
          PetscCall(DMStagVecSetValuesStencil(da, C, 1, to, valTo, INSERT_VALUES));

        }
      }
    }
    PetscCall(VecAssemblyBegin(C));
    PetscCall(VecAssemblyEnd(C));
    PetscCall(DMRestoreLocalVector(da,&ELocal));

    PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EdgeToCellReconstruction_phi(TS ts, Vec E, Vec C, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("EdgeToCellReconstruction_phi",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));


    DM             da;
    Vec            ELocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;



    PetscCall(VecZeroEntries(C));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&ELocal));
    PetscCall(DMGlobalToLocal(da,E,INSERT_VALUES,ELocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[4], to[1];
          PetscScalar valFrom[4], valTo[1];

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = FRONT_LEFT;
          from[0].c = 0;
          from[1].i = er;
          from[1].j = ephi;
          from[1].k = ez;
          from[1].loc = BACK_LEFT;
          from[1].c = 0;
          from[2].i = er;
          from[2].j = ephi;
          from[2].k = ez;
          from[2].loc = FRONT_RIGHT;
          from[2].c = 0;
          from[3].i = er;
          from[3].j = ephi;
          from[3].k = ez;
          from[3].loc = BACK_RIGHT;
          from[3].c = 0;
          PetscCall(DMStagVecGetValuesStencil(da, ELocal, 4, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = ELEMENT;
          to[0].c = 0;
          valTo[0] = 0.25 * (valFrom[0] + valFrom[1] + valFrom[2] + valFrom[3]);
          PetscCall(DMStagVecSetValuesStencil(da, C, 1, to, valTo, INSERT_VALUES));

        }
      }
    }
    PetscCall(VecAssemblyBegin(C));
    PetscCall(VecAssemblyEnd(C));
    PetscCall(DMRestoreLocalVector(da,&ELocal));

    PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode EdgeToCellReconstruction_z(TS ts, Vec E, Vec C, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("EdgeToCellReconstruction_z",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));


    DM             da;
    Vec            ELocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;



    PetscCall(VecZeroEntries(C));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&ELocal));
    PetscCall(DMGlobalToLocal(da,E,INSERT_VALUES,ELocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[4], to[1];
          PetscScalar valFrom[4], valTo[1];

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = UP_LEFT;
          from[0].c = 0;
          from[1].i = er;
          from[1].j = ephi;
          from[1].k = ez;
          from[1].loc = UP_RIGHT;
          from[1].c = 0;
          from[2].i = er;
          from[2].j = ephi;
          from[2].k = ez;
          from[2].loc = DOWN_LEFT;
          from[2].c = 0;
          from[3].i = er;
          from[3].j = ephi;
          from[3].k = ez;
          from[3].loc = DOWN_RIGHT;
          from[3].c = 0;
          PetscCall(DMStagVecGetValuesStencil(da, ELocal, 4, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = ELEMENT;
          to[0].c = 0;
          valTo[0] = 0.25 * (valFrom[0] + valFrom[1] + valFrom[2] + valFrom[3]);
          PetscCall(DMStagVecSetValuesStencil(da, C, 1, to, valTo, INSERT_VALUES));

        }
      }
    }
    PetscCall(VecAssemblyBegin(C));
    PetscCall(VecAssemblyEnd(C));
    PetscCall(DMRestoreLocalVector(da,&ELocal));

    PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

    PetscFunctionReturn(PETSC_SUCCESS);
}

#line 4219

#line 4316

PetscErrorCode FaceToVertexProjection(TS ts, Vec F, Vec V, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("FaceToVertexProjection",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

    User           *user = (User*)ptr;
    DM             da;
    Vec            FLocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;



    PetscCall(VecZeroEntries(V));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&FLocal));
    PetscCall(DMGlobalToLocal(da,F,INSERT_VALUES,FLocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[4], to[1];
          PetscScalar valFrom[4], valTo[1];
          PetscInt nEntries;

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = DOWN;
          from[0].c = 0;
          nEntries = 1;
          if(er!=0 && ez!=0){
            from[1].i = er-1;
            from[1].j = ephi;
            from[1].k = ez-1;
            from[1].loc = DOWN;
            from[1].c = 0;
            from[2].i = er;
            from[2].j = ephi;
            from[2].k = ez-1;
            from[2].loc = DOWN;
            from[2].c = 0;
            from[3].i = er-1;
            from[3].j = ephi;
            from[3].k = ez;
            from[3].loc = DOWN;
            from[3].c = 0;
            nEntries = 4;
          }
          else if(ez!=0){
            from[1].i = er;
            from[1].j = ephi;
            from[1].k = ez-1;
            from[1].loc = DOWN;
            from[1].c = 0;
            nEntries = 2;
          }
          else if(er!=0){
            from[1].i = er-1;
            from[1].j = ephi;
            from[1].k = ez;
            from[1].loc = DOWN;
            from[1].c = 0;
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 1;
          if(nEntries==1){
            valTo[0] = 0.25 * valFrom[0];
          }
          else if(nEntries==2){
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
          }
          else{
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = BACK;
          from[0].c = 0;
          nEntries = 1;
          if((user -> phibtype && !(er == 0)) || ((!(er == 0 || ephi == 0)) && !(user -> phibtype))){
            from[1].i = er-1;
            from[1].j = ephi-1;
            from[1].k = ez;
            from[1].loc = BACK;
            from[1].c = 0;
            from[2].i = er;
            from[2].j = ephi-1;
            from[2].k = ez;
            from[2].loc = BACK;
            from[2].c = 0;
            from[3].i = er-1;
            from[3].j = ephi;
            from[3].k = ez;
            from[3].loc = BACK;
            from[3].c = 0;
            nEntries = 4;
          }
          else if(er!=0){
            from[1].i = er-1;
            from[1].j = ephi;
            from[1].k = ez;
            from[1].loc = BACK;
            from[1].c = 0;
            nEntries = 2;
          }
          else if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
            from[1].i = er;
            from[1].j = ephi-1;
            from[1].k = ez;
            from[1].loc = BACK;
            from[1].c = 0;
            nEntries = 2;
          }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 2;
          if(nEntries==1){
            valTo[0] = 0.25 * valFrom[0];
          }
          else if(nEntries==2){
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
          }
          else{
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = LEFT;
          from[0].c = 0;
          nEntries = 1;
          if ((user -> phibtype && !(ez == 0)) || ((!(ez == 0 || ephi == 0)) && !(user -> phibtype))) {
            from[1].i = er;
            from[1].j = ephi-1;
            from[1].k = ez-1;
            from[1].loc = LEFT;
            from[1].c = 0;
            from[2].i = er;
            from[2].j = ephi-1;
            from[2].k = ez;
            from[2].loc = LEFT;
            from[2].c = 0;
            from[3].i = er;
            from[3].j = ephi;
            from[3].k = ez-1;
            from[3].loc = LEFT;
            from[3].c = 0;
            nEntries = 4;
          }
          else if(ez!=0){
            from[1].i = er;
            from[1].j = ephi;
            from[1].k = ez-1;
            from[1].loc = LEFT;
            from[1].c = 0;
            nEntries = 2;
          }
          else if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
            from[1].i = er;
            from[1].j = ephi-1;
            from[1].k = ez;
            from[1].loc = LEFT;
            from[1].c = 0;
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 0;
          if(nEntries==1){
            valTo[0] = 0.25 * valFrom[0];
          }
          else if(nEntries==2){
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
          }
          else{
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          if(er==N[0]-1){
          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = DOWN;
          from[0].c = 0;
          nEntries = 1;
          if(ez!=0){
            from[1].i = er;
            from[1].j = ephi;
            from[1].k = ez-1;
            from[1].loc = DOWN;
            from[1].c = 0;
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_RIGHT;
          to[0].c = 1;
          if(nEntries==1){
            valTo[0] = 0.25 * valFrom[0];
          }
          else if(nEntries==2){
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
          }
          else{
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = BACK;
          from[0].c = 0;
          nEntries = 1;
          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
            from[1].i = er;
            from[1].j = ephi-1;
            from[1].k = ez;
            from[1].loc = BACK;
            from[1].c = 0;
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_RIGHT;
          to[0].c = 2;
          if(nEntries==1){
            valTo[0] = 0.25 * valFrom[0];
          }
          else if(nEntries==2){
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
          }
          else{
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = RIGHT;
          from[0].c = 0;
          nEntries = 1;
          if ((user -> phibtype && !(ez == 0)) || ((!(ez == 0 || ephi == 0)) && !(user -> phibtype))) {
            from[1].i = er;
            from[1].j = ephi-1;
            from[1].k = ez-1;
            from[1].loc = RIGHT;
            from[1].c = 0;
            from[2].i = er;
            from[2].j = ephi-1;
            from[2].k = ez;
            from[2].loc = RIGHT;
            from[2].c = 0;
            from[3].i = er;
            from[3].j = ephi;
            from[3].k = ez-1;
            from[3].loc = RIGHT;
            from[3].c = 0;
            nEntries = 4;
          }
          else if(ez!=0){
            from[1].i = er;
            from[1].j = ephi;
            from[1].k = ez-1;
            from[1].loc = RIGHT;
            from[1].c = 0;
            nEntries = 2;
          }
          else if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
            from[1].i = er;
            from[1].j = ephi-1;
            from[1].k = ez;
            from[1].loc = RIGHT;
            from[1].c = 0;
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_RIGHT;
          to[0].c = 0;
          if(nEntries==1){
            valTo[0] = 0.25 * valFrom[0];
          }
          else if(nEntries==2){
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
          }
          else{
            valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(ez==N[2]-1){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = DOWN;
            from[0].c = 0;
            nEntries = 1;
            if(er!=0){
            from[1].i = er-1;
            from[1].j = ephi;
            from[1].k = ez;
            from[1].loc = DOWN;
            from[1].c = 0;
            nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_DOWN_LEFT;
            to[0].c = 1;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = LEFT;
            from[0].c = 0;
            nEntries = 1;
            if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
            from[1].i = er;
            from[1].j = ephi-1;
            from[1].k = ez;
            from[1].loc = LEFT;
            from[1].c = 0;
            nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_DOWN_LEFT;
            to[0].c = 0;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = FRONT;
            from[0].c = 0;
            nEntries = 1;
            if ((user -> phibtype && !(er == 0)) || ((!(er == 0 || ephi == 0)) && !(user -> phibtype))) {
              from[1].i = er-1;
              from[1].j = ephi-1;
              from[1].k = ez;
              from[1].loc = FRONT;
              from[1].c = 0;
              from[2].i = er;
              from[2].j = ephi-1;
              from[2].k = ez;
              from[2].loc = FRONT;
              from[2].c = 0;
              from[3].i = er-1;
              from[3].j = ephi;
              from[3].k = ez;
              from[3].loc = FRONT;
              from[3].c = 0;
              nEntries = 4;
            }
            else if(er!=0){
              from[1].i = er-1;
              from[1].j = ephi;
              from[1].k = ez;
              from[1].loc = FRONT;
              from[1].c = 0;
              nEntries = 2;
            }
            else if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
              from[1].i = er;
              from[1].j = ephi-1;
              from[1].k = ez;
              from[1].loc = FRONT;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_DOWN_LEFT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(ephi==N[1]-1 && !user->phibtype){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = UP;
            from[0].c = 0;
            nEntries = 1;
            if(er!=0 && ez!=0){
              from[1].i = er-1;
              from[1].j = ephi;
              from[1].k = ez-1;
              from[1].loc = UP;
              from[1].c = 0;
              from[2].i = er;
              from[2].j = ephi;
              from[2].k = ez-1;
              from[2].loc = UP;
              from[2].c = 0;
              from[3].i = er-1;
              from[3].j = ephi;
              from[3].k = ez;
              from[3].loc = UP;
              from[3].c = 0;
              nEntries = 4;
            }
            else if(ez!=0){
              from[1].i = er;
              from[1].j = ephi;
              from[1].k = ez-1;
              from[1].loc = UP;
              from[1].c = 0;
              nEntries = 2;
            }
            else if(er!=0) {
              from[1].i = er-1;
              from[1].j = ephi;
              from[1].k = ez;
              from[1].loc = UP;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_LEFT;
            to[0].c = 1;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = LEFT;
            from[0].c = 0;
            nEntries = 1;
            if(ez!=0){
              from[1].i = er;
              from[1].j = ephi;
              from[1].k = ez-1;
              from[1].loc = LEFT;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_LEFT;
            to[0].c = 0;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = BACK;
            from[0].c = 0;
            nEntries = 1;
            if(er!=0){
              from[1].i = er-1;
              from[1].j = ephi;
              from[1].k = ez;
              from[1].loc = BACK;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_LEFT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(er==N[0]-1 && ephi==N[1]-1 && !user->phibtype){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = UP;
            from[0].c = 0;
            nEntries = 1;
            if(ez!=0){
              from[1].i = er;
              from[1].j = ephi;
              from[1].k = ez-1;
              from[1].loc = UP;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_RIGHT;
            to[0].c = 1;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = RIGHT;
            from[0].c = 0;
            nEntries = 1;
            if(ez!=0){
              from[1].i = er;
              from[1].j = ephi;
              from[1].k = ez-1;
              from[1].loc = RIGHT;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_RIGHT;
            to[0].c = 0;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = BACK;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_RIGHT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = UP;
            from[0].c = 0;
            nEntries = 1;
            if(er!=0){
              from[1].i = er-1;
              from[1].j = ephi;
              from[1].k = ez;
              from[1].loc = UP;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_LEFT;
            to[0].c = 1;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = FRONT;
            from[0].c = 0;
            nEntries = 1;
            if(er!=0){
              from[1].i = er-1;
              from[1].j = ephi;
              from[1].k = ez;
              from[1].loc = FRONT;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_LEFT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = LEFT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_LEFT;
            to[0].c = 0;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(er==N[0]-1 && ez==N[2]-1){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = DOWN;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_DOWN_RIGHT;
            to[0].c = 1;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = RIGHT;
            from[0].c = 0;
            nEntries = 1;
            if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
              from[1].i = er;
              from[1].j = ephi-1;
              from[1].k = ez;
              from[1].loc = RIGHT;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_DOWN_RIGHT;
            to[0].c = 0;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = FRONT;
            from[0].c = 0;
            nEntries = 1;
            if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
              from[1].i = er;
              from[1].j = ephi-1;
              from[1].k = ez;
              from[1].loc = FRONT;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_DOWN_RIGHT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(er==N[0]-1 && ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = FRONT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_RIGHT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = RIGHT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_RIGHT;
            to[0].c = 0;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = UP;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, FLocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_RIGHT;
            to[0].c = 1;
            if(nEntries==1){
              valTo[0] = 0.25 * valFrom[0];
            }
            else if(nEntries==2){
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0]);
            }
            else{
              valTo[0] = 0.25 * (valFrom[1] + valFrom[0] + valFrom[2] + valFrom[3]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }
        }
      }
    }
    PetscCall(VecAssemblyBegin(V));
    PetscCall(VecAssemblyEnd(V));
    PetscCall(DMRestoreLocalVector(da,&FLocal));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));
    PetscFunctionReturn(PETSC_SUCCESS);
}

#line 6088

#line 6487

PetscErrorCode EdgeToVertexProjection(TS ts, Vec E, Vec V, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("EdgeToVertexProjection",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

    User           *user = (User*)ptr;
    DM             da;
    Vec            ELocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;



    PetscCall(VecZeroEntries(V));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&ELocal));
    PetscCall(DMGlobalToLocal(da,E,INSERT_VALUES,ELocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[2], to[1];
          PetscScalar valFrom[2], valTo[1];
          PetscInt nEntries;

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
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 1;
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = DOWN_LEFT;
          from[0].c = 0;
          nEntries = 1;
          if(ez!=0){
            from[1].i = er;
            from[1].j = ephi;
            from[1].k = ez-1;
            from[1].loc = DOWN_LEFT;
            from[1].c = 0;
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 2;
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = BACK_DOWN;
          from[0].c = 0;
          nEntries = 1;
          if(er!=0){
            from[1].i = er-1;
            from[1].j = ephi;
            from[1].k = ez;
            from[1].loc = BACK_DOWN;
            from[1].c = 0;
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 0;
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          if(er==N[0]-1){
          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = BACK_RIGHT;
          from[0].c = 0;
          nEntries = 1;
          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
            from[1].i = er;
            from[1].j = ephi-1;
            from[1].k = ez;
            from[1].loc = BACK_RIGHT;
            from[1].c = 0;
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_RIGHT;
          to[0].c = 1;
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = DOWN_RIGHT;
          from[0].c = 0;
          nEntries = 1;
          if(ez!=0){
            from[1].i = er;
            from[1].j = ephi;
            from[1].k = ez-1;
            from[1].loc = DOWN_RIGHT;
            from[1].c = 0;
            nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_RIGHT;
          to[0].c = 2;
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = BACK_DOWN;
          from[0].c = 0;
          nEntries = 1;
          PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_RIGHT;
          to[0].c = 0;
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(ez==N[2]-1){
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
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
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
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

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
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
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
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = DOWN_LEFT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
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
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(ephi==N[1]-1 && !user->phibtype){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = BACK_LEFT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_LEFT;
            to[0].c = 1;
            if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
            }
            else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = BACK_UP;
            from[0].c = 0;
            nEntries = 1;
            if(er!=0){
            from[1].i = er-1;
            from[1].j = ephi;
            from[1].k = ez;
            from[1].loc = BACK_UP;
            from[1].c = 0;
            nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_LEFT;
            to[0].c = 0;
            if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
            }
            else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = UP_LEFT;
            from[0].c = 0;
            nEntries = 1;
            if(ez!=0){
            from[1].i = er;
            from[1].j = ephi;
            from[1].k = ez-1;
            from[1].loc = UP_LEFT;
            from[1].c = 0;
            nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_LEFT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(er==N[0]-1 && ephi==N[1]-1 && !user->phibtype){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = BACK_RIGHT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_RIGHT;
            to[0].c = 1;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = BACK_UP;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_RIGHT;
            to[0].c = 0;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = UP_RIGHT;
            from[0].c = 0;
            nEntries = 1;
            if(ez!=0){
              from[1].i = er;
              from[1].j = ephi;
              from[1].k = ez-1;
              from[1].loc = UP_RIGHT;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_RIGHT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = FRONT_LEFT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_LEFT;
            to[0].c = 1;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = UP_LEFT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_LEFT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = FRONT_UP;
            from[0].c = 0;
            nEntries = 1;
            if(er!=0){
              from[1].i = er-1;
              from[1].j = ephi;
              from[1].k = ez;
              from[1].loc = FRONT_UP;
              from[1].c = 0;
              nEntries = 2;
            }
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_LEFT;
            to[0].c = 0;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(er==N[0]-1 && ez==N[2]-1){
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
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
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
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = FRONT_DOWN;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
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
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = DOWN_RIGHT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
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
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }

          if(er==N[0]-1 && ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = UP_RIGHT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_RIGHT;
            to[0].c = 2;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = FRONT_UP;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_RIGHT;
            to[0].c = 0;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));

            from[0].i = er;
            from[0].j = ephi;
            from[0].k = ez;
            from[0].loc = FRONT_RIGHT;
            from[0].c = 0;
            nEntries = 1;
            PetscCall(DMStagVecGetValuesStencil(da, ELocal, nEntries, from, valFrom));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_RIGHT;
            to[0].c = 1;
            if(nEntries==1){
              valTo[0] = 0.5 * valFrom[0];
            }
            else{
              valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
            }
            PetscCall(DMStagVecSetValuesStencil(da, V, 1, to, valTo, INSERT_VALUES));
          }
        }
      }
    }
    PetscCall(VecAssemblyBegin(V));
    PetscCall(VecAssemblyEnd(V));
    PetscCall(DMRestoreLocalVector(da,&ELocal));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));
    PetscFunctionReturn(PETSC_SUCCESS);
}

#line 7714

#line 7919

PetscErrorCode CellToFaceProjection(TS ts, Vec C, Vec F, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("CellToFaceProjection",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

    User           *user = (User*)ptr;
    DM             da;
    Vec            CLocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;



    PetscCall(VecZeroEntries(F));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&CLocal));
    PetscCall(DMGlobalToLocal(da,C,INSERT_VALUES,CLocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil from[2], to[1];
          PetscScalar valFrom[2], valTo[1];
          PetscInt nEntries;

          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = UP;
          to[0].c = 0;

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = ELEMENT;
          from[0].c = 0;
          nEntries = 1;

          if((ephi != N[1]-1 && !(user -> phibtype)) || user -> phibtype){
          from[1].i = er;
          from[1].j = ephi+1;
          from[1].k = ez;
          from[1].loc = ELEMENT;
          from[1].c = 0;

          nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, CLocal, nEntries, from, valFrom));
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, F, 1, to, valTo, INSERT_VALUES));

          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = LEFT;
          to[0].c = 0;

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = ELEMENT;
          from[0].c = 0;
          nEntries = 1;

          if(er!=0){
          from[1].i = er-1;
          from[1].j = ephi;
          from[1].k = ez;
          from[1].loc = ELEMENT;
          from[1].c = 0;

          nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, CLocal, nEntries, from, valFrom));
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, F, 1, to, valTo, INSERT_VALUES));

          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = DOWN;
          to[0].c = 0;

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = ELEMENT;
          from[0].c = 0;
          nEntries = 1;

          if((ephi != 0 && !(user -> phibtype)) || user -> phibtype){
          from[1].i = er;
          from[1].j = ephi-1;
          from[1].k = ez;
          from[1].loc = ELEMENT;
          from[1].c = 0;

          nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, CLocal, nEntries, from, valFrom));
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, F, 1, to, valTo, INSERT_VALUES));

          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = RIGHT;
          to[0].c = 0;

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = ELEMENT;
          from[0].c = 0;
          nEntries = 1;

          if(er != N[0]-1){
          from[1].i = er+1;
          from[1].j = ephi;
          from[1].k = ez;
          from[1].loc = ELEMENT;
          from[1].c = 0;

          nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, CLocal, nEntries, from, valFrom));
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, F, 1, to, valTo, INSERT_VALUES));

          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK;
          to[0].c = 0;

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = ELEMENT;
          from[0].c = 0;
          nEntries = 1;

          if(ez!=0){
          from[1].i = er;
          from[1].j = ephi;
          from[1].k = ez-1;
          from[1].loc = ELEMENT;
          from[1].c = 0;

          nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, CLocal, nEntries, from, valFrom));
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, F, 1, to, valTo, INSERT_VALUES));

          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = FRONT;
          to[0].c = 0;

          from[0].i = er;
          from[0].j = ephi;
          from[0].k = ez;
          from[0].loc = ELEMENT;
          from[0].c = 0;
          nEntries = 1;

          if(ez != N[2]-1){
          from[1].i = er;
          from[1].j = ephi;
          from[1].k = ez+1;
          from[1].loc = ELEMENT;
          from[1].c = 0;

          nEntries = 2;
          }
          PetscCall(DMStagVecGetValuesStencil(da, CLocal, nEntries, from, valFrom));
          if(nEntries==1){
            valTo[0] = 0.5 * valFrom[0];
          }
          else{
            valTo[0] = 0.5 * (valFrom[1] + valFrom[0]);
          }
          PetscCall(DMStagVecSetValuesStencil(da, F, 1, to, valTo, INSERT_VALUES));

        }
      }
    }
    PetscCall(VecAssemblyBegin(F));
    PetscCall(VecAssemblyEnd(F));
    PetscCall(DMRestoreLocalVector(da,&CLocal));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode VertexCrossProduct(TS ts, Vec A, Vec B, Vec C, void *ptr)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("VertexCrossProduct",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

    User           *user = (User*)ptr;
    DM             da;
    Vec            ALocal, BLocal;
    PetscInt startr, startphi, startz, nr, nphi, nz;
    PetscInt N[3], er, ephi, ez;



    PetscCall(VecZeroEntries(C));
    PetscCall(TSGetDM(ts, & da));
    PetscCall(DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]));
    PetscCall(DMStagGetCorners(da, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL));
    PetscCall(DMGetLocalVector(da,&ALocal));
    PetscCall(DMGlobalToLocal(da,A,INSERT_VALUES,ALocal));
    PetscCall(DMGetLocalVector(da,&BLocal));
    PetscCall(DMGlobalToLocal(da,B,INSERT_VALUES,BLocal));

    /* Loop over all local elements */
    for (ez = startz; ez < startz + nz; ++ez) {
      for (ephi = startphi; ephi < startphi + nphi; ++ephi) {
        for (er = startr; er < startr + nr; ++er) {
          DMStagStencil fromA[3], fromB[3], to[3];
          PetscScalar valFromA[3], valFromB[3], valTo[3];
          PetscInt nEntries=3;

          fromA[0].i = er;
          fromA[0].j = ephi;
          fromA[0].k = ez;
          fromA[0].loc = BACK_DOWN_LEFT;
          fromA[0].c = 0;
          fromA[1].i = er;
          fromA[1].j = ephi;
          fromA[1].k = ez;
          fromA[1].loc = BACK_DOWN_LEFT;
          fromA[1].c = 1;
          fromA[2].i = er;
          fromA[2].j = ephi;
          fromA[2].k = ez;
          fromA[2].loc = BACK_DOWN_LEFT;
          fromA[2].c = 2;
          PetscCall(DMStagVecGetValuesStencil(da, ALocal, nEntries, fromA, valFromA));
          fromB[0].i = er;
          fromB[0].j = ephi;
          fromB[0].k = ez;
          fromB[0].loc = BACK_DOWN_LEFT;
          fromB[0].c = 0;
          fromB[1].i = er;
          fromB[1].j = ephi;
          fromB[1].k = ez;
          fromB[1].loc = BACK_DOWN_LEFT;
          fromB[1].c = 1;
          fromB[2].i = er;
          fromB[2].j = ephi;
          fromB[2].k = ez;
          fromB[2].loc = BACK_DOWN_LEFT;
          fromB[2].c = 2;
          PetscCall(DMStagVecGetValuesStencil(da, BLocal, nEntries, fromB, valFromB));
          to[0].i = er;
          to[0].j = ephi;
          to[0].k = ez;
          to[0].loc = BACK_DOWN_LEFT;
          to[0].c = 0;
          valTo[0] = valFromA[1] * valFromB[2] - valFromA[2] * valFromB[1];
          to[1].i = er;
          to[1].j = ephi;
          to[1].k = ez;
          to[1].loc = BACK_DOWN_LEFT;
          to[1].c = 1;
          valTo[1] = valFromA[2] * valFromB[0] - valFromA[0] * valFromB[2];
          to[2].i = er;
          to[2].j = ephi;
          to[2].k = ez;
          to[2].loc = BACK_DOWN_LEFT;
          to[2].c = 2;
          valTo[2] = valFromA[0] * valFromB[1] - valFromA[1] * valFromB[0];
          PetscCall(DMStagVecSetValuesStencil(da, C, 3, to, valTo, INSERT_VALUES));

          if(er==N[0]-1){
            fromA[0].i = er;
            fromA[0].j = ephi;
            fromA[0].k = ez;
            fromA[0].loc = BACK_DOWN_RIGHT;
            fromA[0].c = 0;
            fromA[1].i = er;
            fromA[1].j = ephi;
            fromA[1].k = ez;
            fromA[1].loc = BACK_DOWN_RIGHT;
            fromA[1].c = 1;
            fromA[2].i = er;
            fromA[2].j = ephi;
            fromA[2].k = ez;
            fromA[2].loc = BACK_DOWN_RIGHT;
            fromA[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, ALocal, nEntries, fromA, valFromA));
            fromB[0].i = er;
            fromB[0].j = ephi;
            fromB[0].k = ez;
            fromB[0].loc = BACK_DOWN_RIGHT;
            fromB[0].c = 0;
            fromB[1].i = er;
            fromB[1].j = ephi;
            fromB[1].k = ez;
            fromB[1].loc = BACK_DOWN_RIGHT;
            fromB[1].c = 1;
            fromB[2].i = er;
            fromB[2].j = ephi;
            fromB[2].k = ez;
            fromB[2].loc = BACK_DOWN_RIGHT;
            fromB[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, BLocal, nEntries, fromB, valFromB));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_DOWN_RIGHT;
            to[0].c = 0;
            valTo[0] = valFromA[1] * valFromB[2] - valFromA[2] * valFromB[1];
            to[1].i = er;
            to[1].j = ephi;
            to[1].k = ez;
            to[1].loc = BACK_DOWN_RIGHT;
            to[1].c = 1;
            valTo[1] = valFromA[2] * valFromB[0] - valFromA[0] * valFromB[2];
            to[2].i = er;
            to[2].j = ephi;
            to[2].k = ez;
            to[2].loc = BACK_DOWN_RIGHT;
            to[2].c = 2;
            valTo[2] = valFromA[0] * valFromB[1] - valFromA[1] * valFromB[0];
            PetscCall(DMStagVecSetValuesStencil(da, C, 3, to, valTo, INSERT_VALUES));
          }

          if(ez==N[2]-1){
            fromA[0].i = er;
            fromA[0].j = ephi;
            fromA[0].k = ez;
            fromA[0].loc = FRONT_DOWN_LEFT;
            fromA[0].c = 0;
            fromA[1].i = er;
            fromA[1].j = ephi;
            fromA[1].k = ez;
            fromA[1].loc = FRONT_DOWN_LEFT;
            fromA[1].c = 1;
            fromA[2].i = er;
            fromA[2].j = ephi;
            fromA[2].k = ez;
            fromA[2].loc = FRONT_DOWN_LEFT;
            fromA[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, ALocal, nEntries, fromA, valFromA));
            fromB[0].i = er;
            fromB[0].j = ephi;
            fromB[0].k = ez;
            fromB[0].loc = FRONT_DOWN_LEFT;
            fromB[0].c = 0;
            fromB[1].i = er;
            fromB[1].j = ephi;
            fromB[1].k = ez;
            fromB[1].loc = FRONT_DOWN_LEFT;
            fromB[1].c = 1;
            fromB[2].i = er;
            fromB[2].j = ephi;
            fromB[2].k = ez;
            fromB[2].loc = FRONT_DOWN_LEFT;
            fromB[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, BLocal, nEntries, fromB, valFromB));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_DOWN_LEFT;
            to[0].c = 0;
            valTo[0] = valFromA[1] * valFromB[2] - valFromA[2] * valFromB[1];
            to[1].i = er;
            to[1].j = ephi;
            to[1].k = ez;
            to[1].loc = FRONT_DOWN_LEFT;
            to[1].c = 1;
            valTo[1] = valFromA[2] * valFromB[0] - valFromA[0] * valFromB[2];
            to[2].i = er;
            to[2].j = ephi;
            to[2].k = ez;
            to[2].loc = FRONT_DOWN_LEFT;
            to[2].c = 2;
            valTo[2] = valFromA[0] * valFromB[1] - valFromA[1] * valFromB[0];
            PetscCall(DMStagVecSetValuesStencil(da, C, 3, to, valTo, INSERT_VALUES));
          }

          if(ephi==N[1]-1 && !user->phibtype){
            fromA[0].i = er;
            fromA[0].j = ephi;
            fromA[0].k = ez;
            fromA[0].loc = BACK_UP_LEFT;
            fromA[0].c = 0;
            fromA[1].i = er;
            fromA[1].j = ephi;
            fromA[1].k = ez;
            fromA[1].loc = BACK_UP_LEFT;
            fromA[1].c = 1;
            fromA[2].i = er;
            fromA[2].j = ephi;
            fromA[2].k = ez;
            fromA[2].loc = BACK_UP_LEFT;
            fromA[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, ALocal, nEntries, fromA, valFromA));
            fromB[0].i = er;
            fromB[0].j = ephi;
            fromB[0].k = ez;
            fromB[0].loc = BACK_UP_LEFT;
            fromB[0].c = 0;
            fromB[1].i = er;
            fromB[1].j = ephi;
            fromB[1].k = ez;
            fromB[1].loc = BACK_UP_LEFT;
            fromB[1].c = 1;
            fromB[2].i = er;
            fromB[2].j = ephi;
            fromB[2].k = ez;
            fromB[2].loc = BACK_UP_LEFT;
            fromB[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, BLocal, nEntries, fromB, valFromB));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_LEFT;
            to[0].c = 0;
            valTo[0] = valFromA[1] * valFromB[2] - valFromA[2] * valFromB[1];
            to[1].i = er;
            to[1].j = ephi;
            to[1].k = ez;
            to[1].loc = BACK_UP_LEFT;
            to[1].c = 1;
            valTo[1] = valFromA[2] * valFromB[0] - valFromA[0] * valFromB[2];
            to[2].i = er;
            to[2].j = ephi;
            to[2].k = ez;
            to[2].loc = BACK_UP_LEFT;
            to[2].c = 2;
            valTo[2] = valFromA[0] * valFromB[1] - valFromA[1] * valFromB[0];
            PetscCall(DMStagVecSetValuesStencil(da, C, 3, to, valTo, INSERT_VALUES));
          }

          if(er==N[0]-1 && ephi==N[1]-1 && !user->phibtype){
            fromA[0].i = er;
            fromA[0].j = ephi;
            fromA[0].k = ez;
            fromA[0].loc = BACK_UP_RIGHT;
            fromA[0].c = 0;
            fromA[1].i = er;
            fromA[1].j = ephi;
            fromA[1].k = ez;
            fromA[1].loc = BACK_UP_RIGHT;
            fromA[1].c = 1;
            fromA[2].i = er;
            fromA[2].j = ephi;
            fromA[2].k = ez;
            fromA[2].loc = BACK_UP_RIGHT;
            fromA[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, ALocal, nEntries, fromA, valFromA));
            fromB[0].i = er;
            fromB[0].j = ephi;
            fromB[0].k = ez;
            fromB[0].loc = BACK_UP_RIGHT;
            fromB[0].c = 0;
            fromB[1].i = er;
            fromB[1].j = ephi;
            fromB[1].k = ez;
            fromB[1].loc = BACK_UP_RIGHT;
            fromB[1].c = 1;
            fromB[2].i = er;
            fromB[2].j = ephi;
            fromB[2].k = ez;
            fromB[2].loc = BACK_UP_RIGHT;
            fromB[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, BLocal, nEntries, fromB, valFromB));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = BACK_UP_RIGHT;
            to[0].c = 0;
            valTo[0] = valFromA[1] * valFromB[2] - valFromA[2] * valFromB[1];
            to[1].i = er;
            to[1].j = ephi;
            to[1].k = ez;
            to[1].loc = BACK_UP_RIGHT;
            to[1].c = 1;
            valTo[1] = valFromA[2] * valFromB[0] - valFromA[0] * valFromB[2];
            to[2].i = er;
            to[2].j = ephi;
            to[2].k = ez;
            to[2].loc = BACK_UP_RIGHT;
            to[2].c = 2;
            valTo[2] = valFromA[0] * valFromB[1] - valFromA[1] * valFromB[0];
            PetscCall(DMStagVecSetValuesStencil(da, C, 3, to, valTo, INSERT_VALUES));
          }

          if(ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
            fromA[0].i = er;
            fromA[0].j = ephi;
            fromA[0].k = ez;
            fromA[0].loc = FRONT_UP_LEFT;
            fromA[0].c = 0;
            fromA[1].i = er;
            fromA[1].j = ephi;
            fromA[1].k = ez;
            fromA[1].loc = FRONT_UP_LEFT;
            fromA[1].c = 1;
            fromA[2].i = er;
            fromA[2].j = ephi;
            fromA[2].k = ez;
            fromA[2].loc = FRONT_UP_LEFT;
            fromA[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, ALocal, nEntries, fromA, valFromA));
            fromB[0].i = er;
            fromB[0].j = ephi;
            fromB[0].k = ez;
            fromB[0].loc = FRONT_UP_LEFT;
            fromB[0].c = 0;
            fromB[1].i = er;
            fromB[1].j = ephi;
            fromB[1].k = ez;
            fromB[1].loc = FRONT_UP_LEFT;
            fromB[1].c = 1;
            fromB[2].i = er;
            fromB[2].j = ephi;
            fromB[2].k = ez;
            fromB[2].loc = FRONT_UP_LEFT;
            fromB[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, BLocal, nEntries, fromB, valFromB));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_LEFT;
            to[0].c = 0;
            valTo[0] = valFromA[1] * valFromB[2] - valFromA[2] * valFromB[1];
            to[1].i = er;
            to[1].j = ephi;
            to[1].k = ez;
            to[1].loc = FRONT_UP_LEFT;
            to[1].c = 1;
            valTo[1] = valFromA[2] * valFromB[0] - valFromA[0] * valFromB[2];
            to[2].i = er;
            to[2].j = ephi;
            to[2].k = ez;
            to[2].loc = FRONT_UP_LEFT;
            to[2].c = 2;
            valTo[2] = valFromA[0] * valFromB[1] - valFromA[1] * valFromB[0];
            PetscCall(DMStagVecSetValuesStencil(da, C, 3, to, valTo, INSERT_VALUES));
          }

          if(er==N[0]-1 && ez==N[2]-1){
            fromA[0].i = er;
            fromA[0].j = ephi;
            fromA[0].k = ez;
            fromA[0].loc = FRONT_DOWN_RIGHT;
            fromA[0].c = 0;
            fromA[1].i = er;
            fromA[1].j = ephi;
            fromA[1].k = ez;
            fromA[1].loc = FRONT_DOWN_RIGHT;
            fromA[1].c = 1;
            fromA[2].i = er;
            fromA[2].j = ephi;
            fromA[2].k = ez;
            fromA[2].loc = FRONT_DOWN_RIGHT;
            fromA[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, ALocal, nEntries, fromA, valFromA));
            fromB[0].i = er;
            fromB[0].j = ephi;
            fromB[0].k = ez;
            fromB[0].loc = FRONT_DOWN_RIGHT;
            fromB[0].c = 0;
            fromB[1].i = er;
            fromB[1].j = ephi;
            fromB[1].k = ez;
            fromB[1].loc = FRONT_DOWN_RIGHT;
            fromB[1].c = 1;
            fromB[2].i = er;
            fromB[2].j = ephi;
            fromB[2].k = ez;
            fromB[2].loc = FRONT_DOWN_RIGHT;
            fromB[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, BLocal, nEntries, fromB, valFromB));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_DOWN_RIGHT;
            to[0].c = 0;
            valTo[0] = valFromA[1] * valFromB[2] - valFromA[2] * valFromB[1];
            to[1].i = er;
            to[1].j = ephi;
            to[1].k = ez;
            to[1].loc = FRONT_DOWN_RIGHT;
            to[1].c = 1;
            valTo[1] = valFromA[2] * valFromB[0] - valFromA[0] * valFromB[2];
            to[2].i = er;
            to[2].j = ephi;
            to[2].k = ez;
            to[2].loc = FRONT_DOWN_RIGHT;
            to[2].c = 2;
            valTo[2] = valFromA[0] * valFromB[1] - valFromA[1] * valFromB[0];
            PetscCall(DMStagVecSetValuesStencil(da, C, 3, to, valTo, INSERT_VALUES));
          }

          if(er==N[0]-1 && ez==N[2]-1 && ephi==N[1]-1 && !user->phibtype){
            fromA[0].i = er;
            fromA[0].j = ephi;
            fromA[0].k = ez;
            fromA[0].loc = FRONT_UP_RIGHT;
            fromA[0].c = 0;
            fromA[1].i = er;
            fromA[1].j = ephi;
            fromA[1].k = ez;
            fromA[1].loc = FRONT_UP_RIGHT;
            fromA[1].c = 1;
            fromA[2].i = er;
            fromA[2].j = ephi;
            fromA[2].k = ez;
            fromA[2].loc = FRONT_UP_RIGHT;
            fromA[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, ALocal, nEntries, fromA, valFromA));
            fromB[0].i = er;
            fromB[0].j = ephi;
            fromB[0].k = ez;
            fromB[0].loc = FRONT_UP_RIGHT;
            fromB[0].c = 0;
            fromB[1].i = er;
            fromB[1].j = ephi;
            fromB[1].k = ez;
            fromB[1].loc = FRONT_UP_RIGHT;
            fromB[1].c = 1;
            fromB[2].i = er;
            fromB[2].j = ephi;
            fromB[2].k = ez;
            fromB[2].loc = FRONT_UP_RIGHT;
            fromB[2].c = 2;
            PetscCall(DMStagVecGetValuesStencil(da, BLocal, nEntries, fromB, valFromB));
            to[0].i = er;
            to[0].j = ephi;
            to[0].k = ez;
            to[0].loc = FRONT_UP_RIGHT;
            to[0].c = 0;
            valTo[0] = valFromA[1] * valFromB[2] - valFromA[2] * valFromB[1];
            to[1].i = er;
            to[1].j = ephi;
            to[1].k = ez;
            to[1].loc = FRONT_UP_RIGHT;
            to[1].c = 1;
            valTo[1] = valFromA[2] * valFromB[0] - valFromA[0] * valFromB[2];
            to[2].i = er;
            to[2].j = ephi;
            to[2].k = ez;
            to[2].loc = FRONT_UP_RIGHT;
            to[2].c = 2;
            valTo[2] = valFromA[0] * valFromB[1] - valFromA[1] * valFromB[0];
            PetscCall(DMStagVecSetValuesStencil(da, C, 3, to, valTo, INSERT_VALUES));
          }
        }
      }
    }
    PetscCall(VecAssemblyBegin(C));
    PetscCall(VecAssemblyEnd(C));
    PetscCall(DMRestoreLocalVector(da,&ALocal));
    PetscCall(DMRestoreLocalVector(da,&BLocal));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode getEJArray(TS ts, Vec X, PetscScalar *ge_ER, PetscScalar *ge_EP, PetscScalar *ge_EZ, void *ptr, int code) {
  PetscFunctionBeginUser;
    PetscLogEvent  USER_EVENT;
    PetscClassId   classid;

    PetscCall(PetscClassIdRegister("class name",&classid));
    PetscCall(PetscLogEventRegister("getEJArray",classid,&USER_EVENT));
    PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

    User           *user = (User*)ptr;
    DM             da;
    DM             dmEr, dmEphi,dmEz, daEr, daEphi,daEz;
    PetscInt       startr,startphi,startz,nr,nphi,nz;
    Vec            F, C, vecEr, vecEphi, vecEz, E_r2, E_phi2, E_z2;
    Vec            XLocal;
    PetscInt       N[3],er,ephi,ez;
    const PetscScalar *array;

#line 8647
    PetscCall(TSGetDM(ts,&da));
    PetscCall(VecDuplicate(X, & F));
    PetscCall(VecZeroEntries(F));
    PetscCall(VecDuplicate(X, & C));
    PetscCall(VecZeroEntries(C));
    //NEED TO CREATE VEC C AND RECONSTRUCT FROM EDGES TO CELLS THE THREE COMPONENTS OF E BEFORE CALLING THE SCATTERING
    if (code == 0)
      FormElectricField(ts, X, F, user);
    if (code == 1)
      FormDerivedCurlnomp(ts, X, F, user);
    //DumpEdgeField(ts, 0, X, user);
    //DumpEdgeField(ts, 1, F, user);
    EdgeToCellReconstruction_r(ts,F,C,user);
    //DumpSolution_Cell(ts, 0, C, user);

    DMStagCreateCompatibleDMStag(da,0,0,0,1,&dmEr); /* 1 dof per cell */
    PetscCall(DMSetUp(dmEr));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmEr,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmEr,&E_r2));

    DMStagCreateCompatibleDMStag(da,0,0,0,1,&dmEphi); /* 1 dof per cell */
    PetscCall(DMSetUp(dmEphi));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmEphi,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmEphi,&E_phi2));

    DMStagCreateCompatibleDMStag(da,0,0,0,1,&dmEz); /* 1 dof per cell */
    PetscCall(DMSetUp(dmEz));
    PetscCall(DMStagSetUniformCoordinatesExplicit(dmEz,user->rmin,user->rmax,user->phimin,user->phimax,user->zmin,user->zmax));
    PetscCall(DMCreateGlobalVector(dmEz,&E_z2));

    PetscCall(DMGetLocalVector(da, & XLocal));
    PetscCall(DMGlobalToLocalBegin(da, C, INSERT_VALUES, XLocal));
    PetscCall(DMGlobalToLocalEnd(da, C, INSERT_VALUES, XLocal));

    //PetscPrintf(PETSC_COMM_WORLD,"Before copying E_r values\n");
    PetscCall(DMStagGetCorners(dmEr,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmEr,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[1];
                PetscScalar   valFrom[1];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = ELEMENT;    from[0].c = 0;
                PetscCall(DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmEr,E_r2,1,from,valFrom,INSERT_VALUES));
            }
        }
    }
    PetscCall(VecAssemblyBegin(E_r2));
    PetscCall(VecAssemblyEnd(E_r2));

    DMStagVecSplitToDMDA(dmEr,E_r2,ELEMENT,-1,&daEr,&vecEr); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecEr,"r_component_cell_center_values"));
    //VecView(vecEr, PETSC_VIEWER_STDOUT_WORLD);
    PetscCall(DMRestoreLocalVector(da, & XLocal));


    PetscCall(VecZeroEntries(C));
    EdgeToCellReconstruction_phi(ts,F,C,user);
    //DumpSolution_Cell(ts, 1, C, user);
    PetscCall(DMGetLocalVector(da, & XLocal));
    PetscCall(DMGlobalToLocalBegin(da, C, INSERT_VALUES, XLocal));
    PetscCall(DMGlobalToLocalEnd(da, C, INSERT_VALUES, XLocal));

    //PetscPrintf(PETSC_COMM_WORLD,"Before copying E_phi values\n");
    PetscCall(DMStagGetCorners(dmEphi,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmEphi,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[1];
                PetscScalar   valFrom[1];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = ELEMENT;    from[0].c = 0;
                PetscCall(DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmEphi,E_phi2,1,from,valFrom,INSERT_VALUES));
            }
        }
    }
    PetscCall(VecAssemblyBegin(E_phi2));
    PetscCall(VecAssemblyEnd(E_phi2));

    DMStagVecSplitToDMDA(dmEphi,E_phi2,ELEMENT,-1,&daEphi,&vecEphi); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecEphi,"phi_component_cell_center_values"));
    //VecView(vecEphi, PETSC_VIEWER_STDOUT_WORLD);
    PetscCall(DMRestoreLocalVector(da, & XLocal));


    PetscCall(VecZeroEntries(C));
    EdgeToCellReconstruction_z(ts,F,C,user);
    //DumpSolution_Cell(ts, 2, C, user);
    PetscCall(DMGetLocalVector(da, & XLocal));
    PetscCall(DMGlobalToLocalBegin(da, C, INSERT_VALUES, XLocal));
    PetscCall(DMGlobalToLocalEnd(da, C, INSERT_VALUES, XLocal));

    //PetscPrintf(PETSC_COMM_WORLD,"Before copying E_z values\n");
    PetscCall(DMStagGetCorners(dmEz,&startr,&startphi,&startz,&nr,&nphi,&nz,NULL,NULL,NULL));
    PetscCall(DMStagGetGlobalSizes(dmEz,&N[0],&N[1],&N[2]));
    for (ez = startz; ez<startz+nz; ++ez) {
        for (ephi = startphi; ephi<startphi+nphi; ++ephi) {
            for (er = startr; er<startr+nr; ++er) {
                DMStagStencil from[1];
                PetscScalar   valFrom[1];
                from[0].i = er; from[0].j = ephi; from[0].k = ez; from[0].loc = ELEMENT;    from[0].c = 0;
                PetscCall(DMStagVecGetValuesStencil(da,XLocal,1,from,valFrom));
                PetscCall(DMStagVecSetValuesStencil(dmEz,E_z2,1,from,valFrom,INSERT_VALUES));
            }
        }
    }
    PetscCall(VecAssemblyBegin(E_z2));
    PetscCall(VecAssemblyEnd(E_z2));

    DMStagVecSplitToDMDA(dmEz,E_z2,ELEMENT,-1,&daEz,&vecEz); /* note -3 : pad with zero */
    PetscCall(PetscObjectSetName((PetscObject)vecEz,"z_component_cell_center_values"));
    //VecView(vecEz, PETSC_VIEWER_STDOUT_WORLD);
    PetscCall(DMRestoreLocalVector(da, & XLocal));


    PetscMPIInt rank;

    VecScatter  scat;
    Vec         Xseq, naturalX;

    DMDACreateNaturalVector(daEr,&naturalX);
    DMDAGlobalToNaturalBegin(daEr, vecEr, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daEr, vecEr, INSERT_VALUES, naturalX);

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
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*(user->Nr)*user->Nz);
      PetscCall(VecGetArrayRead(Xseq, &array));
      memcpy(ge_ER, array, user->Nr*user->Nz*user->Nphi*(sizeof(PetscScalar)));
      PetscCall(VecRestoreArrayRead(Xseq, &array));
    }

    PetscCall(VecDestroy(&Xseq));
    VecScatterDestroy(&scat);
    PetscCall(VecDestroy(&naturalX));


    DMDACreateNaturalVector(daEphi,&naturalX);
    DMDAGlobalToNaturalBegin(daEphi, vecEphi, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daEphi, vecEphi, INSERT_VALUES, naturalX);

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
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,(user->Nphi)*user->Nr*user->Nz);
      PetscCall(VecGetArrayRead(Xseq, &array));
      memcpy(ge_EP, array, user->Nphi*user->Nz*user->Nr*(sizeof(PetscScalar)));
      PetscCall(VecRestoreArrayRead(Xseq, &array));
    }

    PetscCall(VecDestroy(&Xseq));
    VecScatterDestroy(&scat);
    PetscCall(VecDestroy(&naturalX));


    DMDACreateNaturalVector(daEz,&naturalX);
    DMDAGlobalToNaturalBegin(daEz, vecEz, INSERT_VALUES, naturalX);
    DMDAGlobalToNaturalEnd(daEz, vecEz, INSERT_VALUES, naturalX);

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
      //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*user->Nr*(user->Nz));
      PetscCall(VecGetArrayRead(Xseq, &array));
      memcpy(ge_EZ, array, (user->Nz)*user->Nphi*user->Nr*(sizeof(PetscScalar)));
      PetscCall(VecRestoreArrayRead(Xseq, &array));
    }

    PetscCall(VecDestroy(&Xseq));
    VecScatterDestroy(&scat);
    PetscCall(VecDestroy(&naturalX));
    PetscCall(VecDestroy(&F));
    PetscCall(VecDestroy(&C));

    PetscCall(DMDestroy(&dmEr));
    PetscCall(DMDestroy(&dmEphi));
    PetscCall(DMDestroy(&dmEz));

    PetscCall(DMDestroy(&daEr));
    PetscCall(DMDestroy(&daEphi));
    PetscCall(DMDestroy(&daEz));

    PetscCall(VecDestroy(&vecEr));
    PetscCall(VecDestroy(&vecEphi));
    PetscCall(VecDestroy(&vecEz));

    PetscCall(VecDestroy(&E_r2));
    PetscCall(VecDestroy(&E_phi2));
    PetscCall(VecDestroy(&E_z2));

  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Recovers global V from Petsc DMStag (vertices)
   [VR Vphi VZ]_[iR, iphi, iZ]

   Ordering of V:
   Component index fastest fastest, then iZ, then iphi, iR slowest:
   { VR_{0,0,0} Vphi_{0,0,0} VZ_{0,0,0} VR_{0,0,1} Vphi_{0,0,1} VZ_{0,0,1} ...  VZ_{0,0,NZ-1}
   VR_{0,1,0} ... VZ_{0,Nphi-1,NZ-1} VR_{1,0,0} ... VZ_{NR-1,Nphi-1,NZ-1}}
*/


PetscErrorCode getVArray(TS ts, Vec X, PetscScalar *gf_V, void *ptr)
{
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
  PetscInt       er,ephi,ez;
  const PetscScalar *array;

#line 8897
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

#line 9123


PetscErrorCode getBArray(TS ts, Vec X, PetscScalar *gf_B, void *ptr, int derivative)
{
  PetscFunctionBeginUser;
  PetscLogEvent  USER_EVENT;
  PetscClassId   classid;

  PetscCall(PetscClassIdRegister("class name",&classid));
  PetscCall(PetscLogEventRegister("UpdateBArray",classid,&USER_EVENT));
  PetscCall(PetscLogEventBegin(USER_EVENT,0,0,0,0));

  User           *user = (User*)ptr;
  DM             da, dmCoord, dmB, daB;
  PetscInt       startr,startphi,startz,nr,nphi,nz;
  Vec            vecB, B, X_Local, coordLocal;
  PetscInt       N[3],er,ephi,ez;
  const PetscScalar *array;

#line 9143
  PetscCall(TSGetDM(ts, & da));

  PetscCall(DMStagCreateCompatibleDMStag(da, 0, 0, 0, 3, & dmB)); /* 3 dofs per element */
  PetscCall(DMSetUp(dmB));
  PetscCall(DMStagSetUniformCoordinatesExplicit(dmB,  user -> rmin / user->L0, user -> rmax / user->L0, user -> phimin, user -> phimax, user -> zmin / user->L0, user -> zmax / user->L0));
  PetscCall(DMCreateGlobalVector(dmB, & B));
  PetscCall(DMGetLocalVector(da, & X_Local));
  PetscCall(DMGlobalToLocal(da, X, INSERT_VALUES, X_Local));

  PetscCall(DMStagGetCorners(dmB, &startr, &startphi, &startz,
                        &nr,     &nphi,     &nz,
                        NULL,    NULL,     NULL));
  PetscCall(DMStagGetGlobalSizes(dmB,&N[0],&N[1],&N[2]));

  PetscCall(DMGetCoordinateDM(dmB, &dmCoord));
  PetscCall(DMGetCoordinatesLocal(dmB, &coordLocal));

  PetscInt nComp = 3; // BR Bphi BZ
  PetscInt nVals = 6; // RBR RBphi RBZ dRBRdR dRBphidphi dRBZdZ
  for (ez = startz; ez < startz + nz; ++ez)
  {
    for (ephi = startphi; ephi < startphi + nphi; ++ephi)
    {
      for (er = startr; er < startr + nr; ++er)
      {
        DMStagStencil from[nVals], to[nComp];
        PetscScalar valFrom[nVals], valTo[nComp];
        PetscScalar  R[nVals];
        for (PetscInt index = 0; index < nVals; ++index)
        {
          from[index].i = er;
          from[index].j = ephi;
          from[index].k = ez;
          from[index].c = 0;
        }
        from[0].loc = LEFT;
        from[1].loc = RIGHT;
        from[2].loc = UP;
        from[3].loc = DOWN;
        from[4].loc = BACK;
        from[5].loc = FRONT;

        PetscCall(DMStagVecGetValuesStencil(da, X_Local, nVals, from, valFrom));
        PetscCall(DMStagVecGetValuesStencil(dmCoord, coordLocal, nVals, from, R));

        PetscScalar volume = PetscAbsReal(PetscSqr(R[1]) - PetscSqr(R[0])) * 0.5 * user->dphi * user->dz;

        for (PetscInt comp = 0; comp < nComp; ++comp)
        {
          to[comp].i = er;
          to[comp].j = ephi;
          to[comp].k = ez;
          to[comp].loc = ELEMENT;
          to[comp].c = comp;

          // if(comp < nComp)
          if(derivative == 0)
            valTo[comp] =  (valFrom[comp*2] * R[comp*2] * surface(er, ephi, ez, from[comp*2].loc, user) + valFrom[comp*2+1] * R[comp*2+1] * surface(er, ephi, ez, from[comp*2+1].loc, user)) /
                    (surface(er, ephi, ez, from[comp*2].loc, user) + surface(er, ephi, ez, from[comp*2+1].loc, user));
          // else
          else if(derivative == 1)
            valTo[comp] = (-valFrom[comp*2] * surface(er, ephi, ez, from[comp*2].loc, user) + valFrom[comp*2+1] * surface(er, ephi, ez, from[comp*2+1].loc, user)) / volume ;
        }
        PetscCall(DMStagVecSetValuesStencil(dmB, B, nComp, to, valTo, INSERT_VALUES));
      }
    }
  }
  PetscCall(VecAssemblyBegin(B));
  PetscCall(VecAssemblyEnd(B));

  PetscCall(DMStagVecSplitToDMDA(dmB, B, ELEMENT, -3, & daB, & vecB)); /* note -3 : pad with zero in 2D case */
  PetscCall(PetscObjectSetName((PetscObject) vecB, "Magneric field"));

  PetscMPIInt rank;

  VecScatter  scat;
  Vec         Xseq, naturalX;

  PetscCall(DMDACreateNaturalVector(daB,&naturalX));
  PetscCall(DMDAGlobalToNaturalBegin(daB, vecB, INSERT_VALUES, naturalX));
  PetscCall(DMDAGlobalToNaturalEnd(daB, vecB, INSERT_VALUES, naturalX));

  PetscCall(VecScatterCreateToAll(naturalX, &scat, &Xseq));
  PetscCall(VecScatterBegin(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterEnd(scat, naturalX, Xseq, INSERT_VALUES, SCATTER_FORWARD));

  if (rank == 0 || 1) {
    PetscInt sizeX;
    PetscCall(VecGetSize(Xseq, &sizeX));
    //PetscPrintf(PETSC_COMM_SELF,"The size of Xseq is %d, and the grid size is %d\n",sizeX,user->Nphi*(user->Nr+1)*user->Nz);
    PetscCall(VecGetArrayRead(Xseq, &array));
    memcpy(gf_B, array, 3*user->Nr*user->Nz*user->Nphi*(sizeof(PetscScalar)));
    PetscCall(VecRestoreArrayRead(Xseq, &array));
  }

  PetscCall(VecDestroy(&Xseq));
  PetscCall(VecScatterDestroy(&scat));
  PetscCall(VecDestroy(&naturalX));

  /* Destroy DMDAs and Vecs */
  PetscCall(VecDestroy( & vecB));
  PetscCall(DMDestroy( & daB));
  PetscCall(VecDestroy( & B));
  PetscCall(DMDestroy( & dmB));


  PetscCall(PetscLogEventEnd(USER_EVENT,0,0,0,0));

  PetscFunctionReturn(PETSC_SUCCESS);
}


#line 9679

#line 9810

#line 10053

#line 10108





