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

PetscScalar betaf(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Face perpendicular to z direction */
  if (loc == BACK || loc == FRONT) {
    if ((ez == 0 && loc == BACK) || (ez == N[2] - 1 && loc == FRONT)) {
      beta = alphafc(er, ephi, ez, user);
    } else if (loc == BACK) {
      beta = alphafc(er, ephi, ez - 1, user) + alphafc(er, ephi, ez, user);
    } else {
      beta = alphafc(er, ephi, ez, user) + alphafc(er, ephi, ez + 1, user);
    }
  }
  /* Face perpendicular to phi direction */
  else if (loc == DOWN || loc == UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && loc == DOWN) || (ephi == N[1] - 1 && loc == UP)) {
        beta = alphafc(er, ephi, ez, user);
      } else if (loc == DOWN) {
        beta = alphafc(er, ephi - 1, ez, user) + alphafc(er, ephi, ez, user);
      } else {
        beta = alphafc(er, ephi, ez, user) + alphafc(er, ephi + 1, ez, user);
      }
    } else {
      if (loc == DOWN) {
        beta = alphafc(er, ephi - 1, ez, user) + alphafc(er, ephi, ez, user);
      } else {
        beta = alphafc(er, ephi, ez, user) + alphafc(er, ephi + 1, ez, user);
      }
    }
  }
  /* Face perpendicular to r direction */
  else if (loc == LEFT || loc == RIGHT) {
    if ((er == 0 && loc == LEFT) || (er == N[0] - 1 && loc == RIGHT)) {
      beta = alphafc(er, ephi, ez, user);
    } else if (loc == LEFT) {
      beta = alphafc(er - 1, ephi, ez, user) + alphafc(er, ephi, ez, user);
    } else {
      beta = alphafc(er, ephi, ez, user) + alphafc(er + 1, ephi, ez, user);
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betaf function");
  }
  return beta;
}

PetscScalar betaf_wmp(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Face perpendicular to z direction */
  if (loc == BACK || loc == FRONT) {
    if ((ez == 0 && loc == BACK) || (ez == N[2] - 1 && loc == FRONT)) {
      beta = alphafc_wmp(er, ephi, ez, user);
    } else if (loc == BACK) {
      beta = alphafc_wmp(er, ephi, ez - 1, user) + alphafc_wmp(er, ephi, ez, user);
    } else {
      beta = alphafc_wmp(er, ephi, ez, user) + alphafc_wmp(er, ephi, ez + 1, user);
    }
  }
  /* Face perpendicular to phi direction */
  else if (loc == DOWN || loc == UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && loc == DOWN) || (ephi == N[1] - 1 && loc == UP)) {
        beta = alphafc_wmp(er, ephi, ez, user);
      } else if (loc == DOWN) {
        beta = alphafc_wmp(er, ephi - 1, ez, user) + alphafc_wmp(er, ephi, ez, user);
      } else {
        beta = alphafc_wmp(er, ephi, ez, user) + alphafc_wmp(er, ephi + 1, ez, user);
      }
    } else {
      if (loc == DOWN) {
        beta = alphafc_wmp(er, ephi - 1, ez, user) + alphafc_wmp(er, ephi, ez, user);
      } else {
        beta = alphafc_wmp(er, ephi, ez, user) + alphafc_wmp(er, ephi + 1, ez, user);
      }
    }
  }
  /* Face perpendicular to r direction */
  else if (loc == LEFT || loc == RIGHT) {
    if ((er == 0 && loc == LEFT) || (er == N[0] - 1 && loc == RIGHT)) {
      beta = alphafc_wmp(er, ephi, ez, user);
    } else if (loc == LEFT) {
      beta = alphafc_wmp(er - 1, ephi, ez, user) + alphafc_wmp(er, ephi, ez, user);
    } else {
      beta = alphafc_wmp(er, ephi, ez, user) + alphafc_wmp(er + 1, ephi, ez, user);
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betaf function");
  }
  return beta;
}

PetscScalar betae(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = alphaec(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaec(er - 1, ephi - 1, ez, user) + alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi - 1, ez, user) + alphaec(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaec(er - 1, ephi + 1, ez, user) + alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi + 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaec(er, ephi + 1, ez, user) + alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi + 1, ez, user) + alphaec(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaec(er - 1, ephi - 1, ez, user) + alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi - 1, ez, user) + alphaec(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaec(er - 1, ephi + 1, ez, user) + alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi + 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaec(er, ephi + 1, ez, user) + alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi + 1, ez, user) + alphaec(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = alphaec(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = alphaec(er - 1, ephi, ez - 1, user) + alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi, ez - 1, user) + alphaec(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = alphaec(er - 1, ephi, ez + 1, user) + alphaec(er - 1, ephi, ez, user) + alphaec(er, ephi, ez + 1, user) + alphaec(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else {
      beta = alphaec(er, ephi, ez + 1, user) + alphaec(er, ephi, ez, user) + alphaec(er + 1, ephi, ez + 1, user) + alphaec(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = alphaec(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaec(er, ephi - 1, ez - 1, user) + alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez - 1, user) + alphaec(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaec(er, ephi - 1, ez + 1, user) + alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez + 1, user) + alphaec(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaec(er, ephi, ez + 1, user) + alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez + 1, user) + alphaec(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaec(er, ephi - 1, ez - 1, user) + alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaec(er, ephi, ez - 1, user) + alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez - 1, user) + alphaec(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaec(er, ephi - 1, ez + 1, user) + alphaec(er, ephi - 1, ez, user) + alphaec(er, ephi, ez + 1, user) + alphaec(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaec(er, ephi, ez + 1, user) + alphaec(er, ephi, ez, user) + alphaec(er, ephi + 1, ez + 1, user) + alphaec(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betae function");
  }
  return beta;
}

PetscScalar betaenores(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = alphaecnores(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecnores(er - 1, ephi - 1, ez, user) + alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi - 1, ez, user) + alphaecnores(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecnores(er - 1, ephi + 1, ez, user) + alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecnores(er, ephi + 1, ez, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi + 1, ez, user) + alphaecnores(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecnores(er - 1, ephi - 1, ez, user) + alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi - 1, ez, user) + alphaecnores(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecnores(er - 1, ephi + 1, ez, user) + alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecnores(er, ephi + 1, ez, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi + 1, ez, user) + alphaecnores(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = alphaecnores(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = alphaecnores(er - 1, ephi, ez - 1, user) + alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi, ez - 1, user) + alphaecnores(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = alphaecnores(er - 1, ephi, ez + 1, user) + alphaecnores(er - 1, ephi, ez, user) + alphaecnores(er, ephi, ez + 1, user) + alphaecnores(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else {
      beta = alphaecnores(er, ephi, ez + 1, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er + 1, ephi, ez + 1, user) + alphaecnores(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = alphaecnores(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecnores(er, ephi - 1, ez - 1, user) + alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez - 1, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecnores(er, ephi - 1, ez + 1, user) + alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez + 1, user) + alphaecnores(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecnores(er, ephi, ez + 1, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez + 1, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecnores(er, ephi - 1, ez - 1, user) + alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecnores(er, ephi, ez - 1, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez - 1, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecnores(er, ephi - 1, ez + 1, user) + alphaecnores(er, ephi - 1, ez, user) + alphaecnores(er, ephi, ez + 1, user) + alphaecnores(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecnores(er, ephi, ez + 1, user) + alphaecnores(er, ephi, ez, user) + alphaecnores(er, ephi + 1, ez + 1, user) + alphaecnores(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betaenores function");
  }
  return beta;
}

PetscScalar betaenomp(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = alphaecnomp(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecnomp(er - 1, ephi - 1, ez, user) + alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi - 1, ez, user) + alphaecnomp(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecnomp(er - 1, ephi + 1, ez, user) + alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecnomp(er, ephi + 1, ez, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi + 1, ez, user) + alphaecnomp(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecnomp(er - 1, ephi - 1, ez, user) + alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi - 1, ez, user) + alphaecnomp(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecnomp(er - 1, ephi + 1, ez, user) + alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecnomp(er, ephi + 1, ez, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi + 1, ez, user) + alphaecnomp(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = alphaecnomp(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = alphaecnomp(er - 1, ephi, ez - 1, user) + alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi, ez - 1, user) + alphaecnomp(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = alphaecnomp(er - 1, ephi, ez + 1, user) + alphaecnomp(er - 1, ephi, ez, user) + alphaecnomp(er, ephi, ez + 1, user) + alphaecnomp(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else {
      beta = alphaecnomp(er, ephi, ez + 1, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er + 1, ephi, ez + 1, user) + alphaecnomp(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = alphaecnomp(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecnomp(er, ephi - 1, ez - 1, user) + alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez - 1, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecnomp(er, ephi - 1, ez + 1, user) + alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez + 1, user) + alphaecnomp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecnomp(er, ephi, ez + 1, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez + 1, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecnomp(er, ephi - 1, ez - 1, user) + alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecnomp(er, ephi, ez - 1, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez - 1, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecnomp(er, ephi - 1, ez + 1, user) + alphaecnomp(er, ephi - 1, ez, user) + alphaecnomp(er, ephi, ez + 1, user) + alphaecnomp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecnomp(er, ephi, ez + 1, user) + alphaecnomp(er, ephi, ez, user) + alphaecnomp(er, ephi + 1, ez + 1, user) + alphaecnomp(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betaenomp function");
  }
  return beta;
}

PetscScalar betaeperp(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = alphaecperp(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecperp(er - 1, ephi - 1, ez, user) + alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi - 1, ez, user) + alphaecperp(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecperp(er - 1, ephi + 1, ez, user) + alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecperp(er, ephi + 1, ez, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi + 1, ez, user) + alphaecperp(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecperp(er - 1, ephi - 1, ez, user) + alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi - 1, ez, user) + alphaecperp(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecperp(er - 1, ephi + 1, ez, user) + alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecperp(er, ephi + 1, ez, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi + 1, ez, user) + alphaecperp(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = alphaecperp(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = alphaecperp(er - 1, ephi, ez - 1, user) + alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi, ez - 1, user) + alphaecperp(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = alphaecperp(er - 1, ephi, ez + 1, user) + alphaecperp(er - 1, ephi, ez, user) + alphaecperp(er, ephi, ez + 1, user) + alphaecperp(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else {
      beta = alphaecperp(er, ephi, ez + 1, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er + 1, ephi, ez + 1, user) + alphaecperp(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = alphaecperp(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecperp(er, ephi - 1, ez - 1, user) + alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez - 1, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecperp(er, ephi - 1, ez + 1, user) + alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez + 1, user) + alphaecperp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecperp(er, ephi, ez + 1, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez + 1, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecperp(er, ephi - 1, ez - 1, user) + alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecperp(er, ephi, ez - 1, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez - 1, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecperp(er, ephi - 1, ez + 1, user) + alphaecperp(er, ephi - 1, ez, user) + alphaecperp(er, ephi, ez + 1, user) + alphaecperp(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecperp(er, ephi, ez + 1, user) + alphaecperp(er, ephi, ez, user) + alphaecperp(er, ephi + 1, ez + 1, user) + alphaecperp(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betaeperp function");
  }
  return beta;
}

PetscScalar betaephi(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = alphaecphi(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecphi(er - 1, ephi - 1, ez, user) + alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi - 1, ez, user) + alphaecphi(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecphi(er - 1, ephi + 1, ez, user) + alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi(er, ephi + 1, ez, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi + 1, ez, user) + alphaecphi(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecphi(er - 1, ephi - 1, ez, user) + alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi - 1, ez, user) + alphaecphi(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecphi(er - 1, ephi + 1, ez, user) + alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi(er, ephi + 1, ez, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi + 1, ez, user) + alphaecphi(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = alphaecphi(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = alphaecphi(er - 1, ephi, ez - 1, user) + alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi, ez - 1, user) + alphaecphi(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = alphaecphi(er - 1, ephi, ez + 1, user) + alphaecphi(er - 1, ephi, ez, user) + alphaecphi(er, ephi, ez + 1, user) + alphaecphi(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else {
      beta = alphaecphi(er, ephi, ez + 1, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er + 1, ephi, ez + 1, user) + alphaecphi(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = alphaecphi(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecphi(er, ephi - 1, ez - 1, user) + alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez - 1, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecphi(er, ephi - 1, ez + 1, user) + alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez + 1, user) + alphaecphi(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi(er, ephi, ez + 1, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez + 1, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecphi(er, ephi - 1, ez - 1, user) + alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecphi(er, ephi, ez - 1, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez - 1, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecphi(er, ephi - 1, ez + 1, user) + alphaecphi(er, ephi - 1, ez, user) + alphaecphi(er, ephi, ez + 1, user) + alphaecphi(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi(er, ephi, ez + 1, user) + alphaecphi(er, ephi, ez, user) + alphaecphi(er, ephi + 1, ez + 1, user) + alphaecphi(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betaephi function");
  }
  return beta;
}

PetscScalar betae2(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = alphaec2(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaec2(er - 1, ephi - 1, ez, user) + alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi - 1, ez, user) + alphaec2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaec2(er - 1, ephi + 1, ez, user) + alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaec2(er, ephi + 1, ez, user) + alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi + 1, ez, user) + alphaec2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaec2(er - 1, ephi - 1, ez, user) + alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi - 1, ez, user) + alphaec2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaec2(er - 1, ephi + 1, ez, user) + alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaec2(er, ephi + 1, ez, user) + alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi + 1, ez, user) + alphaec2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = alphaec2(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = alphaec2(er - 1, ephi, ez - 1, user) + alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi, ez - 1, user) + alphaec2(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = alphaec2(er - 1, ephi, ez + 1, user) + alphaec2(er - 1, ephi, ez, user) + alphaec2(er, ephi, ez + 1, user) + alphaec2(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else {
      beta = alphaec2(er, ephi, ez + 1, user) + alphaec2(er, ephi, ez, user) + alphaec2(er + 1, ephi, ez + 1, user) + alphaec2(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = alphaec2(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaec2(er, ephi - 1, ez - 1, user) + alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez - 1, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaec2(er, ephi - 1, ez + 1, user) + alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez + 1, user) + alphaec2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaec2(er, ephi, ez + 1, user) + alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez + 1, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaec2(er, ephi - 1, ez - 1, user) + alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaec2(er, ephi, ez - 1, user) + alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez - 1, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaec2(er, ephi - 1, ez + 1, user) + alphaec2(er, ephi - 1, ez, user) + alphaec2(er, ephi, ez + 1, user) + alphaec2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaec2(er, ephi, ez + 1, user) + alphaec2(er, ephi, ez, user) + alphaec2(er, ephi + 1, ez + 1, user) + alphaec2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betae2 function");
  }
  return beta;
}

PetscScalar betaephi_isolcell(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = alphaecphi_isolcell(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecphi_isolcell(er - 1, ephi - 1, ez, user) + alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi - 1, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecphi_isolcell(er - 1, ephi + 1, ez, user) + alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi_isolcell(er, ephi + 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi + 1, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecphi_isolcell(er - 1, ephi - 1, ez, user) + alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi - 1, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecphi_isolcell(er - 1, ephi + 1, ez, user) + alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi_isolcell(er, ephi + 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi + 1, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = alphaecphi_isolcell(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = alphaecphi_isolcell(er - 1, ephi, ez - 1, user) + alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez - 1, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = alphaecphi_isolcell(er - 1, ephi, ez + 1, user) + alphaecphi_isolcell(er - 1, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez + 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else {
      beta = alphaecphi_isolcell(er, ephi, ez + 1, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er + 1, ephi, ez + 1, user) + alphaecphi_isolcell(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = alphaecphi_isolcell(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez - 1, user) + alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez - 1, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez + 1, user) + alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez + 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi_isolcell(er, ephi, ez + 1, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez + 1, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez - 1, user) + alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecphi_isolcell(er, ephi, ez - 1, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez - 1, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecphi_isolcell(er, ephi - 1, ez + 1, user) + alphaecphi_isolcell(er, ephi - 1, ez, user) + alphaecphi_isolcell(er, ephi, ez + 1, user) + alphaecphi_isolcell(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi_isolcell(er, ephi, ez + 1, user) + alphaecphi_isolcell(er, ephi, ez, user) + alphaecphi_isolcell(er, ephi + 1, ez + 1, user) + alphaecphi_isolcell(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betaephi_isolcell function");
  }
  return beta;
}

PetscScalar betaeperp2(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = alphaecperp2(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecperp2(er - 1, ephi - 1, ez, user) + alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi - 1, ez, user) + alphaecperp2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecperp2(er - 1, ephi + 1, ez, user) + alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecperp2(er, ephi + 1, ez, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi + 1, ez, user) + alphaecperp2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecperp2(er - 1, ephi - 1, ez, user) + alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi - 1, ez, user) + alphaecperp2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecperp2(er - 1, ephi + 1, ez, user) + alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecperp2(er, ephi + 1, ez, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi + 1, ez, user) + alphaecperp2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = alphaecperp2(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = alphaecperp2(er - 1, ephi, ez - 1, user) + alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi, ez - 1, user) + alphaecperp2(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = alphaecperp2(er - 1, ephi, ez + 1, user) + alphaecperp2(er - 1, ephi, ez, user) + alphaecperp2(er, ephi, ez + 1, user) + alphaecperp2(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else {
      beta = alphaecperp2(er, ephi, ez + 1, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er + 1, ephi, ez + 1, user) + alphaecperp2(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = alphaecperp2(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecperp2(er, ephi - 1, ez - 1, user) + alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez - 1, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecperp2(er, ephi - 1, ez + 1, user) + alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez + 1, user) + alphaecperp2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecperp2(er, ephi, ez + 1, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez + 1, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecperp2(er, ephi - 1, ez - 1, user) + alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecperp2(er, ephi, ez - 1, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez - 1, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecperp2(er, ephi - 1, ez + 1, user) + alphaecperp2(er, ephi - 1, ez, user) + alphaecperp2(er, ephi, ez + 1, user) + alphaecperp2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecperp2(er, ephi, ez + 1, user) + alphaecperp2(er, ephi, ez, user) + alphaecperp2(er, ephi + 1, ez + 1, user) + alphaecperp2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betaeperp2 function");
  }
  return beta;
}

PetscScalar betaephi2(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = alphaecphi2(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecphi2(er - 1, ephi - 1, ez, user) + alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi - 1, ez, user) + alphaecphi2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecphi2(er - 1, ephi + 1, ez, user) + alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi2(er, ephi + 1, ez, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi + 1, ez, user) + alphaecphi2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = alphaecphi2(er - 1, ephi - 1, ez, user) + alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi - 1, ez, user) + alphaecphi2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = alphaecphi2(er - 1, ephi + 1, ez, user) + alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi2(er, ephi + 1, ez, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi + 1, ez, user) + alphaecphi2(er + 1, ephi, ez, user);
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = alphaecphi2(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi, ez + 1, user);
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi, ez, user);
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = alphaecphi2(er - 1, ephi, ez - 1, user) + alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi, ez - 1, user) + alphaecphi2(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = alphaecphi2(er - 1, ephi, ez + 1, user) + alphaecphi2(er - 1, ephi, ez, user) + alphaecphi2(er, ephi, ez + 1, user) + alphaecphi2(er, ephi, ez, user);
    } /* 4 cells sharing the edge */
    else {
      beta = alphaecphi2(er, ephi, ez + 1, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er + 1, ephi, ez + 1, user) + alphaecphi2(er + 1, ephi, ez, user);
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = alphaecphi2(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi, ez + 1, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecphi2(er, ephi - 1, ez - 1, user) + alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez - 1, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecphi2(er, ephi - 1, ez + 1, user) + alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez + 1, user) + alphaecphi2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi2(er, ephi, ez + 1, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez + 1, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez, user);
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = alphaecphi2(er, ephi - 1, ez - 1, user) + alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = alphaecphi2(er, ephi, ez - 1, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez - 1, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = alphaecphi2(er, ephi - 1, ez + 1, user) + alphaecphi2(er, ephi - 1, ez, user) + alphaecphi2(er, ephi, ez + 1, user) + alphaecphi2(er, ephi, ez, user);
      } /* 4 cells sharing the edge */
      else {
        beta = alphaecphi2(er, ephi, ez + 1, user) + alphaecphi2(er, ephi, ez, user) + alphaecphi2(er, ephi + 1, ez + 1, user) + alphaecphi2(er, ephi + 1, ez, user);
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betaephi2 function");
  }
  return beta;
}

PetscScalar betav(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  if (!(user -> phibtype)) {
    //1 Cell connected to the vertex
    if ((loc == BACK_DOWN_LEFT && er == 0 && ephi == 0 && ez == 0) || (loc == BACK_DOWN_RIGHT && er == N[0] - 1 && ephi == 0 && ez == 0) || (loc == BACK_UP_LEFT && er == 0 && ephi == N[1] - 1 && ez == 0) || (loc == BACK_UP_RIGHT && er == N[0] - 1 && ephi == N[1] - 1 && ez == 0) || (loc == FRONT_DOWN_LEFT && er == 0 && ephi == 0 && ez == N[2] - 1) || (loc == FRONT_DOWN_RIGHT && er == N[0] - 1 && ephi == 0 && ez == N[2] - 1) || (loc == FRONT_UP_LEFT && er == 0 && ephi == N[1] - 1 && ez == N[2] - 1) || (loc == FRONT_UP_RIGHT && er == N[0] - 1 && ephi == N[1] - 1 && ez == N[2] - 1)) {
      beta = alphavc(er, ephi, ez, user);
    }

    //2 Cells connected to the vertex
    else if (er == 0 && ephi == 0 && loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user);
    } else if (er == 0 && ephi == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ephi == 0 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user);
    } else if (er == N[0] - 1 && ephi == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi, ez, user);
    } else if (er == 0 && ephi == N[1] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user);
    } else if (er == 0 && ephi == N[1] - 1 && loc == BACK_UP_LEFT) {
      beta = alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ephi == N[1] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user);
    } else if (er == N[0] - 1 && ephi == N[1] - 1 && loc == BACK_UP_RIGHT) {
      beta = alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi, ez, user);
    } else if (er == 0 && ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez, user);
    } else if (er == 0 && ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (er == N[0] - 1 && ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (er == 0 && ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez, user);
    } else if (er == 0 && ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (ephi == 0 && ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == 0 && ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user);
    } else if (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user);
    } else if (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user);
    } else if (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user);
    }

    //4 Cells connected to the vertex
    else if (er == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez - 1, user);
    } else if (er == 0 && loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez + 1, user);
    } else if (er == 0 && loc == BACK_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez - 1, user);
    } else if (er == 0 && loc == FRONT_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez + 1, user);
    } else if (er == N[0] - 1 && loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez - 1, user);
    } else if (er == N[0] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez + 1, user);
    } else if (er == N[0] - 1 && loc == BACK_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez - 1, user);
    } else if (er == N[0] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez + 1, user);
    } else if (ephi == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er - 1, ephi, ez - 1, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == 0 && loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er - 1, ephi, ez + 1, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er + 1, ephi, ez - 1, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == 0 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er + 1, ephi, ez + 1, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && loc == BACK_UP_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er - 1, ephi, ez - 1, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er - 1, ephi, ez + 1, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && loc == BACK_UP_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er + 1, ephi, ez - 1, user) + alphavc(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er + 1, ephi, ez + 1, user) + alphavc(er, ephi, ez, user);
    } else if (ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez, user);
    } else if (ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez, user);
    } else if (ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez, user);
    }

    //8 Cells connected to the vertex
    else if (loc == BACK_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er - 1, ephi - 1, ez, user) + alphavc(er - 1, ephi, ez - 1, user) + alphavc(er, ephi - 1, ez - 1, user) + alphavc(er - 1, ephi - 1, ez - 1, user);
    } else if (loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er + 1, ephi - 1, ez, user) + alphavc(er + 1, ephi, ez - 1, user) + alphavc(er, ephi - 1, ez - 1, user) + alphavc(er + 1, ephi - 1, ez - 1, user);
    } else if (loc == BACK_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er - 1, ephi + 1, ez, user) + alphavc(er - 1, ephi, ez - 1, user) + alphavc(er, ephi + 1, ez - 1, user) + alphavc(er - 1, ephi + 1, ez - 1, user);
    } else if (loc == BACK_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er + 1, ephi + 1, ez, user) + alphavc(er + 1, ephi, ez - 1, user) + alphavc(er, ephi + 1, ez - 1, user) + alphavc(er + 1, ephi + 1, ez - 1, user);
    } else if (loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er - 1, ephi - 1, ez, user) + alphavc(er - 1, ephi, ez + 1, user) + alphavc(er, ephi - 1, ez + 1, user) + alphavc(er - 1, ephi - 1, ez + 1, user);
    } else if (loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er + 1, ephi - 1, ez, user) + alphavc(er + 1, ephi, ez + 1, user) + alphavc(er, ephi - 1, ez + 1, user) + alphavc(er + 1, ephi - 1, ez + 1, user);
    } else if (loc == FRONT_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er - 1, ephi + 1, ez, user) + alphavc(er - 1, ephi, ez + 1, user) + alphavc(er, ephi + 1, ez + 1, user) + alphavc(er - 1, ephi + 1, ez + 1, user);
    } else if (loc == FRONT_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er + 1, ephi + 1, ez, user) + alphavc(er + 1, ephi, ez + 1, user) + alphavc(er, ephi + 1, ez + 1, user) + alphavc(er + 1, ephi + 1, ez + 1, user);
    }
  } else {
    //2 Cells connected to the vertex
    if (er == 0 && ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez, user);
    } else if (er == 0 && ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (er == N[0] - 1 && ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (er == 0 && ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez, user);
    } else if (er == 0 && ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi + 1, ez, user);
    }

    //4 Cells connected to the vertex
    else if (er == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez - 1, user);
    } else if (er == 0 && loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez + 1, user);
    } else if (er == 0 && loc == BACK_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez - 1, user);
    } else if (er == 0 && loc == FRONT_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez + 1, user);
    } else if (er == N[0] - 1 && loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez - 1, user);
    } else if (er == N[0] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez + 1, user);
    } else if (er == N[0] - 1 && loc == BACK_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez - 1, user);
    } else if (er == N[0] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez + 1, user);
    } else if (ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez, user);
    } else if (ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez, user);
    } else if (ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi - 1, ez, user) + alphavc(er, ephi - 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi + 1, ez, user) + alphavc(er, ephi + 1, ez, user);
    }

    //8 Cells connected to the vertex
    else if (loc == BACK_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er - 1, ephi - 1, ez, user) + alphavc(er - 1, ephi, ez - 1, user) + alphavc(er, ephi - 1, ez - 1, user) + alphavc(er - 1, ephi - 1, ez - 1, user);
    } else if (loc == BACK_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er + 1, ephi - 1, ez, user) + alphavc(er + 1, ephi, ez - 1, user) + alphavc(er, ephi - 1, ez - 1, user) + alphavc(er + 1, ephi - 1, ez - 1, user);
    } else if (loc == BACK_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er - 1, ephi + 1, ez, user) + alphavc(er - 1, ephi, ez - 1, user) + alphavc(er, ephi + 1, ez - 1, user) + alphavc(er - 1, ephi + 1, ez - 1, user);
    } else if (loc == BACK_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi, ez - 1, user) + alphavc(er + 1, ephi + 1, ez, user) + alphavc(er + 1, ephi, ez - 1, user) + alphavc(er, ephi + 1, ez - 1, user) + alphavc(er + 1, ephi + 1, ez - 1, user);
    } else if (loc == FRONT_DOWN_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er - 1, ephi - 1, ez, user) + alphavc(er - 1, ephi, ez + 1, user) + alphavc(er, ephi - 1, ez + 1, user) + alphavc(er - 1, ephi - 1, ez + 1, user);
    } else if (loc == FRONT_DOWN_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi - 1, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er + 1, ephi - 1, ez, user) + alphavc(er + 1, ephi, ez + 1, user) + alphavc(er, ephi - 1, ez + 1, user) + alphavc(er + 1, ephi - 1, ez + 1, user);
    } else if (loc == FRONT_UP_LEFT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er - 1, ephi, ez, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er - 1, ephi + 1, ez, user) + alphavc(er - 1, ephi, ez + 1, user) + alphavc(er, ephi + 1, ez + 1, user) + alphavc(er - 1, ephi + 1, ez + 1, user);
    } else if (loc == FRONT_UP_RIGHT) {
      beta = alphavc(er, ephi, ez, user) + alphavc(er + 1, ephi, ez, user) + alphavc(er, ephi + 1, ez, user) + alphavc(er, ephi, ez + 1, user) + alphavc(er + 1, ephi + 1, ez, user) + alphavc(er + 1, ephi, ez + 1, user) + alphavc(er, ephi + 1, ez + 1, user) + alphavc(er + 1, ephi + 1, ez + 1, user);
    }

  }

  if (!(loc == BACK_DOWN_LEFT || loc == BACK_DOWN_RIGHT || loc == BACK_UP_LEFT || loc == BACK_UP_RIGHT || loc == FRONT_DOWN_LEFT || loc == FRONT_DOWN_RIGHT || loc == FRONT_UP_LEFT || loc == FRONT_UP_RIGHT)) {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betav function");
  }
  return beta;
}

PetscScalar betavnomp(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);

  if (!(user -> phibtype)) {
    //1 Cell connected to the vertex
    if ((loc == BACK_DOWN_LEFT && er == 0 && ephi == 0 && ez == 0) || (loc == BACK_DOWN_RIGHT && er == N[0] - 1 && ephi == 0 && ez == 0) || (loc == BACK_UP_LEFT && er == 0 && ephi == N[1] - 1 && ez == 0) || (loc == BACK_UP_RIGHT && er == N[0] - 1 && ephi == N[1] - 1 && ez == 0) || (loc == FRONT_DOWN_LEFT && er == 0 && ephi == 0 && ez == N[2] - 1) || (loc == FRONT_DOWN_RIGHT && er == N[0] - 1 && ephi == 0 && ez == N[2] - 1) || (loc == FRONT_UP_LEFT && er == 0 && ephi == N[1] - 1 && ez == N[2] - 1) || (loc == FRONT_UP_RIGHT && er == N[0] - 1 && ephi == N[1] - 1 && ez == N[2] - 1)) {
      beta = alphavcnomp(er, ephi, ez, user);
    }

    //2 Cells connected to the vertex
    else if (er == 0 && ephi == 0 && loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user);
    } else if (er == 0 && ephi == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ephi == 0 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user);
    } else if (er == N[0] - 1 && ephi == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == 0 && ephi == N[1] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user);
    } else if (er == 0 && ephi == N[1] - 1 && loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ephi == N[1] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user);
    } else if (er == N[0] - 1 && ephi == N[1] - 1 && loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == 0 && ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == 0 && ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (er == N[0] - 1 && ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (er == 0 && ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == 0 && ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (ephi == 0 && ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == 0 && ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user);
    } else if (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user);
    } else if (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user);
    } else if (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user);
    }

    //4 Cells connected to the vertex
    else if (er == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez - 1, user);
    } else if (er == 0 && loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez + 1, user);
    } else if (er == 0 && loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez - 1, user);
    } else if (er == 0 && loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez + 1, user);
    } else if (er == N[0] - 1 && loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez - 1, user);
    } else if (er == N[0] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez + 1, user);
    } else if (er == N[0] - 1 && loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez - 1, user);
    } else if (er == N[0] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez + 1, user);
    } else if (ephi == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er - 1, ephi, ez - 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == 0 && loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er - 1, ephi, ez + 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er + 1, ephi, ez - 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == 0 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er + 1, ephi, ez + 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er - 1, ephi, ez - 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er - 1, ephi, ez + 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er + 1, ephi, ez - 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ephi == N[1] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er + 1, ephi, ez + 1, user) + alphavcnomp(er, ephi, ez, user);
    } else if (ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez, user);
    } else if (ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez, user);
    } else if (ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    }

    //8 Cells connected to the vertex
    else if (loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er - 1, ephi - 1, ez, user) + alphavcnomp(er - 1, ephi, ez - 1, user) + alphavcnomp(er, ephi - 1, ez - 1, user) + alphavcnomp(er - 1, ephi - 1, ez - 1, user);
    } else if (loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er + 1, ephi - 1, ez, user) + alphavcnomp(er + 1, ephi, ez - 1, user) + alphavcnomp(er, ephi - 1, ez - 1, user) + alphavcnomp(er + 1, ephi - 1, ez - 1, user);
    } else if (loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er - 1, ephi + 1, ez, user) + alphavcnomp(er - 1, ephi, ez - 1, user) + alphavcnomp(er, ephi + 1, ez - 1, user) + alphavcnomp(er - 1, ephi + 1, ez - 1, user);
    } else if (loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er + 1, ephi + 1, ez, user) + alphavcnomp(er + 1, ephi, ez - 1, user) + alphavcnomp(er, ephi + 1, ez - 1, user) + alphavcnomp(er + 1, ephi + 1, ez - 1, user);
    } else if (loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er - 1, ephi - 1, ez, user) + alphavcnomp(er - 1, ephi, ez + 1, user) + alphavcnomp(er, ephi - 1, ez + 1, user) + alphavcnomp(er - 1, ephi - 1, ez + 1, user);
    } else if (loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er + 1, ephi - 1, ez, user) + alphavcnomp(er + 1, ephi, ez + 1, user) + alphavcnomp(er, ephi - 1, ez + 1, user) + alphavcnomp(er + 1, ephi - 1, ez + 1, user);
    } else if (loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er - 1, ephi + 1, ez, user) + alphavcnomp(er - 1, ephi, ez + 1, user) + alphavcnomp(er, ephi + 1, ez + 1, user) + alphavcnomp(er - 1, ephi + 1, ez + 1, user);
    } else if (loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er + 1, ephi + 1, ez, user) + alphavcnomp(er + 1, ephi, ez + 1, user) + alphavcnomp(er, ephi + 1, ez + 1, user) + alphavcnomp(er + 1, ephi + 1, ez + 1, user);
    }
  } else {
    //2 Cells connected to the vertex
    if (er == 0 && ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == 0 && ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (er == N[0] - 1 && ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (er == 0 && ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == 0 && ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez, user);
    } else if (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    }

    //4 Cells connected to the vertex
    else if (er == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez - 1, user);
    } else if (er == 0 && loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez + 1, user);
    } else if (er == 0 && loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez - 1, user);
    } else if (er == 0 && loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez + 1, user);
    } else if (er == N[0] - 1 && loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez - 1, user);
    } else if (er == N[0] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez + 1, user);
    } else if (er == N[0] - 1 && loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez - 1, user);
    } else if (er == N[0] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez + 1, user);
    } else if (ez == 0 && loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez, user);
    } else if (ez == 0 && loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (ez == 0 && loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez, user);
    } else if (ez == 0 && loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi - 1, ez, user) + alphavcnomp(er, ephi - 1, ez, user);
    } else if (ez == N[2] - 1 && loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi + 1, ez, user) + alphavcnomp(er, ephi + 1, ez, user);
    }

    //8 Cells connected to the vertex
    else if (loc == BACK_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er - 1, ephi - 1, ez, user) + alphavcnomp(er - 1, ephi, ez - 1, user) + alphavcnomp(er, ephi - 1, ez - 1, user) + alphavcnomp(er - 1, ephi - 1, ez - 1, user);
    } else if (loc == BACK_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er + 1, ephi - 1, ez, user) + alphavcnomp(er + 1, ephi, ez - 1, user) + alphavcnomp(er, ephi - 1, ez - 1, user) + alphavcnomp(er + 1, ephi - 1, ez - 1, user);
    } else if (loc == BACK_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er - 1, ephi + 1, ez, user) + alphavcnomp(er - 1, ephi, ez - 1, user) + alphavcnomp(er, ephi + 1, ez - 1, user) + alphavcnomp(er - 1, ephi + 1, ez - 1, user);
    } else if (loc == BACK_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi, ez - 1, user) + alphavcnomp(er + 1, ephi + 1, ez, user) + alphavcnomp(er + 1, ephi, ez - 1, user) + alphavcnomp(er, ephi + 1, ez - 1, user) + alphavcnomp(er + 1, ephi + 1, ez - 1, user);
    } else if (loc == FRONT_DOWN_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er - 1, ephi - 1, ez, user) + alphavcnomp(er - 1, ephi, ez + 1, user) + alphavcnomp(er, ephi - 1, ez + 1, user) + alphavcnomp(er - 1, ephi - 1, ez + 1, user);
    } else if (loc == FRONT_DOWN_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi - 1, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er + 1, ephi - 1, ez, user) + alphavcnomp(er + 1, ephi, ez + 1, user) + alphavcnomp(er, ephi - 1, ez + 1, user) + alphavcnomp(er + 1, ephi - 1, ez + 1, user);
    } else if (loc == FRONT_UP_LEFT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er - 1, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er - 1, ephi + 1, ez, user) + alphavcnomp(er - 1, ephi, ez + 1, user) + alphavcnomp(er, ephi + 1, ez + 1, user) + alphavcnomp(er - 1, ephi + 1, ez + 1, user);
    } else if (loc == FRONT_UP_RIGHT) {
      beta = alphavcnomp(er, ephi, ez, user) + alphavcnomp(er + 1, ephi, ez, user) + alphavcnomp(er, ephi + 1, ez, user) + alphavcnomp(er, ephi, ez + 1, user) + alphavcnomp(er + 1, ephi + 1, ez, user) + alphavcnomp(er + 1, ephi, ez + 1, user) + alphavcnomp(er, ephi + 1, ez + 1, user) + alphavcnomp(er + 1, ephi + 1, ez + 1, user);
    }

  }

  if (!(loc == BACK_DOWN_LEFT || loc == BACK_DOWN_RIGHT || loc == BACK_UP_LEFT || loc == BACK_UP_RIGHT || loc == FRONT_DOWN_LEFT || loc == FRONT_DOWN_RIGHT || loc == FRONT_UP_LEFT || loc == FRONT_UP_RIGHT)) {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in betavnomp function");
  }
  return beta;
}

PetscScalar alphaec2(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etawall) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user->eta0 / (user->etaout)) / 4.0;
      }

      if (er == N[0] - 1 && ephi == -1 && ez == 0 && user -> debug) {
        PetscPrintf(PETSC_COMM_SELF, "user->dataC[er + ephi*N[0] + ez*N[1]*N[0]] = %15.10e\n", user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]);
        PetscPrintf(PETSC_COMM_SELF, "cellvolume = %15.10e\n", cellvolume);
        PetscPrintf(PETSC_COMM_SELF, "user->mu0 = %15.10e\n", user -> mu0);
        PetscPrintf(PETSC_COMM_SELF, "alphaec2(er==N[0]-1, ephi==-1, ez==0) = %15.10e\n", alpha);
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etawall) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user->eta0 / (user->etaout)) / 4.0;
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etawall) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user->eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user->eta0 / (user->etaout)) / 4.0;
      }
    }
  }

  /* alpha = arrCoord[ez][ephi][er][icp[0]] * cellvolume * (user->mu0 / user->eta) / 4.0; */
  else {
    alpha = cellvolume * (user->eta0 / user -> eta) / 4.0; /* Corrected on 03/01/2021 */
  }

  return alpha;
}

PetscScalar alphavc(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  if (user -> ictype == 9) { // You need to adjust this code to the 3-layer resistivity setup: DONE
    // etawall and etaplasma were swapped here: CORRECTED
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawall) / 8.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 8.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 8.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 8.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 8.0;
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawall) / 8.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 8.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 8.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
      alpha = cellvolume * (user -> mu0 / (user->etaVV)) / 8.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 8.0;
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawall) / 8.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 8.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
          alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 8.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
          alpha = cellvolume * (user -> mu0 / user -> etaVV) / 8.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 8.0;
      }
    }
  } else {
    alpha = cellvolume * (user -> mu0 / user -> eta) / 8.0; /* Corrected on 06/07/2021 */
  }

  return alpha;
}

PetscScalar alphavcnomp(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  alpha = cellvolume / 8.0;

  return alpha;
}

PetscScalar alphaec(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawall) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 4.0;
      }

      if (er == N[0] - 1 && ephi == -1 && ez == 0 && user -> debug) {
        PetscPrintf(PETSC_COMM_SELF, "user->dataC[er + ephi*N[0] + ez*N[1]*N[0]] = %15.10e\n", user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]);
        PetscPrintf(PETSC_COMM_SELF, "cellvolume = %15.10e\n", cellvolume);
        PetscPrintf(PETSC_COMM_SELF, "user->mu0 = %15.10e\n", user -> mu0);
        PetscPrintf(PETSC_COMM_SELF, "alphaec(er==N[0]-1, ephi==-1, ez==0) = %15.10e\n", alpha);
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawall) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 4.0;
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawall) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 4.0;
      }
    }
  }

  /* alpha = arrCoord[ez][ephi][er][icp[0]] * cellvolume * (user->mu0 / user->eta) / 4.0; */
  else {
    alpha = cellvolume * (user -> mu0 / user -> eta) / 4.0; /* Corrected on 03/01/2021 */
  }

  return alpha;
}

PetscScalar alphaecnores(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  alpha = cellvolume * (user -> mu0) / 4.0; /* Corrected on 10/25/2021 */

  return alpha;
}

PetscScalar alphaecnomp(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  alpha = cellvolume / 4.0; /* Corrected on 11/2/2021 */

  return alpha;
}

PetscScalar alphaecperp(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawallperp) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 4.0;
      }

      if (er == N[0] - 1 && ephi == -1 && ez == 0 && user -> debug) {
        PetscPrintf(PETSC_COMM_SELF, "user->dataC[er + ephi*N[0] + ez*N[1]*N[0]] = %15.10e\n", user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]);
        PetscPrintf(PETSC_COMM_SELF, "cellvolume = %15.10e\n", cellvolume);
        PetscPrintf(PETSC_COMM_SELF, "user->mu0 = %15.10e\n", user -> mu0);
        PetscPrintf(PETSC_COMM_SELF, "alphaecperp(er==N[0]-1, ephi==-1, ez==0) = %15.10e\n", alpha);
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawallperp) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 4.0;
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawallperp) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 4.0;
      }
    }
  }

  /* alpha = arrCoord[ez][ephi][er][icp[0]] * cellvolume * (user->mu0 / user->eta) / 4.0; */
  else {
    alpha = cellvolume * (user -> mu0 / user -> eta) / 4.0; /* Corrected on 03/01/2021 */
  }

  return alpha;
}

PetscScalar alphaecphi(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / (user -> etawallphi )) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 4.0;
      }

      if (er == N[0] - 1 && ephi == -1 && ez == 0 && user -> debug) {
        PetscPrintf(PETSC_COMM_SELF, "user->dataC[er + ephi*N[0] + ez*N[1]*N[0]] = %15.10e\n", user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]);
        PetscPrintf(PETSC_COMM_SELF, "cellvolume = %15.10e\n", cellvolume);
        PetscPrintf(PETSC_COMM_SELF, "user->mu0 = %15.10e\n", user -> mu0);
        PetscPrintf(PETSC_COMM_SELF, "alphaecphi(er==N[0]-1, ephi==-1, ez==0) = %15.10e\n", alpha);
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / (user -> etawallphi )) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 4.0;
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / (user -> etawallphi )) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 4.0;
      }
    }
  }

  /* alpha = arrCoord[ez][ephi][er][icp[0]] * cellvolume * (user->mu0 / user->eta) / 4.0; */
  else {
    alpha = cellvolume * (user -> mu0 / user -> eta) / 4.0; /* Corrected on 03/01/2021 */
  }

  return alpha;
}

PetscScalar alphaecperp2(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etawallperp) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> eta0 / (user->etaout)) / 4.0;
      }

      if (er == N[0] - 1 && ephi == -1 && ez == 0 && user -> debug) {
        PetscPrintf(PETSC_COMM_SELF, "user->dataC[er + ephi*N[0] + ez*N[1]*N[0]] = %15.10e\n", user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]);
        PetscPrintf(PETSC_COMM_SELF, "cellvolume = %15.10e\n", cellvolume);
        PetscPrintf(PETSC_COMM_SELF, "user->mu0 = %15.10e\n", user -> mu0);
        PetscPrintf(PETSC_COMM_SELF, "alphaecperp2(er==N[0]-1, ephi==-1, ez==0) = %15.10e\n", alpha);
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etawallperp) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> eta0 / (user->etaout)) / 4.0;
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etawallperp) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> eta0 / (user->etaout)) / 4.0;
      }
    }
  }

  /* alpha = arrCoord[ez][ephi][er][icp[0]] * cellvolume * (user->mu0 / user->eta) / 4.0; */
  else {
    alpha = cellvolume * (user -> eta0 / user -> eta) / 4.0; /* Corrected on 03/01/2021 */
  }

  return alpha;
}

PetscScalar alphaecphi2(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / (user -> etawallphi)) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> eta0 / (user->etaout)) / 4.0;
      }

      if (er == N[0] - 1 && ephi == -1 && ez == 0 && user -> debug) {
        PetscPrintf(PETSC_COMM_SELF, "user->dataC[er + ephi*N[0] + ez*N[1]*N[0]] = %15.10e\n", user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]);
        PetscPrintf(PETSC_COMM_SELF, "cellvolume = %15.10e\n", cellvolume);
        PetscPrintf(PETSC_COMM_SELF, "user->mu0 = %15.10e\n", user -> mu0);
        PetscPrintf(PETSC_COMM_SELF, "alphaecphi2(er==N[0]-1, ephi==-1, ez==0) = %15.10e\n", alpha);
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / (user -> etawallphi)) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> eta0 / (user->etaout)) / 4.0;
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / (user -> etawallphi)) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> eta0 / (user->etaout)) / 4.0;
      }
    }
  }

  /* alpha = arrCoord[ez][ephi][er][icp[0]] * cellvolume * (user->mu0 / user->eta) / 4.0; */
  else {
    alpha = cellvolume * (user -> eta0 / user -> eta) / 4.0; /* Corrected on 03/01/2021 */
  }

  return alpha;
}

PetscScalar alphaecphi_isolcell(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        if((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176)){
        //if(arrCoord[ez][ephi][er][icp[2]] > 0.0){
          alpha = cellvolume * (user -> eta0 / (user -> etawallphi_isol_cell)) / 4.0;
        }
        else{
          alpha = cellvolume * (user -> eta0 / (user -> etawallphi)) / 4.0;
        }
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> eta0 / (user->etaout)) / 4.0;
      }

      if (er == N[0] - 1 && ephi == -1 && ez == 0 && user -> debug) {
        PetscPrintf(PETSC_COMM_SELF, "user->dataC[er + ephi*N[0] + ez*N[1]*N[0]] = %15.10e\n", user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]);
        PetscPrintf(PETSC_COMM_SELF, "cellvolume = %15.10e\n", cellvolume);
        PetscPrintf(PETSC_COMM_SELF, "user->mu0 = %15.10e\n", user -> mu0);
        PetscPrintf(PETSC_COMM_SELF, "alphaecphi(er==N[0]-1, ephi==-1, ez==0) = %15.10e\n", alpha);
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        if((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176)){
        //if(arrCoord[ez][ephi][er][icp[2]] > 0.0){
          alpha = cellvolume * (user -> eta0 / (user -> etawallphi_isol_cell)) / 4.0;
        }
        else{
          alpha = cellvolume * (user -> eta0 / (user -> etawallphi)) / 4.0;
        }
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> eta0 / (user->etaout)) / 4.0;
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaplasma) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        if((er == 28 && ez == 183) || (er == 13 && ez == 176) || (er == 47 && ez == 176)){
        //if(arrCoord[ez][ephi][er][icp[2]] > 0.0){
          alpha = cellvolume * (user -> eta0 / (user -> etawallphi_isol_cell)) / 4.0;
        }
        else{
          alpha = cellvolume * (user -> eta0 / (user -> etawallphi)) / 4.0;
        }
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etasepwal) / 4.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> eta0 / user -> etaVV) / 4.0;
      } else {
        alpha = cellvolume * (user -> eta0 / (user->etaout)) / 4.0;
      }
    }
  }

  /* alpha = arrCoord[ez][ephi][er][icp[0]] * cellvolume * (user->mu0 / user->eta) / 4.0; */
  else {
    alpha = cellvolume * (user -> eta0 / user -> eta) / 4.0; /* Corrected on 03/01/2021 */
  }

  return alpha;
}

PetscScalar alphafc(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  /* alpha = arrCoord[ez][ephi][er][icp[0]] * cellvolume / 4.0; */
  alpha = cellvolume / 2.0; /* Corrected on 03/01/2021 */
  return alpha;
}

PetscScalar alphafc_wmp(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 2.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawall) / 2.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 2.0;
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 2.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 2.0;
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 2.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawall) / 2.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 2.0;
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 2.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 2.0;
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaplasma) / 2.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etawall) / 2.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etasepwal) / 2.0;
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = cellvolume * (user -> mu0 / user -> etaVV) / 2.0;
      } else {
        alpha = cellvolume * (user -> mu0 / (user->etaout)) / 2.0;
      }
    }
  } else {
    alpha = cellvolume * (user -> mu0 / user -> eta) / 2.0;
  }

  return alpha;
}

PetscScalar alphac(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
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
  PetscScalar alpha, cellvolume;

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

  if (ephi == -1 || ephi == N[1] - 1 || ephi == N[1]) {
    cellvolume = user -> dphi * PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  } else {
    cellvolume = PetscAbsReal(arrCoord[ez][ephi][er][icBzp[2]] - arrCoord[ez][ephi][er][icBzm[2]]) *
      PetscAbsReal(arrCoord[ez][ephi][er][icBphip[1]] - arrCoord[ez][ephi][er][icBphim[1]]) * PetscAbsReal(PetscSqr(arrCoord[ez][ephi][er][icBrp[0]]) - PetscSqr(arrCoord[ez][ephi][er][icBrm[0]])) / 2.0; /* INT_c(r dphi dr dz) */
  }

  alpha = cellvolume;
  return alpha;
}

PetscScalar rese(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = resec(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = 0.5 * (resec(er - 1, ephi, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = 0.5 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = 0.5 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er + 1, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = 0.5 * (resec(er - 1, ephi, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er + 1, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = 0.25 * (resec(er - 1, ephi - 1, ez, user) + resec(er - 1, ephi, ez, user) + resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = 0.25 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user) + resec(er + 1, ephi - 1, ez, user) + resec(er + 1, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = 0.25 * (resec(er - 1, ephi + 1, ez, user) + resec(er - 1, ephi, ez, user) + resec(er, ephi + 1, ez, user) + resec(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else {
        beta = 0.25 * (resec(er, ephi + 1, ez, user) + resec(er, ephi, ez, user) + resec(er + 1, ephi + 1, ez, user) + resec(er + 1, ephi, ez, user));
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = 0.5 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = 0.5 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = 0.25 * (resec(er - 1, ephi - 1, ez, user) + resec(er - 1, ephi, ez, user) + resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = 0.25 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user) + resec(er + 1, ephi - 1, ez, user) + resec(er + 1, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = 0.25 * (resec(er - 1, ephi + 1, ez, user) + resec(er - 1, ephi, ez, user) + resec(er, ephi + 1, ez, user) + resec(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else {
        beta = 0.25 * (resec(er, ephi + 1, ez, user) + resec(er, ephi, ez, user) + resec(er + 1, ephi + 1, ez, user) + resec(er + 1, ephi, ez, user));
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = resec(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = 0.5 * (resec(er - 1, ephi, ez, user) + resec(er, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = 0.5 * (resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = 0.5 * (resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = 0.5 * (resec(er, ephi, ez, user) + resec(er + 1, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi, ez + 1, user));
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = 0.5 * (resec(er - 1, ephi, ez, user) + resec(er, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi, ez + 1, user));
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = 0.5 * (resec(er, ephi, ez, user) + resec(er + 1, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = 0.25 * (resec(er - 1, ephi, ez - 1, user) + resec(er - 1, ephi, ez, user) + resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user));
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = 0.25 * (resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user) + resec(er + 1, ephi, ez - 1, user) + resec(er + 1, ephi, ez, user));
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = 0.25 * (resec(er - 1, ephi, ez + 1, user) + resec(er - 1, ephi, ez, user) + resec(er, ephi, ez + 1, user) + resec(er, ephi, ez, user));
    } /* 4 cells sharing the edge */
    else {
      beta = 0.25 * (resec(er, ephi, ez + 1, user) + resec(er, ephi, ez, user) + resec(er + 1, ephi, ez + 1, user) + resec(er + 1, ephi, ez, user));
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = resec(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = 0.5 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = 0.5 * (resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = 0.5 * (resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi, ez + 1, user));
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = 0.5 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi, ez + 1, user));
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = 0.25 * (resec(er, ephi - 1, ez - 1, user) + resec(er, ephi - 1, ez, user) + resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = 0.25 * (resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user) + resec(er, ephi + 1, ez - 1, user) + resec(er, ephi + 1, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = 0.25 * (resec(er, ephi - 1, ez + 1, user) + resec(er, ephi - 1, ez, user) + resec(er, ephi, ez + 1, user) + resec(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else {
        beta = 0.25 * (resec(er, ephi, ez + 1, user) + resec(er, ephi, ez, user) + resec(er, ephi + 1, ez + 1, user) + resec(er, ephi + 1, ez, user));
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = 0.5 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = 0.5 * (resec(er, ephi - 1, ez, user) + resec(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = 0.5 * (resec(er, ephi, ez, user) + resec(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = 0.25 * (resec(er, ephi - 1, ez - 1, user) + resec(er, ephi - 1, ez, user) + resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = 0.25 * (resec(er, ephi, ez - 1, user) + resec(er, ephi, ez, user) + resec(er, ephi + 1, ez - 1, user) + resec(er, ephi + 1, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = 0.25 * (resec(er, ephi - 1, ez + 1, user) + resec(er, ephi - 1, ez, user) + resec(er, ephi, ez + 1, user) + resec(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else {
        beta = 0.25 * (resec(er, ephi, ez + 1, user) + resec(er, ephi, ez, user) + resec(er, ephi + 1, ez + 1, user) + resec(er, ephi + 1, ez, user));
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in rese function");
  }
  return beta;
}

PetscScalar resec(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt N[3];
  PetscScalar alpha;

  DMStagGetGlobalSizes(user -> coorda, & N[0], & N[1], & N[2]);

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = (1.0 / user -> etaplasma);
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = (1.0 / user -> etawall);
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = (1.0 / user -> etasepwal);
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = (1.0 / user -> etaVV);
      } else {
        alpha = (1.0 / (user->etaout));
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = (1.0 / user -> etaplasma);
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = (1.0 / user -> etawall);
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = (1.0 / user -> etasepwal);
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = (1.0 / user -> etaVV);
      } else {
        alpha = (1.0 / (user->etaout));
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = (1.0 / user -> etaplasma);
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = (1.0 / user -> etawall);
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = (1.0 / user -> etasepwal);
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = (1.0 / user -> etaVV);
      } else {
        alpha = (1.0 / (user->etaout));
      }
    }
  } else {
    alpha = (1.0 / user -> eta);
  }

  return alpha;
}

PetscScalar conduc(PetscInt er, PetscInt ephi, PetscInt ez, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt N[3];
  PetscScalar alpha;

  DMStagGetGlobalSizes(user -> coorda, & N[0], & N[1], & N[2]);

  if (user -> ictype == 9) {
    if (ephi > -1 && ephi < N[1]) {
      if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = (user -> etaplasma);
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = (user -> etawall);
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = (user -> etasepwal);
      } else if (fabs(user -> dataC[er + ephi * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = (user -> etaVV);
      } else {
        alpha = ((user->etaout));
      }
    } else if (ephi == -1) {
      if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = (user -> etaplasma);
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = (user -> etawall);
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = (user -> etasepwal);
      } else if (fabs(user -> dataC[er + (N[1] - 1) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = (user -> etaVV);
      } else {
        alpha = ((user->etaout));
      }
    } else { // ephi == N[1]
      if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 1.0) < 1e-12) {
        alpha = (user -> etaplasma);
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]]) < 1e-12) {
        alpha = (user -> etawall);
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] - 2.0) < 1e-12) {
        alpha = (user -> etasepwal);
      } else if (fabs(user -> dataC[er + (0) * N[0] + ez * N[1] * N[0]] + 1.0) < 1e-12) {
        alpha = (user -> etaVV);
      } else {
        alpha = ((user->etaout));
      }
    }
  } else {
    alpha = (user -> eta);
  }

  return alpha;
}

PetscScalar condu(PetscInt er, PetscInt ephi, PetscInt ez, DMStagStencilLocation loc, void * ptr) {
  User * user = (User * ) ptr;
  PetscInt startr, startphi, startz, nr, nphi, nz, N[3];
  PetscScalar beta = 0;
  DM coordDA = user -> coorda;

  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetGlobalSizes(coordDA, & N[0], & N[1], & N[2]);
  /* DEBUG PRINT*/
  /*PetscPrintf(PETSC_COMM_WORLD,"(er, ephi, ez) = (%d,%d,%d)\n",(int)er,(int)ephi,(int)ez);*/
  DMStagGetCorners(coordDA, & startr, & startphi, & startz, & nr, & nphi, & nz, NULL, NULL, NULL);
  /*if (!(startz <= ez && ez<startz+nz && startphi <= ephi && ephi<startphi+nphi && startr <= er && er<startr+nr))  SETERRQ(PetscObjectComm((PetscObject)coordDA),PETSC_ERR_ARG_SIZ,"The cell indices exceed the local range");*/
  /* Edges in z direction */
  if (loc == DOWN_LEFT || loc == DOWN_RIGHT || loc == UP_LEFT || loc == UP_RIGHT) {
    if (!(user -> phibtype)) {
      if ((er == 0 && ephi == 0 && loc == DOWN_LEFT) || (er == N[0] - 1 && ephi == 0 && loc == DOWN_RIGHT) || (er == 0 && ephi == N[1] - 1 && loc == UP_LEFT) || (er == N[0] - 1 && ephi == N[1] - 1 && loc == UP_RIGHT)) {
        beta = conduc(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ephi == 0 && loc == DOWN_LEFT) {
        beta = 0.5 * (conduc(er - 1, ephi, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == DOWN_LEFT) {
        beta = 0.5 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = 0.5 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == DOWN_RIGHT) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er + 1, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_LEFT) {
        beta = 0.5 * (conduc(er - 1, ephi, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == UP_RIGHT) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er + 1, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = 0.25 * (conduc(er - 1, ephi - 1, ez, user) + conduc(er - 1, ephi, ez, user) + conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = 0.25 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user) + conduc(er + 1, ephi - 1, ez, user) + conduc(er + 1, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = 0.25 * (conduc(er - 1, ephi + 1, ez, user) + conduc(er - 1, ephi, ez, user) + conduc(er, ephi + 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else {
        beta = 0.25 * (conduc(er, ephi + 1, ez, user) + conduc(er, ephi, ez, user) + conduc(er + 1, ephi + 1, ez, user) + conduc(er + 1, ephi, ez, user));
      } /* 4 cells sharing the edge */
    } else {
      if (er == 0 && loc == DOWN_LEFT) {
        beta = 0.5 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == DOWN_RIGHT) {
        beta = 0.5 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == 0 && loc == UP_LEFT) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (er == N[0] - 1 && loc == UP_RIGHT) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (loc == DOWN_LEFT) {
        beta = 0.25 * (conduc(er - 1, ephi - 1, ez, user) + conduc(er - 1, ephi, ez, user) + conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == DOWN_RIGHT) {
        beta = 0.25 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user) + conduc(er + 1, ephi - 1, ez, user) + conduc(er + 1, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == UP_LEFT) {
        beta = 0.25 * (conduc(er - 1, ephi + 1, ez, user) + conduc(er - 1, ephi, ez, user) + conduc(er, ephi + 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else {
        beta = 0.25 * (conduc(er, ephi + 1, ez, user) + conduc(er, ephi, ez, user) + conduc(er + 1, ephi + 1, ez, user) + conduc(er + 1, ephi, ez, user));
      } /* 4 cells sharing the edge */
    }
  }
  /* Edges in phi direction */
  else if (loc == BACK_LEFT || loc == BACK_RIGHT || loc == FRONT_LEFT || loc == FRONT_RIGHT) {
    if ((er == 0 && ez == 0 && loc == BACK_LEFT) || (er == N[0] - 1 && ez == 0 && loc == BACK_RIGHT) || (er == 0 && ez == N[2] - 1 && loc == FRONT_LEFT) || (er == N[0] - 1 && ez == N[2] - 1 && loc == FRONT_RIGHT)) {
      beta = conduc(er, ephi, ez, user);
    } /* Only 1 cell connected to the edge */
    else if (ez == 0 && loc == BACK_LEFT) {
      beta = 0.5 * (conduc(er - 1, ephi, ez, user) + conduc(er, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == BACK_LEFT) {
      beta = 0.5 * (conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == BACK_RIGHT) {
      beta = 0.5 * (conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (ez == 0 && loc == BACK_RIGHT) {
      beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er + 1, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (er == 0 && loc == FRONT_LEFT) {
      beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi, ez + 1, user));
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_LEFT) {
      beta = 0.5 * (conduc(er - 1, ephi, ez, user) + conduc(er, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (er == N[0] - 1 && loc == FRONT_RIGHT) {
      beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi, ez + 1, user));
    } /* 2 cells sharing the edge */
    else if (ez == N[2] - 1 && loc == FRONT_RIGHT) {
      beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er + 1, ephi, ez, user));
    } /* 2 cells sharing the edge */
    else if (loc == BACK_LEFT) {
      beta = 0.25 * (conduc(er - 1, ephi, ez - 1, user) + conduc(er - 1, ephi, ez, user) + conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user));
    } /* 4 cells sharing the edge */
    else if (loc == BACK_RIGHT) {
      beta = 0.25 * (conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user) + conduc(er + 1, ephi, ez - 1, user) + conduc(er + 1, ephi, ez, user));
    } /* 4 cells sharing the edge */
    else if (loc == FRONT_LEFT) {
      beta = 0.25 * (conduc(er - 1, ephi, ez + 1, user) + conduc(er - 1, ephi, ez, user) + conduc(er, ephi, ez + 1, user) + conduc(er, ephi, ez, user));
    } /* 4 cells sharing the edge */
    else {
      beta = 0.25 * (conduc(er, ephi, ez + 1, user) + conduc(er, ephi, ez, user) + conduc(er + 1, ephi, ez + 1, user) + conduc(er + 1, ephi, ez, user));
    } /* 4 cells sharing the edge */
  }
  /* Edges in r direction */
  else if (loc == BACK_DOWN || loc == BACK_UP || loc == FRONT_DOWN || loc == FRONT_UP) {
    if (!(user -> phibtype)) {
      if ((ephi == 0 && ez == 0 && loc == BACK_DOWN) || (ephi == N[1] - 1 && ez == 0 && loc == BACK_UP) || (ephi == 0 && ez == N[2] - 1 && loc == FRONT_DOWN) || (ephi == N[1] - 1 && ez == N[2] - 1 && loc == FRONT_UP)) {
        beta = conduc(er, ephi, ez, user);
      } /* Only 1 cell connected to the edge */
      else if (ez == 0 && loc == BACK_DOWN) {
        beta = 0.5 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == BACK_DOWN) {
        beta = 0.5 * (conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == BACK_UP) {
        beta = 0.5 * (conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == 0 && loc == FRONT_DOWN) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi, ez + 1, user));
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = 0.5 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ephi == N[1] - 1 && loc == FRONT_UP) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi, ez + 1, user));
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = 0.25 * (conduc(er, ephi - 1, ez - 1, user) + conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = 0.25 * (conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez - 1, user) + conduc(er, ephi + 1, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = 0.25 * (conduc(er, ephi - 1, ez + 1, user) + conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez + 1, user) + conduc(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else {
        beta = 0.25 * (conduc(er, ephi, ez + 1, user) + conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez + 1, user) + conduc(er, ephi + 1, ez, user));
      } /* 4 cells sharing the edge */
    } else {
      if (ez == 0 && loc == BACK_DOWN) {
        beta = 0.5 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ez == 0 && loc == BACK_UP) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_DOWN) {
        beta = 0.5 * (conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez, user));
      } /* 2 cells sharing the edge */
      else if (ez == N[2] - 1 && loc == FRONT_UP) {
        beta = 0.5 * (conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez, user));
      } /* 2 cells sharing the edge */
      else if (loc == BACK_DOWN) {
        beta = 0.25 * (conduc(er, ephi - 1, ez - 1, user) + conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == BACK_UP) {
        beta = 0.25 * (conduc(er, ephi, ez - 1, user) + conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez - 1, user) + conduc(er, ephi + 1, ez, user));
      } /* 4 cells sharing the edge */
      else if (loc == FRONT_DOWN) {
        beta = 0.25 * (conduc(er, ephi - 1, ez + 1, user) + conduc(er, ephi - 1, ez, user) + conduc(er, ephi, ez + 1, user) + conduc(er, ephi, ez, user));
      } /* 4 cells sharing the edge */
      else {
        beta = 0.25 * (conduc(er, ephi, ez + 1, user) + conduc(er, ephi, ez, user) + conduc(er, ephi + 1, ez + 1, user) + conduc(er, ephi + 1, ez, user));
      } /* 4 cells sharing the edge */
    }
  } else {
    /* DEBUG PRINT*/
    PetscPrintf(PETSC_COMM_WORLD, "Location : %d\n", (int) loc);
    SETERRQ(PetscObjectComm((PetscObject) coordDA), PETSC_ERR_ARG_SIZ, "Incorrect DMStagStencilLocation input in condu function");
  }
  return beta;
}

