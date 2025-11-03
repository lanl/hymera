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

#if !defined(MASS_MATRIX_COEFFICIENTS_H)
#define MASS_MATRIX_COEFFICIENTS_H

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdraw.h>
#include <petscksp.h>
#include <petscdmstag.h>
#include <petscsys.h>
#include <petscvec.h>
#include <mpi.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <petsc/private/dmstagimpl.h>

PetscScalar betaf(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betaf_wmp(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betae(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betaenores(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betaenomp(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betaeperp(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betaephi(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betae2(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betaeperp2(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betaephi2(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betaephi_isolcell(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betav(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar betavnomp(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar alphafc(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphafc_wmp(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphaec(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphaecnores(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphaecnomp(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphaecperp(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphaecphi(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphaec2(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphaecperp2(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphaecphi2(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphaecphi_isolcell(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphavc(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphavcnomp(PetscInt,PetscInt,PetscInt,void*);
PetscScalar alphac(PetscInt,PetscInt,PetscInt,void*);
PetscScalar rese(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar resec(PetscInt,PetscInt,PetscInt,void*);
PetscScalar condu(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscScalar conduc(PetscInt,PetscInt,PetscInt,void*);

#endif /* defined(MASS_MATRIX_COEFFICIENTS_H) */
