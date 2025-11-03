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

#if !defined(MIMETIC_OPERATORS_H)
#define MIMETIC_OPERATORS_H

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

PetscErrorCode FormMaterialPropertiesMatrix(TS,Vec,void*);
PetscErrorCode FormFaceMassMatrix(TS,Vec,void*);
PetscErrorCode FormEdgeMassMatrix(TS,Vec,void*);
PetscErrorCode FormDiscreteDivergence(TS,DM,Mat,Vec,Vec,void*);
PetscErrorCode FormGradDerivedDivergence(TS,Mat,void*);
PetscErrorCode FormDerivedGradDivergence(TS,Mat,void*);
PetscErrorCode FormDerivedGradEtaDivergence(TS,Mat,void*);
PetscErrorCode FormDerivedDivergence(TS,Mat,void*);
PetscErrorCode FormDerivedGradient(TS,Mat,Mat,void*);
PetscErrorCode FormDiscreteGradient(TS,Mat,void*);
PetscErrorCode FormPrimaryCurl(TS,Vec,Vec,void*);
PetscErrorCode FormDerivedCurl(TS,Vec,Vec,void*);
PetscErrorCode FormDerivedCurlnores(TS,Vec,Vec,void*);
PetscErrorCode FormDerivedCurlnomp(TS,Vec,Vec,void*);
PetscErrorCode FormSourceTermPotential(TS,PetscReal,Vec,void*);
PetscErrorCode FormDerivedCurlExt(Vec,Vec,void*);
PetscErrorCode FormDiscreteGradientEP_tilde(TS,Vec,Vec,void*);
PetscErrorCode FormDiscreteGradientEP(TS,Mat,Vec,Vec,void*);
PetscErrorCode FormDiscreteGradientEP_noMat(TS,Vec,Vec,void*);
PetscErrorCode FormDiscreteGradientVectorField(TS,Vec,Vec,Vec,Vec,void*);
PetscErrorCode ApplyDerivedDivergence(TS,Vec,Vec,void*);
PetscErrorCode ApplyDeltastar(TS,Vec,Vec,void*);
PetscErrorCode ApplyDeltastar2(TS,Vec,Vec,void*);
PetscErrorCode ApplyVectorLaplacian(TS,Vec,Vec,void*);
PetscErrorCode FormElectricField(TS,Vec,Vec,void*); /* This routine computes the electric field E and sets it on the edge values of output vector F. It uses tau field and electrostatic potential from input X: E = tau + \nabla (EP). */

#endif /* defined(MIMETIC_OPERATORS_H) */
