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

#if !defined(MONITOR_FUNCTIONS_H)
#define MONITOR_FUNCTIONS_H

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
#include <petscsnes.h>
#include <petsc/private/tsimpl.h>


PetscErrorCode DumpVelocity_Cell(TS,PetscInt,Vec,char*,void*);
PetscErrorCode DumpSolution_Cell(TS,PetscInt,Vec,void*);
PetscErrorCode DumpSolution(TS,PetscInt,Vec,void*);
PetscErrorCode DumpError(TS,PetscInt,Vec,void*);
PetscErrorCode DumpDivergence(TS,DM,PetscInt,Vec,void*);
PetscErrorCode DumpLevelSet(TS,void*);
PetscErrorCode SaveIntermediateSolution(TS,PetscInt,PetscReal,Vec,void*);
PetscErrorCode ComputeCurrent(TS,Vec,void*);
PetscErrorCode Dump1stVertexField(TS,PetscInt,Vec,void*);
PetscErrorCode DumpEdgeField(TS,PetscInt,Vec,void*);
PetscErrorCode multiplybyR(PetscScalar*,const PetscScalar*,const int);
PetscScalar*   createHermiteFD(int*,int*,const PetscScalar*,const PetscScalar*, const int,const int,const int,const int,const int);
PetscErrorCode getHermiteDataFD(const PetscScalar*,const int,const int,const int,PetscScalar*,const int,const int,const int,const int);
PetscErrorCode interpolateHermite(PetscScalar*,PetscScalar*,int,int,void*);
PetscErrorCode interpolate2D(PetscScalar*,PetscScalar*,int,void*);
PetscErrorCode interpolatey(const int,const int,PetscScalar*,PetscScalar*,PetscScalar*,void*);
PetscErrorCode interpolatex(const int,const int,PetscScalar*,PetscScalar*,PetscScalar*,void*);
PetscErrorCode integrateHermite_1D_r(PetscScalar*,PetscScalar*,int,int,const PetscScalar,const PetscScalar);
PetscErrorCode integrateHermite_2D_z(PetscScalar*,PetscScalar*,int,int,const PetscScalar,const PetscScalar);
PetscScalar evalHermite1D(const PetscScalar*,const PetscScalar*,const PetscScalar,const PetscScalar,const PetscScalar*,const int,const int,const int);
PetscScalar evalHermite2D(const PetscScalar*,const PetscScalar*,const PetscScalar, const PetscScalar,const PetscScalar,const PetscScalar,const PetscScalar*,const int,const int,const int,const int,const int,const int);
tLocalCoordinate getLocalCoordinate(const PetscScalar*,const PetscScalar*,const  PetscScalar,const  PetscScalar,const PetscScalar,const PetscScalar);
PetscErrorCode evalPsi(PetscScalar*,const PetscScalar*,const PetscScalar*,int,tHermiteDivFreeFields);
PetscErrorCode DumpPsi_Cell(TS,PetscInt,PetscScalar*,void*);
PetscErrorCode TSAdaptChoose_user(TSAdapt,TS,PetscReal,PetscInt*,PetscReal*,PetscBool*,PetscReal*,PetscReal*,PetscReal*);

#endif /* defined(MONITOR_FUNCTIONS_H) */
