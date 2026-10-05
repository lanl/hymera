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


PetscErrorCode DumpSolution_Cell(TS,PetscInt,Vec,void*);
PetscErrorCode DumpError(TS,PetscInt,Vec,void*);
PetscErrorCode DumpDivergence(TS,DM,PetscInt,Vec,void*);
PetscErrorCode DumpLevelSet(TS,void*);
PetscErrorCode SaveIntermediateSolution(TS,PetscInt,PetscReal,Vec,void*);
PetscErrorCode ComputeCurrent(TS,Vec,void*);

#endif /* defined(MONITOR_FUNCTIONS_H) */
