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

#if !defined(GEOMETRY_H)
#define GEOMETRY_H

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

PetscScalar cyldistance(PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar,PetscScalar);
PetscScalar surface(PetscInt,PetscInt,PetscInt,DMStagStencilLocation,void*);
PetscErrorCode ComputeIsEBoundary(TS,IS*,void*);
PetscErrorCode SaveSolution(TS,Vec,void*);
PetscErrorCode SaveCoordinates(TS,void*);
PetscErrorCode CellToVertexProjectionScalar(TS,Vec,Vec,void*);
PetscErrorCode EdgeToCellReconstruction_r(TS,Vec,Vec,void*);
PetscErrorCode EdgeToCellReconstruction_phi(TS,Vec,Vec,void*);
PetscErrorCode EdgeToCellReconstruction_z(TS,Vec,Vec,void*);
PetscErrorCode VertexToEdgeReconstruction_scalar(TS,Vec,Vec,void*);
PetscErrorCode VertexToEdgeReconstruction(TS,Vec,Vec,void*);
PetscErrorCode VertexToFaceReconstruction(TS,Vec,Vec,void*);
PetscErrorCode FaceToVertexProjection(TS,Vec,Vec,void*);
PetscErrorCode EdgeToVertexProjection(TS,Vec,Vec,void*);
PetscErrorCode CellToFaceProjection(TS,Vec,Vec,void*);
PetscErrorCode VertexCrossProduct(TS,Vec,Vec,Vec,void*);

PetscErrorCode getBArray(TS ts, Vec X, PetscScalar *gf_B, void *ptr, int derivative);
PetscErrorCode getVArray(TS ts, Vec X, PetscScalar *gf_V, void *ptr);
PetscErrorCode getEJArray(TS ts, Vec X, PetscScalar *ge_ER, PetscScalar *ge_EP, PetscScalar *ge_EZ, void *ptr, int code);

#endif /* defined(GEOMETRY_H) */
