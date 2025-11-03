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
PetscErrorCode EBoundaryAdjusters(TS,Mat,Mat,void*);
PetscErrorCode ComputeIsEBoundary(TS,IS*,void*);
PetscErrorCode ComputeIsBBoundary(TS,IS*,void*);
PetscErrorCode ComputeIsCBoundary(TS,IS*,void*);
PetscErrorCode SaveSolution(TS,Vec,void*);
PetscErrorCode SaveCoordinates(TS,void*);
PetscErrorCode ReadDataInVec(DM,Vec,Vec,Vec,Vec,void*);
PetscErrorCode CellToVertexProjectionScalar(TS,Vec,Vec,void*);
PetscErrorCode CellToVertexProjectionVector(TS,Vec,Vec,void*);
PetscErrorCode EdgeToCellReconstruction_r(TS,Vec,Vec,void*);
PetscErrorCode EdgeToCellReconstruction_phi(TS,Vec,Vec,void*);
PetscErrorCode EdgeToCellReconstruction_z(TS,Vec,Vec,void*);
PetscErrorCode VertexToCellReconstruction(TS,Vec,Vec,void*);
PetscErrorCode VertexToEdgeReconstruction_scalar(TS,Vec,Vec,void*);
PetscErrorCode VertexToEdgeReconstruction(TS,Vec,Vec,void*);
PetscErrorCode VertexToEdgeReconstructionMat(TS,Mat,void*);
PetscErrorCode VertexToFaceReconstruction(TS,Vec,Vec,void*);
PetscErrorCode VertexToFaceReconstructionMat(TS,Mat,void*);
PetscErrorCode EdgeToCellReconstructionMat(TS,Mat,void*);
PetscErrorCode FaceToCellReconstructionMat(TS,Mat,void*);
PetscErrorCode FaceToVertexProjection(TS,Vec,Vec,void*);
PetscErrorCode FaceToVertexProjectionMat(TS,Mat,void*);
PetscErrorCode EdgeToVertexProjection_Original(TS,Vec,Vec,void*);
PetscErrorCode EdgeToVertexProjection(TS,Vec,Vec,void*);
PetscErrorCode EdgeToVertexProjectionMat(TS,Mat,void*);
PetscErrorCode CellToFaceProjectionMat(TS,Mat,void*);
PetscErrorCode CellToFaceProjection(TS,Vec,Vec,void*);
PetscErrorCode VertexCrossProduct(TS,Vec,Vec,Vec,void*);
PetscErrorCode Computepsi(TS,PetscInt,Vec,void*);
PetscErrorCode FromPetscVecToArray_EfieldCell(TS,Vec,PetscScalar*,PetscScalar*,PetscScalar*,void*);
PetscErrorCode FromPetscVecToArray(TS,Vec,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,void*);
PetscErrorCode CellCoordArrays(TS,PetscScalar*,PetscScalar*,void*);
PetscErrorCode ScatterTest(TS,void*);

int isInDomain(const double *,const double*,void*);
PetscErrorCode PushParticles(TS,Vec,Vec,void*);
// PetscErrorCode PlotQandPoicare(TS,PetscInt,Vec,void*);

#endif /* defined(GEOMETRY_H) */
