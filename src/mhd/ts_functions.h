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

#if !defined(TS_FUNCTIONS_H)
#define TS_FUNCTIONS_H

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

PetscErrorCode FormIJacobian_BImplicit(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*); /* This routine computes the analytical Jacobian. Not implemented yet */
PetscErrorCode FormIFunction_Vperp_viscosity(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the constraints: [(\nabla x B) x B = - lambda (\nabla^2 V)].e_{R/Z}; V . B = 0 */
PetscErrorCode FormIFunction_InitializeEP(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes EP and tau */
PetscErrorCode FormIFunction_InitializeEP_halo(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes EP and tau for the halo current simulation*/
PetscErrorCode FormIFunction_newequilibrium_Vperp(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes V: [(\nabla x B) x B = -(1/Re) (\nabla^2 V) + n_i * dV/dt].e_{R/Z}; V . B = 0 */
PetscErrorCode FormRHSFunction_BImplicit(TS,PetscReal,Vec,Vec,void*); /* This RHSFunction contains only zero entries */
PetscErrorCode FormInitialSolution(TS,Vec,void*); /* This routine sets up the initial solution X_0 from input data files containing the B field values on face centers and the levelset function on cell centers */
PetscErrorCode FormExactSolution(PetscReal,TS,Vec*,void*); /* This routine sets up the vector which will be used to set the boundary conditions inside FormIFunction_DampingV */
PetscErrorCode Monitor(TS,PetscInt,PetscReal,Vec,void*);
PetscErrorCode FormDummyIJacobian4(TS,Vec,Vec,PetscReal,Mat,Mat,void*); /* When used with PCFieldSplitSetDetectSaddlePoint, this dummy jacobian gives a field splitting where the electrostatic potential Phi is first and n, tau, B and V fields are last */
PetscErrorCode SampleShellPCSetUp(PC); /* This routine sets up a Shell Preconditioner that has the same effect as a 2-field fieldsplit preconditioner where the velocity is split from the remaining unknowns */
PetscErrorCode SampleShellPCApply(PC,Vec,Vec); /* This routine applies a Shell Preconditioner that has the same effect as a 2-field fieldsplit preconditioner where the velocity is split from the remaining unknowns */
PetscErrorCode SampleShellPCDestroy(PC); /* This routine destroys a Shell Preconditioner that has the same effect as a 2-field fieldsplit preconditioner where the velocity is split from the remaining unknowns */
PetscErrorCode ReadInitialData(PetscReal**,PetscInt*,const char*); /* This routine reads data from an input file and stores these in an array. The length of the output array is also computed. */
PetscErrorCode FormInitialSolution_psi(TS,Vec,void*); /* This routine sets up the initial solution X_0 from input data files containing G(psi) values on face centers, psi values on edge centers and the levelset function on cell centers */

PetscErrorCode stag_vec_io(User *user, PetscViewer viewer, Vec X, PetscBool load);

#endif /* defined(TS_FUNCTIONS_H) */
