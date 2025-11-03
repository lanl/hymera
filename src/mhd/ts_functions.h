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
PetscErrorCode FormIFunction_Vperp_viscosity_halo(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the \eta(\nabla x B) in the tau and EP constraints replaced by: \eta_\phi (\nabla x B)_phi + \eta_perp [ (\nabla x B)_perp - (\nabla x B)_\phi]; inside the wall */
PetscErrorCode FormIFunction_Vperp_viscosity_halo_isolcell(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the \eta(\nabla x B) in the tau and EP constraints replaced by: \eta_\phi (\nabla x B)_phi + \eta_perp [ (\nabla x B) - (\nabla x B)_\phi]; inside the wall and \eta_\phi is lower on three separate cells */
PetscErrorCode FormIFunction_Inertia_viscosity(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the constraints: (\nabla x B) x B = n_i m_i (V . \nabla) V - lambda (\nabla^2 V); */
PetscErrorCode FormIFunction_Vperp_viscosity(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the constraints: [(\nabla x B) x B = - lambda (\nabla^2 V)].e_{R/Z}; V . B = 0 */
PetscErrorCode FormIFunction_purediffusionVV(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the constraints: tau = (\nabla x B); dB/dt + \nabla x tau = 0; */
PetscErrorCode FormIFunction_purediffusion(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the constraints: tau = (\nabla x B); dB/dt + \nabla x tau = 0; */
PetscErrorCode FormIFunction_Inertia_V(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the constraints: [(\nabla x B) x B].e_{R/Z} = [n_i(t=0) m_i (V . \nabla) V].e_{R/Z}; [V . \nabla V] . B = 0 */
PetscErrorCode FormIFunction_Inertia_V_ni(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the constraints: (\nabla x B) x B = n_i m_i (V . \nabla) V; */
PetscErrorCode FormIFunction_Inertia(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the constraints: [(\nabla x B) x B].e_{R/Z} = [n_i m_i (V . \nabla) V - lambda (\nabla^2 V)].e_{R/Z}; V . B = 0 */
PetscErrorCode FormIFunction2(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction treats tau and Phi as auxiliary variables : Phi = 0; tau = - V x B + (eta/mu_0) (\nabla x B) */
PetscErrorCode FormIFunction(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction treats tau and Phi as auxiliary variables: Phi = 0; tau = 0 */
PetscErrorCode FormIFunction_BImplicit2(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the following vertex equation: ni dV/dt - Curl(B) x B as a left side */
PetscErrorCode FormIFunction_Initializepsi(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes the poloidal flux function psi: (1/r) Delta^* = mu_0 J0 */
PetscErrorCode FormIFunction_InitializeEP(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes EP and tau */
PetscErrorCode FormIFunction_InitializeEP_halo(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes EP and tau for the halo current simulation*/
PetscErrorCode FormIFunction_InitializeEPV(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes V, EP and tau: [(\nabla x B0) x B0 = - lambda (\nabla^2 V0)].e_{R/Z}; V0 . B0 = 0 */
PetscErrorCode FormIFunction_newequilibrium(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes V: (\nabla x B) x B = (1/Re) (\nabla^2 V) + n_i * dV/dt; */
PetscErrorCode FormIFunction_newequilibrium_Vperp(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes V: [(\nabla x B) x B = -(1/Re) (\nabla^2 V) + n_i * dV/dt].e_{R/Z}; V . B = 0 */
PetscErrorCode FormIFunction_DampingV(TS,PetscReal,Vec,Vec,Vec,void*); /* This IFunction has the following vertex equations: (dampV V - (\nabla x B) x B).e_{R/Z} = 0 ; V.B = 0 */
PetscErrorCode FormRHSFunction_BImplicit(TS,PetscReal,Vec,Vec,void*); /* This RHSFunction contains only zero entries */
PetscErrorCode FormInitialSolution(TS,Vec,void*); /* This routine sets up the initial solution X_0 from input data files containing the B field values on face centers and the levelset function on cell centers */
PetscErrorCode FormExactSolution(PetscReal,TS,Vec*,void*); /* This routine sets up the vector which will be used to set the boundary conditions inside FormIFunction_DampingV */
PetscErrorCode Update_J_RE(TS);
PetscErrorCode Monitor(TS,PetscInt,PetscReal,Vec,void*);
PetscErrorCode FormDummyIJacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*); /* When used with PCFieldSplitSetDetectSaddlePoint, this dummy jacobian gives a field splitting where B field is first and tau, EP, V and n fields are last */
PetscErrorCode FormDummyIJacobian2(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*); /* When used with PCFieldSplitSetDetectSaddlePoint, this dummy jacobian gives a field splitting where tau field is first and EP, B, V and n fields are last */
PetscErrorCode FormDummyIJacobian3(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*); /* When used with PCFieldSplitSetDetectSaddlePoint, this dummy jacobian gives a field splitting where the density n is first and tau, EP, B and V fields are last */
PetscErrorCode FormDummyIJacobian4(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*); /* When used with PCFieldSplitSetDetectSaddlePoint, this dummy jacobian gives a field splitting where the electrostatic potential Phi is first and n, tau, B and V fields are last */
PetscErrorCode SampleShellPCSetUp(PC); /* This routine sets up a Shell Preconditioner that has the same effect as a 2-field fieldsplit preconditioner where the velocity is split from the remaining unknowns */
PetscErrorCode SampleShellPCApply(PC,Vec,Vec); /* This routine applies a Shell Preconditioner that has the same effect as a 2-field fieldsplit preconditioner where the velocity is split from the remaining unknowns */
PetscErrorCode SampleShellPCDestroy(PC); /* This routine destroys a Shell Preconditioner that has the same effect as a 2-field fieldsplit preconditioner where the velocity is split from the remaining unknowns */
PetscErrorCode ReadInitialData(PetscReal**,PetscInt*,const char*); /* This routine reads data from an input file and stores these in an array. The length of the output array is also computed. */
PetscErrorCode ReadALine(const char*,PetscInt,PetscReal*); /* This routine reads an entry from a specific line of an input file and stores it in the output array. */
PetscErrorCode FormInitialSolution_LargeData(TS,Vec,void*); /* This routine sets up the initial solution X_0. It uses ReadALine instead of ReadInitialData */
PetscErrorCode FormExactSolution_LargeData(PetscReal,TS,Vec*,void*); /* This routine sets up the vector which will be used to set the boundary conditions inside FormIFunction_DampingV. It uses ReadALine instead of ReadInitialData  */
PetscErrorCode FormIFunction_InitializeEP_LargeData(TS,PetscReal,Vec,Vec,Vec,void*); /* This is the Ifunction used in the TSSolver that initializes EP and tau. It uses ReadALine instead of ReadInitialData */
PetscErrorCode FormInitialSolution_psi(TS,Vec,void*); /* This routine sets up the initial solution X_0 from input data files containing G(psi) values on face centers, psi values on edge centers and the levelset function on cell centers */
PetscErrorCode FormInitialSolution_psi_fromNphi2(TS,Vec,void*); /* This routine sets up the initial solution X_0 from input data files containing G(psi) values on face centers, psi values on edge centers and the levelset function on cell centers. The data files used correspond to the resolution having only 2 cells in phi direction. Setting up the initial data on Nr x Nphi x Nz mesh from the initial data of Nr x 2 x Nz mesh */
PetscErrorCode FormInitialpsi(TS,Vec,void*); /* This routine sets up the initial psi that forces the current to be nonzero only in the vacuum vessel */

#endif /* defined(TS_FUNCTIONS_H) */
