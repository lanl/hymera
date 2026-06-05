//========================================================================================
// (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
// for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to // the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================

static
const char help[] = "Time-dependent magnetic diffusion PDE in 3d cylindrical coordinates using mimetic finite difference method for a simplified quasi-static perpendicular dynamics model.\n";
/*
  ni_t + div(ni Vi_perp) = 0,
  - (curl(B) x B - dampV V).e_r = 0,
  - (curl(B) x B - dampV V).e_z = 0,
  Vi_perp . B = 0,
  B_t = - curl(tau),
  tau = grad(EP) + eta curl(B) - Vi_perp x B,
  Laplacian(EP) = - Div(eta curl(B) - Vi_perp x B),
  mu0 = 4pi*10^-7,
  eta = 3.617*10^-7 in plasma region,
  eta = 3.617*10^-8 inside the wall region,
  eta = 3.617*10^-22 elsewhere.
  The density n_i (scalar) is defined on cell centers, the magnetic field B (vector) is defined on cell faces whereas the divergence-free component tau (vector) of the electric field E (vector) is defined on cell edges, the electrostatic potential EP (scalar) and the velocity Vi_perp (vector) are defined on vertices.
*/

#include <geometry.h>
#include <mass_matrix_coefficients.h>
#include <mfd_config.h>
#include <mimetic_operators.h>
#include <monitor_functions.h>
#include <ts_functions.h>
#include <petscsys.h>
#include <petscvec.h>
#include <mpi.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <petsc/private/dmstagimpl.h>
#include <fenv.h>

#include "default_petsc_options.h"
#include "geometry.h"
#include "mhd.h"

void view4d_zero(view4d_t v);
void view3d_zero(view3d_t v);

void subview_exclude_d1(view4d_t v4, view3d_t v);

int mhd_PetscInit(int * argc, char *** argv, User** user) {
  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  PetscErrorCode ierr = PetscInitialize( argc, argv, (char * ) 0, help);
  if (ierr) {
    printf("CRITICAL: PetscInitialize returned error, aborting mhd_initialize\n");
    return 1;
  }

  *user = (User*) calloc(1, sizeof(User));
  return 0;
}

int mhd_initialize(User* user) {
  KSP ksp, dummyksp, dummykspB, dummykspn, dummykspEP; /* scalable linear equations solver */
  //char              *prefix[2];
  PC dummypc, dummypcB, dummypcn, dummypcEP; /* preconditioner context */
  PC pc;
  Vec dummyX; /* solution and right-hand side vectors */
  Mat J, Jpre, dummyJ, dummyJn;
  //,Jmf = NULL;       /* jacobian matrix */
  PetscErrorCode ierr = 0;
  SNES snes;
  PetscReal time, ftime;

  PetscBool matrix_free = PETSC_FALSE, matrix_free_FDprec = PETSC_FALSE, user_defined_pc = PETSC_FALSE;
  KSP * subksp, * subsubksp, * subsubsubksp;
  PC subsubsubpc[2] = {NULL,NULL}, subsubpc[2] = {NULL,NULL}, subpc[2] = {NULL,NULL};
  PetscInt n = 1;
  PetscBool removezero = PETSC_FALSE;


  // Avoid command-line options.
  // Some of these probably need to be separated out to separate statements with a default argument.
  PetscCall(default_petsc_options());

  ierr = AppCtxView(PETSC_COMM_WORLD, user);
  if (ierr) {
    printf("CRITICAL: AppCtxView returned error, aborting mhd_initialize\n");
    return 1;
  }

  PetscInt numC = 0;
  char filename[PETSC_MAX_PATH_LEN];

  PetscSNPrintf(filename, sizeof(filename), "%s/veceta_grid%.3Dx%.2Dx%.3D.txt", user->input_folder, user->Nr, user->Nphi, user->Nz);
  ReadInitialData( & (user->dataC), & numC, filename);

  PetscSNPrintf(filename, sizeof(filename), "%s/vecpsi_grid%.3Dx%.2Dx%.3D.txt", user->input_folder, user->Nr, user->Nphi, user->Nz);
  ReadInitialData( & (user->datapsi), & (user->numpsi), filename);

  PetscSNPrintf(filename, sizeof(filename), "%s/vecg_grid%.3Dx%.2Dx%.3D.txt", user->input_folder, user->Nr, user->Nphi, user->Nz);
  ReadInitialData( & (user->datag), & (user->numg), filename);

  PetscPrintf(PETSC_COMM_WORLD, "Read initial data from g and psi input!\n");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create 3D DMStag for the solution, and set up.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  {
    const PetscInt dof0 = 4,
      dof1 = 1,
      dof2 = 1,
      dof3 = 1; /* 1 dof on each edge, face and cell center and 4 dofs on each vertex (3 for vector field V, and 1 for the electrostatic potential EP) */
    const PetscInt stencilWidth = 1;
    PetscInt nr = 0, nphi = 0, nz = 0;

    if (user->phibtype) {
      DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE, user->Nr, user->Nphi, user->Nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, & user->da);
    } else {
      DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, user->Nr, user->Nphi, user->Nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, & user->da);
    }
    DMSetFromOptions(user->da);
    DMSetUp(user->da);

    DMStagGetNumRanks(user->da, & nr, & nphi, & nz);

    if (user->phibtype) {
      DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE, user->Nr, user->Nphi, user->Nz, nr, nphi, nz, 4, 1, 1, 1, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, & (user->coorda));
    } else {
      DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, user->Nr, user->Nphi, user->Nz, nr, nphi, nz, 4, 1, 1, 1, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, & (user->coorda));
    }

    DMSetFromOptions(user->coorda);
    DMSetUp(user->coorda);
    DMStagSetUniformCoordinatesExplicit(user->da, user->rmin/user->L0, user->rmax/user->L0, user->phimin, user->phimax, user->zmin/user->L0, user->zmax/user->L0);
    DMStagSetUniformCoordinatesExplicit(user->coorda, user->rmin/user->L0, user->rmax/user->L0, user->phimin, user->phimax, user->zmin/user->L0, user->zmax/user->L0);

    DM dmCoorda;
    Vec coordaLocal;

    DMGetCoordinateDM(user->coorda, & dmCoorda);
    DMGetCoordinatesLocal(user->coorda, & coordaLocal);
    DMStagVecGetArrayRead(dmCoorda, coordaLocal, &user->arrCoord);

    DMSetApplicationContext(user->da, user);
    DMCreateGlobalVector(user->da, & user->X0);
    DMCreateGlobalVector(user->da, & user->X);
  }
  /* Print out some info */
  {
    PetscInt N[3];
    DMStagGetGlobalSizes(user->da, & N[0], & N[1], & N[2]);
    PetscPrintf(PETSC_COMM_WORLD, "Using a %D x %D x %D mesh\n", N[0], N[1], N[2]);
    PetscPrintf(PETSC_COMM_WORLD, "dr: %g\n", user->dr);
    PetscPrintf(PETSC_COMM_WORLD, "dphi: %g\n", user->dphi);
    PetscPrintf(PETSC_COMM_WORLD, "dz: %g\n", user->dz);
    PetscPrintf(PETSC_COMM_WORLD, "normalized dt: %g\n", user->dt);
    PetscPrintf(PETSC_COMM_WORLD, "non-normalized dt: %g\n", user->dt * user->L0/user->V_A); // Alfven time := L0 / V_A
    PetscPrintf(PETSC_COMM_WORLD, "Characteristic resistive time: %g\n", user->L0*user->L0*user->mu0/user->etaplasma); // tau_eta := mu0 L0^2 / non_normalized_eta
    PetscPrintf(PETSC_COMM_WORLD, "Reynolds parameter for viscosity: %g\n", user->Re);
    PetscPrintf(PETSC_COMM_WORLD, "Lundquist number: %g\n", user->eta0 / user->etaplasma); // tau_eta / Alfven time = mu0 L0 V_A / non_normalized_eta
    PetscPrintf(PETSC_COMM_WORLD, "Resistivity inside the plasma chamber (in Ohm.meter): %g\n",user->etaplasma);
    PetscPrintf(PETSC_COMM_WORLD, "Resistivity outside the vacuum vessel (in Ohm.meter): %g\n", user->etaout );
    PetscPrintf(PETSC_COMM_WORLD, "Resistivity inside the vacuum vessel (in Ohm.meter): %g\n", user->etaVV );
    PetscPrintf(PETSC_COMM_WORLD, "Resistivity inside the blanket module (in Ohm.meter): %g\n", user->etawall );
    //PetscPrintf(PETSC_COMM_WORLD, "CFL value: %g\n", user->dt / PetscMin(user->dr,user->dz));
  }

  VecZeroEntries(user->X);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create timestepping solver context
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSCreate(PETSC_COMM_WORLD, & user->ts);
  TSSetDM(user->ts, user->da);

  if (user->savecoords) {
    SaveCoordinates(user->ts, user);
    return (0);
  }

  //DMSetMatrixPreallocateOnly(user->da,PETSC_TRUE);

  //TSSetProblemType(user->ts, TS_LINEAR);
  TSMonitorSet(user->ts, Monitor, user, NULL); /* Set optional user-defined monitoring routine */

  switch (user->tstype) {
  case 1:
    TSSetType(user->ts, TSEULER); /* Forward Euler method */
    break;
  case 2:
    TSSetType(user->ts, TSBEULER); /* Backward Euler method */
    break;
  case 3:
    TSSetType(user->ts, TSCN); /* Crank-Nicholson method */
    break;
  case 4:
//    TSGetAdapt(user->ts, & adapt);
//    TSAdaptSetType(adapt, TSADAPTNONE);
//    TSSetType(user->ts, TSARKIMEX); /* Additive Runge-Kutta IMEX method */
//    TSARKIMEXSetFullyImplicit(user->ts, PETSC_TRUE);
//    TSARKIMEXSetType(user->ts, TSARKIMEXL2);
//
//    TSSetEquationType(user->ts,TS_EQ_IMPLICIT);
    break;
  case 5:
    //TSSetProblemType(user->ts,TS_NONLINEAR);
    TSSetType(user->ts,TSROSW);
    break;
  case 6:
    TSSetType(user->ts, TSBDF); /* Backward differentiation formula of order 2*/
    TSBDFSetOrder(user->ts,2);
    break;
  case 7:
    TSSetType(user->ts, TSTHETA); /* Implicit Theta method */
    TSThetaSetTheta(user->ts, 0.5); // Default theta is 0.5 but this value can be set through command line using the -ts_theta_theta flag
    break;
  default:
    TSSetType(user->ts, TSEULER); /* Forward Euler method */
    break;
  }

  TSGetSNES(user->ts, & snes);
  if (user->tstype > 1) {
    PetscOptionsGetBool(NULL, NULL, "-snes_mf", & matrix_free, NULL);
    PetscOptionsGetBool(NULL, NULL, "-snes_mf_operator", & matrix_free_FDprec, NULL);
    PetscOptionsGetBool(NULL, NULL, "-removezero", &removezero, NULL);
    if (matrix_free) { //matrix_free mode without any preconditioning matrix
      PetscPrintf(PETSC_COMM_WORLD, "======use matrix-free evaluation and no preconditioning======\n");
    } else if (matrix_free_FDprec) { //matrix_free mode with colored finite difference jacobian for preconditioning
      DMCreateMatrix(user->da, & J);
      TSSetIJacobian(user->ts, J, J, TSComputeIJacobianDefaultColor, NULL);
      PetscPrintf(PETSC_COMM_WORLD, "======use matrix-free evaluation and FD coloring Jacobian for preconditioning======\n");
    } else {
      DMCreateMatrix(user->da, & J);
      DMCreateMatrix(user->da, & Jpre);
      if (user->jtype == 0) {
        TSSetIJacobian(user->ts, J, Jpre, FormIJacobian_BImplicit, user); /* use user provided Jacobian evaluation routine */
        PetscPrintf(PETSC_COMM_WORLD, "======use Analytical Jacobian======\n");
      } else {
        /* use finite difference Jacobian J as preconditioner and '-snes_mf_operator' for Mat*vec */
        /*MatCreateSNESMF(snes,&Jmf);*/
        if (user->jtype == 1) {
          /* slow finite difference J; */
          SNESSetJacobian(snes, J, J, SNESComputeJacobianDefault, PETSC_NULLPTR);
          PetscPrintf(PETSC_COMM_WORLD, "======use FD Jacobian======\n");
        } else if (user->jtype == 2) {
          /* Use coloring to compute finite difference J efficiently */
          TSSetIJacobian(user->ts, J, J, TSComputeIJacobianDefaultColor, PETSC_NULLPTR);
          PetscPrintf(PETSC_COMM_WORLD, "======use FD coloring Jacobian======\n");
        } else {
          SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "This jtype is not supported");
        }
      }
    }
    //TSSetIFunction(user->ts, NULL, FormIFunction_DampingV, user);
    TSSetIFunction(user->ts, NULL, FormIFunction_Vperp_viscosity, user);
    TSSetRHSFunction(user->ts, NULL, FormRHSFunction_BImplicit, user);
  } else {
    TSSetRHSFunction(user->ts, NULL, FormRHSFunction_BImplicit, user);
  }
  SNESSetFromOptions(snes);

  TSSetTime(user->ts, user->itime);
  ftime = user->ftime;
  TSSetMaxTime(user->ts, ftime);
  TSSetExactFinalTime(user->ts, TS_EXACTFINALTIME_STEPOVER);
  TSSetSolution(user->ts, user->X);
  TSSetTimeStep(user->ts, user->dt);

  TSSetFromOptions(user->ts);
  TSSetUp(user->ts);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set index sets
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  IS isEPBVndup, isEPtauVndup, isEPtauBVdup, istauBVndup, isndup, isEPdup, isBdup, istaudup;

  DMCreateGlobalVector(user->da, & dummyX);
  VecCopy(user->X, dummyX);
  DMCreateMatrix(user->da, & dummyJ);

  if (user->tstype > 1 && matrix_free_FDprec) {
    PetscRandom rctx;
    PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rctx));
    PetscCall(PetscRandomSetInterval(rctx, 1.0, 2.0));
    PetscCall(VecSetRandom(user->X, rctx));
    PetscCall(VecSetRandom(dummyX, rctx));
    PetscCall(PetscRandomDestroy(&rctx));
    PetscCall(TSComputeIJacobian(user->ts, 0.0, user->X, dummyX, 2.0, J, J, PETSC_FALSE));
    if (removezero) PetscCall(TSPruneIJacobianColor(user->ts, J, J));
  }

  KSPCreate(PETSC_COMM_WORLD, & dummykspEP);
  KSPSetOptionsPrefix(dummykspEP, "sepKSP_");
  FormDummyIJacobian4(user->ts,dummyX, dummyX, 1.0 / user->dt, dummyJ, dummyJ, user);
  KSPSetOperators(dummykspEP, dummyJ, dummyJ);
  KSPGetPC(dummykspEP, & dummypcEP);
  PCSetType(dummypcEP, PCFIELDSPLIT);
  PCFieldSplitSetDetectSaddlePoint(dummypcEP, PETSC_TRUE);
  PCSetUp(dummypcEP);
  PCFieldSplitSetType(dummypcEP, PC_COMPOSITE_SCHUR);
  PCFieldSplitSetSchurFactType(dummypcEP, PC_FIELDSPLIT_SCHUR_FACT_FULL);
  PCFieldSplitSetSchurPre(dummypcEP, PC_FIELDSPLIT_SCHUR_PRE_SELFP, NULL);
  KSPSetUp(dummykspEP);
  KSPSetTolerances(dummykspEP, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 outer iteration for the dummy solve
  PCFieldSplitGetSubKSP(dummypcEP, & n, & subksp);
  ierr = KSPGetPC(subksp[1], & (subpc[1]));
  if (ierr) {
    printf("CRITICAL: KSPGetPC returned error, aborting mhd_initialize\n");
    return 1;
  }
  KSPSetTolerances(subksp[1], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 inner iteration maximum for the dummy solve
  ierr = KSPSolve(dummykspEP, dummyX, dummyX);
  if (ierr) {
    printf("CRITICAL: KSPSolve returned error, aborting mhd_initialize\n");
    return 1;
  }
  //ISDuplicate(isBdup, & user->isB);
  PCFieldSplitGetISByIndex(dummypcEP, 0, & isEPdup);
  ISDuplicate(isEPdup, & user->isEP);


  VecDestroy( & dummyX);
  KSPDestroy( & dummykspEP);
  MatDestroy( & dummyJ);


  char ** namelist;
  IS * islist, isALL, isALL_V, isBV;
  IS ISV;

  PetscInt len, d = 0;
  ierr = DMCreateFieldDecomposition(user->da, & len, & namelist, & islist, NULL);
  if (ierr) {
    printf("CRITICAL: DMCreateFieldDecomposition returned error, aborting mhd_initialize\n");
    return 1;
  }
  PetscPrintf(PETSC_COMM_WORLD, "The number of subproblems in the field decomposition is: %g\n", (double)(len));
  for (d = 0; d < len; ++d) {
    PetscPrintf(PETSC_COMM_WORLD, "The name of field number %d is: %s.\n", d, namelist[d]);
    //PetscPrintf(PETSC_COMM_WORLD, "The global indices for field number %d are as follows.\n", d);
    //ISView(islist[d],PETSC_VIEWER_STDOUT_SELF);
  }
  PetscBool flagV = PETSC_FALSE, flagE = PETSC_FALSE, flagF = PETSC_FALSE, flagC = PETSC_FALSE;
  ISDifference(islist[0], user->isEP, & user->isV);

  ISDuplicate(islist[1], & user->istau);
  ISDuplicate(islist[2], & user->isB);
  ISDuplicate(islist[3], & user->isni);

  const IS islist2[5] = {user->isV, user->isEP, user->istau, user->isB, user->isni};

  ISConcatenate(PETSC_COMM_WORLD,5,islist2,&isALL);
  ISDifference(isALL, user->isV, & isALL_V);
  ISDestroy( & isALL);

  const IS islist3[5] = {user->isV, user->isB};
  ISConcatenate(PETSC_COMM_WORLD,2,islist3,&isBV);
  ISDestroy( & isBV);

  for (d = 0; d < len; ++d) {
    ISDestroy( & islist[d]);
  }
  PetscFree(islist);
  PetscFree(namelist);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set preconditioner options
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  SNESGetKSP(snes, & ksp);
  ierr = KSPGetPC(ksp, &pc);
  if (ierr) {
    printf("CRITICAL: KSPGetPC, second occurance returned error, aborting mhd_initialize\n");
    return 1;
  }
  /* Set a user-defined "shell" preconditioner if desired */
  PetscOptionsGetBool(NULL,NULL,"-user_defined_pc",&user_defined_pc,NULL);
  if (user_defined_pc) {
    /* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
    PCSetType(pc,PCSHELL);
    PCShellSetContext(pc,user);

    /* Do any setup required for the preconditioner */
    PCShellSetSetUp(pc,SampleShellPCSetUp);

    /* (Required) Set the user-defined routine for applying the preconditioner */
    PCShellSetApply(pc,SampleShellPCApply);

    /* (Optional) Set user-defined function to free objects used by custom preconditioner */
    PCShellSetDestroy(pc,SampleShellPCDestroy);

    /* (Optional) Set a name for the preconditioner, used for PCView() */
    PCShellSetName(pc,"ShellPrec");
  }
  else {/* first level -> split ni from the rest : {ni}, {V Phi tau B}
         second level -> split tau from {V Phi B} : {ni}, {{tau},{V Phi B}}
         third level -> split Phi from {B V} : {ni}, {{tau},{{Phi}, {B V}}}
         fourth level -> split {V} from {B} : {ni}, {{tau},{{Phi}, {{B},{V}}}}
         */
    IS            is[2];
    DMStagStencil stencil0[1], stencil1[10];
    PC            pc_notc, pc_noe;

    const char *name[2] = {"ni", "TEBV"};

    PetscCall(KSPGetPC(ksp,&pc));
    PetscCall(PCSetType(pc,PCFIELDSPLIT));

    // First split is cells
    stencil0[0].loc = DMSTAG_ELEMENT;
    stencil0[0].c = 0;

    // Second split is the rest
    for (PetscInt c=0; c<4; ++c) {
      stencil1[c].loc = DMSTAG_BACK_DOWN_LEFT;
      stencil1[c].c = c;
    }
    stencil1[4].loc = DMSTAG_LEFT;
    stencil1[4].c = 0;
    stencil1[5].loc = DMSTAG_BACK;
    stencil1[5].c = 0;
    stencil1[6].loc = DMSTAG_DOWN;
    stencil1[6].c = 0;
    stencil1[7].loc = DMSTAG_BACK_DOWN;
    stencil1[7].c = 0;
    stencil1[8].loc = DMSTAG_BACK_LEFT;
    stencil1[8].c = 0;
    stencil1[9].loc = DMSTAG_DOWN_LEFT;
    stencil1[9].c = 0;

    PetscCall(DMStagCreateISFromStencils(user->da,1,stencil0,&is[0]));
    PetscCall(DMStagCreateISFromStencils(user->da,10,stencil1,&is[1]));

    for (PetscInt i=0; i<2; ++i) {
      PetscCall(PCFieldSplitSetIS(pc,name[i],is[i]));
    }

    for (PetscInt i=0; i<2; ++i) {
      PetscCall(ISDestroy(&is[i]));
    }

    /* Logic below modifies the PC directly, so this is the last chance to change the solver from the command line */
    PetscCall(KSPSetFromOptions(ksp));

    PetscBool is_fieldsplit;
    /* If the fieldsplit PC wasn't overridden, further split the second split */
    {
      PCType pc_type;


      PetscCall(KSPGetPC(ksp, &pc));
      PetscCall(PCGetType(pc,&pc_type));
      PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
      if (is_fieldsplit) {
        DM            dm_notc;
        KSP           *sub_ksp;

        PetscInt      n_splits;
        DMStagStencil stencil_notc_edges[3], stencil_notc_notedges[7];
        IS            is_notc[2];
        const char    *name_notc[2] = {"tau","EBV"};

        PetscCall(PCSetUp(pc)); // Set up the Fieldsplit PC
        PetscCall(PCFieldSplitGetSubKSP(pc,&n_splits,&sub_ksp));
        PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)user->da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
        PetscCall(KSPGetPC(sub_ksp[1],&pc_notc));
        PetscCall(PetscFree(sub_ksp));

        PetscCall(DMStagCreateCompatibleDMStag(user->da,4,1,1,0,&dm_notc));

        // First split within notc is edges
        stencil_notc_edges[0].loc = DMSTAG_BACK_DOWN;
        stencil_notc_edges[0].c = 0;
        stencil_notc_edges[1].loc = DMSTAG_BACK_LEFT;
        stencil_notc_edges[1].c = 0;
        stencil_notc_edges[2].loc = DMSTAG_DOWN_LEFT;
        stencil_notc_edges[2].c = 0;

        // Second split within notc is faces and vertices
        for (PetscInt c=0; c<3; ++c) {
          stencil_notc_notedges[c].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_notc_notedges[c].c = c;
        }
        stencil_notc_notedges[3].loc = DMSTAG_BACK_DOWN_LEFT;
        stencil_notc_notedges[3].c = 3;
        stencil_notc_notedges[4].loc = DMSTAG_LEFT;
        stencil_notc_notedges[4].c = 0;
        stencil_notc_notedges[5].loc = DMSTAG_BACK;
        stencil_notc_notedges[5].c = 0;
        stencil_notc_notedges[6].loc = DMSTAG_DOWN;
        stencil_notc_notedges[6].c = 0;

        PetscCall(DMStagCreateISFromStencils(dm_notc,3,stencil_notc_edges,&is_notc[0]));
        PetscCall(DMStagCreateISFromStencils(dm_notc,7,stencil_notc_notedges,&is_notc[1]));

        for (PetscInt i=0; i<2; ++i) {
          PetscCall(PCFieldSplitSetIS(pc_notc,name_notc[i],is_notc[i]));
        }

        for (PetscInt i=0; i<2; ++i) {
          PetscCall(ISDestroy(&is_notc[i]));
        }
        PetscCall(DMDestroy(&dm_notc));
      }
    }

    /* If the fieldsplit PC wasn't overridden, further split the second split of the second level */
    if (is_fieldsplit) {
      PCType pc_type;


      PetscCall(PCGetType(pc_notc,&pc_type));
      PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
      if (is_fieldsplit) {
        DM            dm_noe;
        KSP           *sub_ksp;

        PetscInt      n_splits;
        DMStagStencil stencil_noe_EP[1], stencil_noe_notEP[6];
        IS            is_noe[2];
        const char    *name_noe[2] = {"EP", "BV"};

        PetscCall(PCSetUp(pc_notc)); // Set up the Fieldsplit PC
        PetscCall(PCFieldSplitGetSubKSP(pc_notc,&n_splits,&sub_ksp));
        PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)user->da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
        PetscCall(KSPGetPC(sub_ksp[1],&pc_noe));
        PetscCall(PetscFree(sub_ksp));

        PetscCall(DMStagCreateCompatibleDMStag(user->da,4,0,1,0,&dm_noe));

        // First split within notv is 4th dofs on vertices
        stencil_noe_EP[0].loc = DMSTAG_BACK_DOWN_LEFT;
        stencil_noe_EP[0].c = 3;

        // Second split within notv is faces and the first 3 dofs on vertices
        for (PetscInt c=0; c<3; ++c) {
          stencil_noe_notEP[c].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_noe_notEP[c].c = c;
        }
        stencil_noe_notEP[3].loc = DMSTAG_LEFT;
        stencil_noe_notEP[3].c = 0;
        stencil_noe_notEP[4].loc = DMSTAG_BACK;
        stencil_noe_notEP[4].c = 0;
        stencil_noe_notEP[5].loc = DMSTAG_DOWN;
        stencil_noe_notEP[5].c = 0;

        PetscCall(DMStagCreateISFromStencils(dm_noe,1,stencil_noe_EP,&is_noe[0]));
        PetscCall(DMStagCreateISFromStencils(dm_noe,6,stencil_noe_notEP,&is_noe[1]));

        for (PetscInt i=0; i<2; ++i) {
          PetscCall(PCFieldSplitSetIS(pc_noe,name_noe[i],is_noe[i]));
        }

        for (PetscInt i=0; i<2; ++i) {
          PetscCall(ISDestroy(&is_noe[i]));
        }
        PetscCall(DMDestroy(&dm_noe));
      }
    }

    PC            pc_noe_2;

    /* If the fieldsplit PC wasn't overridden, further split the first split of the third level */
    if (is_fieldsplit) {
      PCType pc_type;


      PetscCall(PCGetType(pc_noe,&pc_type));
      PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
      if (is_fieldsplit) {
        DM            dm_notv;
        KSP           *sub_ksp;

        PetscInt      n_splits;
        DMStagStencil stencil_notv_faces[3], stencil_notv_notfaces[3];
        IS            is_notv[2];
        const char    *name_notv[2] = {"V", "B"};

        PetscCall(PCSetUp(pc_noe)); // Set up the Fieldsplit PC
        PetscCall(PCFieldSplitGetSubKSP(pc_noe,&n_splits,&sub_ksp));
        PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)user->da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
        PetscCall(KSPGetPC(sub_ksp[1],&pc_noe_2));
        PetscCall(PetscFree(sub_ksp));

        PetscCall(DMStagCreateCompatibleDMStag(user->da,3,0,1,0,&dm_notv));

        // First split within notv is faces
        stencil_notv_faces[0].loc = DMSTAG_LEFT;
        stencil_notv_faces[0].c = 0;
        stencil_notv_faces[1].loc = DMSTAG_BACK;
        stencil_notv_faces[1].c = 0;
        stencil_notv_faces[2].loc = DMSTAG_DOWN;
        stencil_notv_faces[2].c = 0;

        // Second split within notv is vertices
        for (PetscInt c=0; c<3; ++c) {
          stencil_notv_notfaces[c].loc = DMSTAG_BACK_DOWN_LEFT;
          stencil_notv_notfaces[c].c = c;
        }

        PetscCall(DMStagCreateISFromStencils(dm_notv,3,stencil_notv_notfaces,&is_notv[0]));
        PetscCall(DMStagCreateISFromStencils(dm_notv,3,stencil_notv_faces,&is_notv[1]));


        for (PetscInt i=0; i<2; ++i) {
          PetscCall(PCFieldSplitSetIS(pc_noe_2,name_notv[i],is_notv[i]));
        }

        for (PetscInt i=0; i<2; ++i) {
          PetscCall(ISDestroy(&is_notv[i]));
        }
        PetscCall(DMDestroy(&dm_notv));
      }
    }
   }

  // Decrese reference count of jacobian
  if (user->tstype > 1) {
    MatDestroy( & J);
  }
  FormInitialSolution_psi(user->ts, user->X, user); // Set initial condition
  return 0;
}

int mhd_step(User* user) {
  Vec X;
  TSGetSolution(user->ts, &X);
  // Backup solition state
  VecCopy(X, user->X0);

  // Do one step
  PetscErrorCode ierr = TSStep(user->ts);
  if (ierr) {
    printf("CRITICAL: TSSolve, second occurance returned error, aborting mhd_initialize\n");
    return 1;
  }

  /* update internal time */
  PetscReal t = 0;
  TSGetTime(user->ts, &t);
  t += user->dt;
  TSSetTime(user->ts, t);

  if (user->savesol) {
    SaveSolution(user->ts,user->X,user);
  }

  PetscReal ftime;
  TSGetSolveTime(user->ts, & ftime);

  PetscInt steps;
  TSGetStepNumber(user->ts, & steps);

  TSConvergedReason reason;
  TSGetConvergedReason(user->ts, & reason);
  PetscPrintf(PETSC_COMM_WORLD, "%s at time %g after %D steps\n", TSConvergedReasons[reason], (double) ftime, steps);
  return 0;
}

int mhd_resetState(User* user) {
  Vec X;
  TSGetSolution(user->ts, &X);
  // Backup solition state
  VecCopy(user->X0, X);
  TSSetSolution(user->ts, X);

  /* update internal time */
  PetscReal t = 0;
  TSGetTime(user->ts, &t);
  t -= user->dt;
  TSSetTime(user->ts, t);

  return 0;
}

int mhd_getF(User* user, field_id fid, view3d_t v) {
  Vec X;
  TSGetSolution(user->ts, &X);

  DM da;
  TSGetDM(user->ts,&da);

  size_t NR = user->Nr;
  size_t Nphi = user->Nphi;
  size_t NZ = user->Nz;

  view4d_t v4 = {
      (PetscScalar*) malloc(sizeof(PetscScalar) * (3 * NR * Nphi * NZ)),
      NR, Nphi, NZ, 3,
      3, 3 * NR, 3 * NR * Nphi, 1
  };

  view4d_zero(v4);

  // E field and J field have different ordering, field index is the slowest
  if (fid == fid_E || fid == fid_J) {
    v4.stride0 = 1;
    v4.stride1 = NR;
    v4.stride2 = NR * Nphi;
    v4.stride3 = NR * Nphi * NZ;
  }

  if (fid == fid_B) // compute B
    getBArray(user->ts, X, v4.data, user, 0);
  else if (fid == fid_E) // compute E
    getEJArray(user->ts, X,
        v4.data,
        v4.data + v4.stride3,
        v4.data + 2 * v4.stride3,
    user, 0);
  else if (fid == fid_Jre) {
    view3d_t jre = user->jre;
    for (size_t i = 0; i < jre.dim0; ++i)
    for (size_t j = 0; j < jre.dim1; ++j)
    for (size_t k = 0; k < jre.dim2; ++k)
      v.data[i * v.stride0 + j * v.stride1 + k * v.stride2] =
        jre.data[i * jre.stride0 + j * jre.stride1 + k * jre.stride2];
    free(v4.data);
    return 0;
  }
  else if (fid == fid_J)
    getEJArray(user->ts, X,
        v4.data,
        v4.data + v4.stride3,
        v4.data + 2 * v4.stride3,
    user, 1);
  else if (fid == fid_V)
    getVArray(user->ts, X, v4.data, user);
  else if (fid == fid_GradB)
    getBArray(user->ts, X, v4.data, user, 1);
  else {
    PetscPrintf(PETSC_COMM_WORLD, "ERROR: Wrong field id\n");
    free(v4.data);
    return 1;
  }

  subview_exclude_d1(v4, v);
  free(v4.data);
  return 0;
}

int mhd_destroy(User* user) {
  PetscErrorCode ierr = 0;
  // Free work space.
  ISDestroy( & user->isALL_V);

  {
    DM dmCoorda;
    Vec coordaLocal;

    DMGetCoordinateDM(user->coorda, & dmCoorda);
    DMGetCoordinatesLocal(user->coorda, & coordaLocal);
    DMStagVecRestoreArrayRead(dmCoorda, coordaLocal, user->arrCoord);
  }

  VecDestroy( & user->X);
  VecDestroy( & user->X0);
  TSDestroy( & user->ts);

  DMDestroy( & user->da);

  PetscFinalize();

  return 0;
}

int mhd_savesolution(User* user, const char* filename) {
  PetscViewer viewerX;
  PetscPrintf(PETSC_COMM_WORLD, "Writing X vector into file %s ...\n", filename);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, & viewerX);

  Vec X;
  TSGetSolution(user->ts, &X);
  VecView(X, viewerX);

  PetscViewerDestroy( & viewerX);
  PetscPrintf(PETSC_COMM_WORLD, "Created %s\n", filename);
  return 0;
}

int mhd_loadsolution(User* user, const char* filename) {
  PetscViewer viewerX;
  PetscPrintf(PETSC_COMM_WORLD, "Reading X vector from file %s ...\n", filename);
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_READ, & viewerX);
  Vec X;
  TSGetSolution(user->ts, &X);
  VecLoad(X, viewerX);
  TSSetSolution(user->ts, X);

  PetscViewerDestroy( & viewerX);
  PetscPrintf(PETSC_COMM_WORLD, "Reading from file %s is over.\n", filename);
}

void view4d_zero(view4d_t v) {
  for (size_t i = 0; i < v.dim0 * v.dim1 * v.dim2 * v.dim3; ++i)
        v.data[i] = 0.0;
}

void view3d_zero(view3d_t v) {
  for (size_t i = 0; i < v.dim0 * v.dim1 * v.dim2; ++i)
        v.data[i] = 0.0;
}

void subview_exclude_d1(view4d_t v4, view3d_t  v) {
  for (size_t i = 0; i < v.dim0; ++i)
    for (size_t j = 0; j < v.dim1; ++j)
      for (size_t k = 0; k < v.dim2; ++k)
        v.data[i * v.stride0 + j * v.stride1 + k * v.stride2] =
        v4.data[i * v4.stride0 + j * v4.stride2 + k * v4.stride3];
}






