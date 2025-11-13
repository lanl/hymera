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
#include <hFlux/c_wrapper.h>

#include "mhd_run.h"

int mhd_run(int argc, char ** argv, double* raw_field_ptr) {
  KSP ksp, dummyksp, dummykspB, dummykspn, dummykspEP; /* scalable linear equations solver */
  //char              *prefix[2];
  PC dummypc, dummypcB, dummypcn, dummypcEP; /* preconditioner context */
  PC pc;
  Vec X, dummyX; /* solution and right-hand side vectors */
  Mat J, Jpre, dummyJ, dummyJn;
  //,Jmf = NULL;       /* jacobian matrix */
  PetscInt steps;
  PetscErrorCode ierr = 0;
  DM da;
  SNES snes;
  PetscReal time, ftime;
  User user; /* user-defined work context */
  //TS ts = user.ts; /* time integrator */
  TSConvergedReason reason;
  TSAdapt adapt;
  PetscBool matrix_free = PETSC_FALSE, matrix_free_FDprec = PETSC_FALSE, user_defined_pc = PETSC_FALSE;
  KSP * subksp, * subsubksp, * subsubsubksp;
  PC subsubsubpc[2] = {NULL,NULL}, subsubpc[2] = {NULL,NULL}, subpc[2] = {NULL,NULL};
  PetscInt n = 1;
  PetscBool removezero = PETSC_FALSE;

  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  PetscInitialize( & argc, & argv, (char * ) 0, help);
  if (ierr) return ierr;

  // Avoid command-line options.
  // Some of these probably need to be separated out to separate statements with a default argument.
  PetscOptionsInsertString(NULL, "-ts_adapt_type none -B0 5.3 -dummyKSP_ksp_converged_reason -dummyKSP_ksp_monitor -dummySNES_snes_lag_jacobian 1 -dummySNES_snes_lag_preconditioner 1 -dummySNES_snes_monitor -dummySNES_snes_linesearch_type nleqerr -dummySNES_snes_mf_operator -dummySNES_snes_converged_reason -dummyKSP_ksp_type fgmres -dummyKSP_pc_type lu -dummyKSP_pc_factor_mat_solver_type mumps -dummyKSP_mat_mumps_icntl_6 2 -dummyKSP_mat_mumps_icntl_24 1 -dummyKSP_mat_mumps_cntl_1 1e-5 -dummyKSP_mat_mumps_cntl_3 1e-5 -dummyKSP_mat_mumps_icntl_14 5000 -dummyKSP_pc_fieldsplit_type multiplicative -dummyKSP_fieldsplit_ni_ksp_type fgmres -dummyKSP_fieldsplit_ni_pc_type hypre -dummyKSP_fieldsplit_ni_pc_hypre_type euclid -dummyKSP_fieldsplit_TEBV_ksp_type fgmres -dummyKSP_fieldsplit_TEBV_pc_type fieldsplit -dummyKSP_fieldsplit_TEBV_pc_fieldsplit_type schur -dummyKSP_fieldsplit_TEBV_pc_fieldsplit_schur_precondition selfp -dummyKSP_fieldsplit_TEBV_fieldsplit_tau_ksp_type gmres -dummyKSP_fieldsplit_TEBV_fieldsplit_tau_pc_type bjacobi -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_ksp_type preonly -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_pc_type fieldsplit -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_pc_fieldsplit_type multiplicative -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_EP_ksp_type gmres -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_EP_pc_type hypre -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_ksp_type preonly -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_pc_type lu -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_pc_factor_mat_solver_type mumps -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_icntl_6 2 -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_icntl_24 1 -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_cntl_1 1e-5 -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_cntl_3 1e-5 -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_icntl_14 5000 -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_superlu_dist_replacetinypivot -dummyKSP_ksp_rtol 1e-6 -dummyKSP_ksp_norm_type unpreconditioned -dummySNES_snes_max_funcs 100000000000 -dummySNES_snes_mf_operator -dummySNES_snes_stol 1e-20 -dummySNES_snes_rtol 1e-5 -dummySNES_snes_max_it 20 -dummySNES_snes_ksp_ew -dummySNES_snes_ksp_ew_version 3 -dummySNES_snes_ksp_ew_rtol0 0.2 -dummySNES_snes_ksp_ew_rtolmax 0.9 -dummySNES_snes_ksp_ew_gamma 0.9 -dummySNES_snes_ksp_ew_alpha 1.5 -dummySNES_snes_ksp_ew_alpha2 1.5 -dummySNES_snes_ksp_ew_threshold 0.1 -dummyKSP_fieldsplit_ni_ksp_converged_reason -dummyKSP_fieldsplit_TEBV_ksp_converged_reason -dummyKSP_fieldsplit_TEBV_fieldsplit_tau_ksp_converged_reason -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_EP_ksp_converged_reason -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_ksp_converged_reason -dummyKSP_fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_fieldsplit_V_ksp_converged_reason -ksp_converged_reason -ksp_monitor -snes_monitor -snes_mf_operator -snes_converged_reason -ksp_type fgmres -pc_type lu -pc_factor_mat_solver_type mumps -mat_mumps_icntl_6 2 -mat_mumps_icntl_24 1 -mat_mumps_cntl_1 1e-5 -mat_mumps_cntl_3 1e-5 -mat_mumps_icntl_14 5000 -pc_fieldsplit_type multiplicative -fieldsplit_ni_ksp_type fgmres -fieldsplit_ni_pc_type hypre -fieldsplit_ni_pc_hypre_type euclid -fieldsplit_TEBV_ksp_type fgmres -fieldsplit_TEBV_pc_type fieldsplit -fieldsplit_TEBV_pc_fieldsplit_type schur -fieldsplit_TEBV_pc_fieldsplit_schur_precondition selfp -fieldsplit_TEBV_fieldsplit_tau_ksp_type gmres -fieldsplit_TEBV_fieldsplit_tau_pc_type bjacobi -fieldsplit_TEBV_fieldsplit_EBV_ksp_type preonly -fieldsplit_TEBV_fieldsplit_EBV_pc_type fieldsplit -fieldsplit_TEBV_fieldsplit_EBV_pc_fieldsplit_type multiplicative -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_EP_ksp_type gmres -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_EP_pc_type hypre -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_ksp_type preonly -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_pc_type lu -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_pc_factor_mat_solver_type mumps -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_icntl_6 2 -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_icntl_24 1 -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_cntl_1 1e-5 -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_cntl_3 1e-5 -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_mumps_icntl_14 5000 -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_superlu_dist_replacetinypivot -ksp_rtol 1e-6 -ksp_norm_type unpreconditioned -snes_max_funcs 100000000000 -snes_stol 1e-20 -snes_rtol 1e-4 -jtype 2 -snes_max_it 20 -snes_ksp_ew -snes_ksp_ew_version 3 -snes_ksp_ew_rtol0 0.2 -snes_ksp_ew_rtolmax 0.9 -snes_ksp_ew_gamma 0.9 -snes_ksp_ew_alpha 1.5 -snes_ksp_ew_alpha2 1.5 -snes_ksp_ew_threshold 0.1 -fieldsplit_ni_ksp_converged_reason -fieldsplit_TEBV_ksp_converged_reason -fieldsplit_TEBV_fieldsplit_tau_ksp_converged_reason -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_EP_ksp_converged_reason -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_ksp_converged_reason -fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_fieldsplit_V_ksp_converged_reason -mat_coloring_type lf -snes_lag_jacobian 1 -snes_lag_preconditioner 1 -snes_linesearch_type nleqerr -dampV 1.0 -Re 200.0 -mat_mffd_err 1e-4 -etaout 1.30288e-3 -etawall 4.4e-2 -etaplasma 0.0006101792698927684 -etasepwal 0.0006101792698927684 -etaVV 1.30288e-6 -dump 0 -L0 2.0 -ts_fd_color -ts_fd_color_use_mat -removezero -prestep 1 -dtLA 0.05 -dtCD 50.0 -gamma_BC_fraq 0.5 -RECurrentSeedFraction 0.02 -nParticleSeed 1 -nParticleMax 32000 -DoParticleOrbitConservationTest 0 -delay_kinetic 0 -tempdump 1 -dumpfreq 1");

  /* Initialize user application context */
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Magnetic diffusion options", ""); {
	     // Set all bytes in 'example' to 0
    memset(&user.kctx, 0, sizeof(user.kctx));
    user.mi = 1.67e-27;
    PetscOptionsReal("-mi", "Ion mass", "", user.mi, & user.mi, NULL);
    user.density = 1e+20;
    PetscOptionsReal("-density", "Initial Density", "", user.density, & user.density, NULL);
    user.B0 = 5.0;
    PetscOptionsReal("-B0", "Initial magnetic field magnitude", "", user.B0, & user.B0, NULL);
    user.V_A = 1.1e+7;
    PetscOptionsReal("-V_A", "Alfven speed", "", user.V_A, & user.V_A, NULL);
    user.L0 = 2.0;
    PetscOptionsReal("-L0", "Characteristic length", "", user.L0, & user.L0, NULL);
    user.mu0 = 0.00000125663;
    PetscOptionsReal("-mu0", "Permeability", "", user.mu0, & user.mu0, NULL);
    user.eta0 = (user.L0 * user.V_A * user.mu0);
    user.eta = 1.0 ;
    PetscOptionsReal("-eta", "Resistivity", "", user.eta, & user.eta, NULL);
    user.etawall = 0.044;
    PetscOptionsReal("-etawall", "Resistivity inside the wall", "", user.etawall, & user.etawall, NULL);
    user.etaplasma = 1.71e-5;
    PetscOptionsReal("-etaplasma", "Resistivity inside the plasma chamber", "", user.etaplasma, & user.etaplasma, NULL);
    user.etawallperp = user.etawall;
    PetscOptionsReal("-etawallperp", "Poloidal resistivity in between inner Vacuum Vessel wall and plasma wall", "", user.etawallperp, & user.etawallperp, NULL);
    user.etawallphi = user.etawall;
    PetscOptionsReal("-etawallphi", "Toroidal resistivity between inner Vacuum Vessel wall and plasma wall", "", user.etawallphi, & user.etawallphi, NULL);
    user.etawallphi_isol_cell = user.etawall;
    PetscOptionsReal("-etawallphi_isol_cell", "Lower toroidal resistivity in isolated cells between inner Vacuum Vessel wall and plasma wall", "", user.etawallphi_isol_cell, & user.etawallphi_isol_cell, NULL);
    user.etasepwal = user.etaplasma;
    PetscOptionsReal("-etasepwal", "Resistivity between the separatrix and the wall", "", user.etasepwal, & user.etasepwal, NULL);
    user.etaVV = 1.30288 * 1e-6;
    PetscOptionsReal("-etaVV", "Resistivity inside the Vacuum Vessel", "", user.etaVV, & user.etaVV, NULL);
    user.etaout = 1e-20;
    PetscOptionsReal("-etaout", "Resistivity outside the wall", "", user.etaout, & user.etaout, NULL);
    user.dampV = 0.01;
    PetscOptionsReal("-dampV", "Stabilization coefficient in the velocity constraint", "", user.dampV, & user.dampV, NULL);
    user.Re = 1.0;
    PetscOptionsReal("-Re", "Reynolds parameter in the constraint : - (\nabla x B) x B + (V.\nabla) V - 1/(Re) (\nabla^2 V) = 0", "", user.Re, & user.Re, NULL);
    user.rmin = 3.05;
    PetscOptionsReal("-rmin", "Minimum radius", "", user.rmin, & user.rmin, NULL);
    user.rmax = 9.95;
    PetscOptionsReal("-rmax", "Maximum radius", "", user.rmax, & user.rmax, NULL);
    user.phimin = 0.0;
    PetscOptionsReal("-phimin", "Minimum azimuth", "", user.phimin, & user.phimin, NULL);
    user.phimax = 2.0 * PETSC_PI;
    PetscOptionsReal("-phimax", "Maximum azimuth", "", user.phimax, & user.phimax, NULL);
    user.zmin = -5.95;
    PetscOptionsReal("-zmin", "Minimum height", "", user.zmin, & user.zmin, NULL);
    user.zmax = 5.95;
    PetscOptionsReal("-zmax", "Maximum height", "", user.zmax, & user.zmax, NULL);
    user.dt = 500.0;
    PetscOptionsReal("-dt", "Length of timestep", "", user.dt, & user.dt, NULL);
    user.ictype = 9;
    PetscOptionsInt("-ictype", "Type of initial conditions", "", user.ictype, & user.ictype, NULL);
    user.Nr = 100;
    PetscOptionsInt("-Nr", "Global number of grid points in r direction", "", user.Nr, & user.Nr, NULL);
    user.Nphi = 2;
    PetscOptionsInt("-Nphi", "Global number of grid points in phi direction", "", user.Nphi, & user.Nphi, NULL);
    user.Nz = 200;
    PetscOptionsInt("-Nz", "Global number of grid points in z direction", "", user.Nz, & user.Nz, NULL);
    user.itime = 0.0;
    PetscOptionsReal("-itime", "Initial time", "", user.itime, & user.itime, NULL);
    user.ftime = user.itime + 1.0 * user.dt;
    PetscOptionsReal("-ftime", "Final time", "", user.ftime, & user.ftime, NULL);
    user.phibtype = 1;
    PetscOptionsInt("-phibtype", "Boundary type in phi direction", "", user.phibtype, & user.phibtype, NULL);
    user.dr = (user.rmax - user.rmin) / (user.Nr*user.L0);
    user.dphi = (user.phimax - user.phimin) / user.Nphi;
    user.dz = (user.zmax - user.zmin) / (user.Nz*user.L0);
    user.pred_loop = 0;
    PetscOptionsInt("-pred_loop", "Number of times the predictor step is being executed", "", user.pred_loop, & user.pred_loop, NULL);
    user.tstype = 2;
    PetscOptionsInt("-tstype", "Timestepping method: 1 for Forward Euler, 2 for Backward Euler, 3 for Crank-Nicholson, 4 for Additive Runge-Kutta IMEX", "", user.tstype, & user.tstype, NULL);
    user.jtype = 2;
    PetscOptionsInt("-jtype", "Jacobian type (0: user provide Jacobian, 1: slow finite difference, 2: fd with coloring)", "", user.jtype, & user.jtype, NULL);
    user.adaptdt = 0;
    PetscOptionsInt("-adaptdt", "Adaptive setting of step size (0: deactivated, 1: activated)", "", user.adaptdt, & user.adaptdt, NULL);
    user.debug = 0;
    PetscOptionsInt("-debug", "Debug mode (0: deactivated, 1: activated)", "", user.debug, & user.debug, NULL);
    user.dump = 0;
    PetscOptionsInt("-dump", "Dump mode (0: deactivated, 1: activated)", "", user.dump, & user.dump, NULL);
    user.EnableRelaxation = 0;
    PetscOptionsInt("-EnableRelaxation", "Relaxation mode (0: deactivated, 1: activated)", "", user.EnableRelaxation, & user.EnableRelaxation, NULL);
    user.EnableReadICFromBinary = 1;
    PetscOptionsInt("-EnableReadICFromBinary", "Read from binary mode (0: deactivated, 1: activated)", "", user.EnableReadICFromBinary, & user.EnableReadICFromBinary, NULL);
    user.prestep = 0;
    PetscOptionsInt("-prestep", "Flag for activating the prestep computations to approximate the runaway current contribution ", "", user.prestep, & user.prestep, NULL);
    user.savecoords = 0;
    PetscOptionsInt("-savecoords", "Flag for saving coordinates of cell centers, face centers and edge centers in .m files: if 0 the save option is disabled, otherwise the save option is enabled", "", user.savecoords, & user.savecoords, NULL);
    user.savesol = 0;
    PetscOptionsInt("-savesol", "Flag for saving values of magnetic fields at face centers in .m files: if 0 the save option is disabled, otherwise the save option is enabled. The generated files are: vecC.m (contains cell center coordinates), vecEphi.m (contains edge center coordinates for edges in phi direction,i.e., parallel to e_phi), vecEr.m (contains edge center coordinates for edges in r direction, i.e., parallel to e_r), vecEz.m (contains edge center coordinates for edges in z direction, i.e., parallel to e_z), vecFr.m (contains face center coordinates for faces whose normals are parallel to e_r), vecFz.m (contains face center coordinates for faces whose normals are parallel to e_z), vecFphi.m (contains face center coordinates for faces whose normals are parallel to e_phi)", "", user.savesol, & user.savesol, NULL);
    user.isB = NULL;
    user.isEP = NULL;
    user.istau = NULL;
    user.isV = NULL;
    user.isni = NULL;
    user.isB_boundary = NULL;
    user.isE_boundary = NULL;
    user.isni_boundary = NULL;

    user.delay_kinetic = 3;
    PetscOptionsInt("-delay_kinetic", "Delay kinetic solver by this many timesteps", "", user.delay_kinetic, & user.delay_kinetic, NULL);

    PetscOptionsReal("-dtLA", "Large Angle Collision timestep", "", user.kctx.dt_LA, &user.kctx.dt_LA, NULL);
    PetscOptionsReal("-dtCD", "Current Deposit timestep", "", user.kctx.dt_CD, &user.kctx.dt_CD, NULL);
    PetscOptionsReal("-RECurrentSeedFraction", "Runaway Electron current seed wrt total current", "", user.kctx.RECurrentSeedFraction, &user.kctx.RECurrentSeedFraction, NULL);
    PetscOptionsReal("-gamma_BC_fraq", "fraction between secondary mean", "", user.kctx.g_BC_fraq, &user.kctx.g_BC_fraq, NULL);
    PetscOptionsInt("-nParticleSeed", "Initial number of particles", "", user.kctx.nParticleSeed, &user.kctx.nParticleSeed, NULL);
    PetscOptionsInt("-nParticleMax", "Max number of particles", "", user.kctx.nParticleMax, &user.kctx.nParticleMax, NULL);

    char output_folder[PETSC_MAX_PATH_LEN] = ".";
    PetscOptionsString("-output", "Output folder", "", output_folder, output_folder, PETSC_MAX_PATH_LEN, NULL);

    /* RuKS_chdir(output_folder); */
    chdir(output_folder);

    // Set default location for input data.
    snprintf(user.input_folder, PETSC_MAX_PATH_LEN, "../../inputs");
    PetscOptionsString("-input", "Input folder", "", user.input_folder, user.input_folder, PETSC_MAX_PATH_LEN, NULL);
    strcpy(user.kctx.input_folder, user.input_folder);

    user.kctx.DoParticleOrbitConservationTest = 0;
    PetscOptionsInt("-DoParticleOrbitConservationTest", "Switch it do particle orbit conservation test", "",user.kctx.DoParticleOrbitConservationTest, &user.kctx.DoParticleOrbitConservationTest, NULL);

    user.kctx.DoAdvance = 1;
    PetscOptionsInt("-DoAdvance", "Switch it do advance", "",user.kctx.DoAdvance, &user.kctx.DoAdvance, NULL);

    user.kctx.DoPoincare = 0;
    PetscOptionsInt("-DoPoincare", "Switch it do Poincare", "",user.kctx.DoPoincare, &user.kctx.DoPoincare, NULL);

    user.kctx.DumpFields = 0;
    PetscOptionsInt("-DumpFields", "Switch it output fields", "",user.kctx.DumpFields, &user.kctx.DumpFields, NULL);



    user.kctx.R0 = user.rmin / user.L0;
    user.kctx.Z0 = user.zmin / user.L0;
    user.kctx.dR = user.dr;
    user.kctx.dZ = user.dz;
    user.kctx.NR = user.Nr;
    user.kctx.NZ = user.Nz;
    user.kctx.timeSeed = user.itime;

    user.axis[0] = (user.rmin + user.rmax) / user.L0 / 2.0;
    user.axis[1] = (user.zmin + user.zmax) / user.L0 / 2.0;



    PetscInt numC = 0;
    char filename[PETSC_MAX_PATH_LEN];

    if (user.ictype == 12 || user.ictype == 13 || user.ictype == 10) {

      FILE * pvdfile;

      PetscSNPrintf(filename, sizeof(filename), "%s/veceta_1layer_%.2Dx%.1Dx%.2D.txt", user.input_folder, user.Nr, user.Nphi, user.Nz);
      ReadInitialData( & (user.dataC), & numC, filename);
      }
    if ((user.ictype == 9 || user.ictype == 11) && ! user.savecoords) {
      PetscSNPrintf(filename, sizeof(filename), "%s/veceta_grid%.3Dx%.2Dx%.3D.txt", user.input_folder, user.Nr, user.Nphi, user.Nz);
      ReadInitialData( & (user.dataC), & numC, filename);

      if(0){
        PetscSNPrintf(filename, sizeof(filename), "%s/vecBr.txt", user.input_folder);
        ReadInitialData( & (user.datar), & (user.numr), filename);
        PetscSNPrintf(filename, sizeof(filename), "%s/vecBphi.txt", user.input_folder);
        ReadInitialData( & (user.dataphi), & (user.numphi), filename);
        PetscSNPrintf(filename, sizeof(filename), "%s/vecBz.txt", user.input_folder);
        ReadInitialData( & (user.dataz), & (user.numz), filename);
      }
      else{
        //ReadInitialData( & (user.datapsi), & (user.numpsi), "InitialData/vecpsi.txt");
        PetscSNPrintf(filename, sizeof(filename), "%s/vecpsi_grid%.3Dx%.2Dx%.3D.txt", user.input_folder, user.Nr, user.Nphi, user.Nz);
        ReadInitialData( & (user.datapsi), & (user.numpsi), filename);

        //ReadInitialData( & (user.datag), & (user.numg), "InitialData/vecg.txt");
        PetscSNPrintf(filename, sizeof(filename), "%s/vecg_grid%.3Dx%.2Dx%.3D.txt", user.input_folder, user.Nr, user.Nphi, user.Nz);
        ReadInitialData( & (user.datag), & (user.numg), filename);

        PetscPrintf(PETSC_COMM_WORLD, "Read initial data from g and psi input!\n");
      }
      if (user.debug) {
        /*These prints are just for debugging*/
        PetscPrintf(PETSC_COMM_WORLD, "numr = %d\n", user.numr);
        PetscPrintf(PETSC_COMM_WORLD, "numphi = %d\n", user.numphi);
        PetscPrintf(PETSC_COMM_WORLD, "numz = %d\n", user.numz);
      }
    }
    user.Ebc = 0;
    PetscOptionsInt("-Ebc", "Type of boundary conditions for E field (0: E_tangential(t)=0, 1: E_tangential(t)=E_tangential(t=0))", "", user.Ebc, & user.Ebc, NULL);
    user.tempdump = 0;
    PetscOptionsInt("-tempdump", "Flag for saving intermediate output in .dat files", "", user.tempdump, & user.tempdump, NULL);
    user.dumpfreq = (int) PetscCeilReal((user.ftime) / (10.0 * user.dt));
    PetscOptionsInt("-dumpfreq", "Frequency for saving intermediate output in .dat files", "", user.dumpfreq, & user.dumpfreq, NULL);
    user.testSpGD = 0;
    PetscOptionsInt("-testSpGD", "Flag for testing S the Schur complement matrix vs S + GD the Schur complement matrix augmented by the gradient of derived divergence operator", "", user.testSpGD, & user.testSpGD, NULL);
    user.testSpGDsamerhs = 0;
    PetscOptionsInt("-testSpGDsamerhs", "Flag for setting the right hand side vector while testing S the Schur complement matrix vs S - dt*GD the Schur complement matrix augmented by the gradient of derived divergence operator: 0 for a different right-hand side (S - dt*GD)*E, 1 for the same right-hand side S*E", "", user.testSpGDsamerhs, & user.testSpGDsamerhs, NULL);

    user.oldstep = 0;

    PetscOptionsInt("-oldstep", "Last step in previous simulation when a restart is used", "", user.oldstep, & user.oldstep, NULL);
    user.n_record =0;
    user.n_record_Steady_jRE =0;
    user.CorrectorIdentifier = 1;
  }
  PetscOptionsEnd();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RUNAWAY SOLVER INITIALIZATION
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscPrintf(PETSC_COMM_WORLD, "Initializing Runaway solver...");
  // Creating runaway kinetic solver instance
  /* RuKS_init(&(user.kctx), &(user.runaway_solver) ); */

  // Allocate memory for runaway current density
  // First half for timestep average, second half for delta-tail average
  user.jre_data = (double*) malloc(sizeof(double) * (user.Nr) * (user.Nz) * 3 * 2);
  for (int i = 0; i < user.Nr * user.Nz * 3 * 2; ++i) {
    user.jre_data[i] = 0.0;
  }
  user.jre = user.jre_data;
  user.jreR   = user.jre;
  user.jrephi = user.jre +   (user.Nr) * (user.Nz);
  user.jreZ   = user.jre + 2*(user.Nr) * (user.Nz);


#if 0

  //RuKS_setRECurrentPtr(user.runaway_solver, user.jre_data); */

  if(user.itime > 0 && user.oldstep > 0){
    char filename[PETSC_MAX_PATH_LEN];
    PetscSNPrintf(filename, sizeof(filename), "%s/input.h5", user.input_folder);
    /* RuKS_H5Read(user.runaway_solver, filename); */
  }
  else{
    for (int ii = 0; ii <  (user.Nr) * (user.Nz) * 3 * 2; ++ii){
      user.jre_data[ii] = 0.0;
   }
  }
#endif
  PetscMPIInt rank;
  PetscMPIInt size;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  int nphi = 1, nt = 2, ndims = 3, nfields = 4;
  double* raw_fields = (double*) malloc(user.Nr * user.Nz * nphi * nt * ndims * nfields * sizeof(double));
  user.field_data = raw_fields;
  user.poincare_counter = 0;
  user.field_counter = 0;

  hflux_kokkos_init();
  hflux_init(user.Nr, user.Nz, nfields, nphi, nt, user.rmin / user.L0, user.zmin / user.L0, user.dr, user.dz, &user.field_interpolation);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    RUNAWAY SOLVER INITIALIZATION END
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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

    if (user.phibtype) {
      DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE, user.Nr, user.Nphi, user.Nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, & da);
    } else {
      DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, user.Nr, user.Nphi, user.Nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, dof0, dof1, dof2, dof3, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, & da);
    }
    DMSetFromOptions(da);
    DMSetUp(da);

    DMStagGetNumRanks(da, & nr, & nphi, & nz);

    if (user.phibtype) {
      DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_NONE, user.Nr, user.Nphi, user.Nz, nr, nphi, nz, 4, 1, 1, 1, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, & user.coorda);
    } else {
      DMStagCreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, user.Nr, user.Nphi, user.Nz, nr, nphi, nz, 4, 1, 1, 1, DMSTAG_STENCIL_BOX, stencilWidth, NULL, NULL, NULL, & user.coorda);
    }

    DMSetFromOptions(user.coorda);
    DMSetUp(user.coorda);
    DMStagSetUniformCoordinatesExplicit(da, user.rmin/user.L0, user.rmax/user.L0, user.phimin, user.phimax, user.zmin/user.L0, user.zmax/user.L0);
    DMStagSetUniformCoordinatesExplicit(user.coorda, user.rmin/user.L0, user.rmax/user.L0, user.phimin, user.phimax, user.zmin/user.L0, user.zmax/user.L0);


    DMSetApplicationContext(da, &user);
    DMCreateGlobalVector(da, & user.oldX);
    DMCreateGlobalVector(da, & user.X_star);
  }
  /* Print out some info */
  {
    PetscInt N[3];
    DMStagGetGlobalSizes(da, & N[0], & N[1], & N[2]);
    PetscPrintf(PETSC_COMM_WORLD, "Using a %D x %D x %D mesh\n", N[0], N[1], N[2]);
    PetscPrintf(PETSC_COMM_WORLD, "dr: %g\n", user.dr);
    PetscPrintf(PETSC_COMM_WORLD, "dphi: %g\n", user.dphi);
    PetscPrintf(PETSC_COMM_WORLD, "dz: %g\n", user.dz);
    PetscPrintf(PETSC_COMM_WORLD, "normalized dt: %g\n", user.dt);
    PetscPrintf(PETSC_COMM_WORLD, "non-normalized dt: %g\n", user.dt * user.L0/user.V_A); // Alfven time := L0 / V_A
    PetscPrintf(PETSC_COMM_WORLD, "Characteristic resistive time: %g\n", user.L0*user.L0*user.mu0/user.etaplasma); // tau_eta := mu0 L0^2 / non_normalized_eta
    PetscPrintf(PETSC_COMM_WORLD, "Reynolds parameter for viscosity: %g\n", user.Re);
    PetscPrintf(PETSC_COMM_WORLD, "Lundquist number: %g\n", user.eta0 / user.etaplasma); // tau_eta / Alfven time = mu0 L0 V_A / non_normalized_eta
    PetscPrintf(PETSC_COMM_WORLD, "Resistivity inside the plasma chamber (in Ohm.meter): %g\n",user.etaplasma);
    PetscPrintf(PETSC_COMM_WORLD, "Resistivity outside the vacuum vessel (in Ohm.meter): %g\n", user.etaout );
    PetscPrintf(PETSC_COMM_WORLD, "Resistivity inside the vacuum vessel (in Ohm.meter): %g\n", user.etaVV );
    PetscPrintf(PETSC_COMM_WORLD, "Resistivity inside the blanket module (in Ohm.meter): %g\n", user.etawall );
    //PetscPrintf(PETSC_COMM_WORLD, "CFL value: %g\n", user.dt / PetscMin(user.dr,user.dz));
  }
  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Extract global vectors from DM;
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  DMCreateGlobalVector(da, & X);
  VecZeroEntries(X);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create timestepping solver context
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSCreate(PETSC_COMM_WORLD, & user.ts);
  TS ts = user.ts;
  TSSetDM(ts, da);

  if (user.savecoords) {
    SaveCoordinates(ts, & user);
    return (0);
  }

  //DMSetMatrixPreallocateOnly(da,PETSC_TRUE);

  //TSSetProblemType(ts, TS_LINEAR);
  TSMonitorSet(ts, Monitor, & user, NULL); /* Set optional user-defined monitoring routine */

  switch (user.tstype) {
  case 1:
    TSSetType(ts, TSEULER); /* Forward Euler method */
    break;
  case 2:
    TSSetType(ts, TSBEULER); /* Backward Euler method */
    break;
  case 3:
    TSSetType(ts, TSCN); /* Crank-Nicholson method */
    break;
  case 4:
    TSGetAdapt(ts, & adapt);
    TSAdaptSetType(adapt, TSADAPTNONE);
    TSSetType(ts, TSARKIMEX); /* Additive Runge-Kutta IMEX method */
    TSARKIMEXSetFullyImplicit(ts, PETSC_TRUE);
    TSARKIMEXSetType(ts, TSARKIMEXL2);

    TSSetEquationType(ts,TS_EQ_IMPLICIT);
    break;
  case 5:
    //TSSetProblemType(ts,TS_NONLINEAR);
    TSSetType(ts,TSROSW);
    break;
  case 6:
    TSSetType(ts, TSBDF); /* Backward differentiation formula of order 2*/
    TSBDFSetOrder(ts,2);
    break;
  case 7:
    TSSetType(ts, TSTHETA); /* Implicit Theta method */
    TSThetaSetTheta(ts, 0.5); // Default theta is 0.5 but this value can be set through command line using the -ts_theta_theta flag
    break;
  default:
    TSSetType(ts, TSEULER); /* Forward Euler method */
    break;
  }

  TSGetSNES(ts, & snes);
  if (user.tstype > 1) {
    PetscOptionsGetBool(NULL, NULL, "-snes_mf", & matrix_free, NULL);
    PetscOptionsGetBool(NULL, NULL, "-snes_mf_operator", & matrix_free_FDprec, NULL);
    PetscOptionsGetBool(NULL, NULL, "-removezero", &removezero, NULL);
    if (matrix_free) { //matrix_free mode without any preconditioning matrix
      PetscPrintf(PETSC_COMM_WORLD, "======use matrix-free evaluation and no preconditioning======\n");
    } else if (matrix_free_FDprec) { //matrix_free mode with colored finite difference jacobian for preconditioning
      DMCreateMatrix(da, & J);
      TSSetIJacobian(ts, J, J, TSComputeIJacobianDefaultColor, NULL);
      PetscPrintf(PETSC_COMM_WORLD, "======use matrix-free evaluation and FD coloring Jacobian for preconditioning======\n");
    } else {
      DMCreateMatrix(da, & J);
      DMCreateMatrix(da, & Jpre);
      if (user.jtype == 0) {
        TSSetIJacobian(ts, J, Jpre, FormIJacobian_BImplicit, & user); /* use user provided Jacobian evaluation routine */
        PetscPrintf(PETSC_COMM_WORLD, "======use Analytical Jacobian======\n");
      } else {
        /* use finite difference Jacobian J as preconditioner and '-snes_mf_operator' for Mat*vec */
        /*MatCreateSNESMF(snes,&Jmf);*/
        if (user.jtype == 1) {
          /* slow finite difference J; */
          SNESSetJacobian(snes, J, J, SNESComputeJacobianDefault, PETSC_NULLPTR);
          PetscPrintf(PETSC_COMM_WORLD, "======use FD Jacobian======\n");
        } else if (user.jtype == 2) {
          /* Use coloring to compute finite difference J efficiently */
          TSSetIJacobian(ts, J, J, TSComputeIJacobianDefaultColor, PETSC_NULLPTR);
          PetscPrintf(PETSC_COMM_WORLD, "======use FD coloring Jacobian======\n");
        } else {
          SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "This jtype is not supported");
        }
      }
    }
    //TSSetIFunction(ts, NULL, FormIFunction_DampingV, & user);
    TSSetIFunction(ts, NULL, FormIFunction_Vperp_viscosity, & user);
    TSSetRHSFunction(ts, NULL, FormRHSFunction_BImplicit, & user);
  } else {
    TSSetRHSFunction(ts, NULL, FormRHSFunction_BImplicit, & user);
  }
  SNESSetFromOptions(snes);

  TSSetTime(ts, user.itime);
  ftime = user.ftime;
  TSSetMaxTime(ts, ftime);
  TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER);
  TSSetSolution(ts, X);
  TSSetTimeStep(ts, user.dt);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set runtime options
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if(user.prestep){
    TSSetPreStep(ts, Update_J_RE);
  }
  TSSetFromOptions(ts);
  if(user.adaptdt){
  // Adaptive Time Step Controller
    TSGetAdapt(ts,&adapt);
    TSAdaptSetScaleSolveFailed(adapt, 0.5);
    adapt->ops->choose = TSAdaptChoose_user;
    TSAdaptSetMonitor( adapt, PETSC_TRUE);
    TSSetMaxSNESFailures( ts, -1);
  }
  TSSetUp(ts);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set index sets
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  TSGetTime(ts, & time);

  IS isEPBVndup, isEPtauVndup, isEPtauBVdup, istauBVndup, isndup, isEPdup, isBdup, istaudup;

  DMCreateGlobalVector(da, & dummyX);
  VecCopy(X, dummyX);
  DMCreateMatrix(da, & dummyJ);

  if (user.tstype > 1 && matrix_free_FDprec) {
    PetscRandom rctx;
    PetscCall(PetscRandomCreate(PETSC_COMM_WORLD, &rctx));
    PetscCall(PetscRandomSetInterval(rctx, 1.0, 2.0));
    PetscCall(VecSetRandom(X, rctx));
    PetscCall(VecSetRandom(dummyX, rctx));
    PetscCall(PetscRandomDestroy(&rctx));
    PetscCall(TSComputeIJacobian(ts, 0.0, X, dummyX, 2.0, J, J, PETSC_FALSE));
    if (removezero) PetscCall(TSPruneIJacobianColor(ts, J, J));
  }

  KSPCreate(PETSC_COMM_WORLD, & dummykspEP);
  KSPSetOptionsPrefix(dummykspEP, "sepKSP_");
  FormDummyIJacobian4(ts, time, dummyX, dummyX, 1.0 / user.dt, dummyJ, dummyJ, & user);
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
  CHKERRQ(ierr);
  KSPSetTolerances(subksp[1], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1); //use 1 inner iteration maximum for the dummy solve
  ierr = KSPSolve(dummykspEP, dummyX, dummyX);
  CHKERRQ(ierr);
  //ISDuplicate(isBdup, & user.isB);
  PCFieldSplitGetISByIndex(dummypcEP, 0, & isEPdup);
  ISDuplicate(isEPdup, & user.isEP);


  VecDestroy( & dummyX);
  KSPDestroy( & dummykspEP);
  MatDestroy( & dummyJ);


  char ** namelist;
  IS * islist, isALL, isALL_V, isBV;
  IS ISV;

  PetscInt len, d = 0;
  ierr = DMCreateFieldDecomposition(da, & len, & namelist, & islist, NULL);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "The number of subproblems in the field decomposition is: %g\n", (double)(len));
  for (d = 0; d < len; ++d) {
    PetscPrintf(PETSC_COMM_WORLD, "The name of field number %d is: %s.\n", d, namelist[d]);
    //PetscPrintf(PETSC_COMM_WORLD, "The global indices for field number %d are as follows.\n", d);
    //ISView(islist[d],PETSC_VIEWER_STDOUT_SELF);
  }
  PetscBool flagV = PETSC_FALSE, flagE = PETSC_FALSE, flagF = PETSC_FALSE, flagC = PETSC_FALSE;
  ISDifference(islist[0], user.isEP, & user.isV);

  ISDuplicate(islist[1], & user.istau);
  ISDuplicate(islist[2], & user.isB);
  ISDuplicate(islist[3], & user.isni);

  const IS islist2[5] = {user.isV, user.isEP, user.istau, user.isB, user.isni};

  ISConcatenate(PETSC_COMM_WORLD,5,islist2,&isALL);
  ISDifference(isALL, user.isV, & isALL_V);

  const IS islist3[5] = {user.isV, user.isB};
  ISConcatenate(PETSC_COMM_WORLD,2,islist3,&isBV);

  for (d = 0; d < len; ++d) {
    ISDestroy( & islist[d]);
  }
  PetscFree(islist);
  PetscFree(namelist);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set preconditioner options
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  SNESGetKSP(snes, & ksp);
  ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
  /* Set a user-defined "shell" preconditioner if desired */
  PetscOptionsGetBool(NULL,NULL,"-user_defined_pc",&user_defined_pc,NULL);
  if (user_defined_pc) {
    /* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
    PCSetType(pc,PCSHELL);
    PCShellSetContext(pc,&user);

    /* Do any setup required for the preconditioner */
    PCShellSetSetUp(pc,SampleShellPCSetUp);

    /* (Required) Set the user-defined routine for applying the preconditioner */
    PCShellSetApply(pc,SampleShellPCApply);

    /* (Optional) Set user-defined function to free objects used by custom preconditioner */
    PCShellSetDestroy(pc,SampleShellPCDestroy);

    /* (Optional) Set a name for the preconditioner, used for PCView() */
    PCShellSetName(pc,"ShellPrec");
  }
  else{/* first level -> split ni from the rest : {ni}, {V Phi tau B}
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

   PetscCall(DMStagCreateISFromStencils(da,1,stencil0,&is[0]));
   PetscCall(DMStagCreateISFromStencils(da,10,stencil1,&is[1]));

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
       PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
       PetscCall(KSPGetPC(sub_ksp[1],&pc_notc));
       PetscCall(PetscFree(sub_ksp));

       PetscCall(DMStagCreateCompatibleDMStag(da,4,1,1,0,&dm_notc));

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
       PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
       PetscCall(KSPGetPC(sub_ksp[1],&pc_noe));
       PetscCall(PetscFree(sub_ksp));

       PetscCall(DMStagCreateCompatibleDMStag(da,4,0,1,0,&dm_noe));

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
       PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
       PetscCall(KSPGetPC(sub_ksp[1],&pc_noe_2));
       PetscCall(PetscFree(sub_ksp));

       PetscCall(DMStagCreateCompatibleDMStag(da,3,0,1,0,&dm_notv));

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

   /* If the fieldsplit PC wasn't overridden, further split the first split of the third level */
//   if (0 && is_fieldsplit) {
//     PCType pc_type;
//
//     PetscCall(PCGetType(pc_noe_2,&pc_type));
//     PetscCall(PetscStrcmp(pc_type,PCFIELDSPLIT,&is_fieldsplit));
//     if (is_fieldsplit) {
//       Vec coordLocal;
//       DM dmCoord;
//       PetscScalar ** ** arrCoord;
//
//       Mat subG, G;
//       DM            dm_notv, dm_3;
//       KSP           *sub_ksp;
//       PC            pc_noe_3;
//       PetscInt      n_splits;
//       DMStagStencil stencil_notv_faces[3], stencil_notv_notfaces[3];
//       IS            is_notv[2];
//       const char    *name_notv[2] = {"V", "B"};
//
//          PetscCall(PCSetUp(pc_noe_2)); // Set up the Fieldsplit PC
//          PetscCall(PCFieldSplitGetSubKSP(pc_noe_2,&n_splits,&sub_ksp));
//          PetscAssert(n_splits == 2,PetscObjectComm((PetscObject)da),PETSC_ERR_SUP,"Expected a Fieldsplit PC with two fields");
//          PetscCall(KSPGetPC(sub_ksp[1],&pc_noe_3));
//          PetscCall(PetscFree(sub_ksp));
//
//          PetscCall(DMStagCreateCompatibleDMStag(da,0,0,1,0,&dm_notv));
//
//          DMGetCoordinateDM(dm_notv, & dmCoord);
//          DMGetCoordinatesLocal(dm_notv, & coordLocal);
//          DMStagVecGetArrayRead(dmCoord, coordLocal, & arrCoord);
//          PCSetCoordinates(pc_noe_3, 3, 3, arrCoord);
//
//
//         PetscCall(DMStagCreateCompatibleDMStag(da,3,1,0,0,&dm_3));
//
//
//          FormDiscreteGradient(ts, G, & user);
//         MatCreateSubMatrix(G,user.istau,user.isV,MAT_INITIAL_MATRIX,&subG);
//
//
//          PetscCall(DMDestroy(&dm_notv));
//     }
//   }

  }


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Set initial conditions
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  // FormInitialpsi(ts, X, & user);
  FormInitialSolution_psi(ts, X, & user);

  //Update_J_RE(ts);
  if(0){
    Vec F,C;
    VecDuplicate(X, & F);
    VecZeroEntries(F);
    VecDuplicate(X, & C);
    VecZeroEntries(C);
    //NEED TO CREATE VEC C AND RECONSTRUCT FROM EDGES TO CELLS THE THREE COMPONENTS OF E BEFORE CALLING THE SCATTERING
    //FormElectricField(ts, X, F, &user);
    //DumpEdgeField(ts, 2, X, &user);
    //DumpEdgeField(ts, 3, F, &user);
    //EdgeToCellReconstruction_r(ts,F,C,&user);
    DumpSolution_Cell(ts, 3, X, &user);
    VecDestroy(&C);
    VecDestroy(&F);
  }

  //PetscFinalize();
  //return(0);

  //Could make a copy of X now: CopyX0 := X
  if(0){
    PetscFree(user.datar);
    PetscFree(user.dataphi);
    PetscFree(user.dataz);
  }
  else{
    //PetscFree(user.datag);
    //PetscFree(user.datapsi);
  }


  if (0) {
    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Normalize the system
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    Vec B, E, ni, Vi_perp;
    VecGetSubVector(X, user.isB, & B);
    VecScale(B, 0.2); // B := tilde{B} = (B_0^-1) B
    VecRestoreSubVector(X, user.isB, & B);
    VecGetSubVector(X, user.istau, & E);
    VecScale(E, 0.2); // E := tilde{E} = (E_0^-1) E
    VecScale(E, 1.0 / 11000000); // E := tilde{E} = (E_0^-1) E
    VecRestoreSubVector(X, user.istau, & E);
    VecGetSubVector(X, user.isni, & ni);
    VecScale(ni, 1e-20); // ni := tilde{ni} = (ni_0^-1) ni
    VecRestoreSubVector(X, user.isni, & ni);
    VecGetSubVector(X, user.isV, & Vi_perp);
    VecScale(Vi_perp, 1.0 / 11000000); // Vi_perp := tilde{Vi_perp} = (Vi_perp_0^-1) Vi_perp
    VecRestoreSubVector(X, user.isV, & Vi_perp);
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Solve nonlinear system
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSSolve(ts, X);  CHKERRQ(ierr);


  if (user.savesol) {
    SaveSolution(ts,X,& user);
  }

  TSGetSolveTime(ts, & ftime);
  TSGetStepNumber(ts, & steps);
  TSGetConvergedReason(ts, & reason);
  PetscPrintf(PETSC_COMM_WORLD, "%s at time %g after %D steps\n", TSConvergedReasons[reason], (double) ftime, steps);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Free work space.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ISDestroy( & isALL);
  ISDestroy( & isALL_V);
  ISDestroy( & isBV);


  VecDestroy( & X);
  VecDestroy( & user.oldX);
  VecDestroy( & user.X_star);
  if (user.tstype > 1) {
    if (matrix_free) { //matrix_free mode without any preconditioning matrix
    } else if (matrix_free_FDprec) { //matrix_free mode with colored finite difference jacobian for preconditioning
      MatDestroy( & J);
    } else {
      MatDestroy( & J);
      MatDestroy( & Jpre);
    }
  }
  TSDestroy( & user.ts);

  DMDestroy( & da);
  PetscFinalize();

  hflux_destroy(user.field_interpolation);
  hflux_kokkos_finalize();
  return ierr;
}
