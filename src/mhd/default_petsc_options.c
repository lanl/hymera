#include <petscsys.h>
#include "default_petsc_options.h"

/* Default PETSc options for the MHD solver.
 *
 * These are inserted into the global options database at the start of
 * mhd_initialize. They configure the nonlinear and linear solver stack: a
 * matrix-free SNES whose preconditioner is built from a finite-difference
 * coloured Jacobian, wrapped in a four-level nested PCFIELDSPLIT that peels the
 * saddle-point structure apart as
 *
 *     ni  |  TEBV -> { tau | EBV -> { EP | BV -> { B | V } } }
 *
 * with a MUMPS direct factorization on the innermost BV block. The "dummyKSP_"
 * and "dummySNES_" prefixed copies configure the throwaway solvers used during
 * initial-condition construction in ts_functions.c.
 *
 * Precedence: these are defaults, so they are inserted first and the process
 * command line is re-applied afterwards. PetscInitialize has already consumed
 * argv by the time this runs, so without the re-insertion below every option
 * here would silently override anything the user passed on the command line,
 * and no solver setting could be changed without recompiling.
 */

/* Options shared by the production solver (no prefix) and by the temporary
 * solvers used to build the initial condition (prefixed). Keeping one list and
 * emitting it twice with different prefixes removes the duplication that
 * previously let the two copies drift apart -- they had already diverged:
 * snes_rtol was 1e-4 unprefixed but 1e-5 for dummySNES_.
 */
#define MHD_KSP_OPTS(p) \
  "-" p "ksp_type fgmres " \
  "-" p "ksp_rtol 1e-6 " \
  "-" p "ksp_norm_type unpreconditioned " \
  "-" p "ksp_converged_reason " \
  "-" p "ksp_monitor "

/* MUMPS controls, applied both to the outer LU and to the innermost BV block.
 *   icntl_6  2    column permutation / scaling strategy
 *   icntl_14 5000 working-space growth allowance, in percent
 *   icntl_24 1    detect and handle null pivots
 *   cntl_1/3 1e-5 relative pivoting and null-pivot thresholds
 */
#define MHD_MUMPS_OPTS(p) \
  "-" p "pc_type lu " \
  "-" p "pc_factor_mat_solver_type mumps " \
  "-" p "mat_mumps_icntl_6 2 " \
  "-" p "mat_mumps_icntl_24 1 " \
  "-" p "mat_mumps_cntl_1 1e-5 " \
  "-" p "mat_mumps_cntl_3 1e-5 " \
  "-" p "mat_mumps_icntl_14 5000 "

/* The nested fieldsplit tree. The leaf BV block is solved exactly with MUMPS;
 * the EP block uses hypre; tau uses block Jacobi; the Schur complement of the
 * TEBV block is approximated with selfp.
 */
#define MHD_FIELDSPLIT_OPTS(p) \
  "-" p "pc_fieldsplit_type multiplicative " \
  "-" p "fieldsplit_ni_ksp_type fgmres " \
  "-" p "fieldsplit_ni_pc_type hypre " \
  "-" p "fieldsplit_ni_pc_hypre_type euclid " \
  "-" p "fieldsplit_ni_ksp_converged_reason " \
  "-" p "fieldsplit_TEBV_ksp_type fgmres " \
  "-" p "fieldsplit_TEBV_pc_type fieldsplit " \
  "-" p "fieldsplit_TEBV_pc_fieldsplit_type schur " \
  "-" p "fieldsplit_TEBV_pc_fieldsplit_schur_precondition selfp " \
  "-" p "fieldsplit_TEBV_ksp_converged_reason " \
  "-" p "fieldsplit_TEBV_fieldsplit_tau_ksp_type gmres " \
  "-" p "fieldsplit_TEBV_fieldsplit_tau_pc_type bjacobi " \
  "-" p "fieldsplit_TEBV_fieldsplit_tau_ksp_converged_reason " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_ksp_type preonly " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_pc_type fieldsplit " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_pc_fieldsplit_type multiplicative " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_EP_ksp_type gmres " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_EP_pc_type hypre " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_EP_ksp_converged_reason " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_ksp_type preonly " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_mat_superlu_dist_replacetinypivot " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_ksp_converged_reason " \
  "-" p "fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_fieldsplit_V_ksp_converged_reason " \
  MHD_MUMPS_OPTS(p "fieldsplit_TEBV_fieldsplit_EBV_fieldsplit_BV_")

/* Eisenstat-Walker adaptive linear-solve tolerance, plus the matrix-free SNES
 * setup. nleqerr line search; snes_stol 1e-20 effectively disables the
 * step-length convergence test so only rtol decides.
 */
#define MHD_SNES_OPTS(p, rtol) \
  "-" p "snes_mf_operator " \
  "-" p "snes_monitor " \
  "-" p "snes_converged_reason " \
  "-" p "snes_linesearch_type nleqerr " \
  "-" p "snes_lag_jacobian 1 " \
  "-" p "snes_lag_preconditioner 1 " \
  "-" p "snes_max_funcs 100000000000 " \
  "-" p "snes_stol 1e-20 " \
  "-" p "snes_rtol " rtol " " \
  "-" p "snes_max_it 20 " \
  "-" p "snes_ksp_ew " \
  "-" p "snes_ksp_ew_version 3 " \
  "-" p "snes_ksp_ew_rtol0 0.2 " \
  "-" p "snes_ksp_ew_rtolmax 0.9 " \
  "-" p "snes_ksp_ew_gamma 0.9 " \
  "-" p "snes_ksp_ew_alpha 1.5 " \
  "-" p "snes_ksp_ew_alpha2 1.5 " \
  "-" p "snes_ksp_ew_threshold 0.1 "

PetscErrorCode default_petsc_options(void) {
  PetscFunctionBeginUser;

  /* Fixed timestep: no adaptivity, so the step sequence is deterministic and
   * reproducible run to run. The regression harness relies on this. */
  PetscCall(PetscOptionsInsertString(NULL, "-ts_adapt_type none "));

  /* Solvers used to build the initial condition (ts_functions.c). */
  PetscCall(PetscOptionsInsertString(NULL,
      MHD_KSP_OPTS("dummyKSP_")
      MHD_MUMPS_OPTS("dummyKSP_")
      MHD_FIELDSPLIT_OPTS("dummyKSP_")
      MHD_SNES_OPTS("dummySNES_", "1e-5")));

  /* The production solver. */
  PetscCall(PetscOptionsInsertString(NULL,
      MHD_KSP_OPTS("")
      MHD_MUMPS_OPTS("")
      MHD_FIELDSPLIT_OPTS("")
      MHD_SNES_OPTS("", "1e-4")));

  /* Jacobian: coloured finite differences, reused as the preconditioner for
   * the matrix-free operator. mat_mffd_err sets the differencing increment. */
  PetscCall(PetscOptionsInsertString(NULL,
      "-ts_fd_color "
      "-ts_fd_color_use_mat "
      "-mat_coloring_type lf "
      "-mat_mffd_err 1e-4 "
      "-removezero "));

  /* Re-apply the process command line so that user-supplied options take
   * precedence over every default above. PetscInitialize consumed argv before
   * this function ran, so fetch it back from PETSc. */
  {
    int argc = 0;
    char **argv = NULL;
    PetscCall(PetscGetArgs(&argc, &argv));
    if (argc > 1 && argv) {
      PetscCall(PetscOptionsInsert(NULL, &argc, &argv, NULL));
    }
  }

  PetscFunctionReturn(PETSC_SUCCESS);
}
