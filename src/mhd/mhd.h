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

#ifndef MHD_H_
#define MHD_H_

#include "mfd_config.h"

#ifdef __cplusplus
extern "C" {
#endif

int mhd_PetscInit(int *argc, char *** argv, User** user);
int mhd_initialize(User* mhd_config);
int mhd_step(User* mhd_config);


int mhd_getF(User* mhd_config, field_id fid, view3d_t v);
int mhd_resetState(User* mhd_config);
int mhd_destroy(User* mhd_config);
int mhd_savesolution(User* mhd_config, const char* filename);
int mhd_loadsolution(User* mhd_config, const char* filename);


#ifdef __cplusplus
}
#endif

/* Every function above returns a PETSc error code: 0 on success, nonzero on
 * failure. Ignoring it is not safe. PETSc's own convention is that all library
 * errors are fatal, and the solver's internal error paths use SETERRQ, which
 * records a message and unwinds. If a caller discards the code, execution
 * continues with a partly-constructed or unpopulated state; in practice that has
 * meant a missing input file producing a segfault thousands of lines later
 * instead of a message naming the file.
 *
 * MHD_CHECK evaluates the call once, and on failure prints the failing
 * expression with its source location and aborts the whole MPI job. Aborting is
 * deliberate: these entry points sit at the top of the coupled timestep, there is
 * no recovery path for a failed field solve, and continuing would corrupt the
 * kinetic side's view of the fields.
 *
 * Usage:   MHD_CHECK(mhd_step(cfg));
 */
#ifdef __cplusplus
#include <cstdio>
#include <cstdlib>
extern "C" {
#endif

/* Reports a failed MHD call and aborts. Defined in mhd.c so that both languages
 * get identical behaviour and the MPI abort path lives in one place. */
void mhd_fail(const char *expr, const char *file, int line, int code);

#ifdef __cplusplus
}
#endif

#define MHD_CHECK(call)                                                        \
  do {                                                                         \
    const int mhd_check_rc_ = (call);                                           \
    if (mhd_check_rc_ != 0) {                                                   \
      mhd_fail(#call, __FILE__, __LINE__, mhd_check_rc_);                       \
    }                                                                          \
  } while (0)

#endif
