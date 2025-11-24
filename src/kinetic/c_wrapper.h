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

#ifndef RUNAWAY_C_WRAPPER_H_
#define RUNAWAY_C_WRAPPER_H_

#include "mhd/mfd_config.h"
#ifdef __cplusplus
extern "C" {
#endif

  int  parthenon_init(void ** man, int argc, char *argv[], User* mhd_config);
  int  runaway_init(void * man);
  void runaway_finalize(void* man);

  void runaway_push(void * man);
  void runaway_saveState(void * man);
  void runaway_restoreState(void * man);

#ifdef __cplusplus
}
#endif

#endif




