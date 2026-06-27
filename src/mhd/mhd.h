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
#include "fid.h"


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

#endif
