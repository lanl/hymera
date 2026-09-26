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
//
// ATTIC -- NOT BUILT.
//
// This file is not built and is not part of any CMake target. It holds
// functions that had zero callers anywhere in src/, tests/, or the built
// objects (verified via nm) as of this commit. They were moved out of
// mhd.c verbatim (byte-identical bodies) to shrink that file while keeping
// the code in-tree and greppable for physics review. Neither function was
// declared in mhd.h, so no header edit was needed for this file.
//
// Retained for physics review only; scheduled for deletion in a later,
// separate commit. Known defects: see docs/mhd/physics/open-questions.md.
// Do NOT re-enable (re-add to a CMake target / re-wire callers) without
// first regenerating the regression baselines under tests/regression/.
//

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
#include <petscdmda.h>
#include <petscviewerhdf5.h>
#include <fenv.h>

#include "default_petsc_options.h"
#include "geometry.h"
#include "mhd.h"

PetscErrorCode mhd_save_hdf5(User *user, const char *filename)
{
    PetscViewer viewer;
    Vec X;

    PetscFunctionBeginUser;

    PetscCall(TSGetSolution(user->ts, &X));

    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD,
                                  filename,
                                  FILE_MODE_WRITE,
                                  &viewer));

    PetscCall(stag_vec_io(user, viewer, X, PETSC_FALSE));

    PetscCall(PetscViewerDestroy(&viewer));

    PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode mhd_load_hdf5(User *user, const char *filename)
{
    PetscViewer viewer;
    Vec X;

    PetscFunctionBeginUser;

    PetscCall(TSGetSolution(user->ts, &X));

    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD,
                                  filename,
                                  FILE_MODE_READ,
                                  &viewer));

    PetscCall(stag_vec_io(user, viewer, X, PETSC_TRUE));

    PetscCall(PetscViewerDestroy(&viewer));

    PetscFunctionReturn(PETSC_SUCCESS);
}
