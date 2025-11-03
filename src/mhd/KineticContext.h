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

#ifndef KINETICCONTEXT_H_
#define KINETICCONTEXT_H_
#include <petscsys.h>
typedef struct
{
    int mpiRank, mpiSize;
    double R0, Z0, dR, dZ;
    int NR, NZ;
    double dt_LA;
    double dt_CD;
    double RECurrentSeedFraction;
    double g_BC_fraq;

    int nParticleSeed;
    int nParticleMax;
    double timeSeed;

    int DoParticleOrbitConservationTest;
    int DoAdvance;
    int DoPoincare;
    int DumpFields;

    char input_folder [PETSC_MAX_PATH_LEN];
} KineticContext;
#endif
