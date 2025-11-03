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

#ifndef RUNAWAYKINETICSOLVERWRAPPER_H_
#define RUNAWAYKINETICSOLVERWRAPPER_H_
#include "KineticContext.h"

void RuKS_init(KineticContext *ctx, void** ruks);
void RuKS_setRECurrentPtr(void* ruks, double* jre);
void RuKS_getFieldDataPtr(void* ruks, double** ptr);

void RuKS_advance(void* ruks, double time);
void RuKS_poincare(void* ruks);
void RuKS_reset(void* ruks);
void RuKS_testConservation(void* ruks, double time);

double RuKS_getIre(void* ruks);

void RuKS_H5Write(void* ruks, const char* filename);
void RuKS_H5Read(void* ruks, const char* filename);

void RuKS_chdir(const char* path);
#endif
