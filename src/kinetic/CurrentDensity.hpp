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
#pragma once
#include "util/common.hpp"

template <class Field>
KOKKOS_INLINE_FUNCTION
Real getParticleCurrent(Dim5& X, Real t, Real w, const Field f) {
  Dim3 B, curlB, dBdR, dBdZ, E;
  ERROR_CODE ret = f(X, t, B, curlB, dBdR, dBdZ, E);
  KOKKOS_ASSERT(ret == SUCCESS);

  const Dim5::value_type p = X[0];
  const Dim5::value_type xi = X[1];
  const Dim5::value_type R = X[2];
  const Dim3::value_type b_phi = B[1] / norm_(B) ;

  return -p * xi / gamma_(p) / R / 2.0 / M_PI * w * b_phi;
};

