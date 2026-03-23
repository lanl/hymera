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

// Get current carried by a particle, in dirrection of b_phi
template <class Field>
KOKKOS_INLINE_FUNCTION
Real getParticleCurrent(Dim5& X, Real t, Real w, const Field f) {
  Dim3 B, curlB, dBdR, dBdZ, E, dbdt;
  ERROR_CODE ret = f(X, t, B, curlB, dBdR, dBdZ, E, dbdt);
  KOKKOS_ASSERT(ret == SUCCESS);

  const Dim5::value_type p = X[0];
  const Dim5::value_type xi = X[1];
  const Dim5::value_type R = X[2];
  const Dim3::value_type b_phi = B[1] / norm_(B) ;

  return -p * xi / gamma_(p) / R / 2.0 / M_PI * w * b_phi;
};

// Second order shaping funciton
KOKKOS_INLINE_FUNCTION
Real S2(Real x) {
    if (x < 0.5) return 0.75 - x * x;
    else return (3.0 - 2.0 * x) * (3.0 - 2.0 * x) / 8.0;
}

// Add particle currect contibution wheighted by time interval to a current dencity for averaging
template <class CurrentDensityView, class Field>
KOKKOS_INLINE_FUNCTION
void DepositCurrent(const Dim5& X, const Real t, const Real w, CurrentDensityView jre, const Real time_interval, Field& field) {

  const Dim5::value_type p = X[0];
  const Dim5::value_type xi = X[1];
  const Dim5::value_type R = X[2];

  Real contribution = -p * xi / gamma_(p) / R /
    field.cdg.dR /
    field.cdg.dZ /
    2.0 / M_PI * time_interval * w;

  int i, j;
  int level = field.cdg.indicator(X, i, j);

  Dim2 Xlocd = {};
  field.cdg.getLocalCoordinate(X, i, j, Xlocd);

  if (level < 1) return;

  Dim3 B, curlB, dBdR, dBdZ, E, dbdt;
  ERROR_CODE ret = field(X, t, B, curlB, dBdR, dBdZ, E, dbdt);
  KOKKOS_ASSERT(ret == SUCCESS);


  Real BB = norm_(B);

  for (int ii = -1; ii < 2; ++ii) {
      if(i + ii >= 0 and i + ii < field.data.extent(0)) {
          Real wr = S2(abs(Xlocd[0] - static_cast<Real>(ii)));
          for (int jj = -1; jj < 2; ++jj) {
              if(j + jj >= 0 and j + jj < field.data.extent(1)) {
                  Real wz = S2(abs(Xlocd[1] - static_cast<Real>(jj)));
                  Real weighted_contribution = contribution * wr * wz;
                  for (int kk = 0; kk < 3; ++kk) {
                      Real wcB = weighted_contribution * B[kk] / BB;
                      Kokkos::atomic_add(&(jre(i,j,kk)), wcB);
                  }
              }
          }
      }
  }
}
