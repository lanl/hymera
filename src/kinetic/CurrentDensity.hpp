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
#include "FieldEvaluator.hpp"

// Get current carried by a particle, in dirrection of b_phi
template <class Field>
KOKKOS_INLINE_FUNCTION
Real getParticleCurrent(Dim5& X, Real t, Real w, const Field f) {
  const Real p = X[0];
  const Real xi = X[1];
  const Real R = X[2];
  const Real Z = X[4];

  Kinetic::EvalB B;
  f.eval(B, R, Z, t);

  const Real b_phi = B.B[1] / B.Bmag;
  return -p * xi / gamma_(p) / R / 2.0 / M_PI * w * b_phi;
};

// Second order shaping funciton
KOKKOS_INLINE_FUNCTION
Real S2(Real x) {
    if (x < 0.5) return 0.75 - x * x;
    else return (3.0 - 2.0 * x) * (3.0 - 2.0 * x) / 8.0;
}

// Add particle currect contibution wheighted by time interval to a current dencity for averaging
template <class CurrentDensityView, class Field, class CurrentDensityLocator>
KOKKOS_INLINE_FUNCTION
void DepositCurrent(const Dim5& X, const Real t, const Real w, CurrentDensityView jre, const Real time_interval, const Field field, const CurrentDensityLocator locator) {

  const Dim5::value_type p = X[0];
  const Dim5::value_type xi = X[1];
  const Dim5::value_type R = X[2];
  const Dim5::value_type Z = X[4];

  Real contribution = -p * xi / gamma_(p) / R /
    locator.dR /
    locator.dZ /
    2.0 / M_PI * time_interval * w;

  int i, j;
  Real xiR, xiZ;
  locator.locate(R, Z, i, j, xiR, xiZ);

  Kinetic::EvalB B;
  field.eval(B, R, Z, t);

  for (int ii = -1; ii < 2; ++ii) {
      if(i + ii >= 0 and i + ii < jre.extent(0)) {
          Real wr = S2(Kokkos::abs(xiR - static_cast<Real>(ii)));
          for (int jj = -1; jj < 2; ++jj) {
              if(j + jj >= 0 and j + jj < jre.extent(1)) {
                  Real wz = S2(Kokkos::abs(xiZ - static_cast<Real>(jj)));
                  Real weighted_contribution = contribution * wr * wz;
                  for (int kk = 0; kk < 3; ++kk) {
                      Real wcB = weighted_contribution * B.B[kk] / B.Bmag;
                      Kokkos::atomic_add(&(jre(i,j,kk)), wcB);
                  }
              }
          }
      }
  }
}
