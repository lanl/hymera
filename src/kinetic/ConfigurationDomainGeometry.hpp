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

struct ConfigurationDomainGeometry {

  using IndicatorViewType = Kokkos::View<int**, Kokkos::DefaultExecutionSpace>;

  const Real R0;
  const Real Z0;

  const Real dR;
  const Real dZ;

  const typename IndicatorViewType::value_type outside;
  const IndicatorViewType indicator_view;

  ConfigurationDomainGeometry(Real R0, Real Z0, Real dR, Real dZ, typename IndicatorViewType::value_type outside, IndicatorViewType indicator_view):
    R0(R0), Z0(Z0), dR(dR), dZ(dZ), outside(outside), indicator_view(indicator_view) {};

  KOKKOS_INLINE_FUNCTION
  typename IndicatorViewType::value_type indicator(const Dim5& X, int&i, int&j) const {
    i = static_cast<int> (floor((X[2] - R0) / dR));
    j = static_cast<int> (floor((X[4] - Z0) / dZ));

    if ((i < 0) or (i >= indicator_view.extent(0)) or
        (j < 0) or (j >= indicator_view.extent(1)))
      return outside;
    return indicator_view(i, j);

  };

  KOKKOS_INLINE_FUNCTION
  void getLocalCoordinate(const Dim5& X, int i, int j, Dim2 Xloc) const {
    Xloc[0] = X[2] / dR - static_cast<Real>(i) - 0.5;
    Xloc[1] = X[4] / dZ - static_cast<Real>(j) - 0.5;
  };
};


