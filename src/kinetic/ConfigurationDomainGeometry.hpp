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
#include "hFlux/StructuredLocator.hpp"

// Conventionally
// region = -2 outside wall
// region = -1 inner wall
// region =  0 metal wall
// region =  1 inside separatrix
// region =  2 scrape-off layer
using Region = int;

struct ConfigurationDomainGeometry {

  using IndicatorViewType = Kokkos::View<int**, Kokkos::DefaultExecutionSpace>;

  StructuredLocator indicator_locator;
  StructuredLocator hermite_locator;
  const IndicatorViewType indicator_view;

  ConfigurationDomainGeometry(
      StructuredLocator indicator_locator,
      IndicatorViewType indicator_view):
    indicator_locator(indicator_locator),
    hermite_locator(makeHermiteLocator<7>(indicator_locator)),
    indicator_view(indicator_view) {};

  KOKKOS_INLINE_FUNCTION
  void locate_region(Real R, Real Z, int& iR, int &iZ, Region& region) const {

    indicator_locator.locateCell(R, Z, iR, iZ);
    region = indicator_view(iR, iZ);
  }

  KOKKOS_INLINE_FUNCTION
  void locate(Real R, Real Z, int& iR, int &iZ, Real& xiR, Real& xiZ) const {

    hermite_locator.locate(R, Z, iR, iZ, xiR, xiZ);
  }
};

struct FreeGeometry {
  KOKKOS_INLINE_FUNCTION
  void locate_region(Real R, Real Z, int&iR, int&jZ, Region& region) const {
    region = 1;
  };
};


