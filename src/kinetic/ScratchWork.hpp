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
// The Runge-Kutta work array of one particle in team scratch memory (Simulation/
// push_scratch = 1), instead of the per-thread Kokkos::Array the push keeps in
// registers or local memory.
//
// One team owns a ScratchPad1D of NSlots * 5 * nthreads doubles. Thread t's
// component i of slot s is at ((s * 5 + i) * nthreads + t): neighbouring threads
// touch neighbouring words, so a warp reads a slot component without bank
// conflicts. work[s] is a row handle with row[i]; the solvers (rk4.hpp, rk45.hpp)
// and GuidingCenterEquations::operator() take either form.
#pragma once

#include <parthenon/package.hpp>

#include "util/common.hpp"

namespace Kinetic {

struct ScratchRow {
  Real *base;
  int stride;
  KOKKOS_FORCEINLINE_FUNCTION Real &operator[](const int i) const { return base[i * stride]; }
  static constexpr KOKKOS_FUNCTION int size() { return 5; }
};

template <int NSlots>
struct ScratchWork {
  Real *base;  // this thread's word of slot 0, component 0
  int stride;  // the team's thread count

  KOKKOS_FORCEINLINE_FUNCTION ScratchRow operator[](const int s) const {
    return ScratchRow{base + s * 5 * stride, stride};
  }

  /// Doubles of team scratch for nthreads particles.
  static constexpr std::size_t Count(const int nthreads) {
    return static_cast<std::size_t>(NSlots) * 5 * nthreads;
  }
};

} // namespace Kinetic
