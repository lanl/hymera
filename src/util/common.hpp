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

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>

#include <hFlux/common.hpp>

#include "mfd_config.h"
using RNGPool=Kokkos::Random_XorShift64_Pool<>;

constexpr int _p = 8; // Precision to output doubles

KOKKOS_INLINE_FUNCTION Real ChsPsi(Real x) {
    return 0.5*(erf(x) - (2.0 / sqrt(M_PI)) * x * exp(-x * x)) / x / x;
}

KOKKOS_INLINE_FUNCTION Real gamma_(const Real p) {
  return sqrt(p * p + 1.);
}

KOKKOS_INLINE_FUNCTION Real momentum_(const Real g) {
  return sqrt(g * g - 1.);
}

template<typename T>
KOKKOS_INLINE_FUNCTION
typename T::value_type norm_(const T& v) {
    typename T::value_type ret = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        ret += v[i] * v[i];
    }
    return sqrt(ret);
}

template<typename T>
view3d_t wrap_view(T v) {
    return {
        v.data(),
        v.extent(0),
        v.extent(1),
        v.extent(2),
        v.stride(0),
        v.stride(1),
        v.stride(2)
    };
}

using ExecSpace = Kokkos::DefaultExecutionSpace;
using View3 = Kokkos::View<double**[3], ExecSpace>;
using DualView3 = Kokkos::DualView<double**[3]>;

