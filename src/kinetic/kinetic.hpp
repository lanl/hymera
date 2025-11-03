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

#ifndef _KINETIC_KINETIC_HPP_
#define _KINETIC_KINETIC_HPP_

#include <memory>
#include "Kokkos_Random.hpp"
#include <parthenon/package.hpp>
// #include <interface/swarm_default_names.hpp>

namespace Kinetic {
using namespace parthenon::package::prelude;

typedef Kokkos::Random_XorShift64_Pool<> RNGPool;

#define VARIABLE(ns, varname)                                                            \
  struct varname : public parthenon::variable_names::base_t<false> {                     \
    template <class... Ts>                                                               \
    KOKKOS_INLINE_FUNCTION varname(Ts &&...args)                                         \
        : parthenon::variable_names::base_t<false>(std::forward<Ts>(args)...) {}         \
    static std::string name() { return #ns "." #varname; }                               \
  }

SWARM_VARIABLE(Real, particle, p); // momentum
SWARM_VARIABLE(Real, particle, xi);// pitch
SWARM_VARIABLE(Real, particle, R);
SWARM_VARIABLE(Real, particle, phi);
SWARM_VARIABLE(Real, particle, Z);
SWARM_VARIABLE(Real, particle, weight);
SWARM_VARIABLE(int, particle, will_scatter);
SWARM_VARIABLE(int, particle, secondary_index);

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin);
void ComputeParticleCounts(Mesh *pm);
void reinit(Mesh *pm);

TaskStatus March(Mesh *pm);

} // namespace Kinetic

#endif // _KINETIC_KINETIC_HPP_
