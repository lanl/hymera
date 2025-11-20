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

constexpr bool PartialScreening = true;
constexpr bool EnergyScattering = true;
constexpr bool ModifiedCouLog = true;

namespace Kinetic {


using namespace parthenon;
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
// For collision book keeping
SWARM_VARIABLE(int, particle, will_scatter);
SWARM_VARIABLE(int, particle, secondary_index);

// For save/restore particle state (predictor corrector implementation)
typedef enum STATUS_ENUM {
    PROTECTED = 1,
    ALIVE = 2
} STATUS;
SWARM_VARIABLE(int, particle, status);
SWARM_VARIABLE(Real, particle, saved_p);
SWARM_VARIABLE(Real, particle, saved_xi);
SWARM_VARIABLE(Real, particle, saved_R);
SWARM_VARIABLE(Real, particle, saved_phi);
SWARM_VARIABLE(Real, particle, saved_Z);
SWARM_VARIABLE(Real, particle, saved_w);

// constexpr auto mkParticleDescriptror_r(const std::string swarm_name) = parthenon::MakeSwarmPackDescriptor<
//       swarm_position::x, swarm_position::y, swarm_position::z,
//       p, xi, R, phi, Z, weight,
//       saved_p, saved_xi, saved_R, saved_phi, saved_Z, saved_w>(swarm_name);
// constexpr auto mkParticleDescriptror_i(const std::string swarm_name) = parthenon::MakeSwarmPackDescriptor<
//       will_scatter, secondary_index, status>(std::stringswarm_name);

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin);
void ComputeParticleWeights(Mesh* pm);
void SaveState(Mesh* pm);
void RestoreState(Mesh* pm);
void InitializeDriver(ParthenonManager* man);
void Push(ParthenonManager* man);

} // namespace Kinetic

#endif // _KINETIC_KINETIC_HPP_
