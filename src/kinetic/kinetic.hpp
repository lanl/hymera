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
#include <hFlux/FieldData.hpp>

constexpr bool PartialScreening = true;
constexpr bool EnergyScattering = true;
constexpr bool ModifiedCouLog = true;
constexpr int FIELD_SMOOTHNESS = 2;  // Hermite m
constexpr int FIELD_FD_STENSIL = 7;  // How many points to use for derivative approximation

#include "mhd/mfd_config.h"
#include "kinetic/FieldComponents.h"

namespace Kinetic {


using namespace parthenon;
using namespace parthenon::package::prelude;

using FieldData_t = FieldData<FIELD_SMOOTHNESS, FIELD_FD_STENSIL, Kokkos::DefaultExecutionSpace, FieldComponents::Total>;

typedef Kokkos::Random_XorShift64_Pool<> RNGPool;

#define VARIABLE(ns, varname)                                                            \
  struct varname : public parthenon::variable_names::base_t<false> {                     \
    template <class... Ts>                                                               \
    KOKKOS_INLINE_FUNCTION varname(Ts &&...args)                                         \
        : parthenon::variable_names::base_t<false>(std::forward<Ts>(args)...) {}         \
    static std::string name() { return #ns "." #varname; }                               \
  }

PAR_SWARMVAR(Real, particle, p); // momentum
PAR_SWARMVAR(Real, particle, xi);// pitch
PAR_SWARMVAR(Real, particle, R);
PAR_SWARMVAR(Real, particle, phi);
PAR_SWARMVAR(Real, particle, Z);
PAR_SWARMVAR(Real, particle, weight);
PAR_SWARMVAR(Real, particle, p_phi);
PAR_SWARMVAR(Real, particle, mu);

// For collision book keeping
PAR_SWARMVAR(int, particle, will_scatter);
PAR_SWARMVAR(int, particle, secondary_index);

// For save/restore particle state (predictor corrector implementation)
typedef enum STATUS_ENUM {
    PROTECTED = 1,
    ALIVE = 2,
    DEATH_BY_MOMENTUM=4,
    DEATH_BY_WALL=8
} STATUS;
PAR_SWARMVAR(int, particle, status);
PAR_SWARMVAR(Real, particle, saved_p);
PAR_SWARMVAR(Real, particle, saved_xi);
PAR_SWARMVAR(Real, particle, saved_R);
PAR_SWARMVAR(Real, particle, saved_phi);
PAR_SWARMVAR(Real, particle, saved_Z);
PAR_SWARMVAR(Real, particle, saved_w);

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin, User* mhd_context);
std::shared_ptr<StateDescriptor> InitializeAnalytic(ParameterInput *pin);
void ComputeParticleWeights(Mesh* pm);

TaskStatus MakeOutputs(Outputs* pouts, Mesh* pmesh, ParameterInput* pinput, Real time, int iPR);

void SaveRawFieldData(Mesh * pm, const char* filename);
void LoadRawFieldData(Mesh * pm, const char* filename);

void WorkBeforeOutput(Mesh * pm, ParameterInput * pin, SimTime const & tm, User* mhd_context);
void PlotFieldsTime(Mesh * pm, ParameterInput * pin, SimTime const & tm, User* mhd_context);
void PlotGCETime(Mesh * pm, ParameterInput * pin, SimTime const & tm, User* mhd_context);
void PlotCurrents(Mesh * pm, ParameterInput * pin, SimTime const & tm, User* mhd_context);
void WorkBeforeRestartOutput(Mesh * pm, ParameterInput * pin, OutputParameters * op, User* mhd_context);
void WorkBeforeLoop(Mesh * pm, User* mhd_context);
TaskStatus HijackEField(Mesh * pm);



} // namespace Kinetic

#endif // _KINETIC_KINETIC_HPP_
