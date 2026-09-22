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
#include <memory>
#include <parthenon/driver.hpp>
#include <hFlux/common.hpp>
#include "mhd.h"

namespace Kinetic {

using namespace parthenon;
using namespace parthenon::driver::prelude;

// State of the single tracked particle. Lives on the host in a shared_ptr owned
// by main() so both the driver (which advances it each step) and the
// UserMeshWorkBeforeOutput hook (which writes it out) see the same object.
// X = (p, xi, R, phi, Z); t is the accumulated normalized time (t / tau_c).
struct SingleParticleState {
  Dim5 X = {};
  Real t = 0.0;
  Real p_phi = 0.0;
  Real mu = 0.0;
};

// Single-particle push driver. Per parthenon timestep it performs exactly ONE
// collisionless guiding-center push over the full step followed by ONE
// small-angle collision kick, then recomputes the conserved quantities
// (canonical toroidal momentum p_phi, magnetic moment mu). No swarm, no MHD, no
// current deposit -- the particle is a plain host variable.
class SingleParticleDriver : public EvolutionDriver {
public:
  SingleParticleDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh,
                       std::shared_ptr<SingleParticleState> state);
  TaskCollection MakeTaskCollection(BlockList_t &blocks, SimTime tm) { return TaskCollection(); }
  TaskListStatus Step();

private:
  std::shared_ptr<SingleParticleState> state_;
};

// Fill state->p_phi / state->mu for the current state->X without advancing the
// particle. Used to populate the initial (cycle 0) output row after the field
// grid is loaded.
void ComputeConservedForState(Mesh *pm, std::shared_ptr<SingleParticleState> state);
}
