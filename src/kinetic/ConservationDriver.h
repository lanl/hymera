
//========================================================================================
// (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC // for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================
#pragma once
#include <parthenon/driver.hpp>
#include "mhd.h"

namespace Kinetic {

using namespace parthenon;
using namespace parthenon::driver::prelude;

// Push-only driver used to verify conservation of the guiding-center adiabatic
// invariants (canonical toroidal momentum p_phi, magnetic moment mu). It runs
// the same push sub-cycle as AvalancheDriver but with no avalanche/collision
// tasks, computing the invariants after every push.
class ConservationDriver : public EvolutionDriver {
public:
  ConservationDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh);
  TaskCollection MakeTaskCollection(BlockList_t &blocks, SimTime tm);
  TaskListStatus Step();
};
}
