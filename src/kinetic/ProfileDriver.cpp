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
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <globals.hpp>
#include <parthenon_manager.hpp>
#include <iostream>
#include <iomanip>
#include <limits>
#include <format>

#include <Kokkos_Core.hpp>

using namespace parthenon;
using namespace parthenon::driver::prelude;

#include "kinetic/ProfileDriver.h"

#include "kinetic/kinetic.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/AnalyticField.hpp"
#include "kinetic/ParticleVerificator.hpp"
#include "tasks/Tasks.h"
#include "pgen.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

namespace Kinetic {

ProfileDriver::ProfileDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh)
      : EvolutionDriver(pin, app_in, pmesh) {
  auto pkg = pmesh->packages.Get("Deck");
}

// Push-only task collection: iterate the same nCD x nLA substep grid as the
// avalanche/hybrid drivers, but do nothing except call PushParticles. No
// scatter, secondaries, random-remove, or current-collector tasks -- this
// isolates the guiding-center integrator kernel for profiling.
TaskCollection ProfileDriver::MakeTaskCollection(BlockList_t &blocks, SimTime tm) {
  TaskCollection tc;
  TaskID none(0);

  auto pkg = pmesh->packages.Get("Deck");

  const int nCDperMHDstep = pkg->Param<int>("nCD");
  const int nLAperCD      = pkg->Param<int>("nLA");
  const Real dtLA         = pkg->Param<Real>("dtLA_over_tauC");
  const Real tau_c        = pkg->Param<Real>("tau_c");

  Real dt = tm.dt / tau_c;
  Real dtCD = dt / nCDperMHDstep;

  auto *tl = &tc.AddRegion(1)[0];
  auto dep = none;

  if (Globals::my_rank == 0)
    std::cout << std::format("ProfileDriver task collection: dt = {:e}, dtCD = {:e}, dtLA = {:e}, tau_c = {:e}\n",
        dt, dtCD, dtLA, tau_c);

  for (int iCD = 0; iCD < nCDperMHDstep; ++iCD) {
    for (int iLA = 0; iLA < nLAperCD; ++iLA) {
      Real t0 = iLA * dtLA + iCD * dtCD;
      dep = tl->AddTask(dep, PushParticles, pmesh, t0, dtLA);
    }
  }

  return tc;
}

TaskListStatus ProfileDriver::Step() {
  PARTHENON_INSTRUMENT
  using DriverUtils::ConstructAndExecuteTaskLists;
  TaskListStatus status = ConstructAndExecuteTaskLists<>(this, tm);
  return status;
}

}
