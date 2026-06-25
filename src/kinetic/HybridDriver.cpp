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

#include "HybridDriver.h"

#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/EM_Field.hpp"
#include "pgen.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

namespace Kinetic {

HybridDriver::HybridDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh,
    User* p_mhd_config)
      : EvolutionDriver(pin, app_in, pmesh), p_mhd_config(p_mhd_config) {
  	auto pkg = pmesh->packages.Get("Deck");
}

TaskCollection HybridDriver::MakeTaskCollection(BlockList_t &blocks, SimTime tm) {
  TaskCollection tc;
  TaskID none(0);

  auto pkg = pmesh->packages.Get("Deck");

  const int EnableLargeAngleCollisions = pkg->Param<int>("EnableLargeAngleCollisions");
  const int nPredictorSteps = pkg->Param<int>("nPR");
  const int nCDperMHDstep   = pkg->Param<int>("nCD");
  const int nLAperCD        = pkg->Param<int>("nLA");
  const Real dtLA_over_tauC  = pkg->Param<Real>("dtLA_over_tauC");
  const Real tau_c = pkg->Param<Real>("tau_c");

  Real dt = tm.dt / tau_c;
  Real dtCD = dt / nCDperMHDstep;

  auto dep = none;
  tl = &tc.AddRegion(1)[0];

  AttachFieldPredictor(tl, dep, pmesh, dt);

  for (int iCD = 0; iCD < nCDperMHDstep; ++iCD) {
    for (int iLA = 0; iLA < nLAperCD; ++iLA) {
      Real t0 = iLA * dtLA_over_tauC + iCD * dtCD;
      dep = tl->AddTask(dep, PushParticles, pmesh, t0, dtLA_over_tauC);

      if (EnableLargeAngleCollisions == 1) {
        // these are per block tasklists
        TaskRegion &async_region = tc.AddRegion(blocks.size());
        for (int i = 0; i < blocks.size(); ++i) {
          // required by this MeshData object)
	        auto &pmb = blocks[i];
          auto &tl = async_region[i];
          auto check_scatter = tl.AddTask(none, CheckScatter, pmb.get());
          auto add_secondaries = tl.AddTask(check_scatter, AddSecondaries, pmb.get(), dtLA_over_tauC);
          auto cleanup = tl.AddTask(add_secndaries, AddSecondaries, pmb.get(), dtLA_over_tauC);
        }

        tl = &tc.AddRegion(1)[0];
        dep = none;
      }
    }
    dep = tl->AddTask(dep, RandomRemove, pmesh);
    AttachCurrentCollector(dep, pmesh, dtCD);
    dep = AddTask(dep, UpdateMomentumBoundary, pmesh);

    Real time =  tm.time / tau_c + (iCD+1) * dtCD;
  }

  AttachFieldCorrector(tl, dep, p_mhd_context);
  tl.AddTask(tl, dep, DefragSwarmsMesh, pmesh);

  return tc;
}

void HybridDriver::PreExecute() {
}

void HybridDriver::PostExecute(parthenon::DriverStatus st) {
}

TaskListStatus HybridDriver::Step() {
  PARTHENON_INSTRUMENT
  using DriverUtils::ConstructAndExecuteTaskLists;
  TaskListStatus status = ConstructAndExecuteTaskLists<>(this, tm);
  return status;
}

}
