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

#include "kinetic/AvalancheDriver.h"

#include "kinetic/kinetic.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/AnalyticField.hpp"
#include "kinetic/ParticleVerificator.hpp"
#include "tasks/Tasks.h"
#include "pgen.hpp"
#include "rk45.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

namespace Kinetic {

struct AnalyticParticleVerificator {
  using ResultCode_t = int;
  static constexpr ResultCode_t Success = ParticleVerifyCodes::Success;

  Real p_BC;

  KOKKOS_INLINE_FUNCTION
  ResultCode_t verify(const Dim5 &y) const {
    if (!Kokkos::isfinite(y[0])) {
      return ParticleVerifyCodes::InvalidState;
    }
    if (y[0] < p_BC) {
      return ParticleVerifyCodes::MomentumCutoff;
    }
    return Success;
  }
};
/*
TaskStatus PushParticlesAnalytic(Mesh *pm, Real t0, Real dt) {
  // get mesh data
  auto md = pm->mesh_data.Get();
  auto pkg = pm->packages.Get("Deck");

  const auto h = pkg->Param<Real>("hRK");
  const auto atol = pkg->Param<Real>("atol");
  const auto rtol = pkg->Param<Real>("rtol");
  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");

  const auto gamma_min = pkg->Param<Real>("gamma_min");
  const auto p_BC =      pkg->Param<Real>("p_BC");
  const auto p_RE =      pkg->Param<Real>("p_RE");

  const auto ms = pkg->Param<MollerSource>("MollerSource");
  const auto sa = pkg->Param<
    SmallAngleCollision<PartialScreening, EnergyScattering, ModifiedCouLog>
  >("SmallAngleCollision");

  const Real dtSA_min = sa.getSmallAngleCollisionTimestep(momentum_(1.00002));
  const Real dtSA_max = dt;

  const auto c_aw0  = pkg->Param<Real>("c_aw0");
  const auto ct_a   = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");


  const auto f = pkg->Param<AnalyticField>("Field");
  GuidingCenterEquations<AnalyticField, true, false> gce(f, c_aw0, ct_a, alpha0);
  AnalyticParticleVerificator ver{p_BC};

  Kokkos::Timer timer;

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles"); auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index,
                                         Kinetic::status>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());


  const Real tstart = t0;
  const Real tstop =  t0 + dt;


  int EnableLargeAngleCollisions = pkg->Param<int>("EnableLargeAngleCollisions");
  int EnableSmallAngleCollisions = pkg->Param<int>("EnableSmallAngleCollisions");

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_r.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        // block and particle indices
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm_r.GetContext(b);
        const auto markers_d = pack_swarm_i.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)&&
            (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE) ) {
          Dim5 X;
          X[0] = pack_swarm_r(b, Kinetic::p(), n);

          // Skip particles below momentum threshold
          if (X[0] < p_BC) {
            pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
            if ((pack_swarm_i(b, Kinetic::status(), n) & PROTECTED) == 0)
              swarm_d.MarkParticleForRemoval(n);
            return;
          }

          Real t = tstart;
          X[1] = pack_swarm_r(b, Kinetic::xi(), n);
          X[2] = pack_swarm_r(b, Kinetic::R(), n);
          X[3] = pack_swarm_r(b, Kinetic::phi(), n);
          X[4] = pack_swarm_r(b, Kinetic::Z(), n);
          Real w = pack_swarm_r(b, Kinetic::weight(), n);
          Kokkos::Array<Dim5, 10> work_d;

          bool last_step = false;

          while (last_step == false) {
            Real dtSA =
                sa.getSmallAngleCollisionTimestep(X[0], dtSA_min, dtSA_max);

            if (t + dtSA > tstop) {
              dtSA = tstop - t;
              if (dtSA < 1e-16) {
                break;
              }
              last_step = true;
            }

            auto ret = solve_rk45(gce, ver, X, t, t + dtSA, rtol, atol, h, 1e-9,
                         std::numeric_limits<int>::max(), work_d);
            if (ret != AnalyticParticleVerificator::Success) {
              pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
              if ((pack_swarm_i(b, Kinetic::status(), n) & PROTECTED) == 0)
                swarm_d.MarkParticleForRemoval(n);
              break;
            }
            if (EnableSmallAngleCollisions == 1)
              sa(X[0], X[1], dtSA, rng_pool);
            t += dtSA;
            if (t > tstop)
              break;
          }

        	pack_swarm_r(b, Kinetic::p(), n)   = X[0];
        	pack_swarm_r(b, Kinetic::xi(), n)  = X[1];
        	pack_swarm_r(b, Kinetic::R(), n)   = X[2];
        	pack_swarm_r(b, Kinetic::phi(), n) = X[3];
        	pack_swarm_r(b, Kinetic::Z(), n)   = X[4];

          if (EnableLargeAngleCollisions == 1)
            pack_swarm_i(b, Kinetic::will_scatter(), n) = ms(X[0], w, dt, gamma_min, rng_pool);
        }
      });

  return TaskStatus::complete;

}
*/

AvalancheDriver::AvalancheDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh)
      : EvolutionDriver(pin, app_in, pmesh) {
  	auto pkg = pmesh->packages.Get("Deck");
}

TaskCollection AvalancheDriver::MakeTaskCollection(BlockList_t &blocks, SimTime tm) {

  TaskCollection tc;
  TaskID none(0);

  auto pkg = pmesh->packages.Get("Deck");

  const int EnableLargeAngleCollisions = pkg->Param<int>("EnableLargeAngleCollisions");
  const int nCDperMHDstep   = pkg->Param<int>("nCD");
  const int nLAperCD        = pkg->Param<int>("nLA");
  const Real dtLA  = pkg->Param<Real>("dtLA_over_tauC");
  const Real tau_c = pkg->Param<Real>("tau_c");

  auto * tl = &tc.AddRegion(1)[0];
  auto dep = none;

  Real dt = tm.dt / tau_c;
  Real dtCD = dt / nCDperMHDstep;

  std::cout << "Creating task collection: " << std::format("dt = {:e}\n dtCD = {:e}\n dtLA = {:e}\n tau_c = {:e}\n", dt, dtCD, dtLA, tau_c);
  dep = none;

  for (int iCD = 0; iCD < nCDperMHDstep; ++iCD) {
    for (int iLA = 0; iLA < nLAperCD; ++iLA) {
      Real t0 = iLA * dtLA + iCD * dtCD;
      dep = tl->AddTask(dep, PushParticles, pmesh, t0, dtLA);

      TaskRegion &async_region = tc.AddRegion(blocks.size());
      for (int i = 0; i < blocks.size(); ++i) {
        // required by this MeshData object)
	      auto &pmb = blocks[i];
        auto &tl = async_region[i];
        auto check_scatter = tl.AddTask(none, CheckScatter, pmb.get());
        auto add_secondaries = tl.AddTask(check_scatter, AddSecondaries, pmb.get(), dtLA);
        auto cleanup = tl.AddTask(add_secondaries, CleanupParticles, pmb.get());
      }
      tl = &tc.AddRegion(1)[0];
      dep = none;
    }
  //  dep = tl->AddTask(dep, CollectCurrent, pmesh, iCD, dtCD);
  }

  dep = tl->AddTask(dep, RandomRemove, pmesh); // If there are more alive particles then limit, kill half, doubling the weight
  TaskRegion &async_region = tc.AddRegion(blocks.size());
  for (int i = 0; i < blocks.size(); ++i) {
    // required by this MeshData object)
    auto &pmb = blocks[i];
    auto &tl = async_region[i];
    auto cleanup = tl.AddTask(none, CleanupParticles, pmb.get());
  }
  tl = &tc.AddRegion(1)[0];

  return tc;
}

TaskListStatus AvalancheDriver::Step() {
  PARTHENON_INSTRUMENT
  using DriverUtils::ConstructAndExecuteTaskLists;
  TaskListStatus status = ConstructAndExecuteTaskLists<>(this, tm);
  return status;
}

}
