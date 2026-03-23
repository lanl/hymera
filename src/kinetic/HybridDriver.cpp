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

#include <hFlux/dopri.hpp>

using namespace parthenon;
using namespace parthenon::driver::prelude;

#include "HybridDriver.h"

#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/EM_Field.hpp"
#include "pgen.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

namespace Kinetic {

TaskStatus PushParticles(Mesh *pm, Real t0, Real dt) {
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

  const Real dtSA_min = sa.getSmallAngleCollisionTimestep(momentum_(1.002));
  const Real dtSA_max = dt;

  const auto c_aw0  = pkg->Param<Real>("c_aw0");
  const auto ct_a   = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");

  const auto f = pkg->Param<EM_Field>("Field");
  GuidingCenterEquations<EM_Field, true, false> gce(f, c_aw0, ct_a, alpha0);

  Kokkos::Timer timer;

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index,
                                         Kinetic::status>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  auto jre = pkg->Param<ParArray3D<Real>>("Jre_push_deposit");

  Kokkos::fence();

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
          Real t = tstart;
          X[0] = pack_swarm_r(b, Kinetic::p(), n);
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

            auto ret = solve_dopri5(gce, X, t, t + dtSA, rtol, atol, h, 1e-9,
                         std::numeric_limits<int>::max(), work_d);
            if (ret != ErrorCode::Success) {
              pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
              if ((pack_swarm_i(b, Kinetic::status(), n) & PROTECTED) == 0)
                swarm_d.MarkParticleForRemoval(n);
              break;
            }
            int ii,jj;
            int level = f.cdg.indicator(X, ii,jj);

            if (level < 1 || X[0] < p_BC) {
              pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
              if ((pack_swarm_i(b, Kinetic::status(), n) & PROTECTED) == 0)
                swarm_d.MarkParticleForRemoval(n);
              if (level < 1)
                pack_swarm_i(b, Kinetic::status(), n) |= Kinetic::DEATH_BY_WALL;
              else
                pack_swarm_i(b, Kinetic::status(), n) |= Kinetic::DEATH_BY_MOMENTUM;
              break;
            }

            if (X[0] > p_RE) {
              DepositCurrent(X, t, w, jre, dtSA, f);
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
  Kokkos::fence();

  return TaskStatus::complete;

}

TaskStatus CheckScatter(MeshBlock* pmb) {
  return TaskStatus::complete;

  auto data = pmb->meshblock_data.Get();
  auto swarm = data->GetSwarmData()->Get("particles");
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index,
                                         Kinetic::status>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(data.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(data.get());

  auto swarm_d = swarm->GetDeviceContext();

  Kokkos::parallel_scan(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      KOKKOS_LAMBDA(const int n, int &running_total, const bool final_pass) {
        const int b = 0;
        if (swarm_d.IsActive(n)&& !swarm_d.IsMarkedForRemoval(n) && (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
          if (pack_swarm_i(b, Kinetic::will_scatter(), n) == 1) {
            running_total += 1;
            if (final_pass) {
              pack_swarm_i(b, Kinetic::secondary_index(), n) = running_total;
            }
          } else {
            if (final_pass) {
              pack_swarm_i(b, Kinetic::secondary_index(), n) = 0;
            }
          }
        }
      });
  Kokkos::fence();

	return TaskStatus::complete;
}

TaskStatus CleanupParticles(MeshBlock* pmb) {
  pmb->meshblock_data.Get()
  ->GetSwarmData()->Get("particles")
  ->RemoveMarkedParticles();
	return TaskStatus::complete;
}

TaskStatus AddSecondaries(MeshBlock* pmb, const Real dtLA) {
  return TaskStatus::complete;

  auto pkg = pmb->packages.Get("Deck");
  auto gamma_min = *(pkg->Param<std::shared_ptr<Real>>("gamma_min"));
  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  auto data = pmb->meshblock_data.Get();
  auto swarm = data->GetSwarmData()->Get("particles");
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index,
                                         Kinetic::status>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(data.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(data.get());

  auto swarm_d = swarm->GetDeviceContext();
  int ntot = 0, nalive = 0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      KOKKOS_LAMBDA(const int n, int &nnew, int &nnalive) {
        const int b = 0;
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)&& (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
          nnalive += 1;
          if (pack_swarm_i(b, Kinetic::will_scatter(), n) == 1)
            nnew += 1;
        }
      },
      ntot, nalive);
  Kokkos::fence();
  if (ntot > 0) {
    std::cout << std::format("Adding {} new particles, total alive {}, ratio {}", ntot, nalive, (Real) (ntot + nalive) / (Real) nalive) << std::endl;

    const int oldMaxIndex = pack_swarm_r.GetMaxFlatIndex();
    auto newParticlesContext = swarm->AddEmptyParticles(ntot);
    auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
        swarm_position::x, swarm_position::y, swarm_position::z,
        Kinetic::p, Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z,
        Kinetic::weight>("particles");
    auto desc_swarm_i = parthenon::MakeSwarmPackDescriptor<
        Kinetic::will_scatter, Kinetic::secondary_index, Kinetic::status>("particles");
    pack_swarm_r = desc_swarm_r.GetPack(data.get());
    pack_swarm_i = desc_swarm_i.GetPack(data.get());

    swarm_d = swarm->GetDeviceContext();

    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL, DevExecSpace(), 0,
        newParticlesContext.GetNewParticlesMaxIndex(),
        // new_n ranges from 0 to N_new_particles
        KOKKOS_LAMBDA(const int new_n) {
          // this is the particle index inside the swarm
          const int n = newParticlesContext.GetNewParticleIndex(new_n);
          const int b = 0;
          pack_swarm_i(b, Kinetic::will_scatter(), n) = 0;
          pack_swarm_i(b, Kinetic::secondary_index(), n) = 0;
        });

    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL, DevExecSpace(), 0,
        oldMaxIndex,
        // new_n ranges from 0 to N_new_particles
        KOKKOS_LAMBDA(const int n_primary) {
          const int b = 0;
          // this is the particle index inside the swarm
          if (swarm_d.IsActive(n_primary)&& !swarm_d.IsMarkedForRemoval(n_primary) &&(pack_swarm_i(b, Kinetic::status(), n_primary) & Kinetic::ALIVE))
            if (pack_swarm_i(b, Kinetic::will_scatter(), n_primary) == 1) {
              int new_n =
                  pack_swarm_i(b, Kinetic::secondary_index(), n_primary) - 1;
              const int n = newParticlesContext.GetNewParticleIndex(new_n);
              pack_swarm_r(b, swarm_position::x(), n) =
                  pack_swarm_r(b, swarm_position::x(), n_primary);
              pack_swarm_r(b, swarm_position::y(), n) =
                  pack_swarm_r(b, swarm_position::y(), n_primary);
              pack_swarm_r(b, swarm_position::z(), n) =
                  pack_swarm_r(b, swarm_position::z(), n_primary);
              pack_swarm_r(b, Kinetic::R(), n) =
                  pack_swarm_r(b, Kinetic::R(), n_primary);
              pack_swarm_r(b, Kinetic::phi(), n) =
                  pack_swarm_r(b, Kinetic::phi(), n_primary);
              pack_swarm_r(b, Kinetic::Z(), n) =
                  pack_swarm_r(b, Kinetic::Z(), n_primary);
              Real p = pack_swarm_r(b, Kinetic::p(), n_primary);
              Real xi = pack_swarm_r(b, Kinetic::xi(), n_primary);
              Real w = pack_swarm_r(b, Kinetic::weight(), n_primary);

              LargeAngleCollision(p, xi, w, dtLA, gamma_min, rng_pool);
              if (p > 0.0) {
                pack_swarm_i(b, Kinetic::status(), n) = Kinetic::ALIVE;
                pack_swarm_r(b, Kinetic::p(), n) = p;
                pack_swarm_r(b, Kinetic::xi(), n) = xi;
                pack_swarm_r(b, Kinetic::weight(), n) = w;
              } else {
                pack_swarm_i(b, Kinetic::status(), n) = 0;
                swarm_d.MarkParticleForRemoval(n);
              }
            }
        });

    Kokkos::fence();
  }

	return TaskStatus::complete;
}


TaskStatus CollectCurrent(Mesh *pm, const int iCD, const Real dtCD) {
  auto md = pm->mesh_data.Get();
  auto pkg = pm->packages.Get("Deck");

  auto NR = pkg->Param<int>("NR");
  auto NZ = pkg->Param<int>("NZ");
  auto eta_a3VaB0 = pkg->Param<Real>("eta_a3VaB0");

  auto jre_d = pkg->Param<ParArray3D<Real>>("Jre_push_deposit");
  auto jre_deposit = pkg->Param<ParArrayHost<Real>>("Jre_deposit").KokkosView();
  Kokkos::parallel_for("FillInterpolatedData_plot",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {NR,NZ,3}),
      KOKKOS_LAMBDA(const int i, const int j, const int k) {
        jre_d(i,j,k) *= eta_a3VaB0 / dtCD;
      });
  Kokkos::fence();
  auto jre_h = create_mirror_view_and_copy(Kokkos::HostSpace(),jre_d);
  Kokkos::fence();
  MPI_Allreduce(MPI_IN_PLACE,jre_h.data(),jre_h.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  auto sub = Kokkos::subview(jre_deposit, 0, 0, 0, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iCD);
  Kokkos::deep_copy(sub, jre_h);
  Kokkos::deep_copy(jre_d, jre_h);
  Kokkos::fence();

  auto f = pkg->Param<EM_Field>("Field");
  auto jre_data = Kokkos::subview(f.data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, static_cast<size_t>(fid::Jre), 0);
  Kokkos::deep_copy(jre_data, jre_d);
  Kokkos::fence();
  Kokkos::Array<fid,1> fids = {fid::Jre};
  f.interpolate(fids, 0);
  Kokkos::deep_copy(jre_d, 0.0);
  Kokkos::fence();

  return TaskStatus::complete;
}

TaskStatus MHDStep(User* p_mhd_config) {
  mhd_step(p_mhd_config);
  return TaskStatus::complete;
}

TaskStatus ResetState() {
  return TaskStatus::complete;
}


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

  auto * tl = &tc.AddRegion(1)[0];
  auto dep = none;

  Real dtCD = tm.dt / tau_c / nCDperMHDstep;

  dep = tl->AddTask(dep, Interpolate, pmesh, p_mhd_config);

  for (int iPR = 0; iPR < nPredictorSteps + 1; ++iPR) {
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
            auto cleanup = tl.AddTask(add_secondaries, CleanupParticles, pmb.get());
          }
          tl = &tc.AddRegion(1)[0];
          dep = none;
        }
      }
      dep = tl->AddTask(dep, CollectCurrent, pmesh, iCD, dtCD);
    }
    dep = tl->AddTask(dep, MHDStep, p_mhd_config);
    if (iPR < nPredictorSteps) {
      dep = tl->AddTask(dep, ResetState); // Puts particles back to the start, resets MHD state back to the start
    }
  }

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
