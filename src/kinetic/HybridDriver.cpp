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

TaskStatus PushParticles(Mesh *pm, Real dt) {
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

  auto jre = f->getJreDataSubview();

  const Real tstart = tm.time;
  const Real tstop = tm.time + tm.dt;

  auto field_interpolation = *f;

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
            if (ret != SUCCESS) {
              pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
              if ((pack_swarm_i(b, Kinetic::status(), n) & PROTECTED) == 0)
                swarm_d.MarkParticleForRemoval(n);
              break;
            }
            int ii,jj;
            int level = field_interpolation.cdg.indicator(X, ii,jj);

            if (level < 1 || X[0] < p_BC) {
              pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
              if ((pack_swarm_i(b, Kinetic::status(), n) & PROTECTED) == 0)
                swarm_d.MarkParticleForRemoval(n);
              break;
            }

            if (X[0] > p_RE) {
              DepositCurrent(X, t, w, jre, dtSA, field_interpolation);
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
            pack_swarm_i(b, Kinetic::will_scatter(), n) = ms(X[0], w, tm.dt, gamma_min, rng_pool);
        }
      });
  Kokkos::fence();

  return TaskStatus::complete;

}

TaskStatus CheckScatter(MeshBlock* pmb) {
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


TaskStatus CollectCurrent() {
  return TaskStatus::complete;
}

TaskStatus MHDStep() {
  return TaskStatus::complete;
}

TaskStatus ResetState() {
  return TaskStatus::complete;
}

TaskCollection HybridDriver::MakeTaskCollection(BlockList_t &blocks, SimTime tm) {
  TaskCollection tc;
  TaskID none(0);

  auto pkg = pmesh->packages.Get("Deck");

  int EnableLargeAngleCollisions = pkg->Param<int>("EnableLargeAngleCollisions");
  int nPredictorSteps = pkg->Param<int>("nPR");
  int nCDperMHDstep   = pkg->Param<int>("nCD");
  int nLAperCD        = pkg->Param<int>("nLA");

  auto * tl = &tc.AddRegion(1)[0];
  auto dep = none;

  for (int iPR = 0; iPR < nPredictorSteps; ++iPR) {
    for (int iCD = 0; iCD < nCDperMHDstep; ++iCD) {
      for (int iLA = 0; iLA < nLAperCD; ++iLA) {
        auto dep = tl->AddTask(dep, PushParticles, pmesh, tm);

        if (EnableLargeAngleCollisions == 1) {
          // these are per block tasklists
          TaskRegion &async_region = tc.AddRegion(blocks.size());
          for (int i = 0; i < blocks.size(); ++i) {
            // required by this MeshData object)
	          auto &pmb = blocks[i];
            auto &tl = async_region[i];
            auto check_scatter = tl.AddTask(none, CheckScatter, pmb.get());
            auto add_secondaries = tl.AddTask(check_scatter, AddSecondaries, pmb.get(), tm.dt);
            auto cleanup = tl.AddTask(add_secondaries, CleanupParticles, pmb.get());
          }
          tl = &tc.AddRegion(1)[0];
          dep = none;
        }
      }
      dep = tl->AddTask(dep, CollectCurrent, pmesh, tm);
    }
    dep = tl->AddTask(dep, MHDStep, pmesh, tm);
    if (iPR < nPredictorSteps - 1) {
      dep = tl->AddTask(dep, ResetState); // Puts particles back to the start, resets MHD state back to the start
      // Background field interpolation uses new predicted end state
    }
  }

  return tc;
}

void HybridDriver::PreExecute() {
  auto pkg = pmesh->packages.Get("Deck");
  auto f = pkg->Param<std::shared_ptr<EM_Field>>("Field");

  // Interpolate the jre data thats there and fields if new
  f->interpolate(); // sets jre to electric field
  auto f_d = *f;

  std::shared_ptr<Real> gamma_min = pkg->Param<std::shared_ptr<Real>>("gamma_min");
  std::shared_ptr<Real> p_BC = pkg->Param<std::shared_ptr<Real>>("p_BC");
  std::shared_ptr<Real> p_RE = pkg->Param<std::shared_ptr<Real>>("p_RE");

  const auto time = tm.time;

  Real maxE;
  Kokkos::parallel_reduce("max E",
  Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {f_d.nR_data,f_d.nZ_data}),
  KOKKOS_LAMBDA(int i, int j, Real& Epar) {
    Dim5 X;
    X[2] = f_d.cdg.R0 + f_d.dR * i;
    X[4] = f_d.cdg.Z0 + f_d.dZ * j;

    int ii,jj;
    int level = f_d.cdg.indicator(X, ii,jj);
    if (level > 0) {
      Dim3 B = {}, curlB = {}, dBdR = {}, dBdZ = {}, E = {}, Jre = {}, V = {}, dbdt = {};

      auto ret = f_d(X, time, B, curlB, dBdR, dBdZ, E, Jre, V, dbdt);
      Epar = Kokkos::abs(dot_product(E,B)
           / Kokkos::sqrt(dot_product(B,B)));
    }
    else
      Epar = 0.0;
    },
    Kokkos::Max<double>(maxE)
  );

  Kokkos::fence();

  Real gBalance = Kokkos::min(1.0 + 0.1/maxE, 1.002);;

  *gamma_min = gBalance;

  *p_BC = momentum_(gBalance);
  *p_RE = *p_BC;

  if (Globals::my_rank == 0) {
    std::cout << std::format("Max E field {:.6e}\n", maxE) <<
    std::format("gamma_min, p_BC, p_RE = {:20.14e}, {:20.14e}, {:20.14e}\n", *gamma_min, *p_BC, *p_RE);
  }


  // Zero out locan jre data to start depositing current
  auto jre = f->getJreDataSubview();

  Kokkos::parallel_for(
    PARTHENON_AUTO_LABEL,
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {jre.extent(0), jre.extent(1), jre.extent(2)}),
    // loop over all particles
    KOKKOS_LAMBDA(int i, int j, int k) {
      jre(i,j,k) = 0.0;
    });
  Kokkos::fence();
}

void HybridDriver::PostExecute(parthenon::DriverStatus st) {
  auto pkg = pmesh->packages.Get("Deck");
  auto f = pkg->Param<std::shared_ptr<EM_Field>>("Field");
  auto jre_subview = f->getJreDataSubview();
  auto jre = f->jre_data;
  Kokkos::deep_copy(jre, jre_subview);

  using Host = Kokkos::HostSpace;
  using Unmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
  auto jre_h = Kokkos::create_mirror_view_and_copy(Host(), jre);
  MPI_Allreduce(MPI_IN_PLACE,jre_h.data(),jre_h.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  Kokkos::deep_copy(jre, jre_h);

  auto ts = pkg->Param<std::shared_ptr<int>>("ts");

  const auto filePath = pkg->Param<std::string>("filePath");
  const auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  const auto dt_cd = pkg->Param<Real>("dt_cd");
  const auto dt_mhd = pkg->Param<Real>("dt_mhd");
  std::shared_ptr<Real> gamma_min = pkg->Param<std::shared_ptr<Real>>("gamma_min");
  std::shared_ptr<Real> p_BC = pkg->Param<std::shared_ptr<Real>>("p_BC");
  std::shared_ptr<Real> p_RE = pkg->Param<std::shared_ptr<Real>>("p_RE");

  if (Globals::my_rank == 0) {
    std::cout << std::format("Dumping fields, t = {:.8e}", tm.time) << std::endl;
    dumpToHDF5(*f, *ts, tm.time);
  }

  auto jre_mhd = pkg->Param<Kokkos::View<Real****, Kokkos::LayoutLeft, Host, Unmanaged>>("JreData");
  auto jre_mhd_d = Kokkos::create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), jre_mhd);

  const Real etaec_a3VaB0 = f->etaec_a3VaB0;

  Kokkos::parallel_for(
    PARTHENON_AUTO_LABEL,
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {jre.extent(0), jre.extent(1), jre.extent(2)}),
    // loop over all particles
    KOKKOS_LAMBDA(int i, int j, int k) {
      jre_mhd_d(i,j,k,0) += jre(i,j,k) / dt_mhd * etaec_a3VaB0;
      jre(i,j,k) /= dt_cd;
    });
  Kokkos::fence();
  Kokkos::deep_copy(jre_mhd, jre_mhd_d);

  // Double copy because subviews cannot copy to views throgh host-device interface
  Kokkos::deep_copy(jre_subview, jre);

  if (Globals::my_rank == 0) {
    std::cout << "Interpolating fields!" << std::endl;
  }

  f->interpolate(); // sets jre to electric field
  auto f_d = *f;

  const auto t = tm.time;

  Real maxE;
  Kokkos::parallel_reduce("max E",
  Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {f_d.nR_data,f_d.nZ_data}),
  KOKKOS_LAMBDA(int i, int j, Real& Epar) {
    Dim5 X;
    X[2] = f_d.R0 + f_d.dR * i;
    X[4] = f_d.Z0 + f_d.dZ * j;

    int ii,jj;
    int level = cdg.indicator(X, ii,jj);
    if (level > 0) {
      Dim3 B = {}, curlB = {}, dBdR = {}, dBdZ = {}, E = {}, Jre = {}, V = {}, dbdt = {};

      auto ret = f_d(X, t, B, curlB, dBdR, dBdZ, E, Jre, V, dbdt);
      Epar  = Kokkos::abs(dot_product(E,B)
           / Kokkos::sqrt(dot_product(B,B)));
    }
    else
      Epar = 0.0;
    },
    Kokkos::Max<double>(maxE)
  );

  Kokkos::fence();

  Real gBalance = Kokkos::min(1.0 + 0.1/maxE, 1.002);;

  *gamma_min = gBalance;

  *p_BC = momentum_(gBalance);
  *p_RE = *p_BC;

  if (Globals::my_rank == 0) {
    std::cout << std::format("Max E field {:.6e}\n", maxE) <<
    std::format("gamma_min, p_BC, p_RE = {:20.14e}, {:20.14e}, {:20.14e}\n", *gamma_min, *p_BC, *p_RE);
  }

  auto md = pmesh->mesh_data.Get();
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");

  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  auto field_interpolation = *f;

  const Real p_RE_d = *p_RE;

  const auto c_aw0 = pkg->Param<Real>("c_aw0");
  const auto ct_a = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");
  GuidingCenterEquations<EM_Field, true, false> gce(field_interpolation, c_aw0, ct_a, alpha0);

  Real I_re = 0.0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      // loop over all particles
      KOKKOS_LAMBDA(const int idx, Real &weight) {
        // block and particle indices
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm_r.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n) && (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
          Dim5 X;
          X[0] = pack_swarm_r(b, Kinetic::p(), n);
          X[1] = pack_swarm_r(b, Kinetic::xi(), n);
          X[2] = pack_swarm_r(b, Kinetic::R(), n);
          X[3] = pack_swarm_r(b, Kinetic::phi(), n);
          X[4] = pack_swarm_r(b, Kinetic::Z(), n);
          Real w = pack_swarm_r(b, Kinetic::weight(), n);
          if (X[0] > p_RE_d) {
            weight += getParticleCurrent(X, t, w, field_interpolation);
          }
        }
      },
      I_re);
  MPI_Allreduce(MPI_IN_PLACE,&I_re,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  Real I_re_integral = 0.0;
  Real I_ohmic = 0.0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL,
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {jre.extent(0), jre.extent(1)}),
      // loop over all particles
      KOKKOS_LAMBDA(int i, int j, Real& integral, Real& integral_ohmic) {
        integral += cdg.dR * cdg.dZ * jre(i,j,1);

        Real R = cdg.R0 + i * cdg.dR;
        Real Z = cdg.Z0 + j * cdg.dZ;
        Dim3 B = {}, curlB = {}, dBdR = {}, dBdZ = {}, E = {}, dbdt = {};
        Dim5 X = {0.0,0.0,R,0.0,Z};

        auto ret = field_interpolation(X, t, B, curlB, dBdR, dBdZ, E, dbdt);
        if (ret == SUCCESS) integral_ohmic += cdg.dR * cdg.dZ * curlB[1];
      },
      I_re_integral, I_ohmic);

  Kokkos::fence();

  Real p_phi_total = 0.0;
  Real mu_total = 0.0;
  Real w_total = 0.0;

  int EnableComputeConservedQuantities = pkg->Param<int>("EnableComputeConservedQuantities");

  if (EnableComputeConservedQuantities == 1) {
    Kokkos::View<Real******> psi_hermite_data("psi",
        field_interpolation.hermite_data.extent(0),
        field_interpolation.hermite_data.extent(1),
        field_interpolation.hermite_data.extent(2) + 1,
        field_interpolation.hermite_data.extent(3),
        field_interpolation.hermite_data.extent(6),
        field_interpolation.hermite_data.extent(7));
    Kokkos::parallel_for("psi_compute",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {field_interpolation.nphi_data,field_interpolation.nt}),
    KOKKOS_LAMBDA(int k, int ti){
      auto sbv_hermite_data = Kokkos::subview(field_interpolation.hermite_data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0, Kokkos::ALL, k, ti);
      auto sbv_psi_data = Kokkos::subview(psi_hermite_data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, k, ti);
      computeFlux<2>(sbv_hermite_data, sbv_psi_data, field_interpolation.hR, field_interpolation.hZ);
    });
    Kokkos::parallel_reduce(
        PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
        // loop over all particles
        KOKKOS_LAMBDA(const int idx, Real& p_phi, Real& mu, Real &weight) {
          // block and particle indices
          auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
          const auto swarm_d = pack_swarm_r.GetContext(b);
          if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n) && (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
            Dim5 X;
            X[0] = pack_swarm_r(b, Kinetic::p(), n);
            X[1] = pack_swarm_r(b, Kinetic::xi(), n);
            X[2] = pack_swarm_r(b, Kinetic::R(), n);
            X[3] = pack_swarm_r(b, Kinetic::phi(), n);
            X[4] = pack_swarm_r(b, Kinetic::Z(), n);
            Real w = pack_swarm_r(b, Kinetic::weight(), n);
            weight += w;

            Real my_phi, my_mu;
            Real psi;
            field_interpolation.evalPsi(psi, X, t, psi_hermite_data);
            gce.computeConservedQuantities(X, my_phi, my_mu, t, psi);

            p_phi += my_phi * w;
            mu += my_mu * w;
          }
        },
        p_phi_total, mu_total, w_total);
    MPI_Allreduce(MPI_IN_PLACE,&p_phi_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&mu_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&w_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }


  if (Globals::my_rank == 0) {
    std::ofstream ofs(filePath, std::ios::app);
    ofs << std::format("{:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e}",
        t, I_re * pc::qe * pc::c * .5, I_re_integral * pc::qe * pc::c * .5,
        I_ohmic * 5.3  * 2.0 / pc::mu0,
        p_phi_total, mu_total, w_total) << std::endl;
  }
  *ts += 1;
  Kokkos::fence();
}

}
