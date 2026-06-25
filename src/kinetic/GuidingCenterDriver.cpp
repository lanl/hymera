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
#include "RunawayDriver.h"
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

#include "kinetic/GuidingCenterDriver.h"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "pgen.hpp"
#include "rk4.hpp"


using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

static const bool EnableEfield = false;

namespace Kinetic {

TaskStatus PushParticlesGCE(Mesh *pm, SimTime tm) {
  // get mesh data
  auto md = pm->mesh_data.Get();

  auto pkg = pm->packages.Get("Deck");
  const auto h = pkg->Param<Real>("hRK");
  const auto atol = pkg->Param<Real>("atol");
  const auto rtol = pkg->Param<Real>("rtol");
  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");

  auto f = *(pkg->Param<std::shared_ptr<EM_Field>>("Field"));

  const auto c_aw0 = pkg->Param<Real>("c_aw0");
  const auto ct_a = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");
  GuidingCenterEquations<EM_Field, EnableEfield, false> gce(f, c_aw0, ct_a, alpha0);

  auto RESwarmPackDescriptor = parthenon::MakeSwarmPackDescriptor<p,xi,R,phi,Z,weight,p_phi,mu>("particles");
  auto RESwarmPack = RESwarmPackDescriptor.GetPack(md.get());

  const Real tstart = tm.time;
  const Real tstop = tm.time + tm.dt;

  const auto sa = pkg->Param<SmallAngleCollision<PartialScreening, EnergyScattering, ModifiedCouLog>>("SmallAngleCollision");
  const Real dtSA_min = sa.getSmallAngleCollisionTimestep(momentum_(1.002));
  const Real dtSA_max = tm.dt;

  auto jre = f.getJreDataSubview();

    Kokkos::View<Real******> psi_hermite_data("psi",
        f.hermite_data.extent(0),
        f.hermite_data.extent(1),
        f.hermite_data.extent(2) + 1,
        f.hermite_data.extent(3),
        f.hermite_data.extent(6),
        f.hermite_data.extent(7));
    Kokkos::parallel_for("psi_compute",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {f.nphi_data,f.nt}),
    KOKKOS_LAMBDA(int k, int ti){
      auto sbv_hermite_data = Kokkos::subview(f.hermite_data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0, Kokkos::ALL, k, ti);
      auto sbv_psi_data = Kokkos::subview(psi_hermite_data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, k, ti);
      computeFlux<2>(sbv_hermite_data, sbv_psi_data, f.hR, f.hZ);
    });
  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, RESwarmPack.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        // block and particle indices
        auto [b, n] = RESwarmPack.GetBlockParticleIndices(idx);
        const auto swarm_d = RESwarmPack.GetContext(b);
        const auto markers_d = RESwarmPack.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n) ) {
          Dim5 X;
          Real t = tstart;
          X[0] = RESwarmPack(b, Kinetic::p(), n);
          X[1] = RESwarmPack(b, Kinetic::xi(), n);
          X[2] = RESwarmPack(b, Kinetic::R(), n);
          X[3] = RESwarmPack(b, Kinetic::phi(), n);
          X[4] = RESwarmPack(b, Kinetic::Z(), n);
          Real w = RESwarmPack(b, Kinetic::weight(), n);
          Kokkos::Array<Dim5, 5> work_d;

          bool last_step = false;

          while (last_step == false) {
            Real dtSA = sa.getSmallAngleCollisionTimestep(X[0], dtSA_min, dtSA_max);

            if (t + dtSA > tstop) {
              dtSA = tstop - t;
              if (dtSA < 1e-16) {
                break;
              }
              last_step = true;
            }

 //            auto ret = solve_dopri5(gce, X, t, t + dtSA, rtol, atol, h, 1e-9, std::numeric_limits<int>::max(), work_d);
						auto ret = solve_rk4_fixed(gce, X, t, t+dtSA, h, work_d);
            if (ret != SUCCESS) {
              swarm_d.MarkParticleForRemoval(n);
              break;
            }
            int ii,jj;
            int level = f.cdg.indicator(X, ii,jj);

            DepositCurrent(X, t, w, jre, dtSA, f);
            t += dtSA;
            if (t > tstop)
              break;
          }
          Real my_p_phi, my_mu;
          Real psi;
          f.evalPsi(psi, X, t, psi_hermite_data);
          gce.computeConservedQuantities(X, my_p_phi, my_mu, t, psi);

        	RESwarmPack(b, Kinetic::p(), n)   = X[0];
        	RESwarmPack(b, Kinetic::xi(), n)  = X[1];
        	RESwarmPack(b, Kinetic::R(), n)   = X[2];
        	RESwarmPack(b, Kinetic::phi(), n) = X[3];
        	RESwarmPack(b, Kinetic::Z(), n)   = X[4];
        	RESwarmPack(b, Kinetic::p_phi(), n) = my_p_phi;
        	RESwarmPack(b, Kinetic::mu(), n)   = my_mu;


        }
      });
  Kokkos::fence();

  return TaskStatus::complete;

}

void GuidingCenterDriver::PreExecute() {
  auto pkg = pmesh->packages.Get("Deck");
  auto f = pkg->Param<std::shared_ptr<EM_Field>>("Field");
  // Interpolate the jre data thats there and fields if new
  f->interpolate();

  auto ts = pkg->Param<std::shared_ptr<int>>("ts");
  dumpToHDF5(*f, *ts, tm.time);
  (*ts)++;
}

void GuidingCenterDriver::PostExecute(parthenon::DriverStatus st) {
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


  const auto filePath = pkg->Param<std::string>("filePath");
  const auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  const auto dt_cd = pkg->Param<Real>("dt_cd");
  const auto dt_mhd = pkg->Param<Real>("dt_mhd");

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

  auto md = pmesh->mesh_data.Get();
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");

  auto RESwarmPack = desc_swarm_r.GetPack(md.get());

  auto field_interpolation = *f;

  std::shared_ptr<Real> p_RE = pkg->Param<std::shared_ptr<Real>>("p_RE");
  const Real p_RE_d = *p_RE;

  const auto c_aw0 = pkg->Param<Real>("c_aw0");
  const auto ct_a = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");
  GuidingCenterEquations<EM_Field, EnableEfield, false> gce(field_interpolation, c_aw0, ct_a, alpha0);

  Real I_re = 0.0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, RESwarmPack.GetMaxFlatIndex() + 1,
      // loop over all particles
      KOKKOS_LAMBDA(const int idx, Real &weight) {
        // block and particle indices
        auto [b, n] = RESwarmPack.GetBlockParticleIndices(idx);
        const auto swarm_d = RESwarmPack.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)) {
          Dim5 X;
          X[0] = RESwarmPack(b, Kinetic::p(), n);
          X[1] = RESwarmPack(b, Kinetic::xi(), n);
          X[2] = RESwarmPack(b, Kinetic::R(), n);
          X[3] = RESwarmPack(b, Kinetic::phi(), n);
          X[4] = RESwarmPack(b, Kinetic::Z(), n);
          Real w = RESwarmPack(b, Kinetic::weight(), n);
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
        PARTHENON_AUTO_LABEL, RESwarmPack.GetMaxFlatIndex() + 1,
        // loop over all particles
        KOKKOS_LAMBDA(const int idx, Real& p_phi, Real& mu, Real &weight) {
          // block and particle indices
          auto [b, n] = RESwarmPack.GetBlockParticleIndices(idx);
          const auto swarm_d = RESwarmPack.GetContext(b);
          if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)) {
            Dim5 X;
            X[0] = RESwarmPack(b, Kinetic::p(), n);
            X[1] = RESwarmPack(b, Kinetic::xi(), n);
            X[2] = RESwarmPack(b, Kinetic::R(), n);
            X[3] = RESwarmPack(b, Kinetic::phi(), n);
            X[4] = RESwarmPack(b, Kinetic::Z(), n);
            Real w = RESwarmPack(b, Kinetic::weight(), n);
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
  Kokkos::fence();
}

TaskCollection GuidingCenterDriver::MakeTaskCollection(BlockList_t &blocks, SimTime tm) {
  TaskCollection tc;
  TaskID none(0);

  auto partitions = pmesh->GetDefaultBlockPartitions();
  int num_partitions = partitions.size();

  // note that task within this region that contains one tasklist per pack
  // could still be executed in parallel
  TaskRegion &single_tasklist_per_pack_region = tc.AddRegion(num_partitions);
  for (int i = 0; i < num_partitions; i++) {
    auto &tl = single_tasklist_per_pack_region[i];
    // Initialize the base MeshData for this partition
    // (this automatically initializes the MeshBlockData objects
    // required by this MeshData object)
    auto &mbase = pmesh->mesh_data.Add("base", partitions[i]);

    // add tasks that are per mesh here
    auto push = tl.AddTask(none, PushParticlesGCE, pmesh, tm);
  }
  return tc;
}
}
