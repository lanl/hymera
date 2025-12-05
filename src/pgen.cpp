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

#include <parthenon/package.hpp>
#include <utils/error_checking.hpp>
using namespace parthenon::package::prelude;

#include "kinetic/kinetic.hpp"
#include "kinetic/EM_Field.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "util/common.hpp"

void GenerateParticleRing(parthenon::MeshBlock *pmb, parthenon::ParameterInput *pin) {

  auto pkg = pmb->packages.Get("Deck");

  const Real pmin  = pkg->Param<Real>("pmin");
  const Real pmax  = pkg->Param<Real>("pmax");

  const Real ximin = pkg->Param<Real>("ximin");
  const Real ximax = pkg->Param<Real>("ximax");

  const Real Rc = pkg->Param<Real>("Rcenter");
  const Real Zc = pkg->Param<Real>("Zcenter");
  const Real r_0 = pkg->Param<Real>("r_0");

  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  const int N = pkg->Param<int>("num_particles_per_block");

  // Pull out swarm object
  auto &data = pmb->meshblock_data.Get();
  auto swarm = data->GetSwarmData()->Get("particles");

  // Meshblock geometry
  const IndexRange &ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  const IndexRange &jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  const IndexRange &kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  const int &nx_i = pmb->cellbounds.ncellsi(IndexDomain::interior);
  const int &nx_j = pmb->cellbounds.ncellsj(IndexDomain::interior);
  const int &nx_k = pmb->cellbounds.ncellsk(IndexDomain::interior);
  const Real &dx_i = pmb->coords.Dxf<1>(ib.s);
  const Real &dx_j = pmb->coords.Dxf<2>(jb.s);
  const Real &dx_k = pmb->coords.Dxf<3>(kb.s);
  const Real &minx_i = pmb->coords.Xf<1>(ib.s);
  const Real &minx_j = pmb->coords.Xf<2>(jb.s);
  const Real &minx_k = pmb->coords.Xf<3>(kb.s);

  // Create an accessor to particles, allocate particles
  auto newParticlesContext = swarm->AddEmptyParticles(N);

  // Make a SwarmPack via types to get positions
  static auto desc_swarm =
    parthenon::MakeSwarmPackDescriptor<
    swarm_position::x,
    swarm_position::y,
    swarm_position::z,
    Kinetic::p,
    Kinetic::xi,
    Kinetic::R,
    Kinetic::phi,
    Kinetic::Z,
    Kinetic::weight>("particles");
  static auto desc_markers =
    parthenon::MakeSwarmPackDescriptor<
    Kinetic::will_scatter,
    Kinetic::secondary_index
      >("particles");

  auto pack_swarm = desc_swarm.GetPack(data.get());
  auto pack_markers = desc_markers.GetPack(data.get());

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
      DevExecSpace(), 0,
      newParticlesContext.GetNewParticlesMaxIndex(),
      // new_n ranges from 0 to N_new_particles
      KOKKOS_LAMBDA(const int new_n) {
      const int n = newParticlesContext.GetNewParticleIndex(new_n);
      const int b = 0;
      auto rng_gen = rng_pool.get_state();
      pack_swarm(b, Kinetic::p(), n) = rng_gen.drand(pmin, pmax);
      pack_swarm(b, Kinetic::xi(), n) = rng_gen.drand(ximin, ximax);
      pack_swarm(b, Kinetic::phi(), n) = rng_gen.drand(0, 2.0*M_PI);
      Real theta = rng_gen.drand(-M_PI, M_PI);
      rng_pool.free_state(rng_gen);


      int ind = n % (nx_i * nx_j * nx_k);

      int x_i = ind % nx_i;
      ind /= nx_i;
      int x_j = ind % nx_j;
      ind /= nx_j;
      int x_k = ind;

      pack_swarm(b, swarm_position::x(), n) = minx_i + (x_i+0.5) * dx_i;
      pack_swarm(b, swarm_position::y(), n) = minx_j + (x_j+0.5) * dx_j;
      pack_swarm(b, swarm_position::z(), n) = minx_k + (x_k+0.5) * dx_k;

      pack_swarm(b, Kinetic::R(), n) = Rc + r_0 * cos(theta);
      pack_swarm(b, Kinetic::Z(), n) = Zc + r_0 * sin(theta);
      pack_swarm(b, Kinetic::weight(), n) = 1.0;
      pack_markers(b, Kinetic::will_scatter(), n) = 0;
      pack_markers(b, Kinetic::secondary_index(), n) = 0;

      }
  );
}

void GenerateParticleSquare(parthenon::MeshBlock *pmb, parthenon::ParameterInput *pin) {
  auto &data = pmb->meshblock_data.Get();

  // pull out information/global params from package
  auto pkg = pmb->packages.Get("Kinetic");
  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  const int N = pkg->Param<int>("num_particles_per_block");
  const Real pmin  = pkg->Param<Real>("pmin");
  const Real pmax  = pkg->Param<Real>("pmax");
  const Real ximin = pkg->Param<Real>("ximin");
  const Real ximax = pkg->Param<Real>("ximax");
  const Real Rmin  = pkg->Param<Real>("Rmin");
  const Real Rmax  = pkg->Param<Real>("Rmax");
  const Real Zmin  = pkg->Param<Real>("Zmin");
  const Real Zmax  = pkg->Param<Real>("Zmax");

  // Pull out swarm object
  auto swarm = data->GetSwarmData()->Get("particles");

  // Meshblock geometry
  const IndexRange &ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  const IndexRange &jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  const IndexRange &kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  const int &nx_i = pmb->cellbounds.ncellsi(IndexDomain::interior);
  const int &nx_j = pmb->cellbounds.ncellsj(IndexDomain::interior);
  const int &nx_k = pmb->cellbounds.ncellsk(IndexDomain::interior);
  const Real &dx_i = pmb->coords.Dxf<1>(ib.s);
  const Real &dx_j = pmb->coords.Dxf<2>(jb.s);
  const Real &dx_k = pmb->coords.Dxf<3>(kb.s);
  const Real &minx_i = pmb->coords.Xf<1>(ib.s);
  const Real &minx_j = pmb->coords.Xf<2>(jb.s);
  const Real &minx_k = pmb->coords.Xf<3>(kb.s);


  // Create an accessor to particles, allocate particles
  auto newParticlesContext = swarm->AddEmptyParticles(N);

  // Make a SwarmPack via types to get positions
  static auto desc_swarm =
    parthenon::MakeSwarmPackDescriptor<
    swarm_position::x,
    swarm_position::y,
    swarm_position::z,
    Kinetic::p,
    Kinetic::xi,
    Kinetic::R,
    Kinetic::phi,
    Kinetic::Z,
    Kinetic::weight>("particles");
  auto pack_swarm = desc_swarm.GetPack(data.get());
  auto swarm_d = swarm->GetDeviceContext();

  // loop over new particles created
  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
      DevExecSpace(), 0,
      newParticlesContext.GetNewParticlesMaxIndex(),
      // new_n ranges from 0 to N_new_particles
      KOKKOS_LAMBDA(const int new_n) {
      // this is the particle index inside the swarm
      const int n = newParticlesContext.GetNewParticleIndex(new_n);
      auto rng_gen = rng_pool.get_state();

      // Normally b would be free-floating and set by pack.GetBlockparticleIndices
      // but since we're on a single meshblock for this loop, it's just 0
      // because block index = 0
      const int b = 0;
      //auto [b, n] = pack_swarm.GetBlockparticleIndices(idx);

      // Find the cell index in 1D
      int ind = n % (nx_i * nx_j * nx_k);

      int x_i = ind % nx_i;
      ind /= nx_i;
      int x_j = ind % nx_j;
      ind /= nx_j;
      int x_k = ind;

      pack_swarm(b, swarm_position::x(), n) = minx_i + (x_i+0.5) * dx_i;
      pack_swarm(b, swarm_position::y(), n) = minx_j + (x_j+0.5) * dx_j;
      pack_swarm(b, swarm_position::z(), n) = minx_k + (x_k+0.5) * dx_k;
      // randomly sample particle positions

      // set canonical momentum
      pack_swarm(b, Kinetic::p(), n) = rng_gen.drand(pmin, pmax);
      pack_swarm(b, Kinetic::xi(), n) = rng_gen.drand(ximin, ximax);
      // randomly sample particle positions (TODO: r, phi, z)
      pack_swarm(b, Kinetic::R(), n) = rng_gen.drand(Rmin, Rmax);
      pack_swarm(b, Kinetic::phi(), n) = rng_gen.drand(0, 2.0*M_PI);
      pack_swarm(b, Kinetic::Z(), n) = rng_gen.drand(Zmin, Zmax);

      // set weights to 1
      pack_swarm(b, Kinetic::weight(), n) = 1.0;

      // release random number generator
      rng_pool.free_state(rng_gen);
      });
}

void GenerateParticleCurrentDensity(parthenon::MeshBlock *pmb, parthenon::ParameterInput *pin) {
  auto &data = pmb->meshblock_data.Get();

  // pull out information/global params from package
  auto pkg = pmb->packages.Get("Deck");
  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  int N = pkg->Param<int>("num_particles_per_block");
  const Real pmin  = pkg->Param<Real>("pmin");
  const Real pmax  = pkg->Param<Real>("pmax");
  const Real ximin = pkg->Param<Real>("ximin");
  const Real ximax = pkg->Param<Real>("ximax");
  const Real Rmin  = pkg->Param<Real>("Rmin");
  const Real Rmax  = pkg->Param<Real>("Rmax");
  const Real Zmin  = pkg->Param<Real>("Zmin");
  const Real Zmax  = pkg->Param<Real>("Zmax");
  const Real Rc  = pkg->Param<Real>("Rc");
  const Real Zc  = pkg->Param<Real>("Zc");

  const auto f = pkg->Param<std::shared_ptr<EM_Field>>("Field");
  auto field_interpolation = *f;
  field_interpolation.t_a = 0.0;
  field_interpolation.t_b = 1.0;;
  const auto cdg = field_interpolation.cdg;
  const auto seed_current = pkg->Param<Real>("seed_current");
  const auto p_RE = pkg->Param<Real>("p_RE");

  KOKKOS_ASSERT(p_RE < pmax);

  // Pull out swarm object
  auto swarm = data->GetSwarmData()->Get("particles");

  // Meshblock geometry
  const IndexRange &ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
  const IndexRange &jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
  const IndexRange &kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
  const int &nx_i = pmb->cellbounds.ncellsi(IndexDomain::interior);
  const int &nx_j = pmb->cellbounds.ncellsj(IndexDomain::interior);
  const int &nx_k = pmb->cellbounds.ncellsk(IndexDomain::interior);
  const Real &dx_i = pmb->coords.Dxf<1>(ib.s);
  const Real &dx_j = pmb->coords.Dxf<2>(jb.s);
  const Real &dx_k = pmb->coords.Dxf<3>(kb.s);
  const Real &minx_i = pmb->coords.Xf<1>(ib.s);
  const Real &minx_j = pmb->coords.Xf<2>(jb.s);
  const Real &minx_k = pmb->coords.Xf<3>(kb.s);


  // Create an accessor to particles, allocate particles
  // Dirty fix for darwin GH

  //if (parthenon::Globals::my_rank == 0) N = 0;
  auto newParticlesContext = swarm->AddEmptyParticles(N);

  // Make a SwarmPack via types to get positions
  static auto desc_swarm =
    parthenon::MakeSwarmPackDescriptor<
    swarm_position::x,
    swarm_position::y,
    swarm_position::z,
    Kinetic::p,
    Kinetic::xi,
    Kinetic::R,
    Kinetic::phi,
    Kinetic::Z,
    Kinetic::weight>("particles");
  static auto desc_markers =
    parthenon::MakeSwarmPackDescriptor<
    Kinetic::status
      >("particles");

  auto pack_swarm = desc_swarm.GetPack(data.get());
  auto pack_status = desc_markers.GetPack(data.get());

  auto indicator_h = create_mirror_view_and_copy(Kokkos::HostSpace(), cdg.indicator_view);
  FILE* fs = std::fopen("indicator.txt", "w");
  std::fprintf(fs, "%le, %le, %le, %le\n", cdg.R0, cdg.Z0, cdg.dR, cdg.dZ);
  for (int i = 0; i < indicator_h.extent(0); ++i) {
    for (int j = 0; j < indicator_h.extent(1); ++j) {
      std::fprintf(fs, "%d ", indicator_h(i,j));
    }
    std::fprintf(fs, "\n");
  }
  std::fclose(fs);

  // loop over new particles created
  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
      DevExecSpace(), 0,
      newParticlesContext.GetNewParticlesMaxIndex(),
      // new_n ranges from 0 to N_new_particles
      KOKKOS_LAMBDA(const int new_n) {

      // this is the particle index inside the swarm
      const int n = newParticlesContext.GetNewParticleIndex(new_n);

      // Normally b would be free-floating and set by pack.GetBlockparticleIndices
      // but since we're on a single meshblock for this loop, it's just 0
      // because block index = 0
      const int b = 0;
      //auto [b, n] = pack_swarm.GetBlockparticleIndices(idx);

      // Find the cell index in 1D
      int ind = n % (nx_i * nx_j * nx_k);

      int x_i = ind % nx_i;
      ind /= nx_i;
      int x_j = ind % nx_j;
      ind /= nx_j; int x_k = ind;

      pack_swarm(b, swarm_position::x(), n) = minx_i + (x_i+0.5) * dx_i;
      pack_swarm(b, swarm_position::y(), n) = minx_j + (x_j+0.5) * dx_j;
      pack_swarm(b, swarm_position::z(), n) = minx_k + (x_k+0.5) * dx_k;
      pack_swarm(b, Kinetic::R(), n) = Rc;
      pack_swarm(b, Kinetic::Z(), n) = Zc;
      // randomly sample particle positions

      Dim5 X;
      Real t = 0.0;
      X[2] = Rc;
      X[4] = Zc;
      int i, j;
      int level = cdg.indicator(X, i, j);
      KOKKOS_ASSERT(level == 2);


      Dim3 B = {}, dBdR = {}, dBdZ = {}, curlB = {}, E = {}, dbdt = {};
      Dim3 B_center = {}, curlB_center = {};
      ERROR_CODE status = field_interpolation(X, t, B_center, curlB_center, dBdR, dBdZ, E, dbdt);
      KOKKOS_ASSERT(status == SUCCESS);

      // Generate particles:
      // Generate within level according to curl

      for (;;) {
        // Real theta = theta_dist(RNGs[thread_id]);
        auto rng_gen = rng_pool.get_state();

        Real randNum = rng_gen.drand(abs(curlB_center[1]));
        X[0] =  rng_gen.drand(pmin, pmax);
        X[1] =  rng_gen.drand(ximin, ximax);
        X[2] = rng_gen.drand(Rmin, Rmax);
        X[4] = rng_gen.drand(Zmin, Zmax);

        rng_pool.free_state(rng_gen);

        int i, j;
        int level = cdg.indicator(X, i, j);
        if (level != 2)
          continue;
        status = field_interpolation(X, t, B, curlB, dBdR, dBdZ, E, dbdt);
        KOKKOS_ASSERT(status == SUCCESS);

        if (randNum < abs(curlB[1])) break;
      }

      pack_swarm(b, Kinetic::p(), n)   = X[0];
      pack_swarm(b, Kinetic::xi(), n)  = X[1];
      pack_swarm(b, Kinetic::R(), n)   = X[2];
      pack_swarm(b, Kinetic::phi(), n) = X[3];
      pack_swarm(b, Kinetic::Z(), n)   = X[4];

      // set weights to 1
      pack_swarm(b, Kinetic::weight(), n) = 1.0;
        pack_status(b, Kinetic::status(), n) = Kinetic::ALIVE | Kinetic::PROTECTED;
      });

}
