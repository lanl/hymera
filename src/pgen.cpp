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

using namespace parthenon;
using namespace parthenon::package::prelude;

#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/FieldEvaluator.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "util/common.hpp"
#include "kinetic/Log.hpp"
#include <hFlux/Evaluator.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>


namespace {

// Interior cell geometry of one mesh block.
struct BlockCells {
  int nx_i, nx_j, nx_k;
  Real dx_i, dx_j, dx_k;
  Real minx_i, minx_j, minx_k;

  explicit BlockCells(parthenon::MeshBlock *pmb) {
    const IndexRange &ib = pmb->cellbounds.GetBoundsI(IndexDomain::interior);
    const IndexRange &jb = pmb->cellbounds.GetBoundsJ(IndexDomain::interior);
    const IndexRange &kb = pmb->cellbounds.GetBoundsK(IndexDomain::interior);
    nx_i = pmb->cellbounds.ncellsi(IndexDomain::interior);
    nx_j = pmb->cellbounds.ncellsj(IndexDomain::interior);
    nx_k = pmb->cellbounds.ncellsk(IndexDomain::interior);
    dx_i = pmb->coords.Dxf<1>(ib.s);
    dx_j = pmb->coords.Dxf<2>(jb.s);
    dx_k = pmb->coords.Dxf<3>(kb.s);
    minx_i = pmb->coords.Xf<1>(ib.s);
    minx_j = pmb->coords.Xf<2>(jb.s);
    minx_k = pmb->coords.Xf<3>(kb.s);
  }
};

// Sets the mesh position (swarm x, y, z) of new particle n: the centre of cell
// n mod (cells in the block). The mesh only load-balances particles; their
// physical position is (R, phi, Z). This is the one place a new particle's mesh
// position is chosen, so mapping momentum onto the mesh coordinate replaces it.
template <typename Pack>
KOKKOS_INLINE_FUNCTION void PlaceInCell(const Pack &pack, const int b, const int n,
                                        const BlockCells &c) {
  int ind = n % (c.nx_i * c.nx_j * c.nx_k);

  int x_i = ind % c.nx_i;
  ind /= c.nx_i;
  int x_j = ind % c.nx_j;
  ind /= c.nx_j; int x_k = ind;

  pack(b, swarm_position::x(), n) = c.minx_i + (x_i+0.5) * c.dx_i;
  pack(b, swarm_position::y(), n) = c.minx_j + (x_j+0.5) * c.dx_j;
  pack(b, swarm_position::z(), n) = c.minx_k + (x_k+0.5) * c.dx_k;
}

// Seeding parameters shared by the generators, from the Deck package.
struct SeedParams {
  Kinetic::RNGPool rng_pool;
  int N;                        ///< particles per block
  Real pmin, pmax, ximin, ximax;
  Real Rmin, Rmax, Zmin, Zmax;
  Real Rc, Zc;                  ///< magnetic axis

  explicit SeedParams(const std::shared_ptr<StateDescriptor> &pkg)
      : rng_pool(pkg->Param<Kinetic::RNGPool>("rng_pool")),
        N(pkg->Param<int>("num_particles_per_block")),
        pmin(pkg->Param<Real>("pmin")), pmax(pkg->Param<Real>("pmax")),
        ximin(pkg->Param<Real>("ximin")), ximax(pkg->Param<Real>("ximax")),
        Rmin(pkg->Param<Real>("Rmin")), Rmax(pkg->Param<Real>("Rmax")),
        Zmin(pkg->Param<Real>("Zmin")), Zmax(pkg->Param<Real>("Zmax")),
        Rc(pkg->Param<Real>("Rc")), Zc(pkg->Param<Real>("Zc")) {}
};

// Cells of a fine (R, Z) grid inside the separatrix (region 1), with |J_phi| at
// the cell centre, from the host copy of the field. The seed density is
// proportional to |J_phi| dR dZ, as in the rejection sampler.
struct SeedCells {
  std::vector<Real> R, Z, J;
  Real dR = 0.0, dZ = 0.0, Jmax = 0.0;
};

template <class FieldH, class IndH>
SeedCells PlasmaCells(const ConfigurationDomainGeometry &cdg, const FieldH &f_h, const IndH &ind,
                      Real Rmin, Real Rmax, Real Zmin, Real Zmax, int nR, int nZ) {
  SeedCells c;
  c.dR = (Rmax - Rmin) / nR;
  c.dZ = (Zmax - Zmin) / nZ;
  const auto &l = cdg.indicator_locator;
  for (int i = 0; i < nR; ++i)
    for (int j = 0; j < nZ; ++j) {
      const Real R = Rmin + (i + 0.5) * c.dR, Z = Zmin + (j + 0.5) * c.dZ;
      int iR, jZ;
      Real xiR, xiZ;
      cdg.locate(R, Z, iR, jZ, xiR, xiZ);
      if (iR < 0 || iR >= cdg.hermite_locator.nR || jZ < 0 || jZ >= cdg.hermite_locator.nZ)
        continue;
      if (R < l.R0 || R >= l.R0 + l.nR * l.dR || Z < l.Z0 || Z >= l.Z0 + l.nZ * l.dZ) continue;
      l.locateCell(R, Z, iR, jZ);
      if (ind(iR, jZ) != 1) continue;
      Kinetic::EvalGCE ev;
      f_h.eval(ev, R, Z, 0.0);
      const Real J = std::abs(ev.J[1]);
      c.R.push_back(R);
      c.Z.push_back(Z);
      c.J.push_back(J);
      c.Jmax = std::max(c.Jmax, J);
    }
  return c;
}

} // namespace

void GenerateParticleCurrentDensity(parthenon::MeshBlock *pmb, parthenon::ParameterInput *pin) {
  auto &data = pmb->meshblock_data.Get();

  // pull out information/global params from package
  auto pkg = pmb->packages.Get("Deck");
  const SeedParams sp(pkg);
  auto rng_pool = sp.rng_pool;
  const int N = sp.N;
  const Real pmin = sp.pmin, pmax = sp.pmax, ximin = sp.ximin, ximax = sp.ximax;
  const Real Rmin = sp.Rmin, Rmax = sp.Rmax, Zmin = sp.Zmin, Zmax = sp.Zmax;
  const Real Rc = sp.Rc, Zc = sp.Zc;

  auto field_data = pkg->Param<Kinetic::FieldData_t>("FieldData");
  const auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  Kinetic::FieldEvaluator f{cdg.hermite_locator, field_data.hermite_data.view_device()};
  Kinetic::FieldEvaluator f_h{cdg.hermite_locator, field_data.hermite_data.view_host()};
  const auto seed_current = pkg->Param<Real>("seed_current");

  field_data.hermite_data.sync_host();
  field_data.hermite_data.sync_device();

  // Plasma cells of a fine grid: the envelope max |J_phi| inside the separatrix
  // for the rejection sampler, and the cells of the stratified seed. Built once.
  // (J on the axis was the envelope: where |J| exceeds it, every point was kept,
  // so the seed followed min(|J|, |J_axis|) rather than |J|.)
  static SeedCells seed_cells;
  if (seed_cells.R.empty()) {
    auto ind = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cdg.indicator_view);
    seed_cells = PlasmaCells(cdg, f_h, ind, Rmin, Rmax, Zmin, Zmax, 1000, 1000);
    Kinetic::Log::Info("Seed: {} plasma cells, max |J_phi| {:g}", seed_cells.R.size(), seed_cells.Jmax);
  }
  const Real Jc = seed_cells.Jmax;
  const int stratified = pkg->Param<int>("seed_stratified");

  // Pull out swarm object
  auto swarm = data->GetSwarmData()->Get("particles");

  // Meshblock geometry
  const BlockCells cells(pmb);

  auto newParticlesContext = swarm->AddEmptyParticles(N);

  // ParticleSeed/stratified = 1: positions from one global systematic sample of
  // |J_phi| dR dZ over the plasma cells ordered by psi, then by the poloidal angle,
  // so every flux-surface band holds its share of the seed to within one particle.
  // Mesh block g takes samples g, g + nb, g + 2 nb, ... (nb blocks in total), which
  // spreads each band over all blocks. The list is the same on every rank.
  Kokkos::View<Real *[2]> pos("seed_pos", stratified ? N : 0);
  if (stratified) {
    static std::vector<int> order;
    static std::vector<Real> cdf;
    const auto &c = seed_cells;
    if (order.empty()) {
      Evaluator pev{field_data.hermite_locator};
      field_data.psi_data.sync_host();
      auto psi_h = field_data.psi_data.view_host();
      std::vector<Real> psi(c.R.size()), th(c.R.size());
      for (size_t k = 0; k < c.R.size(); ++k) {
        pev.evalPsi(psi[k], c.R[k], c.Z[k], psi_h);
        th[k] = std::atan2(c.Z[k] - Zc, c.R[k] - Rc);
      }
      // Bands of psi (|psi| grows away from the axis) of equal seed share; within
      // a band the poloidal angle orders the cells.
      order.resize(c.R.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](int a, int b) { return std::abs(psi[a]) < std::abs(psi[b]); });
      const int nband = pkg->Param<int>("seed_bands");
      Real total = 0.0;
      for (int k : order) total += c.J[k];
      std::vector<int> band(c.R.size());
      Real acc = 0.0;
      for (int k : order) {
        band[k] = std::min(nband - 1, static_cast<int>(nband * acc / total));
        acc += c.J[k];
      }
      std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
        return band[a] != band[b] ? band[a] < band[b] : th[a] < th[b];
      });
      cdf.resize(order.size() + 1, 0.0);
      for (size_t k = 0; k < order.size(); ++k) cdf[k + 1] = cdf[k] + c.J[order[k]] / total;
    }
    const int nb = pmb->pmy_mesh->nbtotal;
    const long long Ntot = static_cast<long long>(N) * nb;
    // One offset for the whole mesh (the same rng_seed on every rank), and a
    // jitter within the cell per particle.
    std::mt19937_64 gen(pkg->Param<int>("seed_rng_base"));
    std::uniform_real_distribution<Real> U(0.0, 1.0);
    const Real u0 = U(gen);
    std::mt19937_64 jit(pkg->Param<int>("seed_rng_base") + 7919 * (pmb->gid + 1));
    auto pos_h = Kokkos::create_mirror_view(pos);
    for (int m = 0; m < N; ++m) {
      const long long q = static_cast<long long>(m) * nb + pmb->gid;
      const Real u = (q + u0) / Ntot;
      const size_t k = std::upper_bound(cdf.begin(), cdf.end(), u) - cdf.begin() - 1;
      const int cell = order[std::min(k, order.size() - 1)];
      pos_h(m, 0) = c.R[cell] + (U(jit) - 0.5) * c.dR;
      pos_h(m, 1) = c.Z[cell] + (U(jit) - 0.5) * c.dZ;
    }
    Kokkos::deep_copy(pos, pos_h);
  }

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
    Kinetic::weight,
    Kinetic::p_phi,
    Kinetic::mu>("particles");
  static auto desc_markers =
    parthenon::MakeSwarmPackDescriptor<
    Kinetic::status
      >("particles");

  auto pack_swarm = desc_swarm.GetPack(data.get());
  auto pack_status = desc_markers.GetPack(data.get());




  // loop over new particles created
  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
      DevExecSpace(), 0,
      newParticlesContext.GetNewParticlesMaxIndex(),
      // new_n ranges from 0 to N_new_particles
      KOKKOS_LAMBDA(const int new_n) {

      // this is the particle index inside the swarm
      const int n = newParticlesContext.GetNewParticleIndex(new_n);

      // The pack spans this one mesh block, so the block index is 0.
      const int b = 0;

      PlaceInCell(pack_swarm, b, n, cells);

      Dim5 X;
      Real t = 0.0;

      // Generate particles:
      // Generate within level according to curl


      if (stratified) {
        auto rng_gen = rng_pool.get_state();
        X[0] = rng_gen.drand(pmin, pmax);
        X[1] = rng_gen.drand(ximin, ximax);
        rng_pool.free_state(rng_gen);
        X[2] = pos(new_n, 0);
        X[4] = pos(new_n, 1);
      } else
      for (;;) {
        // Real theta = theta_dist(RNGs[thread_id]);
        auto rng_gen = rng_pool.get_state();

        Real randNum = rng_gen.drand(abs(Jc));
        X[0] = rng_gen.drand(pmin, pmax);
        X[1] = rng_gen.drand(ximin, ximax);
        X[2] = rng_gen.drand(Rmin, Rmax);
        X[4] = rng_gen.drand(Zmin, Zmax);

        rng_pool.free_state(rng_gen);

        int i, j;
        Real xiR, xiZ;
        cdg.locate(X[2], X[4], i, j, xiR, xiZ);
        if (i < 0 || i >= cdg.hermite_locator.nR || j < 0 || j >= cdg.hermite_locator.nZ)
          continue;

        Region level;
        cdg.locate_region(X[2], X[4], i, j, level);
        // Only keep particles within separatrix
        if (level != 1)
          continue;
        Kinetic::EvalGCE ev;
        f.eval(ev, X[2], X[4], t);
        if (randNum < abs(ev.J[1])) break;
      }

      Real my_phi = 0.0, my_mu = 0.0;

      pack_swarm(b, Kinetic::p(), n)   = X[0];
      pack_swarm(b, Kinetic::xi(), n)  = X[1];
      pack_swarm(b, Kinetic::R(), n)   = X[2];
      pack_swarm(b, Kinetic::phi(), n) = X[3];
      pack_swarm(b, Kinetic::Z(), n)   = X[4];
      pack_swarm(b, Kinetic::p_phi(), n) = my_phi;
      pack_swarm(b, Kinetic::mu(), n)   = my_mu;

      // set weights to 1
      pack_swarm(b, Kinetic::weight(), n) = 1.0;
      pack_status(b, Kinetic::status(), n) = Kinetic::ALIVE | Kinetic::PROTECTED;
      });

}

void GenerateParticlePoint(parthenon::MeshBlock *pmb, parthenon::ParameterInput *pin) {

  auto &data = pmb->meshblock_data.Get();

  // pull out information/global params from package
  auto pkg = pmb->packages.Get("Deck");
  const SeedParams sp(pkg);
  auto rng_pool = sp.rng_pool;
  const int N = sp.N;
  const Real pmin = sp.pmin, pmax = sp.pmax, ximin = sp.ximin, ximax = sp.ximax;
  const Real Rseed  = pkg->Param<Real>("Rseed");
  const Real Zseed  = pkg->Param<Real>("Zseed");

  // Not used here; read so that run scripts passing ParticleSeed/r_0 keep working
  // and the parameter dump keeps recording it.
  [[maybe_unused]] const Real r_0 = pin->GetOrAddReal("ParticleSeed", "r_0", 0.0);


  // Pull out swarm object
  auto swarm = data->GetSwarmData()->Get("particles");

  // Meshblock geometry
  const BlockCells cells(pmb);

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


  Kinetic::Log::Debug("Block {}: placing {} particles at (R, Z) = ({}, {})", pmb->gid, N, Rseed, Zseed);

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
      DevExecSpace(), 0,
      newParticlesContext.GetNewParticlesMaxIndex(),
      // new_n ranges from 0 to N_new_particles
      KOKKOS_LAMBDA(const int new_n) {

      // this is the particle index inside the swarm
      const int n = newParticlesContext.GetNewParticleIndex(new_n);

      // The pack spans this one mesh block, so the block index is 0.
      const int b = 0;

      PlaceInCell(pack_swarm, b, n, cells);

      Dim5 X;

      auto rng_gen = rng_pool.get_state();
      X[0] = rng_gen.drand(pmin, pmax);
      X[1] = rng_gen.drand(ximin, ximax);
      rng_pool.free_state(rng_gen);

      X[2] = Rseed;
      X[4] = Zseed;

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

// Ring seed of the fig11a avalanche benchmark: p and xi uniform as in
// GenerateParticlePoint, and the poloidal angle theta uniform on the circle of
// radius ParticleSeed/r_0 around the magnetic axis (Rc, Zc). The draws are made
// in the order p, xi, theta, as the generator that produced the reference
// results did. With seed_on_surface the particles are moved later, along the
// ray at angle theta, to a flux surface (MoveSeedsToFluxSurface).
void GenerateParticleRing(parthenon::MeshBlock *pmb, parthenon::ParameterInput *pin) {
  auto &data = pmb->meshblock_data.Get();

  auto pkg = pmb->packages.Get("Deck");
  const SeedParams sp(pkg);
  auto rng_pool = sp.rng_pool;
  const int N = sp.N;
  const Real pmin = sp.pmin, pmax = sp.pmax, ximin = sp.ximin, ximax = sp.ximax;
  const Real Rc = sp.Rc, Zc = sp.Zc;
  const Real r_0 = pkg->Param<Real>("seed_r0");

  auto swarm = data->GetSwarmData()->Get("particles");
  const BlockCells cells(pmb);
  auto newParticlesContext = swarm->AddEmptyParticles(N);

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

  Kinetic::Log::Debug("Block {}: placing {} particles on the ring r_0 = {} around ({}, {})",
                      pmb->gid, N, r_0, Rc, Zc);

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
      DevExecSpace(), 0,
      newParticlesContext.GetNewParticlesMaxIndex(),
      KOKKOS_LAMBDA(const int new_n) {
      const int n = newParticlesContext.GetNewParticleIndex(new_n);
      const int b = 0;

      PlaceInCell(pack_swarm, b, n, cells);

      auto rng_gen = rng_pool.get_state();
      const Real p = rng_gen.drand(pmin, pmax);
      const Real xi = rng_gen.drand(ximin, ximax);
      const Real theta = rng_gen.drand(0.0, 2 * M_PI);
      rng_pool.free_state(rng_gen);

      pack_swarm(b, Kinetic::p(), n)   = p;
      pack_swarm(b, Kinetic::xi(), n)  = xi;
      pack_swarm(b, Kinetic::R(), n)   = Rc + r_0 * Kokkos::cos(theta);
      pack_swarm(b, Kinetic::phi(), n) = 0.0;
      pack_swarm(b, Kinetic::Z(), n)   = Zc + r_0 * Kokkos::sin(theta);

      pack_swarm(b, Kinetic::weight(), n) = 1.0;
      pack_status(b, Kinetic::status(), n) = Kinetic::ALIVE | Kinetic::PROTECTED;
   });
}
