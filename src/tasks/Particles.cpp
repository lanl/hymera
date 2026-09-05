#include "Tasks.h"
#include <fstream>
#include <format>
#include <globals.hpp>
#include "kinetic/ParticleVerificator.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/FieldEvaluator.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/rk45.hpp"
#include "kinetic/rk4.hpp"

// Number of scratch Dim5 slots each stepper's work array needs. DP45 uses
// slots 0..9 (rk45.hpp); RK4 uses slots 0..4 (rk4.hpp). Sizing the per-thread
// stack array to exactly what the chosen stepper touches keeps the GPU
// local-memory / register footprint minimal (see plan notes on register burden).
template <int METHOD> struct WorkSize;
template <> struct WorkSize<INTEGRATOR_DOPRI5> { static constexpr int value = 10; };
template <> struct WorkSize<INTEGRATOR_RK4>    { static constexpr int value = 5; };

// Advance a single particle's state X from t to t+dtSA with the compile-time
// selected stepper. Returns the verificator result code.
template <int METHOD, class GCE, class Ver, class Work>
KOKKOS_INLINE_FUNCTION typename Ver::ResultCode_t
StepParticle(const GCE &gce, const Ver &ver, Dim5 &X, const Real t,
             const Real dtSA, const Real rtol, const Real atol, const Real h,
             Work &work_d) {
  if constexpr (METHOD == INTEGRATOR_RK4) {
    return solve_rk4_fixed(gce, ver, X, t, t + dtSA, h, work_d);
  } else {
    return solve_rk45(gce, ver, X, t, t + dtSA, rtol, atol, h, 1e-9,
                      std::numeric_limits<int>::max(), work_d);
  }
}

template <int METHOD>
static TaskStatus PushParticlesImpl(Mesh *pm, Real t0, Real dt) {

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

//  const Real dtSA_min = Kokkos::max(h, sa.getSmallAngleCollisionTimestep(momentum_(1.002)));
  const Real dtSA_min = sa.getSmallAngleCollisionTimestep(momentum_(1.002));
  const Real dtSA_max = dt;

  const auto c_aw0  = pkg->Param<Real>("c_aw0");
  const auto ct_a   = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");


  auto data = pkg->Param<Kinetic::FieldData_t>("FieldData");
  auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  Kinetic::FieldEvaluator f{cdg.hermite_locator, data.hermite_data.view_device()};
  GuidingCenterEquations<decltype(f), true, false> gce(f, c_aw0, ct_a, alpha0);
  ParticleVerificator ver{cdg, p_BC};

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

  auto jre = pkg->Param<View3>("Jre");

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
          Kokkos::Array<Dim5, WorkSize<METHOD>::value> work_d;

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

            auto ret = StepParticle<METHOD>(gce, ver, X, t, dtSA, rtol, atol, h,
                                            work_d);
            if (ret != ParticleVerificator::Success) {

              pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
              swarm_d.MarkParticleForRemoval(n);
//              if (ret == ParticleVerifyCodes::MomentumCutoff) {
//                pack_swarm_i(b, Kinetic::status(), n) |= Kinetic::DEATH_BY_MOMENTUM;
//              } else if (ret == ParticleVerifyCodes::WallImpact) {
//                pack_swarm_i(b, Kinetic::status(), n) |= Kinetic::DEATH_BY_WALL;
//              }
              break;
            }

            if (X[0] > p_RE) {
              DepositCurrent(X, t, w, jre, dtSA, f, cdg.indicator_locator);
            }
            if (EnableSmallAngleCollisions == 1) {
              sa(X[0], X[1], dtSA, rng_pool);
            }

            t += dtSA;
            if (t > tstop)
              break;
          }

        	pack_swarm_r(b, Kinetic::p(), n)   = X[0];
        	pack_swarm_r(b, Kinetic::xi(), n)  = X[1];
        	pack_swarm_r(b, Kinetic::R(), n)   = X[2];
        	pack_swarm_r(b, Kinetic::phi(), n) = X[3];
        	pack_swarm_r(b, Kinetic::Z(), n)   = X[4];

          if (EnableLargeAngleCollisions == 1 &&
              (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
            pack_swarm_i(b, Kinetic::will_scatter(), n) = ms(X[0], w, dt, gamma_min, rng_pool);
          }
        }
      });

  return TaskStatus::complete;

}

// Dispatch to the compile-time specialized push based on the runtime
// Simulation/integrator param. Branching here (host side) rather than inside the
// kernel means each stepper compiles to its own kernel with its own minimal work
// array — the DP45 kernel's register cost is not paid on RK4 runs and vice versa.
TaskStatus PushParticles(Mesh *pm, Real t0, Real dt) {
  const int integrator = pm->packages.Get("Deck")->Param<int>("integrator");
  if (integrator == INTEGRATOR_RK4)
    return PushParticlesImpl<INTEGRATOR_RK4>(pm, t0, dt);
  return PushParticlesImpl<INTEGRATOR_DOPRI5>(pm, t0, dt);
}

// Compute the two guiding-center adiabatic invariants (canonical toroidal
// momentum p_phi and magnetic moment mu) for every alive particle, store them
// in the particle.p_phi / particle.mu swarm variables (for per-particle .phdf
// output) and log a weight-averaged aggregate to a text file for a quick
// conservation check. Gated by the "EnableComputeConservedQuantities" param.
TaskStatus ComputeConservedQuantities(Mesh *pm, Real t) {

  auto pkg = pm->packages.Get("Deck");
  const int EnableComputeConservedQuantities =
      pkg->Param<int>("EnableComputeConservedQuantities");
  if (EnableComputeConservedQuantities != 1) return TaskStatus::complete;

  auto md = pm->mesh_data.Get();

  const auto c_aw0  = pkg->Param<Real>("c_aw0");
  const auto ct_a   = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");

  auto data = pkg->Param<Kinetic::FieldData_t>("FieldData");
  auto cdg  = pkg->Param<ConfigurationDomainGeometry>("CDG");
  Kinetic::FieldEvaluator f{cdg.hermite_locator, data.hermite_data.view_device()};
  GuidingCenterEquations<decltype(f), true, false> gce(f, c_aw0, ct_a, alpha0);

  // Psi (poloidal flux) is stored separately from B/J/E; evaluate it with the
  // hFlux Taylor evaluator on the psi_data grid (filled during field init).
  Evaluator psi_ev{data.hermite_locator};
  auto psi_view = data.psi_data.view_device();

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      Kinetic::p, Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z,
      Kinetic::weight, Kinetic::p_phi, Kinetic::mu>("particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  Real p_phi_total = 0.0;
  Real mu_total = 0.0;
  Real w_total = 0.0;

  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      KOKKOS_LAMBDA(const int idx, Real &lp, Real &lmu, Real &lw) {
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm_r.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n) &&
            (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
          const Real p  = pack_swarm_r(b, Kinetic::p(), n);
          const Real xi = pack_swarm_r(b, Kinetic::xi(), n);
          const Real R  = pack_swarm_r(b, Kinetic::R(), n);
          const Real Z  = pack_swarm_r(b, Kinetic::Z(), n);
          const Real w  = pack_swarm_r(b, Kinetic::weight(), n);

          Real Psi = 0.0;
          psi_ev.evalPsi(Psi, R, Z, psi_view);

          Real p_phi = 0.0, mu = 0.0;
          gce.computeConservedQuantities(p_phi, mu, p, xi, R, Z, t, Psi);

          pack_swarm_r(b, Kinetic::p_phi(), n) = p_phi;
          pack_swarm_r(b, Kinetic::mu(), n)    = mu;

          lp  += w * p_phi;
          lmu += w * mu;
          lw  += w;
        }
      },
      p_phi_total, mu_total, w_total);
  Kokkos::fence();

  MPI_Allreduce(MPI_IN_PLACE, &p_phi_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &mu_total,    1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &w_total,     1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (Globals::my_rank == 0) {
    const Real inv_w = (w_total > 0.0) ? 1.0 / w_total : 0.0;
    std::ofstream ofs(pkg->Param<std::string>("conservation_log"), std::ios::app);
    ofs << std::format("{:20.14e} {:20.14e} {:20.14e} {:20.14e}",
                       t, p_phi_total * inv_w, mu_total * inv_w, w_total)
        << std::endl;
  }

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
  auto gamma_min = pkg->Param<Real>("gamma_min");
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

  // ntot must be final
  Kokkos::fence();
  if (ntot > 0) {
    //std::cout << std::format("Adding {} new particles, total alive {}, ratio {}", ntot, nalive, (Real) (ntot + nalive) / (Real) nalive) << std::endl;

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
  }

	return TaskStatus::complete;
}

TaskStatus RandomRemove(Mesh* pm) {
  std::cout << "Random remove start" << std::endl;
  auto pkg = pm->packages.Get("Deck");
  int num_particles = 0;

  auto desc_swarm = parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");
  auto md = pm->mesh_data.Get();
  auto pack_swarm = desc_swarm.GetPack(md.get());

  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm.GetMaxFlatIndex() + 1,
      // loop over all particles
      KOKKOS_LAMBDA(const int idx, int& num) {
        // block and particle indices
        auto [b, n] = pack_swarm.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n) && (pack_swarm(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
          num++;
        }
      },
      num_particles);
  MPI_Allreduce(MPI_IN_PLACE,&num_particles,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  std::cout << "Number of particles " << num_particles << std::endl;
  auto num_particles_max = pkg->Param<int>("num_particles_max");
  if (num_particles_max > num_particles) return TaskStatus::complete;

  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  auto desc_swarm_i = parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<Kinetic::weight>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_i.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        auto [b_i, n_i] = pack_swarm_i.GetBlockParticleIndices(idx);
        auto [b_r, n_r] = pack_swarm_r.GetBlockParticleIndices(idx);
        if (pack_swarm_i(b_i, Kinetic::status(), n_i) & Kinetic::ALIVE) {
          auto rng_gen = rng_pool.get_state();
          auto isKilled = (rng_gen.urand(2));
          rng_pool.free_state(rng_gen);
          // block and particle indices

          if (isKilled == 1) {
              pack_swarm_i(b_i, Kinetic::status(), n_i) &= ~Kinetic::ALIVE;
              pack_swarm_i(b_i, Kinetic::status(), n_i) &= ~Kinetic::PROTECTED;
              const auto swarm = pack_swarm_i.GetContext(b_i);
              swarm.MarkParticleForRemoval(n_i);
          } else {
              pack_swarm_r(b_r, Kinetic::weight(), n_r) *= 2.0;
          }
        }
      });

  return TaskStatus::complete;
}

