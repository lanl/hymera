#include "Tasks.h"
#include <fstream>
#include <format>
#include <globals.hpp>
#include "kinetic/ParticleVerificator.hpp"
#include "kinetic/Log.hpp"
#include "kinetic/DiagnosticsFile.hpp"
#include "kinetic/PushTiming.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/FieldEvaluator.hpp"
#include "kinetic/AnalyticFieldEvaluator.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/rk45.hpp"
#include "kinetic/rk4.hpp"
#include "kinetic/ScratchWork.hpp"
#include "kinetic/FeedbackGain.hpp"

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

// Runs body(idx) for idx in [0, max_index], over all particles of the pack.
// Simulation/parallel_push = 0 (default) keeps Parthenon's DEFAULT_LOOP_PATTERN,
// which on a host build is one serial loop: random numbers are drawn in particle
// order, so collision runs are reproducible bit for bit and independent of the
// thread count. 1 runs a Kokkos parallel loop over all threads: much faster on
// CPU, but with collisions the random draws per particle depend on the thread
// schedule. On CUDA both settings give the same parallel kernel.
template <typename Body>
static void ParticleLoop(const std::shared_ptr<StateDescriptor> &pkg, const int max_index,
                         const Body &body) {
  if (pkg->Param<int>("parallel_push") == 1)
    parthenon::par_for(parthenon::loop_pattern_flatrange_tag, PARTHENON_AUTO_LABEL,
                       DevExecSpace(), 0, max_index, body);
  else
    parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                       DevExecSpace(), 0, max_index, body);
}

// The push of every alive particle over [t0, t0 + dt] in field f, checked by
// ver after each step. Deposit = false skips the current deposit (no grid).
//
// Scratch = false (default): one thread per particle slot, with the Runge-Kutta
// work array a per-thread Kokkos::Array (registers, spilling to local memory).
// Scratch = true (Simulation/push_scratch = 1): a team kernel; each team pushes
// push_team_size consecutive slots and holds their work arrays in Parthenon team
// scratch (ScratchPad1D, level 0: shared memory on a GPU). The push body is the
// same code in both, so the arithmetic is identical.
template <int METHOD, bool Deposit, bool Scratch, class Field, class Ver>
static TaskStatus PushParticlesImpl(Mesh *pm, Real t0, Real dt, const Field &f,
                                    const Ver &ver) {

  // get mesh data
  auto md = pm->mesh_data.Get();
  auto pkg = pm->packages.Get("Deck");

  const auto h = pkg->Param<Real>("hRK");
  const auto atol = pkg->Param<Real>("atol");
  const auto rtol = pkg->Param<Real>("rtol");
  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");

  // gamma_min_dyn: set from the maximum E each current-deposition step when
  // Hybrid/update_gamma_min = 1 (UpdateMomentumBoundary); else the deck value.
  const Real gamma_min = pkg->AllParams().hasKey("gamma_min_dyn")
                             ? pkg->Param<Real>("gamma_min_dyn")
                             : pkg->Param<Real>("gamma_min");
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


  auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  GuidingCenterEquations<Field, true, false> gce(f, c_aw0, ct_a, alpha0);

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

  // Hybrid/feedback_gain = 1: also deposit the parallel-field response of the
  // runaway current (FeedbackGain.hpp), with the same shape and time weights.
  const bool resp_on = Deposit && pkg->AllParams().hasKey("Jre_resp");
  Kokkos::View<Real **> resp;
  if (resp_on) resp = pkg->Param<Kokkos::View<Real **>>("Jre_resp");

  using RegisterWork = Kokkos::Array<Dim5, WorkSize<METHOD>::value>;
  using Work = std::conditional_t<Scratch, Kinetic::ScratchWork<WorkSize<METHOD>::value>,
                                  RegisterWork>;

  auto push = KOKKOS_LAMBDA(const int idx, Work &work_d) {
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
          // Captured here, not first inside if constexpr (not allowed by nvcc
          // for extended lambdas).
          const auto &jre_c = jre;
          const auto &loc_c = cdg.indicator_locator;
          const Real p_RE_c = p_RE;
          const bool resp_on_c = resp_on;
          const auto &resp_c = resp;

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
            if (ret != Ver::Success) {

              pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
              swarm_d.MarkParticleForRemoval(n);
//              if (ret == ParticleVerifyCodes::MomentumCutoff) {
//                pack_swarm_i(b, Kinetic::status(), n) |= Kinetic::DEATH_BY_MOMENTUM;
//              } else if (ret == ParticleVerifyCodes::WallImpact) {
//                pack_swarm_i(b, Kinetic::status(), n) |= Kinetic::DEATH_BY_WALL;
//              }
              break;
            }

            if constexpr (Deposit) {
              if (X[0] > p_RE_c) {
                DepositCurrent(X, t, w, jre_c, dtSA, gce.field, loc_c);
                if (resp_on_c) Kinetic::DepositResponse(X, w, resp_c, dtSA, loc_c);
              }
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
      };
  // Simulation/push_timing = 1: time the push kernel alone (fenced before and
  // after, so only the timing is affected) for /timeseries/push.
  const bool timed = pkg->AllParams().hasKey("push_timing");
  if (timed) Kokkos::fence();
  Kokkos::Timer push_timer;
  const int max_index = pack_swarm_r.GetMaxFlatIndex();
  if constexpr (Scratch) {
    using SW = Kinetic::ScratchWork<WorkSize<METHOD>::value>;
    const int T = pkg->Param<int>("push_team_size");
    const int nchunks = (max_index + T) / T;
    const std::size_t bytes = parthenon::ScratchPad1D<Real>::shmem_size(SW::Count(T));
    // The team's threads share its T slots (TeamThreadRange), so the layout of
    // the work arrays depends on T only, not on the team size. On a GPU the team
    // is T threads, one per slot; on the host Kokkos picks it (the OpenMP team
    // size is bounded by the thread count).
    constexpr bool host = Kokkos::SpaceAccessibility<DevExecSpace,
                                                     Kokkos::HostSpace>::accessible;
    auto policy = host ? parthenon::team_policy(DevExecSpace(), nchunks, Kokkos::AUTO)
                       : parthenon::team_policy(DevExecSpace(), nchunks, T);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes));
    Kokkos::parallel_for(
        "PushParticlesScratch", policy,
        KOKKOS_LAMBDA(parthenon::team_mbr_t member) {
          parthenon::ScratchPad1D<Real> scratch(member.team_scratch(0), SW::Count(T));
          const int first = member.league_rank() * T;
          Kokkos::parallel_for(Kokkos::TeamThreadRange(member, T), [&](const int i) {
            const int idx = first + i;
            if (idx > max_index) return;
            SW work{scratch.data() + i, T};
            push(idx, work);
          });
        });
  } else {
    ParticleLoop(pkg, max_index, KOKKOS_LAMBDA(const int idx) {
      RegisterWork work_d;
      push(idx, work_d);
    });
  }
  if (timed) {
    Kokkos::fence();
    Kinetic::PushTiming::Add(push_timer.seconds(), pack_swarm_r.GetMaxFlatIndex() + 1);
  }

  return TaskStatus::complete;

}

// Dispatch to the compile-time specialized push based on the runtime
// Simulation/integrator param. Branching here (host side) rather than inside the
// kernel means each stepper compiles to its own kernel with its own minimal work
// array — the DP45 kernel's register cost is not paid on RK4 runs and vice versa.
//
// With the Deck param "analytic_field" (set by `avalanche` with
// Simulation/field_source=analytic) the push uses the exact
// AnalyticFieldEvaluator, the wall-free OpenDomain check and no current deposit.
template <int METHOD>
static TaskStatus PushParticlesIn(Mesh *pm, Real t0, Real dt) {
  auto pkg = pm->packages.Get("Deck");
  const Real p_BC = pkg->Param<Real>("p_BC");
  const bool scratch = pkg->AllParams().hasKey("push_scratch");
  if (pkg->AllParams().hasKey("analytic_field")) {
    const auto af = pkg->Param<Kinetic::AnalyticFieldEvaluator>("analytic_field");
    const Kinetic::OpenDomain ver{p_BC};
    return scratch ? PushParticlesImpl<METHOD, false, true>(pm, t0, dt, af, ver)
                   : PushParticlesImpl<METHOD, false, false>(pm, t0, dt, af, ver);
  }
  auto data = pkg->Param<Kinetic::FieldData_t>("FieldData");
  auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  Kinetic::FieldEvaluator f{cdg.hermite_locator, data.hermite_data.view_device()};
  const ParticleVerificator ver{cdg, p_BC};
  return scratch ? PushParticlesImpl<METHOD, true, true>(pm, t0, dt, f, ver)
                 : PushParticlesImpl<METHOD, true, false>(pm, t0, dt, f, ver);
}

TaskStatus PushParticles(Mesh *pm, Real t0, Real dt) {
  const int integrator = pm->packages.Get("Deck")->Param<int>("integrator");
  if (integrator == INTEGRATOR_RK4)
    return PushParticlesIn<INTEGRATOR_RK4>(pm, t0, dt);
  return PushParticlesIn<INTEGRATOR_DOPRI5>(pm, t0, dt);
}

// Compute the two guiding-center adiabatic invariants (canonical toroidal
// momentum p_phi and magnetic moment mu) for every alive particle, store them
// in the particle.p_phi / particle.mu swarm variables (for per-particle .phdf
// output) and log a weight-averaged aggregate to a text file for a quick
// conservation check. Gated by the "EnableComputeConservedQuantities" param.
TaskStatus ComputeConservedQuantities(Mesh *pm, Real t, Real step_start) {

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

  const Real inv_w = (w_total > 0.0) ? 1.0 / w_total : 0.0;
  if (auto *diag = DiagnosticsFile::Get()) {
    // Run time [s]: the step start plus the time within the step.
    const Real time = step_start + t * pkg->Param<Real>("tau_c");
    diag->Append("conservation", {time, p_phi_total * inv_w, mu_total * inv_w, w_total});
  }
  if (Globals::my_rank == 0 && pkg->Param<int>("legacy_dat") == 1) {
    // Legacy text file: the time column is the time within the step [tau_c].
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
  const Real gamma_min = pkg->AllParams().hasKey("gamma_min_dyn")
                             ? pkg->Param<Real>("gamma_min_dyn")
                             : pkg->Param<Real>("gamma_min");
  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  // la_fast (set by avalanche/hybrid with <..>/la_sampler = fast): the
  // secondary (gamma, xi) from LargeAngleCollisionFast, same distribution with
  // efficient proposals; otherwise the reference rejection sampler.
  const bool la_fast = pkg->AllParams().hasKey("la_fast");
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

              if (la_fast)
                LargeAngleCollisionFast(p, xi, w, dtLA, gamma_min, rng_pool);
              else
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

  Log::Debug("RandomRemove: {} alive particles", num_particles);
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

