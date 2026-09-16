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
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <globals.hpp>
#include <parthenon_manager.hpp>
#include <iostream>
#include <limits>
#include <format>

#include <Kokkos_Core.hpp>

using namespace parthenon;
using namespace parthenon::driver::prelude;

#include "kinetic/SingleParticleDriver.h"

#include "kinetic/kinetic.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/FieldEvaluator.hpp"
#include "kinetic/ParticleVerificator.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/rk45.hpp"
#include "kinetic/rk4.hpp"
#include "tasks/Tasks.h"

namespace Kinetic {

SingleParticleDriver::SingleParticleDriver(ParameterInput *pin, ApplicationInput *app_in,
                                           Mesh *pmesh,
                                           std::shared_ptr<SingleParticleState> state)
    : EvolutionDriver(pin, app_in, pmesh), state_(state) {}

// Advance the single particle by one full parthenon step: one collisionless
// guiding-center push (t -> t+dt, no small-angle collisions inside the
// integrator) followed by one small-angle collision kick over the whole step.
// Then recompute (p_phi, mu). Runs a 1-element device kernel so the exact same
// device-side code (StepParticle, small-angle operator, conserved-quantity
// formula) as the swarm-based PushParticles path is exercised.
TaskListStatus SingleParticleDriver::Step() {
  PARTHENON_INSTRUMENT
  auto pkg = pmesh->packages.Get("Deck");

  const auto h    = pkg->Param<Real>("hRK");
  const auto atol = pkg->Param<Real>("atol");
  const auto rtol = pkg->Param<Real>("rtol");
  const auto p_BC = pkg->Param<Real>("p_BC");
  const auto c_aw0  = pkg->Param<Real>("c_aw0");
  const auto ct_a   = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");
  const auto tau_c  = pkg->Param<Real>("tau_c");
  const int integrator = pkg->Param<int>("integrator");
  const int EnableSmallAngleCollisions = pkg->Param<int>("EnableSmallAngleCollisions");
  (void)tau_c; // this driver works directly in tau_c normalization

  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  auto sa = pkg->Param<
      SmallAngleCollision<PartialScreening, EnergyScattering, ModifiedCouLog>>(
      "SmallAngleCollision");

  auto data = pkg->Param<Kinetic::FieldData_t>("FieldData");
  auto cdg  = pkg->Param<ConfigurationDomainGeometry>("CDG");
  Kinetic::FieldEvaluator f{cdg.hermite_locator, data.hermite_data.view_device()};
  GuidingCenterEquations<decltype(f), true, false> gce(f, c_aw0, ct_a, alpha0);
  ParticleVerificator ver{cdg, p_BC};

  // Psi (poloidal flux) evaluator on the separate psi grid (as in
  // ComputeConservedQuantities).
  Evaluator psi_ev{data.hermite_locator};
  auto psi_view = data.psi_data.view_device();

  // This driver works directly in tau_c normalization: the parthenon timestep
  // (dt_force in the input) IS a single collisionless step in units of tau_c.
  // One push + one small-angle kick are taken over this whole step -- no
  // MHD/Alfven-time conversion (unlike the hybrid driver) and no inner subcycle.
  const Real dt = tm.dt;
  const Real tstart = state_->t;

  // Device-resident scratch: X (5), plus outputs p_phi, mu, and the RK work array.
  Kokkos::View<Real *> Xd("Xd", 5);
  Kokkos::View<Real *> outd("outd", 2); // [0]=p_phi, [1]=mu
  auto Xh = Kokkos::create_mirror_view(Xd);
  for (int i = 0; i < 5; ++i) Xh(i) = state_->X[i];
  Kokkos::deep_copy(Xd, Xh);

  const int method = integrator;

  Kokkos::parallel_for(
      "single_particle_step", 1, KOKKOS_LAMBDA(const int) {
        Dim5 X;
        for (int i = 0; i < 5; ++i) X[i] = Xd(i);

        // One collisionless push over the full step.
        int ret = ParticleVerificator::Success;
        if (method == INTEGRATOR_RK4) {
          Kokkos::Array<Dim5, 5> work;
          ret = solve_rk4_fixed(gce, ver, X, tstart, tstart + dt, h, work);
        } else {
          Kokkos::Array<Dim5, 10> work;
          ret = solve_rk45(gce, ver, X, tstart, tstart + dt, rtol, atol, h, 1e-9,
                     std::numeric_limits<int>::max(), work);
        }
        if (ret != ParticleVerificator::Success)
          Kokkos::printf("single_particle: push verify failed (code %d) at "
                         "p=%g xi=%g R=%g Z=%g -- particle left valid domain\n",
                         ret, X[0], X[1], X[2], X[4]);

        // One small-angle collision kick over the full step.
        if (EnableSmallAngleCollisions == 1) sa(X[0], X[1], dt, rng_pool);

        // Conserved quantities at the end of the step.
        Real Psi = 0.0;
        psi_ev.evalPsi(Psi, X[2], X[4], psi_view);
        Real p_phi = 0.0, mu = 0.0;
        gce.computeConservedQuantities(p_phi, mu, X[0], X[1], X[2], X[4],
                                       tstart + dt, Psi);

        for (int i = 0; i < 5; ++i) Xd(i) = X[i];
        outd(0) = p_phi;
        outd(1) = mu;
      });
  Kokkos::fence();

  Kokkos::deep_copy(Xh, Xd);
  auto outh = Kokkos::create_mirror_view(outd);
  Kokkos::deep_copy(outh, outd);

  for (int i = 0; i < 5; ++i) state_->X[i] = Xh(i);
  state_->p_phi = outh(0);
  state_->mu = outh(1);
  state_->t = tstart + dt;

  return TaskListStatus::complete;
}

// Evaluate (p_phi, mu) for state->X at its current time without stepping. Runs a
// 1-element device kernel, mirroring the conserved-quantity block of Step().
void ComputeConservedForState(Mesh *pm, std::shared_ptr<SingleParticleState> state) {
  auto pkg = pm->packages.Get("Deck");
  const auto c_aw0  = pkg->Param<Real>("c_aw0");
  const auto ct_a   = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");

  auto data = pkg->Param<Kinetic::FieldData_t>("FieldData");
  auto cdg  = pkg->Param<ConfigurationDomainGeometry>("CDG");
  Kinetic::FieldEvaluator f{cdg.hermite_locator, data.hermite_data.view_device()};
  GuidingCenterEquations<decltype(f), true, false> gce(f, c_aw0, ct_a, alpha0);
  Evaluator psi_ev{data.hermite_locator};
  auto psi_view = data.psi_data.view_device();

  Kokkos::View<Real *> Xd("Xd", 5);
  Kokkos::View<Real *> outd("outd", 2);
  auto Xh = Kokkos::create_mirror_view(Xd);
  for (int i = 0; i < 5; ++i) Xh(i) = state->X[i];
  Kokkos::deep_copy(Xd, Xh);
  const Real t = state->t;

  Kokkos::parallel_for(
      "single_particle_conserved", 1, KOKKOS_LAMBDA(const int) {
        Real Psi = 0.0;
        psi_ev.evalPsi(Psi, Xd(2), Xd(4), psi_view);
        Real p_phi = 0.0, mu = 0.0;
        gce.computeConservedQuantities(p_phi, mu, Xd(0), Xd(1), Xd(2), Xd(4), t, Psi);
        outd(0) = p_phi;
        outd(1) = mu;
      });
  Kokkos::fence();

  auto outh = Kokkos::create_mirror_view(outd);
  Kokkos::deep_copy(outh, outd);
  state->p_phi = outh(0);
  state->mu = outh(1);
}

}
