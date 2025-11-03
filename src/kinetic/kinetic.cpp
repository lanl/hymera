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

#include <memory>
#include <numeric>
#include <typeinfo>  //for 'typeid' to work
#include <parthenon/package.hpp>

#include "numerical/rk4.hpp"
#include "../old_hybrid/src/kinetic/GuidingCenterEquations.hpp"
#include "AnalyticField.hpp"

#include "kinetic.hpp"

namespace Kinetic {

auto &GetCoords(std::shared_ptr<MeshBlock> &pmb) { return pmb->coords; }
auto &GetCoords(MeshBlock *pmb) { return pmb->coords; }
auto &GetCoords(Mesh *pm) { return pm->block_list[0]->coords; }

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {
  auto pkg = std::make_shared<StateDescriptor>("Kinetic");

  const Real final_time = pin->GetOrAddReal("Pusher", "final_time", 5.0);
  pkg->AddParam("final_time", finalTime);
  const Real timeStep = pin->GetOrAddReal("Pusher", "timeStep", 1.e-6);
  pkg->AddParam("timeStep", timeStep);

  const Real dtLA    = pin->GetOrAddReal("Collisions", "dtLA", 1.e-5);
  const Real c_vTe   = pin->GetOrAddReal("Collisions", "c_vTe", 355.6647481959091);
  const Real Coulog0 = pin->GetOrAddReal("Collisions", "Coulog0", 6.387682959207553);
  const Real Zeff    = pin->GetOrAddReal("Collisions", "Zeff", 1.0);
  const Real Z0      = pin->GetOrAddReal("Collisions", "Z0", 10.0);
  const Real ZI      = pin->GetOrAddReal("Collisions", "ZI", 1.0);
  const Real fI      = pin->GetOrAddReal("Collisions", "fI", 100.0);
  const Real NeI     = pin->GetOrAddReal("Callisions", "NeI", 9.0);
  pkg->AddParam("dtLA", dtLA);
  pkg->AddParam("c_vTe", c_vTe);
  pkg->AddParam("Coulog0", Coulog0);
  pkg->AddParam("Zeff", Zeff);
  pkg->AddParam("Z0", Z0);
  pkg->AddParam("ZI", ZI);
  pkg->AddParam("fI", fI);
  pkg->AddParam("NeI", NeI);

  const Real q0  = pin->GetOrAddReal("AnalyticField", "q0", 2.1);
  const Real q2  = pin->GetOrAddReal("AnalyticField", "q2", 2.0);
  const Real R_a = pin->GetOrAddReal("AnalyticField", "R_a", 3.0);
  const Real E_0 = pin->GetOrAddReal("AnalyticField", "E_0", 70.0);
  const Real c_aw0  = pin->GetOrAddReal("GCE", "c_aw0", 0.00016080273811573234);
  const Real ct_a   = pin->GetOrAddReal("GCE", "ct_a", 77666.15939805658);
  const Real alpha0 = pin->GetOrAddReal("GCE", "alpha0", 0.0028213404017273245);
  pkg->AddParam("GuidingCenterEquations",
    GuidingCenterEquations(
      AnalyticField(q0, q2, R_a, E_0),
      c_aw0, ct_a, alpha0
    )
  );

  const int npart = pin->GetOrAddInteger("ParticleSeed", "num_particles_per_block", 16);
  pkg->AddParam("num_particles_per_block", npart);
  // Initialize random number generator pool
  int rng_seed = pin->GetOrAddInteger("ParticleSeed", "rng_seed", 1234);
  pkg->AddParam("rng_seed", rng_seed);
  RNGPool rng_pool(rng_seed);
  pkg->AddParam("rng_pool", rng_pool);

  const Real pmin = pin->GetOrAddInteger("ParticleSeed", "pmin", 0.001);
  pkg->AddParam("pmin", pmin);
  const Real pmax = pin->GetOrAddInteger("ParticleSeed", "pmax", 20.0);
  pkg->AddParam("pmax", pmax);
  const Real ximin = pin->GetOrAddInteger("ParticleSeed", "ximin", -1.0);
  pkg->AddParam("ximin", ximin);
  const Real ximax = pin->GetOrAddInteger("ParticleSeed", "ximax",  1.0);
  pkg->AddParam("ximax", ximax);
  const Real Rmin = pin->GetOrAddInteger("ParticleSeed", "Rmin", 2.0);
  pkg->AddParam("Rmin", Rmin);
  const Real Rmax = pin->GetOrAddInteger("ParticleSeed", "Rmax", 4.0);
  pkg->AddParam("Rmax", Rmax);
  const Real Zmin = pin->GetOrAddInteger("ParticleSeed", "Zmin", -2.0);
  pkg->AddParam("Zmin", Zmin);
  const Real Zmax = pin->GetOrAddInteger("ParticleSeed", "Zmax",  2.0);
  pkg->AddParam("Zmax", Zmax);

  Metadata swarm_metadata({Metadata::Provides, Metadata::None});
  pkg->AddSwarm("particles", swarm_metadata);

  Metadata real_swarmvalue_metadata({Metadata::Real});
  Metadata int_swarmvalue_metadata({Metadata::Integer});
  pkg->AddSwarmValue(p::name(), "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(xi::name(),  "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(R::name(),      "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(phi::name(),    "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(Z::name(),      "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(weight::name(), "particles", real_swarmvalue_metadata);

  Metadata int_swarmvalue_metadata({Metadata::Integer});
  pkg->AddSwarmValue(will_scatter::name(), "particles", int_swarmvalue_metadata);
  pkg->AddSwarmValue(secondary_index::name(), "particles", int_swarmvalue_metadata);

  return pkg;
}


//advance particles
TaskStatus March(Mesh *pm) {
    std::cout<<"In March()\n";
    // get mesh data
    auto md = pm->mesh_data.Get();

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
    auto pack_swarm = desc_swarm.GetPack(md.get());

    // Pull out final time from params
    auto pkg = pm->packages.Get("Kinetic");
    auto tend = pkg->Param<Real>("final_time");
    auto dt = pkg->Param<Real>("timeStep");

    auto GCE_system = pkg->Param<GuidingCenterEquations<AnalyticField>>("GuidingCenterEquations");
    printf("dt=%e, tend=%e\n",dt, tend);
    std::cout<<"Device="<<DevExecSpace().name()<<std::endl;

    Kokkos::Timer timer;

    Real totalMuConservationError = 0.0;
    Real totalPphiConservationError = 0.0;
    Real totalWeight = 0.0;

    Kokkos::parallel_reduce(PARTHENON_AUTO_LABEL,
        pack_swarm.GetMaxFlatIndex(),
        // loop over all particles
        KOKKOS_LAMBDA(const int idx, Real& muConservationError, Real& pphiConservationError, Real& weight) {
            // block and particle indices
            auto [b, n] = pack_swarm.GetBlockParticleIndices(idx);
            const auto swarm_d = pack_swarm.GetContext(b);
            if (swarm_d.IsActive(n)) {
                Dim5 X;
                X[0] = pack_swarm(b, Kinetic::p(), n);
                X[1] = pack_swarm(b, Kinetic::xi(), n);
                X[2] = pack_swarm(b, Kinetic::R(), n);
                X[3] = pack_swarm(b, Kinetic::phi(), n);
                X[4] = pack_swarm(b, Kinetic::Z(), n);

	            Real p_phi0,mu0;
	            GCE_system.computeConservedQuantities(X, p_phi0,  mu0, 0.0);
	            rk4(X, 0.0, tend, dt, GCE_system);
	            Real p_phiend,muend;
	            GCE_system.computeConservedQuantities(X, p_phiend,  muend, 0.0);

                pphiConservationError += abs(p_phiend - p_phi0) * pack_swarm(b, Kinetic::weight(), n);
                muConservationError += abs(muend - mu0) * pack_swarm(b, Kinetic::weight(), n);
                weight += pack_swarm(b, Kinetic::weight(), n);
            }
        },
        totalMuConservationError, totalPphiConservationError, totalWeight);

    Real elapsed_time = timer.seconds();
    std::cout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;
    std::cout << "mu error: "   << totalMuConservationError / totalWeight << std::endl;
    std::cout << "pphi error: " << totalPphiConservationError / totalWeight << std::endl;
    std::cout << "total particles: " << totalWeight << std::endl;
    return TaskStatus::complete;
}

} // namespace Kinetic
