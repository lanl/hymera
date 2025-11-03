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
#include <parthenon_manager.hpp>
using namespace parthenon::driver::prelude;

#include "kinetic/AnalyticField.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/kinetic.hpp"
#include "numerical/rk4.hpp"
#include "pgen.hpp"

//advance particles
TaskStatus ComputeConservationError(Mesh *pm, const Real dt, std::array<Real,2>& err) {
    // get mesh data
    auto md = pm->mesh_data.Get();

    static auto desc_swarm =
    parthenon::MakeSwarmPackDescriptor<
                                     swarm_position::x,
                                     swarm_position::y,
                                     swarm_position::z,
                                       Kinetic::canmom,
                                       Kinetic::pitch,
                                       Kinetic::R,
                                       Kinetic::phi,
                                       Kinetic::Z,
                                       Kinetic::weight>("markers");
    auto pack_swarm = desc_swarm.GetPack(md.get());

    Real tend = 1e-4;

    const AnalyticField af(3.0, 2.1, 2.0, 5.0);
    const Real ct_a = 77666.15939805658;
    const Real c_aw0 = 0.00016080273811573234;
    const Real alpha0 = 0.0;

    const GuidingCenterEquations GCE_system(ct_a, c_aw0, alpha0, af);
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
                X[0] = pack_swarm(b, Kinetic::canmom(), n);
                X[1] = pack_swarm(b, Kinetic::pitch(), n);
                X[2] = pack_swarm(b, Kinetic::R(), n);
                X[3] = pack_swarm(b, Kinetic::phi(), n);
                X[4] = pack_swarm(b, Kinetic::Z(), n);

	            Real p_phi0,mu0;
	            GCE_system.computeConservedQuantities(X, p_phi0,  mu0, 0.0);
	            rk4(X, 0.0, tend, dt, GCE_system);
	            Real p_phiend,muend;
	            GCE_system.computeConservedQuantities(X, p_phiend,  muend, tend);

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

    err[0] = totalMuConservationError / totalWeight;
    err[1] = totalPphiConservationError / totalWeight;

    return TaskStatus::complete;
}

int main(int argc, char* argv[]) {

    parthenon::ParthenonManager pman;

    // Set up kokkos and read pin
    auto manager_status = pman.ParthenonInitEnv(argc, argv);
    if (manager_status == ParthenonStatus::complete) {
      pman.ParthenonFinalize();
      return 0;
    }
    if (manager_status == ParthenonStatus::error) {
      pman.ParthenonFinalize();
      return 1;
    }

     // Redefine parthenon defaults
    pman.app_input->ProcessPackages = [](std::unique_ptr<ParameterInput> &pin) {
      Packages_t packages;
      packages.Add(Kinetic::Initialize(pin.get()));
      return packages;
    };
    pman.app_input->ProblemGenerator = GenerateParticleSquare;

    // call ParthenonInit to set up the mesh
    // scope so that the mesh object, kokkos views, etc, all get cleaned
    // up before kokkos::finalize
    std::array<Real, 2> err1 = {}, err2 = {};
    Real dt1 = 2e-5, dt2 = 1e-5;

    pman.ParthenonInitPackagesAndMesh();
    {
        ComputeConservationError(pman.pmesh.get(), dt1, err1);
        ComputeConservationError(pman.pmesh.get(), dt2, err2);
    }
    pman.ParthenonFinalize();
    // MPI and Kokkos can no longer be used

    Real error_ratio = err1[0]/err2[0];
    Real dt_ratio4 = pow(dt1/dt2,4);
    Real dt_ratio5 = pow(dt1/dt2,5);

    if (error_ratio < dt_ratio4)  {
        std::cerr << "mu convergence order is too low";
        return 1;
    }
    if (error_ratio > dt_ratio5) {
        std::cerr << "mu convergence order is too high";
        return 1;
    }

    error_ratio = err1[1]/err2[1];
    if (error_ratio < dt_ratio4)  {
        std::cerr << "pphi convergence order is too low";
        return 1;
    }
    if (error_ratio > dt_ratio5) {
        std::cerr << "pphi convergence order is too high";
        return 1;
    }

    return (0);
}
