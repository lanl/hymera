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

#include <iostream>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>

#include <Kokkos_Core.hpp>

#include <hFlux/dopri.hpp>

using namespace parthenon;

#include "kinetic/kinetic.hpp"
#include "kinetic/RunawayDriver.h"
#include "pgen.hpp"
using namespace Kinetic;

int main(int argc, char *argv[]) {

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
    packages.Add(Initialize(pin.get()));
    return packages;
  };
  pman.app_input->ProblemGenerator = GenerateParticleCurrentDensity;

  // call ParthenonInit to set up the mesh
  // scope so that the mesh object, kokkos views, etc, all get cleaned
  // up before kokkos::finalize

  pman.ParthenonInitPackagesAndMesh();
  {
    ComputeParticleWeights(pman.pmesh.get());
    RunawayDriver driver(pman.pinput.get(), pman.app_input.get(), pman.pmesh.get());
		auto driver_status = driver.Execute();

    // Separate dumping, current deposit for mhd and actuall avalanching phase.
  }
  pman.ParthenonFinalize();
  return (0);
}
