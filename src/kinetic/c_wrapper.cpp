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
#include <parthenon_manager.hpp>

#include "kinetic/AnalyticField.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/EM_Field.hpp"
#include "pgen.hpp"

constexpr bool PartialScreening = true;
constexpr bool EnergyScattering = false;
constexpr bool ModifiedCouLog = true;

using namespace parthenon;
using namespace parthenon::driver::prelude;


void runaway_init(void ** man, int argc, char *argv[]) {

  ParthenonManager* pman = new ParthenonManager();

  // Set up kokkos and read pin
  auto manager_status = pman->ParthenonInitEnv(argc, argv);
  if (manager_status == ParthenonStatus::complete) {
    pman->ParthenonFinalize();
    return 0;
  }
  if (manager_status == ParthenonStatus::error) {
    pman->ParthenonFinalize();
    return 1;
  }

  // Redefine parthenon defaults
  pman->app_input->ProcessPackages = [](std::unique_ptr<ParameterInput> &pin) {
    Packages_t packages;
    packages.Add(Initialize(pin.get()));
    return packages;
  };
  pman->app_input->ProblemGenerator = GenerateParticleCurrentDensity;

  pman->ParthenonInitPackagesAndMesh();
  ComputeParticleWeights(pman->pmesh.get());

  *man = (void*) pman;
}

void runaway_finalize(void* man) {
  ParthenonManager* pman = (ParthenonManager*) man;
  pman->ParthenonFinalize();
  delete pman;
}

void runaway_push(void * man) {
  ParthenonManager* pman = (ParthenonManager*) man;
  RunawayDriver driver(pman->pinput.get(), pman->app_input.get(), pman->pmesh.get());
	auto driver_status = driver.Execute();
}

void runaway_reset(void * man) {
}

