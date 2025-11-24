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
using namespace parthenon::driver::prelude;

#include "kinetic/c_wrapper.h"
#include "kinetic/kinetic.hpp"
#include "pgen.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;


int parthenon_init(void ** man, int argc, char *argv[], User* mhd_config) {

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

  Kinetic::InitializeMHDConfig(pman->pinput.get(), mhd_config);

  *man = (void*) pman;

  return 0;
}

int runaway_init(void * man) {
  ParthenonManager* pman = (ParthenonManager*) man;
  // Redefine parthenon defaults
  pman->app_input->ProcessPackages = [](std::unique_ptr<ParameterInput> &pin) {
    Packages_t packages;
    packages.Add(Kinetic::Initialize(pin.get()));
    return packages;
  };
  pman->app_input->ProblemGenerator = GenerateParticleCurrentDensity;

  pman->ParthenonInitPackagesAndMesh();
  Kinetic::ComputeParticleWeights(pman->pmesh.get());
  Kinetic::InitializeDriver(pman);

  return 0;
}

void runaway_finalize(void* man) {
  ParthenonManager* pman = (ParthenonManager*) man;
  pman->ParthenonFinalize();
  delete pman;
}

void runaway_push(void * man) {
  ParthenonManager* pman = (ParthenonManager*) man;
  Kinetic::Push(pman);
}

void runaway_saveState(void * man) {
  ParthenonManager* pman = (ParthenonManager*) man;
  Kinetic::SaveState(pman->pmesh.get());
}

void runaway_restoreState(void * man) {
  ParthenonManager* pman = (ParthenonManager*) man;
  Kinetic::RestoreState(pman->pmesh.get());
}


