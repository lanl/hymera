#include <iostream>
#include <format>
#include <fstream>
#include <sstream>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>

#include <Kokkos_Core.hpp>
#include <EM_Field.hpp>
#include "HybridDriver.h"
#include "pgen.hpp"
#include "kinetic.hpp"
#include "mhd.h"

using namespace parthenon;
using namespace parthenon::driver::prelude;


int main(int argc, char *argv[]) {
  User * p_mhd_config;
  mhd_PetscInit(argc, argv, &p_mhd_config);

  ParthenonManager pman;
  auto manager_status = pman.ParthenonInitEnv(argc, argv);

  if (manager_status == ParthenonStatus::complete) {
    pman.ParthenonFinalize();
    return 0;
  }
  if (manager_status == ParthenonStatus::error) {
    pman.ParthenonFinalize();
    return 1;
  }

  pman.app_input->ProcessPackages = [=](std::unique_ptr<ParameterInput> &pin) {
    Packages_t packages;
    packages.Add(Kinetic::Initialize(pin.get(), p_mhd_config));
    return packages;
  };
  pman.app_input->ProblemGenerator = GenerateParticleCurrentDensity;
  pman.app_input->UserWorkBeforeLoop = [=](Mesh * pm, ParameterInput * pin, SimTime const & tm) {
    Kinetic::WorkBeforeLoop(pm, p_mhd_config);
  };

  pman.app_input->UserWorkBeforeRestartOutput = [=](Mesh * pm, ParameterInput * pin, SimTime const & tm, OutputParameters* op) {
    Kinetic::WorkBeforeRestartOutput(pm, pin, op, p_mhd_config);
  };

  pman.app_input->UserMeshWorkBeforeOutput = [=](Mesh * pm, ParameterInput * pin, SimTime const & tm) {
    Kinetic::WorkBeforeOutput(pm, pin, tm, p_mhd_config);
  };

  pman.ParthenonInitPackagesAndMesh();
  Kinetic::ComputeParticleWeights(pman.pmesh.get());

  Kinetic::HybridDriver driver(pman.pinput.get(), pman.app_input.get(), pman.pmesh.get(),
      p_mhd_config);
  driver.Execute();

  mhd_destroy(p_mhd_config);
  pman.ParthenonFinalize();
  return 0;
}
