#include <iostream>
#include <format>
#include <fstream>
#include <sstream>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>

#include <Kokkos_Core.hpp>
#include "AvalancheDriver.h"
#include "pgen.hpp"
#include "kinetic.hpp"

using namespace parthenon;
using namespace parthenon::driver::prelude;


int main(int argc, char *argv[]) {
  User * p_mhd_config;

  ParthenonManager pman;
  auto manager_status = pman.ParthenonInitEnv(argc, argv);
  mhd_PetscInit(&argc, &argv, &p_mhd_config);

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
  pman.app_input->ProblemGenerator = GenerateParticlePoint;

  pman.app_input->UserWorkBeforeLoop = [=](Mesh * pm, ParameterInput * pin, SimTime const & tm) {
    Kinetic::WorkBeforeLoop(pm, p_mhd_config);
    Kinetic::HijackEField(pm);
  };

  pman.ParthenonInitPackagesAndMesh();

  Kinetic::AvalancheDriver driver(pman.pinput.get(), pman.app_input.get(), pman.pmesh.get());
  driver.Execute();
  pman.ParthenonFinalize();
  return 0;
}
