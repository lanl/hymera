#include <iostream>
#include <format>
#include <fstream>
#include <sstream>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>

#include <Kokkos_Core.hpp>
#include <EM_Field.hpp>
#include "pgen.hpp"
#include "kinetic.hpp"
#include "mhd.h"

using namespace parthenon;
using namespace parthenon::driver::prelude;

void writefile(int i, Real* B, Real* V, int N) {
  if (Globals::my_rank != 0) return;
  std::string filename = std::format("field_{:05d}.txt", i);
  std::ofstream outfile(filename);
  for (int j = 0; j < N; ++j) {
    outfile << std::format("{:.15e} {:.15e}\n", B[j], V[j]);
  }
  outfile.close();


  const Real Rmin = 1.525; ///< Minimum R [-]
  const Real Rmax = 4.975; ///< Maximum R [-]
  const Real Zmin = -2.975;///<  Minimum Z [-]
  const Real Zmax =  2.975; ///< Maximum Z [-]

  const int NR    = 100;
  const int NZ    = 200;

  const Real dR = (Rmax - Rmin) / (Real) NR;
  const Real dZ = (Zmax - Zmin) / (Real) NZ;

  const Real RminCellCenter = Rmin + .5 * dR;
  const Real ZminCellCenter = Zmin + .5 * dZ;

  const std::string configurationdomain_file = "../inputs/AxisSymmetricGeometry.dat";
  ConfigurationDomainGeometry::IndicatorViewType indicator("indicator", NR, NZ);
  std::ifstream ifs(configurationdomain_file);
  auto indicator_h = Kokkos::create_mirror_view(indicator);
  for (int i = 0; i < NR; ++i) {
    for (int j = 0; j < NZ; ++j) {
      ifs >> indicator_h(i,j);
    }
  }

  Kokkos::deep_copy(indicator, indicator_h);
  ConfigurationDomainGeometry cdg(RminCellCenter, ZminCellCenter, dR, dZ, -3, indicator);
  EM_Field f(NR, NZ, 1, 2, RminCellCenter, ZminCellCenter, dR, dZ, 0.0, 0.0, 0.0, cdg);
  auto field_data = f.getDataRef();
  using Host = Kokkos::HostSpace;
  using Unmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
  Kokkos::View<Real******, Kokkos::LayoutLeft, Host> field_data_h("field_data_h", NR, NZ, 4, 3, 1, 2);
  Kokkos::View<Real***, Kokkos::LayoutLeft, Host, Unmanaged> Bh(B, NR, NZ, 3);
  auto Bsub = Kokkos::subview(field_data_h, Kokkos::ALL, Kokkos::ALL, 0, Kokkos::ALL, 0, 0);
  Kokkos::deep_copy(Bsub, Bh);
  auto Vsub = Kokkos::subview(field_data_h, Kokkos::ALL, Kokkos::ALL, 1, Kokkos::ALL, 0, 0);
  Kokkos::View<Real***, Kokkos::LayoutLeft, Host, Unmanaged> Vh(B, NR, NZ, 3);
  Kokkos::deep_copy(Vsub, Vh);
  Kokkos::deep_copy(field_data, field_data_h);
  f.interpolate();
  dumpToHDF5(f, i, 0.0);

}

int main(int argc, char *argv[]) {
  User mhd_config;
  mhd_PetscInit(argc, argv);

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

 User * p_mhd_config = &mhd_config;
  pman.app_input->ProcessPackages = [=](std::unique_ptr<ParameterInput> &pin) {
    Packages_t packages;
    packages.Add(Kinetic::Initialize(pin.get(), p_mhd_config));
    return packages;
  };
  pman.app_input->ProblemGenerator = GenerateParticleCurrentDensity;
  pman.ParthenonInitPackagesAndMesh();
  Kinetic::ComputeParticleWeights(pman->pmesh.get());

  HybridDriver driver(pman.pinput.get(), pman.app_input.get(), pman.pmesh.get());
  driber.Execute();


  mhd_destroy(&mhd_config);
  pman.ParthenonFinalize();
  return 0;
}
