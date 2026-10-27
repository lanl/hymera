#include <iostream>
#include <format>
#include <fstream>
#include <sstream>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>

#include <Kokkos_Core.hpp>
#include "HybridDriver.h"
#include "pgen.hpp"
#include "kinetic.hpp"
#include "AppSetup.hpp"
#include "Log.hpp"
#include "tasks/Tasks.h"
#include "mhd.h"

using namespace parthenon;
using namespace parthenon::driver::prelude;


int main(int argc, char *argv[]) {
  User * p_mhd_config;

  ParthenonManager pman;
  if (auto rc = Kinetic::InitApp(pman, argc, argv, &p_mhd_config)) return *rc;

  // <Hybrid> options, read only here so other executables' parameter dumps are
  // unchanged; both off by default.
  //   midplane_profile = 1   write /profiles/midplane_c<cycle> at each output
  //   jre_map = 1            write /profiles/jre_c<cycle>, the 2-D runaway current passed to MHD
  //   update_gamma_min = 1   secondary threshold gamma_min = 2 (gamma_BC - 1) + 1
  //                          following the momentum boundary set by max |E_par|
  // Read only when the block is in the deck, so runs without it keep their
  // parameter dump (and the regression fingerprint) unchanged.
  const bool hb = pman.pinput->DoesBlockExist("Hybrid");
  const int midplane = hb ? pman.pinput->GetOrAddInteger("Hybrid", "midplane_profile", 0) : 0;
  // Read only when present, so that the parameter dump of other runs is unchanged.
  const int jre_map = hb && pman.pinput->DoesParameterExist("Hybrid", "jre_map")
                          ? pman.pinput->GetInteger("Hybrid", "jre_map") : 0;
  const int dyn_gmin = hb ? pman.pinput->GetOrAddInteger("Hybrid", "update_gamma_min", 0) : 0;
  const bool la_fast = hb && pman.pinput->GetOrAddString("Hybrid", "la_sampler", "reference") == "fast";
  const Real gmin_floor =
      dyn_gmin == 1 ? pman.pinput->GetOrAddReal("Hybrid", "gamma_min_floor", 1.02) : 1.02;
  // ParticleSeed/stratified = 1: seed on flux-surface bands (GenerateParticleCurrentDensity).
  // Read only when present, so that the parameter dump of other runs is unchanged.
  const int stratified = pman.pinput->DoesParameterExist("ParticleSeed", "stratified")
                             ? pman.pinput->GetInteger("ParticleSeed", "stratified") : 0;
  const int seed_bands = stratified == 1
                             ? pman.pinput->GetOrAddInteger("ParticleSeed", "bands", 200) : 0;
  const int seed_rng_base = pman.pinput->GetOrAddInteger("ParticleSeed", "rng_seed", 1234);
  pman.app_input->ProcessPackages = [=](std::unique_ptr<ParameterInput> &pin) {
    Packages_t packages;
    auto pkg = Kinetic::Initialize(pin.get(), p_mhd_config);
    if (dyn_gmin == 1) {
      pkg->AddParam("gamma_min_dyn", pkg->Param<Real>("gamma_min"), Params::Mutability::Mutable);
      pkg->AddParam("gamma_min_floor", gmin_floor);
    }
    if (la_fast) pkg->AddParam("la_fast", 1);
    pkg->AddParam("seed_stratified", stratified);
    if (stratified == 1) {
      pkg->AddParam("seed_bands", seed_bands);
      pkg->AddParam("seed_rng_base", seed_rng_base);
    }
    if (midplane == 1 || jre_map == 1) {
      const auto j = pkg->Param<DualView3>("Jre_mhd");
      pkg->AddParam("Jre_last", View3("Jre_last", j.extent(0), j.extent(1)));
    }
    packages.Add(pkg);
    return packages;
  };
  pman.app_input->ProblemGenerator = GenerateParticleCurrentDensity;
  pman.app_input->UserWorkBeforeLoop = [=](Mesh * pm, ParameterInput * pin, SimTime const & tm) {
    Kinetic::WorkBeforeLoop(pm, p_mhd_config);
    // With a dynamic threshold, set it (and p_BC, p_RE) from the initial field
    // before the first push rather than after the first current deposition.
    if (dyn_gmin == 1) {
      // t = 0: the start of the step, the time origin of the in-step fields E + Et t.
      UpdateMomentumBoundary(pm, 0.0, 0.0);
      Kinetic::Log::Info("Momentum boundary from max |E_par|: p_BC = {:g}, gamma_min = {:g}",
                         pm->packages.Get("Deck")->Param<Real>("p_BC"),
                         pm->packages.Get("Deck")->Param<Real>("gamma_min_dyn"));
    }
  };

  pman.app_input->UserWorkBeforeRestartOutput = [=](Mesh * pm, ParameterInput * pin, SimTime const & tm, OutputParameters* op) {
    Kinetic::WorkBeforeRestartOutput(pm, pin, op, p_mhd_config);
  };

  pman.app_input->UserMeshWorkBeforeOutput = [=](Mesh * pm, ParameterInput * pin, SimTime const & tm) {
    Kinetic::WorkBeforeOutput(pm, pin, tm, p_mhd_config);
    if (midplane == 1) Kinetic::WriteMidplaneProfile(pm, tm);
    if (jre_map == 1) Kinetic::WriteJreMap(pm, tm);
  };

  Kinetic::AttachDiagnostics(pman);
  pman.ParthenonInitPackagesAndMesh();

  // The interpolated Hermite field grid is fully populated by Kinetic::Initialize
  // (MHD init + interpolation), which runs during ParthenonInitPackagesAndMesh --
  // no driver step is needed. Dump it here for the `profile` executable and exit
  // before any time integration.
  if (pman.pinput->GetOrAddInteger("Simulation", "dump_raw_fields", 0) == 1) {
    const std::string fname = pman.pinput->GetOrAddString("Simulation", "raw_fields_file", "fields.h5");
    Kinetic::SaveRawFieldData(pman.pmesh.get(), fname.c_str());
    Kinetic::FinishDiagnostics();
    MHD_CHECK(mhd_destroy(p_mhd_config));
    pman.ParthenonFinalize();
    return 0;
  }

  Kinetic::ComputeParticleWeights(pman.pmesh.get());

  Kinetic::HybridDriver driver(pman.pinput.get(), pman.app_input.get(), pman.pmesh.get(),
      p_mhd_config);
  driver.Execute();
  Kinetic::FinishDiagnostics();

  MHD_CHECK(mhd_destroy(p_mhd_config));
  pman.ParthenonFinalize();
  return 0;
}
