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
//
// `single_particle` executable: advance ONE guiding-center particle with one
// collisionless push + one small-angle collision kick per parthenon timestep,
// tracking the conserved quantities (p_phi, mu). No swarm, no MHD -- the field
// grid is loaded from an HDF5 dump (fields.h5) exactly as the conservation and
// profile drivers do. The particle is a plain host variable seeded from the
// <SingleParticle> input block. A text file (one row per parthenon output
// cadence) records t, p, xi, R, phi, Z, p_phi, mu.

#include <iostream>
#include <fstream>
#include <format>
#include <memory>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>
#include <globals.hpp>

#include <Kokkos_Core.hpp>
#include "SingleParticleDriver.h"
#include "kinetic.hpp"

using namespace parthenon;
using namespace parthenon::driver::prelude;

int main(int argc, char *argv[]) {

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

  // Shared particle state, owned here so the driver and output hook see the
  // same object. Seeded from the <SingleParticle> input block.
  auto state = std::make_shared<Kinetic::SingleParticleState>();

  pman.app_input->ProcessPackages = [=](std::unique_ptr<ParameterInput> &pin) {
    Packages_t packages;
    packages.Add(Kinetic::Initialize(pin.get(), NULL));
    return packages;
  };

  // Particle generation does nothing: the single particle is a host variable,
  // not a swarm particle.
  pman.app_input->ProblemGenerator = [](MeshBlock *pmb, ParameterInput *pin) {};

  const std::string output_file =
      pman.pinput->GetOrAddString("SingleParticle", "output_file", "single_particle.dat");

  pman.app_input->UserWorkBeforeLoop =
      [=](Mesh *pm, ParameterInput *pin, SimTime const &tm) {
        // Seed the particle state from input.
        state->X[0] = pin->GetReal("SingleParticle", "p");
        state->X[1] = pin->GetReal("SingleParticle", "xi");
        state->X[2] = pin->GetReal("SingleParticle", "R");
        state->X[3] = pin->GetOrAddReal("SingleParticle", "phi", 0.0);
        state->X[4] = pin->GetReal("SingleParticle", "Z");
        state->t = 0.0;

        // Load the pre-interpolated field grid (no MHD).
        const std::string load_fields =
            pin->GetOrAddString("Simulation", "load_fields", "fields.h5");
        Kinetic::LoadRawFieldData(pm, load_fields.c_str());
        Kinetic::HijackEField(pm);

        // Populate the initial conserved quantities and start the output file.
        Kinetic::ComputeConservedForState(pm, state);
        if (Globals::my_rank == 0) {
          std::ofstream ofs(output_file, std::ios::trunc);
          ofs << "# t p xi R phi Z p_phi mu" << std::endl;
        }
      };

  // Append one row per parthenon output trigger (cycle 0, each output dt, final).
  pman.app_input->UserMeshWorkBeforeOutput =
      [=](Mesh *pm, ParameterInput *pin, SimTime const &tm) {
        if (Globals::my_rank != 0) return;
        std::ofstream ofs(output_file, std::ios::app);
        ofs << std::format(
                   "{:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} "
                   "{:20.14e} {:20.14e}",
                   state->t, state->X[0], state->X[1], state->X[2], state->X[3],
                   state->X[4], state->p_phi, state->mu)
            << std::endl;
      };

  pman.ParthenonInitPackagesAndMesh();

  Kinetic::SingleParticleDriver driver(pman.pinput.get(), pman.app_input.get(),
                                       pman.pmesh.get(), state);
  driver.Execute();
  pman.ParthenonFinalize();
  return 0;
}
