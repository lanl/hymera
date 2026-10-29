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

// Start-up steps shared by the kinetic executables' main functions.
#pragma once

#include <memory>
#include <optional>
#include <string>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>
#include <globals.hpp>

#include "kinetic/kinetic.hpp"
#include "kinetic/DiagnosticsFile.hpp"
#include "kinetic/PushTiming.hpp"
#include "kinetic/FeedbackGain.hpp"
#include "mhd/mhd.h"

namespace Kinetic {

/// Initialises Parthenon and, if mhd_context is given, PETSc (in that order).
/// Returns an exit code when the run is already over (e.g. `-h`, or a bad
/// command line), or nothing when main should go on to set up the run.
inline std::optional<int> InitApp(parthenon::ParthenonManager &pman, int &argc,
                                  char **&argv, User **mhd_context = nullptr) {
  using parthenon::ParthenonStatus;
  auto manager_status = pman.ParthenonInitEnv(argc, argv);
  if (mhd_context != nullptr) MHD_CHECK(mhd_PetscInit(&argc, &argv, mhd_context));

  if (manager_status == ParthenonStatus::complete) {
    pman.ParthenonFinalize();
    return 0;
  }
  if (manager_status == ParthenonStatus::error) {
    pman.ParthenonFinalize();
    return 1;
  }
  return std::nullopt;
}

/// Registers the kinetic "Deck" package. With mhd_context == nullptr the MHD
/// field fill is skipped and the field grid is expected from LoadFieldsFromInput.
inline void RegisterKineticPackage(parthenon::ParthenonManager &pman, User *mhd_context) {
  pman.app_input->ProcessPackages = [=](std::unique_ptr<parthenon::ParameterInput> &pin) {
    parthenon::Packages_t packages;
    packages.Add(Initialize(pin.get(), mhd_context));
    return packages;
  };
}

/// Loads the pre-interpolated field grid named by Simulation/load_fields
/// (default fields.h5), as written by `hybrid` with Simulation/dump_raw_fields=1.
inline void LoadFieldsFromInput(parthenon::Mesh *pm, parthenon::ParameterInput *pin) {
  const std::string load_fields = pin->GetOrAddString("Simulation", "load_fields", "fields.h5");
  LoadRawFieldData(pm, load_fields.c_str());
}

/// Connects the run's diagnostics file to the Parthenon hooks. Call after main
/// has set its own hooks and before ParthenonInitPackagesAndMesh: each hook main
/// set still runs, first, unchanged.
///   UserWorkBeforeLoop        on restart, drop rows after the restart time
///   UserMeshWorkBeforeOutput  flush, so the file is complete at every output
/// Parthenon writes the final output after UserWorkAfterLoop, so the file is
/// closed by FinishDiagnostics, after the driver has returned.
inline void AttachDiagnostics(parthenon::ParthenonManager &pman) {
  using parthenon::Mesh;
  using parthenon::ParameterInput;
  using parthenon::SimTime;
  auto &in = *pman.app_input;

  auto before_loop = in.UserWorkBeforeLoop;
  in.UserWorkBeforeLoop = [before_loop](Mesh *pm, ParameterInput *pin, SimTime &tm) {
    if (before_loop) before_loop(pm, pin, tm);
    if (auto *d = DiagnosticsFile::Get(); d && parthenon::Globals::is_restart)
      d->TruncateAfter(tm.time);
  };

  auto before_output = in.UserMeshWorkBeforeOutput;
  in.UserMeshWorkBeforeOutput = [before_output](Mesh *pm, ParameterInput *pin,
                                                SimTime const &tm) {
    if (before_output) before_output(pm, pin, tm);
    PushTiming::Flush(tm.time);
    FeedbackGain::Flush(tm.time);
    FlushDiagnostics();
  };

}

/// Writes the remaining rows and closes the diagnostics file. Call after
/// driver.Execute().
inline void FinishDiagnostics() { DiagnosticsFile::Close(); }

} // namespace Kinetic
