
//========================================================================================
// (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC // for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <globals.hpp>
#include <parthenon_manager.hpp>
#include <iostream>
#include <iomanip>
#include <limits>
#include <format>

#include <Kokkos_Core.hpp>

using namespace parthenon;
using namespace parthenon::driver::prelude;

#include "MHDDriver.h"

#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "tasks/Tasks.h"
#include "pgen.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

MHDDriver::MHDDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh,
    User* p_mhd_config)
      : EvolutionDriver(pin, app_in, pmesh), p_mhd_config(p_mhd_config) {
  	auto pkg = pmesh->packages.Get("Deck");
}

TaskCollection MHDDriver::MakeTaskCollection(BlockList_t &blocks, SimTime tm) {
  TaskCollection tc;
  TaskID none(0);

  auto pkg = pmesh->packages.Get("Deck");
  FieldData_t data = pkg->Param<FieldData_t>("FieldData");
  View3 E_base = pkg->Param<View3>("E_base");
  const Real eta_norm = pkg->Param<Real>("eta_norm");
  const Real En       = pkg->Param<Real>("En");

  const int EnableLargeAngleCollisions = pkg->Param<int>("EnableLargeAngleCollisions");
  const Real tau_c = pkg->Param<Real>("tau_c");


  Real dt = tm.dt / tau_c;

  auto dep = none;
  auto *tl = &tc.AddRegion(1)[0];

  dep = tl->AddTask(dep, CommunicateBJV, data, p_mhd_config, FieldComponents::B);
  dep = tl->AddTask(dep, ComputeBaseElectricField, E_base, data, En, eta_norm, FieldComponents::B);

  AttachFieldCorrector(tl, dep, pmesh, p_mhd_config);

  dep = tl->AddTask(dep, CommunicateBJV, data, p_mhd_config, FieldComponents::Bt);
  dep = tl->AddTask(dep, ComputeBaseElectricField_in_place, data, En, eta_norm, FieldComponents::Bt);
  dep = tl->AddTask(dep, InterpolateTime, data, dt);
  dep = tl->AddTask(dep, InterpolateHermiteBJE, data, FieldComponents::B);
  dep = tl->AddTask(dep, InterpolateHermiteBJE, data, FieldComponents::Bt);

  return tc;
}

void MHDDriver::PreExecute() {
}

void MHDDriver::PostExecute(parthenon::DriverStatus st) {
}

TaskListStatus MHDDriver::Step() {
  PARTHENON_INSTRUMENT
  using DriverUtils::ConstructAndExecuteTaskLists;
  TaskListStatus status = ConstructAndExecuteTaskLists<>(this, tm);
  return status;
}

