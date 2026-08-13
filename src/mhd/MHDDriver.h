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
#pragma once
#include <parthenon/driver.hpp>
#include "mhd.h"

using namespace parthenon;
using namespace parthenon::driver::prelude;

class MHDDriver : public EvolutionDriver {
public:

  MHDDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh, User* p_mhd_config);

  void PreExecute();
  void PostExecute(parthenon::DriverStatus st);

  TaskCollection MakeTaskCollection(BlockList_t &blocks, SimTime tm);
  TaskListStatus Step();

private:
  User* p_mhd_config;
};

