#pragma once
#include <parthenon/driver.hpp>

namespace Kinetic {

using namespace parthenon;
using namespace parthenon::driver::prelude;

class RunawayDriver : public EvolutionDriver {
public:

  RunawayDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh)
      : EvolutionDriver(pin, app_in, pmesh) {
  }

  TaskListStatus Step() {
    PARTHENON_INSTRUMENT
    using DriverUtils::ConstructAndExecuteTaskLists;
    TaskListStatus status = ConstructAndExecuteTaskLists<>(this, tm);
    return status;
  }

  void PreExecute();
  void PostExecute(parthenon::DriverStatus st);

  TaskCollection MakeTaskCollection(BlockList_t &blocks, SimTime tm);

};
}

