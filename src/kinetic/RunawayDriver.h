#pragma once
#include <parthenon/driver.hpp>

namespace Kinetic {

using namespace parthenon;
using namespace parthenon::driver::prelude;

class RunawayDriver : public EvolutionDriver {
public:
	Real gamma_min;

  RunawayDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh)
      : EvolutionDriver(pin, app_in, pmesh) {

  	auto pkg = pmesh->packages.Get("Deck");
  	gamma_min = pkg->Param<Real>("gamma_min");

  }

  TaskListStatus Step() {
    PARTHENON_INSTRUMENT
    using DriverUtils::ConstructAndExecuteTaskLists;
    TaskListStatus status = ConstructAndExecuteTaskLists<>(this, tm);
    return status;
  }

  void PostExecute(parthenon::DriverStatus st);

  TaskCollection MakeTaskCollection(BlockList_t &blocks, SimTime tm);

};
}

