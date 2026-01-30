#include "mhd/mfd_config.h"
#include "kinetic/kinetic.hpp"
#include "kinetic/c_wrapper.h"
#include "GuidingCenterDriver.h"

int main(int argc, char *argv[]) {

  void* man;
	User mfd_config;

	parthenon_init(&man, argc, argv, &mfd_config);
	Kinetic::LoadRawFieldData(&mfd_config, "raw_field.h5");
	runaway_init(man, &mfd_config);

	parthenon::ParthenonManager * pman = (parthenon::ParthenonManager*) man;
  Kinetic::GuidingCenterDriver driver(pman->pinput.get(), pman->app_input.get(), pman->pmesh.get());
  driver.tm.tlim = 100.0;
  driver.Execute();

	return 0;
}

