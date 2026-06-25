#include "kinetic.hpp"

TaskStatus AdvanceBackgroundFields(User* p_mhd_config, DualView3 Jre) {

  Jre.host_sync();
  mhd_step(p_mhd_config);
  return TaskStatus::complete;
}

TaskStatus ResetBackgroundFields(User *p_mhd_config) {
  mhd_resetState(p_mhd_config);
  return TaskStatus::complete;
}
