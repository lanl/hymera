#include "mhd/mhd.h"
#include "Tasks.h"

TaskStatus AdvanceBackgroundFields(User* p_mhd_config, DualView3 Jre) {

  Jre.sync_host();
  MHD_CHECK(mhd_step(p_mhd_config));
  Kokkos::deep_copy(Jre.view_host(), 0.0);
  return TaskStatus::complete;
}

TaskStatus ResetBackgroundFields(User *p_mhd_config) {
  MHD_CHECK(mhd_resetState(p_mhd_config));
  return TaskStatus::complete;
}
