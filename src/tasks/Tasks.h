#pragma once
#include <parthenon/package.hpp>

#include "common.hpp"
#include "mhd/mfd_config.h"
using namespace parthenon;
using namespace parthenon::package::prelude;

TaskStatus AdvanceBackgroundFields(User* p_mhd_config, DualView3 Jre);
TaskStatus ResetBackgroundFields(User *p_mhd_config);
TaskStatus ComputeBaseElectricField(View3 E_base, FieldData_t data_V, const Real En, const Real eta_norm, const int component0);
TaskStatus ComputeBaseElectricField_in_place(FieldData_t data_V, const Real En, const Real eta_norm, const int component0);
TaskStatus CommunicateBJV(FieldData_t data, User *p_mhd_config, int component0);
TaskStatus InterpolateTime(FieldData_t data, const Real dt);
TaskStatus InterpolateHermiteBJE(FieldData_t data, const int component0);
TaskStatus InterpolateHermiteE(FieldData_t data, const int component0);
TaskStatus UpdateMomentumBoundary(Mesh* pm, const Real time_0, const Real time_1);

void AttachFieldPredictor(TaskList* tl, TaskID& dep, Mesh* pm, User* p_mhd_config, const Real dtField);
void AttachFieldCorrector(TaskList* tl, TaskID& dep, User* p_mhd_config);
void AttachCurrentCollector(TaskList* tl, TaskID& dep, Mesh* pm, const Real dtCD);

TaskStatus RandomRemove(Mesh* pm);
TaskStatus PushParticles(Mesh *pm, Real t0, Real dt);
TaskStatus CheckScatter(MeshBlock* pmb);
TaskStatus CleanupParticles(MeshBlock* pmb);
TaskStatus AddSecondaries(MeshBlock* pmb, const Real dtLA);
