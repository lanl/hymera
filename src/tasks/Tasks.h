#pragma once
#include <parthenon/package.hpp>

#include "kinetic.hpp"
#include "util/common.hpp"

#include "mfd_config.h"
using namespace parthenon;
using namespace parthenon::package::prelude;
using namespace Kinetic;

TaskStatus AdvanceBackgroundFields(User* p_mhd_config, DualView3 Jre);
TaskStatus ResetBackgroundFields(User *p_mhd_config);
TaskStatus ComputeBaseElectricField(View3 E_base, FieldData_t data_V, const Real En, const Real eta_norm, const int component0);
TaskStatus ComputeBaseElectricField_in_place(FieldData_t data_V, const Real En, const Real eta_norm, const int component0);
TaskStatus ComputeAdjustedElectricField(FieldData_t data, View3 E_base, View3 Jre, const Real En);
TaskStatus CommunicateBJV(FieldData_t data, User *p_mhd_config, int component0);
TaskStatus InterpolateTime(FieldData_t data, const Real dt);
TaskStatus InterpolateHermiteBJE(FieldData_t data, const int component0);
TaskStatus InterpolateHermiteE(FieldData_t data, const int component0);
TaskStatus UpdateMomentumBoundary(Mesh* pm, const Real time_0, const Real time_1);
TaskStatus ZeroView(View3 view);

void AttachFieldPredictor(TaskList* tl, TaskID& dep, Mesh* pm, User* p_mhd_config, const Real dtField);
void AttachFieldCorrector(TaskList* tl, TaskID& dep, Mesh* pm, User* p_mhd_config);
void AttachCurrentCollector(TaskList* tl, TaskID& dep, Mesh* pm, const Real dtCD);

TaskStatus RandomRemove(Mesh* pm);
TaskStatus PushParticles(Mesh *pm, Real t0, Real dt);
TaskStatus ComputeConservedQuantities(Mesh *pm, Real t);
TaskStatus CheckScatter(MeshBlock* pmb);
TaskStatus CleanupParticles(MeshBlock* pmb);
TaskStatus AddSecondaries(MeshBlock* pmb, const Real dtLA);
