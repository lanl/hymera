#include "Tasks.h"
#include "ViewManipulation.h"
#include "kinetic/FeedbackGain.hpp"

void AttachFieldPredictor(TaskList* tl, TaskID& dep,
    Mesh* pm,
    User* p_mhd_config, const Real dtField) {

  auto pkg = pm->packages.Get("Deck");

  FieldData_t data = pkg->Param<FieldData_t>("FieldData");
  View3 E_base = pkg->Param<View3>("E_base");
  View3 Jre =  pkg->Param<View3>("Jre");
  const Real eta_norm = pkg->Param<Real>("eta_norm");
  const Real En       = pkg->Param<Real>("En");
  auto Jre_mhd = pkg->Param<DualView3>("Jre_mhd");

  View3 Jre_seed = pkg->Param<View3>("Jre_seed");

  dep = tl->AddTask(dep, CommunicateBJV, data, p_mhd_config, FieldComponents::B);
  dep = tl->AddTask(dep, ComputeBaseElectricField, E_base, data, En, eta_norm, FieldComponents::B);
  // CommunicateBJV leaves V in the E slot. Put E there too (the same arithmetic as
  // E_base), so that InterpolateTime takes Et = (E_end - E_start) / dt, as in
  // SetupFieldGrid. Before, it took (E_end - V) / dt, about E_end / dt.
  dep = tl->AddTask(dep, ComputeBaseElectricField_in_place, data, En, eta_norm, FieldComponents::B);

  dep = tl->AddTask(dep, AdvanceBackgroundFields, p_mhd_config, Jre_mhd);

  dep = tl->AddTask(dep, CommunicateBJV, data, p_mhd_config, FieldComponents::Bt);
  dep = tl->AddTask(dep, ComputeBaseElectricField_in_place, data, En, eta_norm, FieldComponents::Bt);

  dep = tl->AddTask(dep, InterpolateTime, data, dtField);

  // The first current-deposition interval has no current of its own yet: correct
  // its E with the previous step's runaway current (the time average MHD received),
  // not with zero. Jre itself stays zero for the deposition.
  dep = tl->AddTask(dep, ComputeAdjustedElectricField, data, E_base, Jre_seed, En);
  dep = tl->AddTask(dep, InterpolateHermiteBJE, data, FieldComponents::B);

  dep = tl->AddTask(dep, InterpolateHermiteBJE, data, FieldComponents::Bt);

  dep = tl->AddTask(dep, ResetBackgroundFields, p_mhd_config);
}

void AttachFieldCorrector(TaskList* tl, TaskID& dep,
    Mesh* pm, User* p_mhd_config) {

  auto pkg = pm->packages.Get("Deck");
  auto Jre_mhd = pkg->Param<DualView3>("Jre_mhd");
  // Jre_mhd holds the sum of the nCD current-deposition averages of this MHD step
  // (AttachCurrentCollector); the MHD solver expects their time average over the
  // step, in the same units (eta J in B_0 V_A).
  const int nCD = pkg->Param<int>("nCD");
  dep = tl->AddTask(dep, [Jre_mhd, nCD]() mutable {
    Jre_mhd.sync_device();
    auto j = Jre_mhd.view_device();
    const Real scl = 1.0 / nCD;
    Kokkos::parallel_for("Average Jre_mhd",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {j.extent_int(0), j.extent_int(1)}),
        KOKKOS_LAMBDA(const int i, const int k) {
          j(i, k, 0) *= scl;
          j(i, k, 1) *= scl;
          j(i, k, 2) *= scl;
        });
    Jre_mhd.modify_device();
    return TaskStatus::complete;
  });
  // Keep this step's runaway current, as passed to MHD, for the first deposition
  // interval of the next step (AttachFieldPredictor): the corrector zeroes Jre_mhd.
  auto seed = pkg->Param<View3>("Jre_seed");
  dep = tl->AddTask(dep, [Jre_mhd, seed]() mutable {
    Jre_mhd.sync_device();
    Kokkos::deep_copy(seed, Jre_mhd.view_device());
    return TaskStatus::complete;
  });
  // Keep this step's runaway current, as passed to MHD, for the output diagnostics
  // (Jre_last, set by hybrid with Hybrid/midplane_profile = 1): the corrector zeroes Jre_mhd.
  if (pkg->AllParams().hasKey("Jre_last")) {
    auto last = pkg->Param<View3>("Jre_last");
    dep = tl->AddTask(dep, [Jre_mhd, last]() mutable {
      Jre_mhd.sync_device();
      Kokkos::deep_copy(last, Jre_mhd.view_device());
      return TaskStatus::complete;
    });
  }
  dep = tl->AddTask(dep, AdvanceBackgroundFields, p_mhd_config, Jre_mhd);
}

void AttachCurrentCollector(TaskList* tl, TaskID& dep,
    Mesh* pm, const Real dtCD) {

  auto pkg = pm->packages.Get("Deck");

  FieldData_t data = pkg->Param<FieldData_t>("FieldData");
  View3 E_base = pkg->Param<View3>("E_base");
  View3 Jre =  pkg->Param<View3>("Jre");
  auto Jre_mhd = pkg->Param<DualView3>("Jre_mhd");
  const Real eta_a3VaB0 = pkg->Param<Real>("eta_a3VaB0");
  const Real En       = pkg->Param<Real>("En");

  dep = tl->AddTask(dep, ScaleView,  Jre, eta_a3VaB0 / dtCD);
  dep = tl->AddTask(dep, ReduceView, Jre);
  dep = tl->AddTask(dep, AccumulateView, Jre_mhd, Jre);
  // Hybrid/feedback_gain = 1: the loop gain of this interval's E correction.
  if (pkg->AllParams().hasKey("Jre_resp"))
    dep = tl->AddTask(dep, [pm]() {
      Kinetic::FeedbackGain::Collect(pm);
      return TaskStatus::complete;
    });
  dep = tl->AddTask(dep, ComputeAdjustedElectricField, data, E_base, Jre, En);
  dep = tl->AddTask(dep, InterpolateHermiteE, data, FieldComponents::E);
  dep = tl->AddTask(dep, ZeroView, Jre);
}
