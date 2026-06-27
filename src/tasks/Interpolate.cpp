#include <hFlux/Interpolator.hpp>
#include "fid.h"
#include "mhd.h"
#include "Tasks.h"

TaskStatus CommunicateBJV(FieldData_t data, User *p_mhd_config, int component0) {

  auto field_h = data.data.view_host();

  auto B_h = Kokkos::subview(field_h, Kokkos::ALL, Kokkos::ALL,
      Kokkos::make_pair(component0 + FieldComponents::B, component0 + FieldComponents::B+3));
  auto J_h = Kokkos::subview(field_h, Kokkos::ALL, Kokkos::ALL,
      Kokkos::make_pair(component0 + FieldComponents::J, component0 + FieldComponents::J+3));
  // Use E memory to store V
  auto V_h = Kokkos::subview(field_h, Kokkos::ALL, Kokkos::ALL,
      Kokkos::make_pair(component0 + FieldComponents::E, component0 + FieldComponents::E+3));

  mhd_getF(p_mhd_config, fid_B, wrap_view(B_h));
  mhd_getF(p_mhd_config, fid_V, wrap_view(V_h));
  mhd_getF(p_mhd_config, fid_J, wrap_view(J_h));
  data.data.modify_host();

  return TaskStatus::complete;
}

TaskStatus InterpolateTime(FieldData_t data, const Real dt) {

  auto field_d = data.data.view_device();

  const int NR = field_d.extent_int(0);
  const int NZ = field_d.extent_int(1);


  data.data.sync_device();
  Kokkos::parallel_for(
    "Interpolate time",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {NR, NZ, 9}),

    KOKKOS_LAMBDA(const int i, const int j, const int k) {
        field_d(i, j, FieldComponents::Bt + k) = (
            field_d(i, j, FieldComponents::Bt + k)
            - field_d(i, j, FieldComponents::B + k)
            ) / dt;
    });
  data.data.modify_device();

  return TaskStatus::complete;
}

TaskStatus InterpolateHermiteBJE(FieldData_t data, const int component0) {

  Interpolator<FIELD_SMOOTHNESS, FIELD_FD_STENSIL> itrp;

  data.data.sync_device();
  itrp.interpolateRange<9>(data.fd_locator,
      data.hermite_locator,
      data.data.view_device(),
      data.hermite_data.view_device(), component0);
  itrp.cleanDivergence(data.hermite_locator,
      data.hermite_data.view_device(), component0);
  data.hermite_data.modify_device();

  return TaskStatus::complete;
}

TaskStatus InterpolateHermiteE(FieldData_t data, const int component0) {

  Interpolator<FIELD_SMOOTHNESS, FIELD_FD_STENSIL> itrp;

  data.hermite_data.sync_device();
  itrp.interpolateRange<3>(data.fd_locator,
      data.hermite_locator,
      data.data.view_device(),
      data.hermite_data.view_device(), component0);
  data.hermite_data.modify_device();

  return TaskStatus::complete;
}
