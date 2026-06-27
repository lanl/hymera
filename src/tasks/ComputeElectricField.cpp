#include "Tasks.h"

TaskStatus ComputeBaseElectricField(View3 E_base, FieldData_t data_V, const Real En, const Real eta_norm,
    const int component0) {

  auto field_d = data_V.data.view_device();
  auto B_d = Kokkos::subview(field_d, Kokkos::ALL, Kokkos::ALL,
      Kokkos::make_pair(component0 + FieldComponents::B, component0 + FieldComponents::B+3));
  auto J_d = Kokkos::subview(field_d, Kokkos::ALL, Kokkos::ALL,
      Kokkos::make_pair(component0 + FieldComponents::J, component0 + FieldComponents::J+3));
  auto V_d = Kokkos::subview(field_d, Kokkos::ALL, Kokkos::ALL,
      Kokkos::make_pair(component0 + FieldComponents::E, component0 + FieldComponents::E+3));

  const int NR = field_d.extent_int(0);
  const int NZ = field_d.extent_int(1);

  data_V.data.sync_device();
  Kokkos::parallel_for(
    "Pack MHD fields",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {NR, NZ}),
    KOKKOS_LAMBDA(const int i, const int j) {
      const Real V0 = V_d(i, j, 0);
      const Real V1 = V_d(i, j, 1);
      const Real V2 = V_d(i, j, 2);

      const Real B0 = B_d(i, j, 0);
      const Real B1 = B_d(i, j, 1);
      const Real B2 = B_d(i, j, 2);

      const Real J0 = J_d(i, j, 0);
      const Real J1 = J_d(i, j, 1);
      const Real J2 = J_d(i, j, 2);

      const Real VxB0 = V1 * B2 - V2 * B1;
      const Real VxB1 = V2 * B0 - V0 * B2;
      const Real VxB2 = V0 * B1 - V1 * B0;

      E_base(i, j, 0) = En * (-VxB0 + eta_norm * J0);
      E_base(i, j, 1) = En * (-VxB1 + eta_norm * J1);
      E_base(i, j, 2) = En * (-VxB2 + eta_norm * J2);
    });
  return TaskStatus::complete;
}

TaskStatus ComputeBaseElectricField_in_place(FieldData_t data_V, const Real En, const Real eta_norm,
    const int component0) {

  auto field_d = data_V.data.view_device();
  auto B_d = Kokkos::subview(field_d, Kokkos::ALL, Kokkos::ALL,
      Kokkos::make_pair(component0 + FieldComponents::B, component0 + FieldComponents::B+3));
  auto J_d = Kokkos::subview(field_d, Kokkos::ALL, Kokkos::ALL,
      Kokkos::make_pair(component0 + FieldComponents::J, component0 + FieldComponents::J+3));
  auto V_d = Kokkos::subview(field_d, Kokkos::ALL, Kokkos::ALL,
      Kokkos::make_pair(component0 + FieldComponents::E, component0 + FieldComponents::E+3));

  const int NR = field_d.extent_int(0);
  const int NZ = field_d.extent_int(1);

  data_V.data.sync_device();
  Kokkos::parallel_for(
    "Pack MHD fields",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {NR, NZ}),
    KOKKOS_LAMBDA(const int i, const int j) {
      const Real V0 = V_d(i, j, 0);
      const Real V1 = V_d(i, j, 1);
      const Real V2 = V_d(i, j, 2);

      const Real B0 = B_d(i, j, 0);
      const Real B1 = B_d(i, j, 1);
      const Real B2 = B_d(i, j, 2);

      const Real J0 = J_d(i, j, 0);
      const Real J1 = J_d(i, j, 1);
      const Real J2 = J_d(i, j, 2);

      const Real VxB0 = V1 * B2 - V2 * B1;
      const Real VxB1 = V2 * B0 - V0 * B2;
      const Real VxB2 = V0 * B1 - V1 * B0;

      V_d(i, j, 0) = En * (-VxB0 + eta_norm * J0);
      V_d(i, j, 1) = En * (-VxB1 + eta_norm * J1);
      V_d(i, j, 2) = En * (-VxB2 + eta_norm * J2);
    });
  data_V.data.modify_device();

  return TaskStatus::complete;
}

TaskStatus ComputeAdjustedElectricField(FieldData_t data, View3 E_base, View3 Jre, const Real En) {

  auto field_d = data.data.view_device();

  const int NR = field_d.extent_int(0);
  const int NZ = field_d.extent_int(1);

  Kokkos::parallel_for(
    "Pack MHD fields",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {NR, NZ}),
    KOKKOS_LAMBDA(const int i, const int j) {
      field_d(i,j,FieldComponents::E + 0) = E_base(i,j,0) - En * Jre(i,j,0);
      field_d(i,j,FieldComponents::E + 1) = E_base(i,j,1) - En * Jre(i,j,1);
      field_d(i,j,FieldComponents::E + 2) = E_base(i,j,2) - En * Jre(i,j,2);
    });
  data.data.modify_device();

  return TaskStatus::complete;
}
