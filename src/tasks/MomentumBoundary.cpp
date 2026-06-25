#include "Tasks.h"
#include "kinetic/FieldEvaluator.hpp"

TaskStatus UpdateMomentumBoundary(
    Mesh* pm,
    const Real time_0, const Real time_1) {
  auto pkg = pm->packages.Get("Deck");

  auto hpd_R = pkg->Param<ParArray1D<Real>>("Hermite_Field_Plot_data_R");
  auto hpd_Z = pkg->Param<ParArray1D<Real>>("Hermite_Field_Plot_data_Z");
  const int NR_plot = hpd_R.size();
  const int NZ_plot = hpd_Z.size();
  auto data = pkg->Param<FieldData_t>("FieldData");
  auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  FieldEvaluator fev{cdg.hermite_locator, data.hermite_data.view_device()};

  Real maxE = 0.01;
  Kokkos::parallel_reduce(
      "max E",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {NR_plot, NZ_plot}),
      KOKKOS_LAMBDA(const int i, const int j, Real& Epar) {
        const Real R = hpd_R(i);
        const Real Z = hpd_Z(j);

        EvalBE be;
        fev.eval(be, R, Z, time_0);

        Real value = Kokkos::abs(dot_product(be.B, be.E) / be.Bmag);
        if (value > Epar) {
          Epar = value;
        }

        fev.eval(be, R, Z, time_1);

        value = Kokkos::abs(dot_product(be.B, be.E) / be.Bmag);
        if (value > Epar) {
          Epar = value;
        }
      },
      Kokkos::Max<Real>(maxE));

  const Real p_BC_update = momentum_(1.0 + 0.1 / maxE);

  pkg->UpdateParam("p_BC", p_BC_update);
  pkg->UpdateParam("p_RE", p_BC_update);

  return TaskStatus::complete;
}
