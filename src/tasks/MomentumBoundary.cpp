#include "Tasks.h"

#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/FieldEvaluator.hpp"

TaskStatus UpdateMomentumBoundary(
    Mesh* pm,
    const Real time_0, const Real time_1) {
  auto pkg = pm->packages.Get("Deck");

  auto data = pkg->Param<FieldData_t>("FieldData");
  auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  FieldEvaluator fev{cdg.hermite_locator, data.hermite_data.view_device()};

  // Uniform NR_maxE x NZ_maxE sample points just inside the Hermite grid. The
  // expressions match the ones that built the field-plot grid this search used
  // to read, so the sampled points are bit-identical at the default sizes.
  const int NR_maxE = pkg->Param<int>("NR_maxE");
  const int NZ_maxE = pkg->Param<int>("NZ_maxE");
  const auto &hl = data.hermite_locator;
  const Real R0_s = hl.R0 + 1e-10;
  const Real Z0_s = hl.Z0 + 1e-10;
  const Real dR_s = (hl.nR * hl.dR - 2e-10) / static_cast<Real>(NR_maxE);
  const Real dZ_s = (hl.nZ * hl.dZ - 2e-10) / static_cast<Real>(NZ_maxE);

  // The largest |E_par| where particles can be: regions >= 1 (inside the
  // separatrix and the scrape-off layer), as in ParticleVerificator. Before, the
  // search covered the whole grid, and once the plasma field fell the maximum came
  // from the vacuum vessel (tag -1), which no particle sees: p_BC fell with it.
  Real maxE = 0.01;
  Kokkos::parallel_reduce(
      "max E",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {NR_maxE, NZ_maxE}),
      KOKKOS_LAMBDA(const int i, const int j, Real& Epar) {
        const Real R = R0_s + i * dR_s;
        const Real Z = Z0_s + j * dZ_s;

        int iR, iZ;
        Region region;
        cdg.locate_region(R, Z, iR, iZ, region);
        if (region < 1) return;

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

  // gamma_BC - 1 = factor E_c / max |E_par|. The critical energy is about
  // E_c / (2 E), so factor = 0.1 keeps particles down to 0.2 of it: a larger factor
  // is cheaper (fewer slow electrons) with less margin below critical.
  // Hybrid/p_BC_factor sets it; read only when present (Deck param p_BC_factor).
  const Real factor = pkg->AllParams().hasKey("p_BC_factor")
                          ? pkg->Param<Real>("p_BC_factor") : 0.1;
  const Real p_BC_update = momentum_(1.0 + factor / maxE);

  pkg->UpdateParam("p_BC", p_BC_update);
  pkg->UpdateParam("p_RE", p_BC_update);
  // Hybrid/update_gamma_min = 1: the secondary threshold follows the boundary,
  // gamma_min = 2 (gamma_BC - 1) + 1, the relation of the deck default, so that
  // secondaries are not born just above p_BC, but not below Hybrid/gamma_min_floor
  // (default 1.02). The reference large-angle sampler's acceptance falls as
  // (gamma_min - 1) and stalls below about 1.02; the fast sampler (la_sampler =
  // fast) does not, and the floor then only cuts the secondary source.
  if (pkg->AllParams().hasKey("gamma_min_dyn")) {
    const Real floor = pkg->Param<Real>("gamma_min_floor");
    pkg->UpdateParam("gamma_min_dyn", Kokkos::max(floor, 2.0 * (factor / maxE) + 1.0));
  }

  return TaskStatus::complete;
}
