//========================================================================================
// (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
//========================================================================================
#include "kinetic/FeedbackGain.hpp"

#include <mpi.h>

#include <parthenon/package.hpp>

#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/DiagnosticsFile.hpp"

namespace Kinetic {

namespace {
double g_amax = 0.0, g_R = 0.0, g_Z = 0.0;
std::int64_t g_intervals = 0;
bool g_declared = false;
} // namespace

bool FeedbackGain::enabled = false;

int FeedbackGain::Collect(parthenon::Mesh *pm) {
  auto pkg = pm->packages.Get("Deck");
  auto resp = pkg->Param<Kokkos::View<Real **>>("Jre_resp");
  const Real scale = pkg->Param<Real>("En") * pkg->Param<Real>("eta_a3VaB0");
  const auto loc = pkg->Param<ConfigurationDomainGeometry>("CDG").indicator_locator;

  Kokkos::fence();
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), resp);
  MPI_Allreduce(MPI_IN_PLACE, h.data(), h.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (std::size_t i = 0; i < h.extent(0); ++i)
    for (std::size_t j = 0; j < h.extent(1); ++j) {
      const double a = scale * h(i, j);
      if (a > g_amax) {
        g_amax = a;
        g_R = loc.R0 + (i + 0.5) * loc.dR;
        g_Z = loc.Z0 + (j + 0.5) * loc.dZ;
      }
    }
  g_intervals += 1;
  Kokkos::deep_copy(resp, 0.0);
  return 0;
}

void FeedbackGain::Flush(double time) {
  if (!enabled) return;
  if (auto *diag = DiagnosticsFile::Get()) {
    if (!g_declared) {
      diag->DeclareSeries("feedback", {"t", "a_max", "R", "Z", "intervals"},
                          {"s", "1", "a", "a", "1"});
      g_declared = true;
    }
    diag->Append("feedback", {time, g_amax, g_R, g_Z, static_cast<double>(g_intervals)});
  }
  g_amax = 0.0;
  g_R = g_Z = 0.0;
  g_intervals = 0;
}

} // namespace Kinetic
