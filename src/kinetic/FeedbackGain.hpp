//========================================================================================
// (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
//========================================================================================
// Gain of the explicit runaway-current feedback within an MHD step
// (Hybrid/feedback_gain = 1, a diagnostic only).
//
// Each current-deposition interval of length dtCD corrects E with the runaway
// current of the interval before (AttachCurrentCollector), so a perturbation of E
// comes back as -eta dJ_re one interval later. A deposited particle's -v_par
// responds to a parallel field dE (E_c units, time in tau_c) as
//   d(-v_par)/dt = dE (1 + p_perp^2) / gamma^3,
// so the deposited E correction (En * Jre, E_c units) changes at the rate
// kappa dE with kappa = En eta_a3VaB0 sum_i w_i (1 + p_perp,i^2) / gamma_i^3 S_i
// / (2 pi R dR dZ). The loop gain per interval is a = kappa dtCD. With the
// response moment deposited like the current (time weights dtSA, summed over the
// interval), a = En eta_a3VaB0 * deposit in each cell. The lagged update is
// stable for a < 2 and rings for a > 0.34 (hymera_timestepping.tex, Defect 2).
//
// This counts the free acceleration of the deposited particles only: collisions
// reduce the response of slow electrons, and the avalanche adds a term of order
// eta J_re dGamma/dE dtCD, so a is an estimate, not a bound.
//
// Written to /timeseries/feedback at each output: the largest a of any cell and
// interval since the previous output, where it was, and the number of intervals.
#pragma once

#include <iostream>
#include <memory>

#include <Kokkos_DualView.hpp>

#include "util/common.hpp"
#include "kinetic/CurrentDensity.hpp"

namespace parthenon { class Mesh; }

namespace Kinetic {

/// Adds the field response of one particle over dtSA, with the shape weights of
/// DepositCurrent.
template <class View, class Locator>
KOKKOS_INLINE_FUNCTION void DepositResponse(const Dim5 &X, const Real w, const View &resp,
                                            const Real dtSA, const Locator &locator) {
  const Real p = X[0];
  const Real xi = X[1];
  const Real R = X[2];
  const Real Z = X[4];
  const Real g = gamma_(p);
  const Real pperp2 = p * p * (1.0 - xi * xi);
  const Real contribution =
      (1.0 + pperp2) / (g * g * g) / R / locator.dR / locator.dZ / 2.0 / M_PI * dtSA * w;

  int i, j;
  Real xiR, xiZ;
  locator.locate(R, Z, i, j, xiR, xiZ);
  for (int ii = -1; ii < 2; ++ii) {
    if (i + ii >= 0 and i + ii < resp.extent(0)) {
      const Real wr = S2(Kokkos::abs(xiR - static_cast<Real>(ii)));
      for (int jj = -1; jj < 2; ++jj) {
        if (j + jj >= 0 and j + jj < resp.extent(1)) {
          const Real wz = S2(Kokkos::abs(xiZ - static_cast<Real>(jj)));
          Kokkos::atomic_add(&resp(i + ii, j + jj), contribution * wr * wz);
        }
      }
    }
  }
}

struct FeedbackGain {
  /// Sums the interval's response over ranks, records its largest a and zeroes
  /// the deposit. A task of AttachCurrentCollector.
  static int Collect(parthenon::Mesh *pm);
  /// Appends (t, a_max, R, Z, intervals) to /timeseries/feedback and resets.
  /// No-op unless the diagnostic is on.
  static void Flush(double time);
  static bool enabled;
};

} // namespace Kinetic
