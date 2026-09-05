//========================================================================================
// (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
// for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================

#ifndef RK45_HPP
#define RK45_HPP

#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <type_traits>
#include "util/common.hpp"

/**
 * A single Dormand–Prince 4(5) step.
 *
 * @param f    The ODE right-hand side, f(t, y, dydx), which computes dydx =
 * f(t,y).
 * @param t    Current independent variable.
 * @param y    Current dependent variable vector.
 * @param h    Step size (on output, can be adapted by calling code if step is
 * rejected).
 * @param yout On return, yout holds y + one accepted step of size h.
 * @param yerr On return, yerr is the estimate of the local error in each
 * component.
 */
template <class System>
KOKKOS_INLINE_FUNCTION void rk45_step(
    const System &f, const double t, const typename System::value_type &y, double h,
    typename System::value_type &yout, typename System::value_type &yerr,
    Kokkos::Array<typename System::value_type, 10> &k) {
  // Coefficients for Dormand–Prince 4(5):
  constexpr double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0,
                      c5 = 8.0 / 9.0, c6 = 1.0,
                      c7 = 1.0; // same as c6, used for clarity

  constexpr double a21 = 1.0 / 5.0, a31 = 3.0 / 40.0, a32 = 9.0 / 40.0,
                      a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0,
                      a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0,
                      a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0,
                      a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0,
                      a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0,
                      a65 = -5103.0 / 18656.0, a71 = 35.0 / 384.0, a72 = 0.0,
                      a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0,
                      a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;

  // 5th-order weights (b) and embedded 4th-order weights (b*) for error
  // estimate
  constexpr double b1 = 35.0 / 384.0, b2 = 0.0, b3 = 500.0 / 1113.0,
                      b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0,
                      b6 = 11.0 / 84.0,
                      b7 = 0.0; // 5th order

  constexpr double b1s = 5179.0 / 57600.0, b2s = 0.0, b3s = 7571.0 / 16695.0,
                      b4s = 393.0 / 640.0, b5s = -92097.0 / 339200.0,
                      b6s = 187.0 / 2100.0,
                      b7s = 1.0 / 40.0; // 4th order

  const size_t n = y.size();

  // 1) k[1] = f(t, y)
  f(t, y, k[1]);

  // 2) k[2] = f(t + c2*h, y + h*(a21*k[1]))
  for (int i = 0; i < n; ++i)
    k[0][i] = y[i] + h * a21 * k[1][i];
  f(t + c2 * h, k[0], k[2]);

  // 3) k[3]
  for (int i = 0; i < n; ++i)
    k[0][i] = y[i] + h * (a31 * k[1][i] + a32 * k[2][i]);
  f(t + c3 * h, k[0], k[3]);

  // 4) k[4]
  for (int i = 0; i < n; ++i)
    k[0][i] = y[i] + h * (a41 * k[1][i] + a42 * k[2][i] + a43 * k[3][i]);
  f(t + c4 * h, k[0], k[4]);

  // 5) k[5]
  for (int i = 0; i < n; ++i)
    k[0][i] = y[i] + h * (a51 * k[1][i] + a52 * k[2][i] + a53 * k[3][i] +
                          a54 * k[4][i]);
  f(t + c5 * h, k[0], k[5]);

  // 6) k[6]
  for (int i = 0; i < n; ++i)
    k[0][i] = y[i] + h * (a61 * k[1][i] + a62 * k[2][i] + a63 * k[3][i] +
                          a64 * k[4][i] + a65 * k[5][i]);
  f(t + c6 * h, k[0], k[6]);

  // 7) k[7]
  for (int i = 0; i < n; ++i)
    k[0][i] = y[i] + h * (a71 * k[1][i] + a72 * k[2][i] + a73 * k[3][i] +
                          a74 * k[4][i] + a75 * k[5][i] + a76 * k[6][i]);
  f(t + c7 * h, k[0], k[7]);

  // Compute 5th-order update and 4th-order (embedded) update
  for (int i = 0; i < n; ++i) {
    double dy5 = b1 * k[1][i] + b2 * k[2][i] + b3 * k[3][i] + b4 * k[4][i] +
                 b5 * k[5][i] + b6 * k[6][i] + b7 * k[7][i];
    double dy4 = b1s * k[1][i] + b2s * k[2][i] + b3s * k[3][i] + b4s * k[4][i] +
                 b5s * k[5][i] + b6s * k[6][i] + b7s * k[7][i];

    yout[i] = y[i] + h * dy5;  // 5th order solution
    yerr[i] = h * (dy5 - dy4); // local error estimate in each component
  }
}


/**
 * Adaptive solver that uses rk45_step() to integrate from t0 to tf.
 *
 * @param f      ODE RHS: f(t, y, dydx).
 * @param y      On input, initial y vector; on return, final y vector at tf (or
 * last step).
 * @param t0     Initial t.
 * @param tf   Final t.
 * @param rtol   Relative tolerance(s); scalar or length = n.
 * @param atol   Absolute tolerance(s); scalar or length = n.
 * @param itol   0 => rtol,atol are scalars; 1 => they are vectors.
 * @param h      Suggested initial step size (will be adapted).
 * @param hmin   Minimum step size allowed.
 * @param nmax   Max number of steps.
 */

template <class System, class Verificator, typename T, bool VERBOSE = false>
KOKKOS_INLINE_FUNCTION typename Verificator::ResultCode_t
solve_rk45(const System &f,
             const Verificator &v, typename System::value_type &y, const double t0,
             const double tf, const T rtol, const T atol, double h, const double hmin, const int nmax,
             Kokkos::Array<typename System::value_type, 10> &work) {
  const size_t n = work[0].size();

  double t = t0;
  auto &ytemp = work[8];
  auto &yerr = work[9];

  // Helper lambdas to fetch tolerances

  auto get_rtol = [&](int i) -> double {
    if constexpr (std::is_same_v<typename System::value_type, T>)
      return rtol[i];
    else if constexpr (std::is_same_v<T, double>)
      return rtol;
    else
      return 0.0;
  };
  auto get_atol = [&](int i) -> double {
    if constexpr (std::is_same_v<typename System::value_type, T>)
      return atol[i];
    else if constexpr (std::is_same_v<T, double>)
      return atol;
    else
      return 0.0;
  };
 //  nmax = 1000;
  auto st = v.verify(y);
  if (st != Verificator::Success) return st;

  for (int stepCount = 0; stepCount < nmax; ++stepCount) {

    if ((t + 1.0e-14) >= tf) {
      // We are close enough or past the end
      break;
    }
    if (Kokkos::fabs(h) < hmin) {
      h = hmin;
    }
    if (t + h > tf) {
      // Don’t overshoot tf
      h = tf - t;
    }

    // 1) Take a trial step
    rk45_step(f, t, y, h, ytemp, yerr, work);

    // 2) Estimate error norm
    double err = 0.0;
    for (int i = 0; i < n; ++i) {
      double sc = get_atol(i) + get_rtol(i) * Kokkos::fabs(ytemp[i]);
      double e = yerr[i] / sc;
      err += e * e;
    }
    err = Kokkos::sqrt(err / n);

    // 3) Compare with 1.0 to accept/reject
    const double safety = 0.9;
    const double p = 0.2; // exponent = 1/(order+1) = 1/5
    if (err <= 1.0) {
      // Accept the step
      t += h;
      y = ytemp;

      st = v.verify(y);
      if (st != Verificator::Success) return st;


      // Step-size update
      double factor = safety * Kokkos::pow(err + 1.0e-10, -p);
      factor = Kokkos::min(5.0, Kokkos::max(0.2, factor)); // clamp
      h *= factor;
      if constexpr (VERBOSE) {
        printf("%20.14le ", t);
        for (auto x : y)
          printf("%20.14le ", x);
        printf("%3c\n", 'A');
      }
    } else {
      // Reject, reduce h and retry
      double factor = safety * Kokkos::pow(err, -p);
      factor = Kokkos::min(5.0, Kokkos::max(0.2, factor));
      h *= factor;
      if constexpr (VERBOSE) {
        printf("%20.14le ", t);
        for (auto x : y)
          printf("%20.14le ", x);
        printf("%3c\n", 'R');
      }
    }
  }
  return Verificator::Success;
}


template <class System, class Verificator, bool VERBOSE = false>
KOKKOS_INLINE_FUNCTION typename Verificator::ResultCode_t solve_rk45_fixed(
    const System &f,
    const Verificator& v, typename System::value_type &y, const double t0, const double tf,
    double h,
    Kokkos::Array<typename System::value_type, 10> &work)
{

  const size_t n = y.size();
  double t = t0;
  int nmax = static_cast<int>(Kokkos::ceil((tf - t0) / h));
  auto &ytemp = work[8];
  auto &yerr = work[9];

  for (int stepCount = 0; stepCount < nmax; ++stepCount) {
    auto st = v.verify(y);
    if (st != Verificator::Success) return st;

    if (stepCount == nmax - 1)
      h = tf - t;

    rk45_step(f, t, y, h, ytemp, yerr, work);
    y = ytemp;

    st = v.verify(y);
    if (st != Verificator::Success) return st;

    // Accept the step
    t += h;

    if constexpr (VERBOSE) {
      printf("%20.14le ", t);
      for (auto x : y)
        printf("%20.14le ", x);
      printf("A\n");
    }
  }

  return Verificator::Success;
}

#endif // RK45_HPP
