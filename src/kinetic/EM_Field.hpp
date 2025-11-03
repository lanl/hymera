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
#pragma once

#include <hFlux/FieldInterpolation.hpp>
#include "ConfigurationDomainGeometry.hpp"

struct EM_Field: public FieldInterpolation<2,5> {
  const Real E_n, eta_mu0aVa, etaec_a3VaB0;
  const ConfigurationDomainGeometry cdg;

  EM_Field(int nR_data, int nZ_data, int nphi_data, int nt,
          Real R0, Real Z0, Real dR, Real dZ, Real E_n, Real eta_mu0aVa, Real etaec_a3VaB0, ConfigurationDomainGeometry cdg): FieldInterpolation<2,5>(nR_data, nZ_data, 4, nphi_data, nt, R0, Z0, dR, dZ), E_n(E_n), eta_mu0aVa(eta_mu0aVa), etaec_a3VaB0(etaec_a3VaB0), cdg(cdg) {};

  KOKKOS_INLINE_FUNCTION
  ERROR_CODE operator() (const Dim5& X, const Real t, Dim3& B, Dim3& curlB, Dim3& dBdR, Dim3& dBdZ, Dim3& E) const {
    int ii,jj;
    int level = cdg.indicator(X, ii, jj);
    if (level != 2) return ERROR_CODE::WALL_IMPACT;

    Real r =  X[2] - hR0;
    Real z =  X[4] - hZ0;
    ii = static_cast<int> (floor(r / hR));
    jj = static_cast<int> (floor(z / hZ));

    r = r/hR - ii - 0.5;
    z = z/hZ - jj - 0.5;

    KOKKOS_ASSERT(std::abs(r) <= 0.5);
    KOKKOS_ASSERT(std::abs(z) <= 0.5);
    KOKKOS_ASSERT(hermite_data.extent(0) > ii && ii >= 0);
    if (!(hermite_data.extent(1) > jj && jj >= 0)) {
        printf("%le %le\n", X[2], X[4]);
        KOKKOS_ASSERT(false);
    }
    B = {};
    dBdR = {};
    dBdZ = {};
    E = {};
    curlB = {};
    Dim2 Xloc = {r, z};

    auto sbv = Kokkos::subview(hermite_data, ii, jj, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0, Kokkos::ALL);

    Dim3 V = {}, J_re = {};

    Real sclr = 1.0;
    for (int i = 0; i < sbv.extent(0); ++i) {
      Real sclz = 1.0;
      for (int j = 0; j < sbv.extent(1); ++j) {
        Real mon = sclr * sclz;
        for (int k = 0; k < 3; ++k) {
          B[k] += sbv(i, j, 0, k, 0) * mon;
          if (i + 1 < sbv.extent(0))
            dBdR[k] += static_cast<Real>(i + 1) *
                       sbv(i + 1, j, 0, k, 0) * mon;
          if (j + 1 < sbv.extent(1))
            dBdZ[k] += static_cast<Real>(j + 1) *
                       sbv(i, j + 1, 0, k, 0) * mon;

          V[k] += sbv(i, j, 1, k, 0) * mon;
          J_re[k] += sbv(i, j, 2, k, 0) * mon;
        }
        sclz *= Xloc[1];
      }
      sclr *= Xloc[0];
    }

    const Real &XR = X[2];


      for (int i = 0; i < 3; ++i) B[i] /= XR;

    curlB[2] = dBdR[1] / XR / hR;

    for (int k = 0; k < 3; ++k) {
      dBdR[k] = (dBdR[k] / hR - B[k]) / XR;
      dBdZ[k] /= XR * hZ;
    }

    curlB[0] = -dBdZ[1];
    curlB[1] = dBdZ[0] - dBdR[2];

    cross_product(B, V, E);
    // E = - VxB + \eta / mu0 (\nabla x B - muJre)
    for (int k = 0; k < 3; ++k)
      E[k] = E_n * (E[k] + eta_mu0aVa * curlB[k] - etaec_a3VaB0 * J_re[k]);

    return ERROR_CODE::SUCCESS;
  };

  KOKKOS_INLINE_FUNCTION
  ERROR_CODE operator() (const Dim5& X, const Real t, Dim3& B, Dim3& curlB, Dim3& dBdR, Dim3& dBdZ, Dim3& E, Dim3& J_re, Dim3 V) const {
    int ii,jj;
    int level = cdg.indicator(X, ii, jj);
    if (level != 2) return ERROR_CODE::WALL_IMPACT;

    Real r =  X[2] - hR0;
    Real z =  X[4] - hZ0;
    ii = static_cast<int> (floor(r / hR));
    jj = static_cast<int> (floor(z / hZ));

    r = r/hR - ii - 0.5;
    z = z/hZ - jj - 0.5;

    KOKKOS_ASSERT(std::abs(r) <= 0.5);
    KOKKOS_ASSERT(std::abs(z) <= 0.5);
    KOKKOS_ASSERT(hermite_data.extent(0) > ii && ii >= 0);
    if (!(hermite_data.extent(1) > jj && jj >= 0)) {
        printf("%le %le\n", X[2], X[4]);
        KOKKOS_ASSERT(false);
    }
    B = {};
    dBdR = {};
    dBdZ = {};
    E = {};
    V = {};
    J_re = {};
    curlB = {};
    Dim2 Xloc = {r, z};

    auto sbv = Kokkos::subview(hermite_data, ii, jj, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0, Kokkos::ALL);

    Real sclr = 1.0;
    for (int i = 0; i < sbv.extent(0); ++i) {
      Real sclz = 1.0;
      for (int j = 0; j < sbv.extent(1); ++j) {
        Real mon = sclr * sclz;
        for (int k = 0; k < 3; ++k) {
          B[k] += sbv(i, j, 0, k, 0) * mon;
          if (i + 1 < sbv.extent(0))
            dBdR[k] += static_cast<Real>(i + 1) *
                       sbv(i + 1, j, 0, k, 0) * mon;
          if (j + 1 < sbv.extent(1))
            dBdZ[k] += static_cast<Real>(j + 1) *
                       sbv(i, j + 1, 0, k, 0) * mon;

          V[k] += sbv(i, j, 1, k, 0) * mon;
          J_re[k] += sbv(i, j, 2, k, 0) * mon;
        }
        sclz *= Xloc[1];
      }
      sclr *= Xloc[0];
    }

    const Real &XR = X[2];


      for (int i = 0; i < 3; ++i) B[i] /= XR;

    curlB[2] = dBdR[1] / XR / hR;

    for (int k = 0; k < 3; ++k) {
      dBdR[k] = (dBdR[k] / hR - B[k]) / XR;
      dBdZ[k] /= XR * hZ;
    }

    curlB[0] = -dBdZ[1];
    curlB[1] = dBdZ[0] - dBdR[2];

    cross_product(B, V, E);
    // E = - VxB + \eta / mu0 (\nabla x B - muJre)
    for (int k = 0; k < 3; ++k)
      E[k] = E_n * (E[k] + eta_mu0aVa * curlB[k] - etaec_a3VaB0 * J_re[k]);

    return ERROR_CODE::SUCCESS;
  };

  auto getJreDataSubview() const {
    return Kokkos::subview(data, Kokkos::ALL, Kokkos::ALL, 2, Kokkos::ALL, 0, 0);
  }

  // TODO: Add method to only reinterpolate current

};

void dumpToHDF5(EM_Field f, int i, const Real t = 0.0);
// TODO: make const, depends on getCorners from hFlux to get const
