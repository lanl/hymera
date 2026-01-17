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

  using DataViewType = decltype(data);
  using JreViewType = Kokkos::View<Real***, DataViewType::array_layout, DataViewType::device_type>;
  JreViewType jre_data;

  Real t_a, t_b;

  EM_Field(int nR_data, int nZ_data, int nphi_data, int nt,
          Real R0, Real Z0, Real dR, Real dZ, Real E_n, Real eta_mu0aVa, Real etaec_a3VaB0, ConfigurationDomainGeometry cdg): FieldInterpolation<2,5>(nR_data, nZ_data, 4, nphi_data, nt, R0, Z0, dR, dZ), E_n(E_n), eta_mu0aVa(eta_mu0aVa), etaec_a3VaB0(etaec_a3VaB0), cdg(cdg), jre_data("jre_data", nR_data, nZ_data, 3), t_a(0.0), t_b(1.0) {};

  KOKKOS_INLINE_FUNCTION
  ERROR_CODE operator() (const Dim5& X, const Real t, Dim3& B, Dim3& curlB, Dim3& dBdR, Dim3& dBdZ, Dim3& E, Dim3& dbdt) const {
    int ii,jj;
    int level = cdg.indicator(X, ii, jj);
    if (level < 1) return ERROR_CODE::WALL_IMPACT;

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
    dbdt = {};
    Dim2 Xloc = {r, z};

    auto sbv = Kokkos::subview(hermite_data, ii, jj, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0, Kokkos::ALL);

    Dim3 V = {}, J_re = {}, dBdt = {};

    Real dt = t_b - t_a;
    Kokkos::Array<Real, 2> tcof = {(t_b - t) / dt, (t - t_a) / dt};
    Kokkos::Array<Real, 2> dtcof = {-1.0 / dt, 1.0 / dt};

    Real sclr = 1.0;
    for (int i = 0; i < sbv.extent(0); ++i) {
      Real sclz = 1.0;
      for (int j = 0; j < sbv.extent(1); ++j) {
        for (int k = 0; k < 3; ++k) {
          Real mon = sclr * sclz;
          J_re[k] += sbv(i, j, 2, k, 0) * mon;
          for (int it = 0; it < sbv.extent(4); ++it) {
            mon = sclr * sclz * tcof[it];

            B[k] += sbv(i, j, 0, k, it) * mon;
            if (i + 1 < sbv.extent(0))
              dBdR[k] += static_cast<Real>(i + 1) *
                         sbv(i + 1, j, 0, k, it) * mon;
            if (j + 1 < sbv.extent(1))
              dBdZ[k] += static_cast<Real>(j + 1) *
                         sbv(i, j + 1, 0, k, it) * mon;
            V[k] += sbv(i, j, 1, k, it) * mon;

            mon = sclr * sclz * dtcof[it];
            dBdt[k] += sbv(i, j, 0, k, it) * mon;

          }
        }
        sclz *= Xloc[1];
      }
      sclr *= Xloc[0];
    }

    const Real &XR = X[2];

    Real BB = 0.0;
    Real BBprime = 0.0;
    for (int i = 0; i < 3; ++i) {
      dBdt[i] /= XR;
      B[i] /= XR;
      BB += B[i]*B[i];
      BBprime += B[i] * dBdt[i];
    }

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

    for (int k = 0; k < 3; ++k) {
      dbdt[k] = (dBdt[k] - B[k] * BBprime / BB)  / sqrt(BB);
    }

    return ERROR_CODE::SUCCESS;
  };

  KOKKOS_INLINE_FUNCTION
  ERROR_CODE operator() (const Dim5& X, const Real t, Dim3& B, Dim3& curlB, Dim3& dBdR, Dim3& dBdZ, Dim3& E, Dim3& J_re, Dim3& V, Dim3& dbdt) const {
    int ii,jj;
    int level = cdg.indicator(X, ii, jj);
    if (level < 1) return ERROR_CODE::WALL_IMPACT;

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
    dbdt = {};
    Dim2 Xloc = {r, z};

    Dim3 dBdt = {};

    auto sbv = Kokkos::subview(hermite_data, ii, jj, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0, Kokkos::ALL);

    Real dt = t_b - t_a;
    Kokkos::Array<Real, 2> tcof = {(t_b - t) / dt, (t - t_a) / dt};
    Kokkos::Array<Real, 2> dtcof = {-1.0 / dt, 1.0 / dt};

    Real sclr = 1.0;
    for (int i = 0; i < sbv.extent(0); ++i) {
      Real sclz = 1.0;
      for (int j = 0; j < sbv.extent(1); ++j) {
        for (int k = 0; k < 3; ++k) {
          Real mon = sclr * sclz;
          J_re[k] += sbv(i, j, 2, k, 0) * mon;
          for (int it = 0; it < sbv.extent(4); ++it) {
            mon = sclr * sclz * tcof[it];
            B[k] += sbv(i, j, 0, k, it) * mon;
            if (i + 1 < sbv.extent(0))
              dBdR[k] += static_cast<Real>(i + 1) *
                         sbv(i + 1, j, 0, k, it) * mon;
            if (j + 1 < sbv.extent(1))
              dBdZ[k] += static_cast<Real>(j + 1) *
                         sbv(i, j + 1, 0, k, it) * mon;
            V[k] += sbv(i, j, 1, k, it) * mon;

            mon = sclr * sclz * dtcof[it];
            dBdt[k] += sbv(i, j, 0, k, it) * mon;
          }
        }
        sclz *= Xloc[1];
      }
      sclr *= Xloc[0];
    }

    const Real &XR = X[2];


    Real BB = 0.0;
    Real BBprime = 0.0;
    for (int i = 0; i < 3; ++i) {
      dBdt[i] /= XR;
      B[i] /= XR;
      BB += B[i]*B[i];
      BBprime += B[i] * dBdt[i];
    }

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

    for (int k = 0; k < 3; ++k) {
      dbdt[k] = (dBdt[k] - B[k] * BBprime / BB)  / sqrt(BB);
    }

    return ERROR_CODE::SUCCESS;
  };

  template<class PsiViewType>
  KOKKOS_INLINE_FUNCTION
  ERROR_CODE evalPsi(Real& val, Dim5 X, const Real t, PsiViewType hermite_data) const {
    Real r =  X[2] - hR0;
    Real z =  X[4] - hZ0;
    int ii = static_cast<int> (floor(r / hR));
    int jj = static_cast<int> (floor(z / hZ));

    r = r/hR - ii - 0.5;
    z = z/hZ - jj - 0.5;

    KOKKOS_ASSERT(std::abs(r) <= 0.5);
    KOKKOS_ASSERT(std::abs(z) <= 0.5);
    KOKKOS_ASSERT(hermite_data.extent(0) > ii && ii >= 0);
    KOKKOS_ASSERT(hermite_data.extent(1) > jj && jj >= 0);

    auto sbv = Kokkos::subview(hermite_data, ii, jj, Kokkos::ALL, Kokkos::ALL, 0, Kokkos::ALL);
    Kokkos::Array<Real, 2> tcof = {(t_b - t) / (t_b - t_a), (t - t_a) / (t_b - t_a)};


    val = 0.0;
    Real sclr = 1.0;
    for (int i = 0; i < sbv.extent(0); ++i) {
      Real sclz = 1.0;
      for (int j = 0; j < sbv.extent(1); ++j) {
        for (int it = 0; it < sbv.extent(3); ++it) {
          Real mon = sclr * sclz * tcof[it];
          val += mon * sbv(i, j, it);
          sclz *= z;
        }
      }
      sclr *= r;
    }

    return ERROR_CODE::SUCCESS;
  }

  auto getJreDataSubview() const {
    return Kokkos::subview(data, Kokkos::ALL, Kokkos::ALL, 2, Kokkos::ALL, 0, 0);
  }

  // TODO: Add method to only reinterpolate current

};

void dumpToHDF5(EM_Field f, int i, const Real t = 0.0);
// TODO: make const, depends on getCorners from hFlux to get const
