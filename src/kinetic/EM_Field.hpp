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

#include <parthenon/package.hpp>
using namespace parthenon;
using namespace parthenon::package::prelude;

enum class fid: std::size_t {
  B = 0,
  E = 1,
  Jre = 2,
  J = 3,
  V = 4,
  GradB = 5,
  Count = 6
};

struct EM_Field {
  static const int m = 2;
  static const int swidth = 5;
  static constexpr std::size_t nfields =
    static_cast<std::size_t>(fid::Count);

  const int nR_data, nZ_data;
  const static int nt = 2; // B, V, J_re
  const static int ndims = 3;

  const Real R0, Z0;
  const Real dR, dZ;

  const int nR_hermite_data, nZ_hermite_data;
  const Real hR0, hZ0;
  const Real hR, hZ;

  const Real E_n;
  const ConfigurationDomainGeometry cdg;

  Kokkos::View<Real**[ndims][nfields][nt], Kokkos::LayoutRight> data;
  Kokkos::View<Real*******, Kokkos::LayoutLeft> hermite_data;

  EM_Field(int nR_data, int nZ_data,
    Real R0, Real Z0, Real dR, Real dZ, Real E_n, ConfigurationDomainGeometry cdg):
    nR_data(nR_data), nZ_data(nZ_data),
    R0(R0), Z0(Z0), dR(dR), dZ(dZ),
    nR_hermite_data((nR_data-1) / (swidth-1) - 1), nZ_hermite_data((nZ_data-1) / (swidth-1) - 1),
    hR0(R0 + (swidth-1)/2*dR), hZ0(Z0 + (swidth-1)/2 *dZ),
    hR(dR * (swidth - 1)), hZ(dZ * (swidth - 1)),
    data("data", nR_data, nZ_data),
    hermite_data("hermite_data", 2*m+3, 2*m+3, nt, ndims, nfields, nR_hermite_data, nZ_hermite_data),
    E_n(E_n), cdg(cdg)
  {};

  void interpolate(std::span<fid> field_indeces, size_t i_t) {
		 for (fid i_f : field_indeces)
		   for (size_t i_d = 0; i_d < ndims; ++i_d) {
         size_t i_ff = static_cast<size_t>(i_f);
         auto hh = Kokkos::subview(hermite_data, Kokkos::ALL, Kokkos::ALL, i_t, i_d, i_ff, Kokkos::ALL, Kokkos::ALL);
		     auto dd = Kokkos::subview(data, Kokkos::ALL, Kokkos::ALL, i_d, i_ff, i_t);
         compute_derivatives_grid<m,swidth>(dd, hh, hR / dR, hZ / dZ);
         interpolate_grid<m>(hh);
		   }
  };

  void cleanDiv(fid i_f, size_t i_t) {
     size_t i_ff = static_cast<size_t>(i_f);
		 auto hh = Kokkos::subview(hermite_data, Kokkos::ALL, Kokkos::ALL, i_t, Kokkos::ALL, i_ff, Kokkos::ALL, Kokkos::ALL);
     cleanDivergence<m>(hh, hR, hZ);
  };

  Kokkos::Array<Real, 4> getCorners() {
    return {hR0, hR0 + nR_hermite_data * hR, hZ0, hZ0 + nZ_hermite_data * hZ};
  };

  KOKKOS_INLINE_FUNCTION
  ErrorCode operator() (const Dim5& X, const Real t, Dim3& B, Dim3& curlB, Dim3& dBdR, Dim3& dBdZ, Dim3& E, Dim3& dbdt) const {

    Real r =  X[2] - hR0;
    Real z =  X[4] - hZ0;
    int ii,jj;
    int level = cdg.indicator(X, ii,jj);
    if (level < 1) return ErrorCode::WallImpact;

    ii = static_cast<int> (floor(r / hR));
    jj = static_cast<int> (floor(z / hZ));

    r = r/hR - ii - 0.5;
    z = z/hZ - jj - 0.5;

    KOKKOS_ASSERT(std::abs(r) <= 0.5);
    KOKKOS_ASSERT(std::abs(z) <= 0.5);
    KOKKOS_ASSERT(hermite_data.extent(5) > ii && ii >= 0);
    KOKKOS_ASSERT(hermite_data.extent(6) > jj && jj >= 0);

    B = {};
    dBdR = {};
    dBdZ = {};
    E = {};
    curlB = {};
    dbdt = {};

    auto sbv = Kokkos::subview(hermite_data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, ii, jj);

    Dim3 J_re = {}, dBdt = {};

    const size_t Pr = hermite_data.extent(0);
    const size_t Pz = hermite_data.extent(1);

    Real sclr = 1.0;
    for (int i = 0; i < Pr; ++i) {
      Real sclz = 1.0;
      for (int j = 0; j < Pz; ++j) {
        for (int i_d = 0; i_d < 3; ++i_d) {
          Real mon = sclr * sclz;
          dBdt[i_d] += sbv(i, j, 1, i_d, static_cast<size_t>(fid::B)  ) * mon;
          J_re[i_d] += sbv(i, j, 0, i_d, static_cast<size_t>(fid::Jre)) * mon;

          B[i_d] += sbv(i, j, 0, i_d, static_cast<size_t>(fid::B)) * mon;
          B[i_d] += sbv(i, j, 1, i_d, static_cast<size_t>(fid::B)) * mon * t;

          if (i + 1 < Pr) {
            dBdR[i_d] += static_cast<Real>(i + 1) * sbv(i + 1, j, 0, i_d, static_cast<size_t>(fid::B)) * mon;
            dBdR[i_d] += static_cast<Real>(i + 1) * sbv(i + 1, j, 1, i_d, static_cast<size_t>(fid::B)) * mon * t;
          }

          if (j + 1 < Pz) {
            dBdZ[i_d] += static_cast<Real>(j + 1) * sbv(i, j + 1, 0, i_d, static_cast<size_t>(fid::B)) * mon;
            dBdZ[i_d] += static_cast<Real>(j + 1) * sbv(i, j + 1, 1, i_d, static_cast<size_t>(fid::B)) * mon * t;
          }

          curlB[i_d] += sbv(i, j, 0, i_d, static_cast<size_t>(fid::J)) * mon;
          curlB[i_d] += sbv(i, j, 1, i_d, static_cast<size_t>(fid::J)) * mon * t;

          E[i_d] += sbv(i, j, 0, i_d, static_cast<size_t>(fid::E)) * mon;
          E[i_d] += sbv(i, j, 1, i_d, static_cast<size_t>(fid::E)) * mon * t;
        }
        sclz *= z;
      }
      sclr *= r;
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


    for (int k = 0; k < 3; ++k) {
      dBdR[k] = (dBdR[k] / hR - B[k]) / XR;
      dBdZ[k] /= XR * hZ;
    }

    for (int k = 0; k < 3; ++k)
      E[k] = E[k] - J_re[k];

    for (int k = 0; k < 3; ++k) {
      dbdt[k] = (dBdt[k] - B[k] * BBprime / BB)  / sqrt(BB);
    }

    return ErrorCode::Success;
  };

  template<class View>
  KOKKOS_INLINE_FUNCTION
  ErrorCode eval_all(const Real R, const Real Z, const Real t, View& vals) const {

    Real r =  R - hR0;
    Real z =  Z - hZ0;

    int ii = static_cast<int> (floor(r / hR));
    int jj = static_cast<int> (floor(z / hZ));

    r = r/hR - ii - 0.5;
    z = z/hZ - jj - 0.5;

    KOKKOS_ASSERT(std::abs(r) <= 0.5);
    KOKKOS_ASSERT(std::abs(z) <= 0.5);
    KOKKOS_ASSERT(hermite_data.extent(0) > ii && ii >= 0);

    auto sbv = Kokkos::subview(hermite_data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, ii, jj);

    for (int i_f = 0; i_f < static_cast<size_t>(fid::Count); ++i_f) {
      for (int k = 0; k < 3; ++k) {
        vals(k, i_f) = 0.0;
        Real sclr = 1.0;
        for (int i = 0; i < sbv.extent(0); ++i) {
          Real sclz = 1.0;
          for (int j = 0; j < sbv.extent(1); ++j) {
            Real mon = sclr * sclz;
            vals(k, i_f) += sbv(i, j, 0, k, i_f) * mon;
            vals(k, i_f) += sbv(i, j, 1, k, i_f) * mon * t;
            sclz *= z;
          }
          sclr *= r;
        }
      }
    }
    return ErrorCode::Success;
  };

  template<class PsiViewType>
  KOKKOS_INLINE_FUNCTION
  ErrorCode evalPsi(const Real R, const Real Z, const Real t, PsiViewType hermite_data, Real& val) const {
    Real r =  R - hR0;
    Real z =  Z - hZ0;

    int ii = static_cast<int> (floor(r / hR));
    int jj = static_cast<int> (floor(z / hZ));

    r = r/hR - ii - 0.5;
    z = z/hZ - jj - 0.5;

    KOKKOS_ASSERT(std::abs(r) <= 0.5);
    KOKKOS_ASSERT(std::abs(z) <= 0.5);
    KOKKOS_ASSERT(hermite_data.extent(0) > ii && ii >= 0);
    KOKKOS_ASSERT(hermite_data.extent(1) > jj && jj >= 0);

    auto sbv = Kokkos::subview(hermite_data, ii, jj, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);

    val = 0.0;
    Real sclr = 1.0;
    for (int i = 0; i < sbv.extent(0); ++i) {
      Real sclz = 1.0;
      for (int j = 0; j < sbv.extent(1); ++j) {
        Real mon = sclr * sclz;
        val += mon * sbv(i, j, 0);
        val += mon * sbv(i, j, 1) * t;
        sclz *= z;
      }
      sclr *= r;
    }

    return ErrorCode::Success;
  }
};

