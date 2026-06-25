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
#include <hFlux/hFlux.hpp>
#include "FieldComponents.h"

namespace Kinetic {

// Output field structs by the use case
struct EvalB {
  Dim3 B = {};
  Real Bmag = {};
  Real Bsq = {};
};

// Output field structs by the use case
struct EvalBE {
  Dim3 B = {};
  Dim3 E = {};
  Real Bmag = {};
  Real Bsq = {};
};


struct EvalGCE {
  Dim3 B = {};
  Dim3 dBdR = {};
  Dim3 dBdZ = {};
  Dim3 J = {};
  Dim3 E = {};
  Dim3 dbdt = {};
  Real Bsq;
  Real Bmag;
};


// TODO: wrap in concept ... requires ... clause to shorten the code

template<class CoeffView>
KOKKOS_INLINE_FUNCTION
void evalTaylorTimed(EvalB& out, Real x, Real y, Real t, CoeffView pcofs) {
  const int Px = pcofs.extent_int(0);
  const int Py = pcofs.extent_int(1);

  out = {};

  Real sclx = 1.0;
  for (int i = 0; i < Px; ++i) {
    Real mon = sclx;
    for (int j = 0; j < Py; ++j) {
       for (int i_d = 0; i_d < 3; ++i_d) {
          out.B[i_d] += mon * (pcofs(i, j, FieldComponents::B  + i_d) +
                               pcofs(i, j, FieldComponents::Bt + i_d) * t);
       }
       mon *= y;
    }
    sclx *= x;
  }
}

template<class CoeffView>
KOKKOS_INLINE_FUNCTION
void evalTaylorTimed(EvalBE& out, Real x, Real y, Real t, CoeffView pcofs) {
  const int Px = pcofs.extent_int(0);
  const int Py = pcofs.extent_int(1);

  out = {};

  Real sclx = 1.0;
  for (int i = 0; i < Px; ++i) {
    Real mon = sclx;
    for (int j = 0; j < Py; ++j) {
       for (int i_d = 0; i_d < 3; ++i_d) {
          out.B[i_d] += mon * (pcofs(i, j, FieldComponents::B  + i_d) +
                               pcofs(i, j, FieldComponents::Bt + i_d) * t);
          out.E[i_d] += mon * (pcofs(i, j, FieldComponents::E  + i_d) +
                               pcofs(i, j, FieldComponents::Et + i_d) * t);
       }
       mon *= y;
    }
    sclx *= x;
  }
}

template<class CoeffView>
KOKKOS_INLINE_FUNCTION
void evalTaylorTimed(EvalGCE& out, Real x, Real y, Real t, CoeffView pcofs) {
  const int Px = pcofs.extent_int(0);
  const int Py = pcofs.extent_int(1);

  out = {};

  Real sclx = 1.0;
  for (int i = 0; i < Px; ++i) {
    Real mon = sclx;
    for (int j = 0; j < Py; ++j) {
       for (int i_d = 0; i_d < 3; ++i_d) {
          out.B[i_d] += mon * (pcofs(i, j, FieldComponents::B  + i_d) +
                               pcofs(i, j, FieldComponents::Bt + i_d) * t);

          out.dbdt[i_d] += pcofs(i, j, FieldComponents::Bt + i_d) * mon;

          if (i < Px - 1) {
            out.dBdR[i_d] += (i + 1) * mon * (
              pcofs(i + 1, j, FieldComponents::B  + i_d) +
              pcofs(i + 1, j, FieldComponents::Bt + i_d) * t);
          }

          if (j < Py - 1) {
            out.dBdZ[i_d] += (j + 1) * mon * (
              pcofs(i, j + 1, FieldComponents::B  + i_d) +
              pcofs(i, j + 1, FieldComponents::Bt + i_d) * t);
          }

          out.J[i_d] += mon * (pcofs(i, j, FieldComponents::J  + i_d) +
                               pcofs(i, j, FieldComponents::Jt + i_d) * t);

          out.E[i_d] += mon * (pcofs(i, j, FieldComponents::E  + i_d) +
                               pcofs(i, j, FieldComponents::Et + i_d) * t);
        }
        mon *= y;
      }
      sclx *= x;
    }
}



template<class CoeffView>
struct FieldEvaluator {

  StructuredLocator locator;
	CoeffView coeffs;

	KOKKOS_INLINE_FUNCTION
	void eval(EvalB& out, const Real R, const Real Z, const Real t) const {

		int iR = 0;
		int iZ = 0;

    Real xiR = 0.0;
    Real xiZ = 0.0;

    locator.locate(R, Z, iR, iZ, xiR, xiZ);
    KOKKOS_ASSERT(coeffs.extent(2) >= FieldComponents::Total);
    KOKKOS_ASSERT(coeffs.extent(3) > iR && iR >= 0);
    KOKKOS_ASSERT(coeffs.extent(4) > iZ && iZ >= 0);

    auto cell = Kokkos::subview(coeffs, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iR, iZ);

    evalTaylorTimed(out, xiR, xiZ, t, cell);
    out.Bsq = 0.0;
    for (int i = 0; i < 3; ++i) {
      out.B[i] /= R;
      out.Bsq += out.B[i] * out.B[i];
    }
    out.Bmag = Kokkos::sqrt(out.Bsq);
	}

  KOKKOS_INLINE_FUNCTION
  void eval(EvalBE& out, const Real R, const Real Z, const Real t) const {

    int iR = 0;
    int iZ = 0;

    Real xiR = 0.0;
    Real xiZ = 0.0;

    locator.locate(R, Z, iR, iZ, xiR, xiZ);
    KOKKOS_ASSERT(coeffs.extent(2) >= FieldComponents::Total);
    KOKKOS_ASSERT(coeffs.extent(3) > iR && iR >= 0);
    KOKKOS_ASSERT(coeffs.extent(4) > iZ && iZ >= 0);

    auto cell = Kokkos::subview(coeffs, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iR, iZ);

    evalTaylorTimed(out, xiR, xiZ, t, cell);
    out.Bsq = 0.0;
    for (int i = 0; i < 3; ++i) {
      out.B[i] /= R;
      out.Bsq += out.B[i] * out.B[i];
    }
    out.Bmag = Kokkos::sqrt(out.Bsq);
  }



	KOKKOS_INLINE_FUNCTION
	void eval(EvalGCE& out, const Real R, const Real Z, const Real t) const {

		int iR = 0;
		int iZ = 0;

    Real xiR = 0.0;
    Real xiZ = 0.0;

    locator.locate(R, Z, iR, iZ, xiR, xiZ);
    KOKKOS_ASSERT(coeffs.extent(2) >= FieldComponents::Total);
    KOKKOS_ASSERT(coeffs.extent(3) > iR && iR >= 0);
    KOKKOS_ASSERT(coeffs.extent(4) > iZ && iZ >= 0);

    auto cell = Kokkos::subview(coeffs, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, iR, iZ);

    evalTaylorTimed(out, xiR, xiZ, t, cell);

    out.Bsq = 0.0;
    Real BBt = 0.0;
    for (int i = 0; i < 3; ++i) {
      out.B[i] /= R;
      out.dbdt[i] /= R;

      out.Bsq  += out.B[i] * out.B[i];
      BBt += out.B[i] * out.dbdt[i];

      out.dBdR[i] = (out.dBdR[i] / locator.dR - out.B[i]) / R;
      out.dBdZ[i] /= R * locator.dZ;
    }

    out.Bmag = Kokkos::sqrt(out.Bsq);

    for (int i = 0; i < 3; ++i) {
      out.dbdt[i] = (out.dbdt[i] - out.B[i] * BBt / out.Bsq)  / out.Bmag;
    }
  }
};
} // namespace kinetic
