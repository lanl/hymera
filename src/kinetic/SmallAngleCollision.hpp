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
#include "util/common.hpp"
#include "hFlux/common.hpp"

#include <limits>

struct CollisionCoefficients {
  Real psi, CB, CF, CouLogee_ratio;
};

template <bool PartialScreening = false, bool EnergyDiffusion = true, bool ModifiedCouLog = true>
struct SmallAngleCollision {
  const Real c_vTe;
  const Real Zeff;
  const int NSA;
  const Real CouLog0;
  const Real k;

  // Partial screening parameters, set to 0.0 by defeault
  const Real aI;
  const Real FineStructure;
  const Real Z0;
  const Real ZI;
  const Real NeI;
  const Real II;
  const Real fI;

  SmallAngleCollision(const Real c_vTe, const Real Zeff, const int NSA, const Real CouLog0,
    const Real k,
    const Real aI = 0.0,
    const Real FineStructure = 0.0,
    const Real Z0 = 0.0,
    const Real ZI = 0.0,
    const Real NeI = 0.0,
    const Real II = 0.0,
    const Real fI = 0.0):
    c_vTe(c_vTe), Zeff(Zeff), NSA(NSA), CouLog0(CouLog0), k(k),
    aI(aI), FineStructure(FineStructure),
    Z0(Z0), ZI(ZI),
    NeI(NeI), II(II),
    fI(fI) {}

  KOKKOS_INLINE_FUNCTION
  CollisionCoefficients getCollisionCoefficients(const Real p) const {
    CollisionCoefficients ret = {};

    Real g = gamma_(p);
    Real x = c_vTe * p / g;
    ret.psi = ChsPsi(x);


    Real CouLogee_ratio = 1.0;
    Real CouLogei_ratio = 1.0;
    Real CouLogee = 1.0;
    Real CouLogei = 1.0;

    if constexpr (ModifiedCouLog) {
      CouLogee = CouLog0 + 1.0 / k * log(1.0 + pow(2.0 * (g - 1) * c_vTe * c_vTe, k / 2.0));
      CouLogei = CouLog0 + 1.0 / k * log(1.0 + pow(2.0 * p * c_vTe, k));
      CouLogee_ratio = CouLogee / CouLog0;
      CouLogei_ratio = CouLogei / CouLog0;
    }

    ret.CB = (g / p) * (CouLogei_ratio * Zeff + CouLogee_ratio * (erf(x) - ret.psi + 0.5 / pow(c_vTe, 4) * x * x));
    ret.CF = 2.0 * c_vTe * c_vTe * ret.psi * CouLogee_ratio;
    ret.CouLogee_ratio = CouLogee_ratio;

    if constexpr (PartialScreening) {
      Real yI1d5 = pow(2.0 * aI * p / FineStructure, 1.5);
      Real gI = (2.0 / 3.0) * (Z0 * Z0 - ZI * ZI) * log(yI1d5 + 1.0) -
                (2.0 / 3.0) * NeI * NeI * yI1d5 / (yI1d5 + 1.0);
      ret.CB += (g / p) / CouLog0 * (fI / (1.0 + ZI * fI)) * gI;

      Real hI = p * sqrt(g - 1) / II;
      ret.CF += 2.0 * c_vTe * c_vTe * ret.psi / CouLog0 * (fI / (1.0 + ZI * fI)) * NeI *
                     ((1.0 / k) * log(1.0 + pow(hI, k)) - (p * p) / (g * g));
    }
    return ret;
  }

  template <typename RNGPoolType>
  KOKKOS_INLINE_FUNCTION void operator()(Real &p, Real &xi, const Real dtSA,
                                         RNGPoolType rng_pool) const {
    Real g = gamma_(p); // gamma_
    auto cc = getCollisionCoefficients(p);

    // Pitch Angle Scattering
    Real nuDdtSA = cc.CB / (p * p) * dtSA; // pitch-angle scattering frequency

    // caps pitch-angle scattering rate. Needed due to singularity as v->0.
    if (nuDdtSA > 0.5)
      nuDdtSA = 0.5;

    Real sigma = sqrt((1.0 - xi * xi) * nuDdtSA);

    //  Update pitch of electron
    auto rng_gen = rng_pool.get_state();
    int pm1 = 2 * static_cast<int>(rng_gen.urand(2)) - 1;
    int pm2 = 2 * static_cast<int>(rng_gen.urand(2)) - 1;
    rng_pool.free_state(rng_gen);

    xi = xi * (1.0 - nuDdtSA) + pm1 * sigma;

    // Reflect pitch about |xi| = 1 boundary, if particle has |xi| > 1
    if (xi > 1.0)
      xi = 1.0 - (xi - 1.0);
    else if (xi < -1.0)
      xi = -1.0 - (xi + 1.0);

    // Energy Scattering
    Real dg = 1.E-7;
    Real gp1 = g + dg;
    Real pp1 = momentum_(gp1);
    Real xp1 = c_vTe * pp1 / gp1;

    Real ChandFuncp1 = ChsPsi(xp1);
    Real ChandFuncm1;

    if (g - dg > 1.0) {
      Real gm1 = g - dg;
      Real pm1 = momentum_(gm1);
      Real xm1 = c_vTe * pm1 / gm1;
      ChandFuncm1 = ChsPsi(xm1);
    } else {
      ChandFuncm1 = cc.psi;
    }
    Real dChandFuncdg;
    if (g - dg > 1.0) {
      dChandFuncdg = (ChandFuncp1 - ChandFuncm1) / (2.0 * dg);
    } else {
      dChandFuncdg = (ChandFuncp1 - ChandFuncm1) / dg;
    }

    Real gConv = cc.CouLogee_ratio * (2.0 * (cc.psi / p) + p / g * dChandFuncdg);
    sigma = sqrt(2.0 * cc.CouLogee_ratio * p / g * cc.psi * dtSA);

    Real kick = gConv * dtSA + pm2 * sigma;

    // Energy Drag
    Real drag = cc.CF * p / g * dtSA;
    if constexpr (!EnergyDiffusion) kick = 0.0;
    g += kick - drag;
    if (g < 1.0)
      g = 2.0 - g + 4.E-14;

    p = momentum_(g);
  }

  /**
  Compute timesteps based on particle energy (momemtum)
  */
  KOKKOS_INLINE_FUNCTION Real getSmallAngleCollisionTimestep(
      const Real p, const Real dtSA_min = 0.0,
      const Real dtSA_max = std::numeric_limits<Real>::max()) const {

    auto cc = getCollisionCoefficients(p);
    Real dtSA = p * p / (cc.CB * NSA);

    if (dtSA > dtSA_max)
      return dtSA_max;
    else if (dtSA < dtSA_min)
      return dtSA_min;
    else
      return dtSA;
  }

  void print() {
    if constexpr (PartialScreening) std::cout << "PartialScreening = on\n";
    else std::cout << "PartialScreening = off\n";
    std::cout << "c_vTe" <<  c_vTe << std::endl;
    std::cout << "Zeff" <<  Zeff  << std::endl;
    std::cout << "NSA"  <<  NSA   << std::endl;
    std::cout << "k            " <<k             << std::endl;
    std::cout << "CouLog0      " <<CouLog0       << std::endl;

    // Partial screening parameters, set to 0.0 by defeault
    std::cout << "aI           " <<aI            << std::endl;
    std::cout << "FineStructure" <<FineStructure << std::endl;
    std::cout << "Z0           " <<Z0            << std::endl;
    std::cout << "ZI           " <<ZI            << std::endl;
    std::cout << "NeI          " <<NeI           << std::endl;
    std::cout << "II           " <<II            << std::endl;
    std::cout << "fI           " <<fI            << std::endl;

  }
};

