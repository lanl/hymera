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

/**  Total moller cross section, which is finite when \f[ \gamma >
  2\gamma_{\min}-1 \f]
  *  \f[
  *     \sigma_M = \frac{2\pi}{(\gamma^2-1)} \left[
                        \frac{1}{2}(\gamma+1)-\gamma_{\min}-\gamma^2\left(
                            \frac{1}{\gamma-\gamma_{\min}}-\frac{1}{\gamma_{\min}-1}
                        \right)
                        + \frac{2\gamma-1}{\gamma-1}\ln\left(
                            \frac{\gamma_{\min}-1}{\gamma-\gamma_{\min}}
                        \right)
                    \right]
  *  \f]
*/

/// Fully ionized \f$ D_{nRa}\f$
/// \f[ D_{n\text{RA}} =
/// \frac{n_e}{4\pi\text{Coulog0}}\frac{dt_{\text{LA}}}{w}\frac p\gamma
/// \sigma_M \f]

/// Partial screening \f$ D_{nRa}\f$
/// \f[ D_{\text{nRA}} =
/// \frac{n_e}{4\pi\text{Coulog0}}\frac{dt_{\text{LA}}}{w}\frac p\gamma
/// \sigma_M \left(1 + N_{eI} f_I \frac{Z_I f_I}{1+Z_I f_I}\right) \f] where
/// \f$N_{eI} = Z_0 - Z_I \f$ is a nubmer of bound electrons, \f$ f_I =
/// 133.33\f$ is a fraction of impuritity density, normalized to deuterium
/// density.
/// \f$Z_0\f$ is a charge of primary runaway electron, \f$Z_I\f$ is a charge
/// of impurity.

struct MollerSource {

  const Real Coulog0;
  const Real PartialScreeningCoefficient;

  MollerSource(const Real Coulog0, const Real PartialScreeningCoefficient = 1.0):
    Coulog0(Coulog0), PartialScreeningCoefficient(PartialScreeningCoefficient) {};

  template <typename RNGPoolType>
  KOKKOS_INLINE_FUNCTION int operator()(const Real p, const Real w,
                                      const Real dtLA,
                                      const Real gSecondaryMin, RNGPoolType rng_pool) const {
  __const__ Real neNorm = 1.0;
  __const__ Real SecWeightFactor = 1.0;

    Real g = gamma_(p);

    if (g < (2.0 * gSecondaryMin - 1.0)) {
      return 0;
    }

    Real sigmaM =
        (2.0 * M_PI / (g * g - 1.0)) *
        (0.5 * (g + 1.0) - gSecondaryMin -
         g * g * (1.0 / (g - gSecondaryMin) - 1.0 / (gSecondaryMin - 1.0)) +
         ((2.0 * g - 1.0) / (g - 1.0)) *
             log((gSecondaryMin - 1.0) / (g - gSecondaryMin)));

    Real DnRA = neNorm / (4.0 * M_PI * Coulog0) * dtLA * w * p / g * sigmaM *
                PartialScreeningCoefficient;

    Real Prob = (1.0 / SecWeightFactor) * DnRA / w;

    auto rng_gen = rng_pool.get_state();
    Real RandomNumber = rng_gen.drand(0.0, 1.0);
    rng_pool.free_state(rng_gen);

    if (RandomNumber < Prob)
      return 1;

    return 0;
  };
};

template <typename RNGPoolType>
KOKKOS_INLINE_FUNCTION void LargeAngleCollision(Real &p, Real &xi, const Real w,
                        const Real dtLA, const Real gSecondaryMin,
                        RNGPoolType rng_pool) {

  __const__ int MaxTrys = 1E7;
  Real g = gamma_(p);

  Real gSecondary;
  Real xiSecondary;

  // Defines normalization for sampling. This should be the maximum value
  // of the source term. However, since the function "Pi" is singular,
  // "MaxPi" is used to set an effective limit.

  Real x = (g - 1.0) * (g - 1.0) / ((gSecondaryMin - 1.0) * (g - gSecondaryMin));
  Real H = (g > 2.0 * gSecondaryMin - 1.0) ? 1.0 : 0.0; // Heavside function
  Real dsigmaMdg =
      2.0 * M_PI * g * g / ((g - 1.0) * (g - 1.0) * (g - 1.0) * (g + 1.0)) *
      (x * x - 3.0 * x + ((g - 1.0) / g) * ((g - 1.0) / g) * (1.0 + x)) *
      H;           //*0.5*(1.0+tanh((g-(2.0*gSecondaryMin-1.0))/1.0E-12));
  Real Dxi = 0.02; // value used to remove integrable singularity in Pi function
                   // of secondary runaway electron source term
  Real MaxPi =
      (1.0 / M_PI) *
      (1.0 / Dxi); // max value of probability distribution used to determine if
                   // secondary runaway electron will be added

  Real FuncMaxGammaDist = (sqrt(g * g - 1.0) / g) * dsigmaMdg;
  Real FuncMaxPitchDist = MaxPi;
  Real gMaxTry =
      (g + 1.0) /
      2.0; // Max gamma to try, which is set by primary runaway electron's gamma
  int AcceptValue = 0;
  int NumberTrys = 0;
  while (AcceptValue == 0 && NumberTrys < MaxTrys) {
    auto rng_gen = rng_pool.get_state();
    Real TryG = rng_gen.drand(gSecondaryMin, gMaxTry);
    Real TryFunction = rng_gen.drand(FuncMaxGammaDist);
    ;
    rng_pool.free_state(rng_gen);

    Real x = (g - 1.0) * (g - 1.0) / ((TryG - 1.0) * (g - TryG));
    Real H1 = 0.5 * (1 + tanh((g - TryG - (gSecondaryMin - 1.0)) / 1.0E-12));
    Real H2 = 0.5 * (1 + tanh((g - (2.0 * TryG - 1.0)) / 1.0E-12));

    Real TrialFuncVec =
        (2 * M_PI * (sqrt(g * g - 1.0) / g) * g * g /
         ((g - 1.0) * (g - 1.0) * (g - 1.0) * (g + 1.0)) *
         (x * x - 3.0 * x + ((g - 1.0) / g) * ((g - 1.0) / g) * (1.0 + x)) *
         H1 * H2);

    if (TryFunction < TrialFuncVec) {
      gSecondary = TryG;
      AcceptValue = 1;
    }
    NumberTrys += 1;
  }
  if (NumberTrys == MaxTrys) {
    p = -1.0;
    return;
  }
  AcceptValue = 0;
  NumberTrys = 0;

  while (AcceptValue == 0 && NumberTrys < MaxTrys) {
    Real TryG = gSecondary;
    Real xi1 = xi * sqrt((g + 1.0) * (TryG - 1.0) / ((g - 1.0) * (TryG + 1.0)));
    Real xi2 =
        sqrt(2 * (g - TryG) * (1 - xi * xi) / ((g - 1.0) * (TryG + 1.0)));

    auto rng_gen = rng_pool.get_state();
    Real TryXi = rng_gen.drand(-1.0, 1.0);
    Real TryFunction = rng_gen.drand(FuncMaxPitchDist);
    ;
    rng_pool.free_state(rng_gen);

    // Real x = (g - 1.0) * (g - 1.0) / ((TryG - 1.0) * (g - TryG));
    Kokkos::complex<Real> z_arg =
        xi2 * xi2 - (TryXi - xi1) * (TryXi - xi1) + Dxi * Dxi;
    Real PiFunc = Kokkos::real((1.0 / M_PI) / sqrt(z_arg) * 0.5 *
                               (1 + tanh((xi2 - abs(TryXi - xi1)) / 1E-12)));

    // std::complex<Real> foo = (xi2 * xi2 - (TryXi - xi1) * (TryXi - xi1) + Dxi
    // * Dxi) * 0.5 * (1 + tanh((xi2 - abs(TryXi - xi1)) / 1E-12)); Real PiFunc
    // = ((1.0 / M_PI) / sqrt(foo));
    Real TrialFuncVec = std::real(PiFunc);

    if (TryFunction < TrialFuncVec) {
      xiSecondary = TryXi;
      AcceptValue = 1;
    }
    NumberTrys += 1;
  }

  if (NumberTrys == MaxTrys) {
    p = -1.0;
    return;
  }

  p = momentum_(gSecondary);
  xi = xiSecondary;
}

template <typename RNGPoolType>
KOKKOS_INLINE_FUNCTION void LargeAngleCollisionRAMc(Real &p, Real &xi, const Real w,
                        const Real dtLA, const Real gSecondaryMin,
                        RNGPoolType rng_pool) {

    Real g = sqrt(p * p + 1.0);

    // Real gMin = std::max(BC::gCrit, BC::gCut);
    Real gMin = gSecondaryMin;

    Real neNorm = 1.0;
    Real SecWeightFactor = 1.0;


    Real Dxi = 0.02; // value used to remove integrable singularity in Pi function of secondary source term.
    Real MaxPi = (1.0/M_PI) * (1.0/Dxi); // Maximum value of probablity distribution used to generate secondaries
    Real x; // parameter used to define Moller cross section
    Real dsigmadg; // differential Moller cross section
    Real heaviside; // a heaviside function
    Real heaviside1; // a heaviside function
    Real heaviside2; // a heaviside function

    Real  funcMax;
    Real  gMaxTry;
    Real  TryG; // trial value of gamma
    Real  TryXi; // trial value of pitch
    Real  TryZ; // trial value of secondary source
    Real  xi1; // pitch value used in "Pi" function
    Real  xi2; // second pitch value used in "Pi" function
    Real  PiFunc; // "Pi" function used in trial distribution
    Real  TrialFunc; // trial distribution function
    Real  AcceptG; // stores "accepted" gamma for secondary
    Real  AcceptXi; // stores "acdepted" pitch for secondary
    Real  NumTrys; // keeps track of the number of trys made in acceptance-rejection method
    Real  MaxTrys = 10000000; // defines the maximum number of trys so system does not enter an infinite loop
    Real  AcceptValue; // flag to indicate that an accepatable value was obtained.

    //
    // if (1){
    x = (g-1.0)*(g-1.0)/((gMin-1.0)*(g-gMin));
    heaviside = (Real)( g > (2.0*gMin-1.0) );
    dsigmadg = 2.0*M_PI*g*g/((g-1.0)*(g-1.0)*(g-1.0)*(g+1.0)) * ( x*x - 3.0*x + ((g-1.0)/g)*((g-1.0)/g)*(1.0+x) ) * heaviside;
    funcMax = MaxPi * w * (sqrt(g*g-1.0)/g)*dsigmadg;

    // Maximum value of gamma to try. Set by the energy of the primary
    gMaxTry = (g+1.0)/2.0;

    NumTrys = 0;
    AcceptValue = 0;
    while(AcceptValue == 0 && NumTrys < MaxTrys) {
        // Generates "trial" gamma, xi and function values
        auto rng_gen = rng_pool.get_state();
        Real RandomNumber = rng_gen.drand(0.0,1.0);
        TryG = gMin + (gMaxTry-gMin) * RandomNumber;
        RandomNumber = rng_gen.drand(0.0,1.0);
        TryXi = -1.0 + 2.0 * RandomNumber;
        RandomNumber = rng_gen.drand(0.0,1.0);
        rng_pool.free_state(rng_gen);
        // RandomNumber = (double)rand() / (Real)RAND_MAX;
        TryZ = funcMax * RandomNumber;

        // Defines the gamma and xi dependence of the source term. This function
        // is used as the probability function
        x = (g-1.0)*(g-1.0) / ((TryG-1.0)*(g-TryG));
        xi1 = xi*sqrt((g+1.0)*(TryG-1.0)/((g-1.0)*(TryG+1.0)));
        xi2 = sqrt(2.0*(g-TryG)*(1.0-xi*xi)/((g-1.0)*(TryG+1.0)));
        heaviside = (Real)( xi2-abs(TryXi-xi1) > 0.0 );
        // heaviside = ( xi2-abs(TryXi-xi1) > 0.0 );
        PiFunc = (1.0/M_PI)/sqrt(xi2*xi2-(TryXi-xi1)*(TryXi-xi1)+Dxi*Dxi) * heaviside;
        // LOG("xi1 = " + std::to_string(xi1));
        heaviside1 = (Real)( g-TryG-(gMin-1.0) > 0.0 );
        heaviside2 = (Real)( g-(2.0*TryG-1.0) > 0.0 ); //check
                                                       // heaviside1 = ( g-TryG-(gMin-1.0) > 0.0 );
                                                       // heaviside2 = ( g-(2.0*TryG-1.0) > 0.0 );
        TrialFunc = 2.0*M_PI*w*(sqrt(g*g-1.0)/g)*g*g/((g-1.0)*(g-1.0)*(g-1.0)*(g+1.0)) * ( x*x - 3*x + ((g-1.0)/g)*((g-1.0)/g)*(1.0+x) ) * PiFunc * heaviside1 * heaviside2;

        // If TryZ is within TrialFunc, accept values, i.e. use these values of gamma and pitch for secondary
        if ( TryZ < TrialFunc ) {
            AcceptG = TryG;
            AcceptXi = TryXi;
            AcceptValue = 1;
            NumTrys += 1;
        }
    }
    if (NumTrys == MaxTrys) {
        return;
    }
    p = sqrt(AcceptG * AcceptG - 1.0);
    xi = AcceptXi;
    return;
}

