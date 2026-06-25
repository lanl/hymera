#pragma once
#include "ConfigurationDomainGeometry.hpp"

struct ParticleVerifyCodes {
  static constexpr int Success = 0;
  static constexpr int InvalidState = 1;
  static constexpr int WallImpact = 2;
  static constexpr int MomentumCutoff = 3;
};

struct ParticleVerificator {
  using ResultCode_t = int;
  static constexpr ResultCode_t Success = ParticleVerifyCodes::Success;

  ConfigurationDomainGeometry cdg;
  Real p_BC;

  KOKKOS_INLINE_FUNCTION
  ResultCode_t verify(const Dim5 &y) const {
    const Real p = y[0];
    const Real R = y[2];
    const Real Z = y[4];

    if (!Kokkos::isfinite(p) || !Kokkos::isfinite(R) || !Kokkos::isfinite(Z)) {
      return ParticleVerifyCodes::InvalidState;
    }

    if (p < p_BC) {
      return ParticleVerifyCodes::MomentumCutoff;
    }

    const Real Rmin = cdg.indicator_locator.R0;
    const Real Zmin = cdg.indicator_locator.Z0;
    const Real Rmax = Rmin + cdg.indicator_locator.nR * cdg.indicator_locator.dR;
    const Real Zmax = Zmin + cdg.indicator_locator.nZ * cdg.indicator_locator.dZ;
    if (R < Rmin || R >= Rmax || Z < Zmin || Z >= Zmax) {
      return ParticleVerifyCodes::WallImpact;
    }

    int i, j;
    Region region;
    cdg.locate_region(R, Z, i ,j, region);

    if (region < 1) {
      return ParticleVerifyCodes::WallImpact;
    }

    return Success;
  }
};
