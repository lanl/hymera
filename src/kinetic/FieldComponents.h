#pragma once

#include <Kokkos_DualView.hpp>
#include "util/common.hpp"

struct FieldComponents {
  static constexpr int B     = 0,
                       J     = 3,
                       E     = 6,
                       Bt    = 9,
                       Jt    = 12,
                       Et    = 15,
                       Total = 18;
};
