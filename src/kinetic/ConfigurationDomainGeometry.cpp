#include <parthenon/package.hpp>
#include "ConfigurationDomainGeometry.hpp"

using namespace parthenon;

void InitializeGeometry(ParameterInput *pin, std::stared_ptr<StateDescriptor> pkg) {

  ConfigurationDomainGeometry::IndicatorViewType indicator("indicator", NR, NZ);
  std::ifstream ifs(configurationdomain_file);
  auto indicator_h = Kokkos::create_mirror_view(indicator);
  for (int i = 0; i < NR; ++i) {
    for (int j = 0; j < NZ; ++j) {
      ifs >> indicator_h(i,j);
    }
  }
}

