# ========================================================================================
#  (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
#
#  This program was produced under U.S. Government contract 89233218CNA000001 for Los
#  Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC
#  for the U.S. Department of Energy/National Nuclear Security Administration. All rights
#  in the program are reserved by Triad National Security, LLC, and the U.S. Department
#  of Energy/National Nuclear Security Administration. The Government is granted for
#  itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
#  license in this material to reproduce, prepare derivative works, distribute copies to
#  the public, perform publicly and display publicly, and to permit others to do so.
# ========================================================================================

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
)
include("Utils")

add_parthenon_test(avalanche_test
  SRCS
    ../src/util/common.hpp
    ../src/kinetic/AnalyticField.hpp
    ../src/kinetic/GuidingCenterEquations.hpp
    ../src/kinetic/SmallAngleCollision.hpp
    ../src/kinetic/LargeAngleCollision.hpp
    ../src/pgen.cpp
    ../src/pgen.hpp
    ../src/kinetic/EM_Field.hpp
    ../src/kinetic/CurrentDensity.hpp
    kinetic/avalanche_test.cpp
  INPUT
    avalanche_test.input
)

add_parthenon_test(current_replacement
  SRCS
    ../src/util/common.hpp
    ../src/kinetic/AnalyticField.hpp
    ../src/kinetic/GuidingCenterEquations.hpp
    ../src/kinetic/SmallAngleCollision.hpp
    ../src/kinetic/LargeAngleCollision.hpp
    ../src/kinetic/RunawayDriver.cpp
    ../src/kinetic/kinetic.cpp
    ../src/kinetic/EM_Field.cpp
    ../src/pgen.cpp
    ../src/pgen.hpp
    ../src/kinetic/EM_Field.hpp
    ../src/kinetic/CurrentDensity.hpp
    kinetic/current_replacement.cpp
    ../src/kinetic/EM_Field.cpp
    ../src/kinetic/c_wrapper.cpp
  INPUT
    avalanche_test.input
)

find_package(HDF5 COMPONENTS CXX REQUIRED)
