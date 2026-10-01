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

# add_parthenon_test(avalanche_test
#   SRCS
#     ../src/util/common.hpp
#     ../src/kinetic/AnalyticField.hpp
#     ../src/kinetic/GuidingCenterEquations.hpp
#     ../src/kinetic/SmallAngleCollision.hpp
#     ../src/kinetic/LargeAngleCollision.hpp
#     ../src/pgen.cpp
#     ../src/pgen.hpp
#     ../src/kinetic/EM_Field.hpp
#     ../src/kinetic/CurrentDensity.hpp
#     kinetic/avalanche_test.cpp
#   INPUT
#     avalanche_test.input
# )

# add_parthenon_test(current_replacement
#   SRCS
#     ../src/util/common.hpp
#     ../src/kinetic/AnalyticField.hpp
#     ../src/kinetic/GuidingCenterEquations.hpp
#     ../src/kinetic/SmallAngleCollision.hpp
#     ../src/kinetic/LargeAngleCollision.hpp
#     ../src/kinetic/RunawayDriver.cpp
#     ../src/kinetic/kinetic.cpp
#     ../src/kinetic/EM_Field.cpp
#     ../src/pgen.cpp
#     ../src/pgen.hpp
#     ../src/kinetic/EM_Field.hpp
#     ../src/kinetic/CurrentDensity.hpp
#     kinetic/current_replacement.cpp
#     ../src/kinetic/EM_Field.cpp
#     ../src/kinetic/c_wrapper.cpp
#   INPUT
#     avalanche_test.input
# )

find_package(HDF5 COMPONENTS CXX REQUIRED)

# --------------------------------------------------------------------------------------
# MHD regression harness, tier 0: the discrete mass-matrix coefficients.
#
# Evaluates every live alpha-/beta-/rese/condu function in mass_matrix_coefficients.c
# over the full index/location/material space (~9.8M evaluations, a few seconds) and
# compares a SHA-256 of the exact hexadecimal output against a recorded value. Links
# mhd_core directly; needs no input files, no solver and no time step.
#
# This is the gate for deduplicating that file: any refactor of it must leave the
# hash unchanged. On failure, regenerate the full output on both revisions and diff
# them -- tests/regression/run_t0.sh does this.
# --------------------------------------------------------------------------------------
enable_language(C)
add_executable(t0_coefficients ${CMAKE_CURRENT_SOURCE_DIR}/regression/t0_coefficients.c)
target_link_libraries(t0_coefficients PRIVATE mhd_core MPI::MPI_C)
add_test(NAME mhd_t0_coefficients
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regression/run_t0.sh
                 $<TARGET_FILE:t0_coefficients>)
set_tests_properties(mhd_t0_coefficients PROPERTIES LABELS "tier0" TIMEOUT 120)

# Tier 0, second test: the mimetic operators that consume the coefficients above.
# Applies each live Vec-to-Vec operator to a fixed input on the production DOF
# layout and hashes the exact output. Gate for operator-level refactors.
add_executable(t0_operators ${CMAKE_CURRENT_SOURCE_DIR}/regression/t0_operators.c)
target_link_libraries(t0_operators PRIVATE mhd_core MPI::MPI_C m)
add_test(NAME mhd_t0_operators
         COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/regression/run_t0.sh
                 $<TARGET_FILE:t0_operators> --baseline t0_operators)
set_tests_properties(mhd_t0_operators PROPERTIES LABELS "tier0" TIMEOUT 300)
