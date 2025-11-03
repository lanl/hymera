#========================================================================================
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
#========================================================================================

# ------------------------------------------------------------
# Define a helper to add a Parthenon‐based CTest target
# ------------------------------------------------------------
# Usage:
#   add_parthenon_test(<test-name>
#     SRCS    <list of source files...>
#     INPUT   <input-file-basename>   # looked for under ../inputs/<name>
#   )
#
# The INPUT argument should be just the filename (e.g. "conservation_test.input");
# the function will prefix it with ${CMAKE_CURRENT_SOURCE_DIR}/../inputs/
#
function(add_parthenon_test test_name)
  cmake_parse_arguments(
    PT   # prefix for options/args
    ""   # no boolean options
    "INPUT"   # one-value args
    "SRCS"    # multi-value args
    ${ARGN}
  )

  if(NOT PT_SRCS)
    message(FATAL_ERROR "add_parthenon_test(${test_name}): missing SRCS")
  endif()
  if(NOT PT_INPUT)
    message(FATAL_ERROR "add_parthenon_test(${test_name}): missing INPUT")
  endif()

  # 1) Add the executable
  add_executable(${test_name} ${PT_SRCS})

  # 2) Include dirs
  target_include_directories(${test_name} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
  )

  # 3) C++17
  target_compile_features(${test_name} PRIVATE cxx_std_17)

  # 4) Link against Parthenon
  target_link_libraries(${test_name} PRIVATE HDF5::HDF5 Parthenon::parthenon hflux::hflux)

  # 5) Register the test with CTest, passing the input file
  add_test(
    NAME ${test_name}
    COMMAND ${test_name}
            -i  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../inputs/${PT_INPUT}>
  )
endfunction()
