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

set(SRC_LIST
  mhd/main.c
  mhd/ts_functions.c
  mhd/mass_matrix_coefficients.c
  mhd/mimetic_operators.c
  mhd/geometry.c
  mhd/monitor_functions.c
  mhd/mhd_run.c
)

add_executable(mhd ${SRC_LIST})
target_include_directories(mhd PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/mhd
  )

target_link_libraries(mhd PRIVATE PkgConfig::PETSC)
target_link_libraries(mhd PRIVATE HDF5::HDF5 Parthenon::parthenon hflux::hflux)
