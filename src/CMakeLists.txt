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
add_library(mhd_core
  mhd/ts_functions.c
  mhd/mass_matrix_coefficients.c
  mhd/mimetic_operators.c
  mhd/geometry.c
  mhd/monitor_functions.c
  mhd/mfd_config.c
  mhd/mhd.c
  mhd/default_petsc_options.c
)

target_include_directories(mhd_core PRIVATE ${MPI_C_INCLUDE_DIRS})
target_link_libraries(mhd_core PRIVATE MPI::MPI_C)
target_link_libraries(mhd_core PUBLIC PkgConfig::PETSC)
target_include_directories(mhd_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/mhd
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/kinetic
  )

add_library(infrastructure
  kinetic/EM_Field.cpp
  kinetic/kinetic.cpp
)

target_link_libraries(infrastructure PUBLIC PkgConfig::PETSC)
target_link_libraries(infrastructure PUBLIC HDF5::HDF5 Parthenon::parthenon hflux::hflux)
target_include_directories(infrastructure PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/mhd
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/kinetic
  )
target_compile_features(infrastructure PUBLIC cxx_std_20)

add_library(particles_core
  pgen.cpp
)
target_link_libraries(particles_core PUBLIC infrastructure)

set(HYBRID_SRC_LIST
  kinetic/HybridDriver.cpp
  kinetic/hybrid.cpp
)

add_executable(hybrid ${HYBRID_SRC_LIST})
target_link_libraries(hybrid PRIVATE particles_core mhd_core infrastructure)

add_executable(gce kinetic/gce_run.cpp)
target_link_libraries(gce PRIVATE particles_core)

add_executable(setup kinetic/setup.cpp)
target_link_libraries(setup PRIVATE mhd_core particles_core infrastructure)

