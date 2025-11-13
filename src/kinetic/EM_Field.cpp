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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <hdf5.h>
#include <Kokkos_DualView.hpp>
#include "EM_Field.hpp"

void dumpToHDF5(EM_Field f, const int i_file, const Real t) {
  const size_t Nplot = 400;

  // Create a rectangular grid withing the borders
  auto corners = f.getCorners();
  const Real Rmin = corners[0] + 1e-10;
  const Real Rmax = corners[1] - 1e-10;
  const Real Zmin = corners[2] + 1e-10;
  const Real Zmax = corners[3] - 1e-10;

  const Real dR = (Rmax - Rmin) / (Nplot-1);
  const Real dZ = (Zmax - Zmin) / (Nplot-1);

  Kokkos::DualView<Real****> dv("dv", Nplot, Nplot, 8, 3);

  dv.modify<Kokkos::DefaultExecutionSpace>();
  auto dview = dv.d_view;

  Kokkos::View<Real******> psi_hermite_data("psi",
      f.hermite_data.extent(0),
      f.hermite_data.extent(1),
      f.hermite_data.extent(2) + 1,
      f.hermite_data.extent(3),
      f.hermite_data.extent(6),
      f.hermite_data.extent(7));
  Kokkos::parallel_for("psi_compute",
  Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {f.nphi_data,f.nt}),
  KOKKOS_LAMBDA(int k, int ti){
    auto sbv_hermite_data = Kokkos::subview(f.hermite_data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, 0, Kokkos::ALL, k, ti);
    auto sbv_psi_data = Kokkos::subview(psi_hermite_data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, k, ti);
    computeFlux<2>(sbv_hermite_data, sbv_psi_data, f.hR, f.hZ);
  });

  Kokkos::parallel_for("eval fields",
  Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {Nplot,Nplot}),
  KOKKOS_LAMBDA(int i, int j){
    Real R = Rmin + dR * i;
    Real Z = Zmin + dZ * j;
    Dim3 B = {}, curlB = {}, dBdR = {}, dBdZ = {}, E = {}, Jre = {}, V = {};
    Dim5 X = {0.0,0.0,R,0.0,Z};

    auto ret = f(X, t, B, curlB, dBdR, dBdZ, E, Jre, V);
    KOKKOS_ASSERT((ret == SUCCESS) || (ret == WALL_IMPACT));
    Real psi;
    ret = f.evalPsi(psi, X, psi_hermite_data);
    KOKKOS_ASSERT(ret == SUCCESS);

    dview(i,j,0,0) = R;
    dview(i,j,0,1) = psi;
    dview(i,j,0,2) = Z;
    for (int dim = 0; dim < 3; ++dim) {
      dview(i,j,1,dim) = B[dim];
      dview(i,j,2,dim) = curlB[dim] * f.eta_mu0aVa * f.E_n;
      dview(i,j,3,dim) = dBdR[dim];
      dview(i,j,4,dim) = dBdZ[dim];
      dview(i,j,5,dim) = E[dim];
      dview(i,j,6,dim) = Jre[dim] * f.etaec_a3VaB0 * f.E_n;
      dview(i,j,7,dim) = V[dim];
    }
  });
  Kokkos::fence();

  dv.sync<Kokkos::HostSpace>();

  auto hview = dv.h_view;
  {
    std::ostringstream oss;
    oss << "fields_" << std::setw(9) << std::setfill('0') << i_file;
    oss << ".h5";
    auto filename = oss.str();
    std::cout << "Writing fields into " << filename << std::endl;

    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims[4] = {
      3, 8,
        static_cast<hsize_t>(Nplot),
        static_cast<hsize_t>(Nplot)
    };
    hid_t space = H5Screate_simple(4, dims, nullptr);
    hid_t dset  = H5Dcreate(file, "dv", H5T_NATIVE_DOUBLE, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write contiguous host data
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
             H5P_DEFAULT, dv.h_view.data());

    H5Dclose(dset);
    H5Sclose(space);
    H5Fclose(file);
  }

}

//TODO: using ndim = 3
