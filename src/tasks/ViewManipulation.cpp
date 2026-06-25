#include "ViewManipulation.h"

TaskStatus ScaleView(View3 view, const Real scl) {

  const int NR = view.extent_int(0);
  const int NZ = view.extent_int(1);
  Kokkos::parallel_for("ScaleView",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {NR,NZ}),
      KOKKOS_LAMBDA(const int i, const int j, const int k) {
        view(i,j,0) *= scl;
        view(i,j,1) *= scl;
        view(i,j,2) *= scl;
      });

  return TaskStatus::complete;
}

TaskStatus ReduceView(View3 view) {

  // Communicate Runaway current accross MPI Ranks, TODO exploit CUDA aware MPI
  Kokkos::fence();
  auto Jre_host = create_mirror_view_and_copy(Kokkos::HostSpace(), Jre);
  MPI_Allreduce(MPI_IN_PLACE, Jre_host.data(), Jre_host.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Kokkos::deep_copy(Jre, Jre_host);
}

TaskStatus AccumulateView(DualView3 acc_view, View3 view) {

  const int NR = view.extent_int(0);
  const int NZ = view.extent_int(1);

  auto acc_view_d = acc_view.view_device();

  acc_view.device_sync();
  Kokkos::parallel_for("Accumulate View",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {NR,NZ}),
      KOKKOS_LAMBDA(const int i, const int j, const int k) {
        acc_view_d(i,j,0) += view(i,j,0);
        acc_view_d(i,j,1) += view(i,j,1);
        acc_view_d(i,j,2) += view(i,j,2);
      });

  acc_view.device_modify();

  return TaskStatus::complete;
}

