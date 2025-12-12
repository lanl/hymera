#include <mpi.h>
#include <Kokkos_Core.hpp>

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>
using namespace parthenon;
using namespace parthenon::driver::prelude;

int detect_num_gpus() {
    // Option 1: SLURM_GPUS_ON_NODE (typical on Slurm systems)
    if (const char* env = std::getenv("SLURM_GPUS_ON_NODE")) {
        int n = std::atoi(env);
        if (n > 0) return n;
    }

    // Option 2: CUDA_VISIBLE_DEVICES list length, if present
    if (const char* env = std::getenv("CUDA_VISIBLE_DEVICES")) {
        // e.g. "0,1,3" -> 3 GPUs visible
        int count = 0;
        bool in_token = false;
        for (const char* p = env; *p; ++p) {
            if (*p == ',' || *p == ' ') {
                if (in_token) {
                    ++count;
                    in_token = false;
                }
            } else {
                in_token = true;
            }
        }
        if (in_token) ++count;
        if (count > 0) return count;
    }

    // Fallback (you can also abort here if you require GPUs)
    return 1;
}

int MPI_Comm_get_GPU(const MPI_Comm primary_comm, MPI_Comm& gpu_comm) {

  // 1. Create node-local communicator: all ranks on same node
  MPI_Comm node_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD,
                      MPI_COMM_TYPE_SHARED,
                      0, MPI_INFO_NULL,
                      &node_comm);

  int local_rank = -1, local_size = 0;
  MPI_Comm_rank(node_comm, &local_rank);
  MPI_Comm_size(node_comm, &local_size);

  int gpu_id = Kokkos::device_id();
  MPI_Comm_split(node_comm,
                 gpu_id,      // color: all ranks with same gpu_id together
                 local_rank,  // key: ordering
                 &gpu_comm);

//  MPI_Comm_free(&node_comm);
  return 0;
}

int main(int argc, char** argv) {
   MPI_Init(&argc, &argv);
   ParthenonManager pman;

   // Set up kokkos and read pin
   auto manager_status = pman.ParthenonInitEnv(argc, argv);
   if (manager_status == ParthenonStatus::complete) {
     pman.ParthenonFinalize();
     return 0;
   }
   if (manager_status == ParthenonStatus::error) {
     pman.ParthenonFinalize();
     return 1;
   }

   MPI_Comm gpu_comm;
   MPI_Comm_get_GPU(gpu_comm);

   int gpu_comm_rank = -1, gpu_comm_size = 0;
   MPI_Comm_rank(gpu_comm, &gpu_comm_rank);
   MPI_Comm_size(gpu_comm, &gpu_comm_size);

   char hostname[256];
   std::memset(hostname, 0, sizeof(hostname));
   gethostname(hostname, sizeof(hostname) - 1);
   std::cout << std::format("Hostname: {}, GPU Master: {}\n",
       hostname, (gpu_comm_rank == 0)?"Yes":"No");

   MPI_Comm_free(&gpu_comm);
   pman.ParthenonFinalize();
   return 0;
}
