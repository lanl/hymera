#include "RunawayDriver.h"

#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>

#include <Kokkos_Core.hpp>

#include <hFlux/dopri.hpp>

using namespace parthenon;
using namespace parthenon::driver::prelude;

#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/EM_Field.hpp"
#include "pgen.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

namespace Kinetic {

KOKKOS_INLINE_FUNCTION Real S2(Real x) {
    if (x < 0.5) return 0.75 - x * x;
    else return (3.0 - 2.0 * x) * (3.0 - 2.0 * x) / 8.0;
}

template <class CurrentDensityView, class CDG, class Field>
KOKKOS_INLINE_FUNCTION
void DepositCurrent(const Dim5& X, const Real t, const Real w, CurrentDensityView jre, const Real time_interval, Field& field, CDG cdg) {

  const Dim5::value_type p = X[0];
  const Dim5::value_type xi = X[1];
  const Dim5::value_type R = X[2];

  Real contribution = -p * xi / gamma_(p) / R / cdg.dR / cdg.dZ / 2.0 / M_PI * time_interval * w;

  Dim3 B, curlB, dBdR, dBdZ, E;
  ERROR_CODE ret = field(X, t, B, curlB, dBdR, dBdZ, E);
  KOKKOS_ASSERT(ret == SUCCESS);

  int i, j;
  int level = cdg.indicator(X, i, j);

  Dim2 Xlocd = {};
  cdg.getLocalCoordinate(X, i, j, Xlocd);

  if (level < 1) return;

  Real BB = norm_(B);

  for (int ii = -1; ii < 2; ++ii) {
      if(i + ii >= 0 and i + ii < jre.extent(0)) {
          Real wr = S2(abs(Xlocd[0] - static_cast<Real>(ii)));
          for (int jj = -1; jj < 2; ++jj) {
              if(j + jj >= 0 and j + jj < jre.extent(1)) {
                  Real wz = S2(abs(Xlocd[1] - static_cast<Real>(jj)));
                  Real weighted_contribution = contribution * wr * wz;
                  for (int kk = 0; kk < 3; ++kk) {
                      Real wcB = weighted_contribution * B[kk] / BB;
                      Kokkos::atomic_add(&(jre(i,j,kk)), wcB);
                  }
              }
          }
      }
  }
}

TaskStatus PushParticles(Mesh *pm, const Real dtLA) {
  // get mesh data
  auto md = pm->mesh_data.Get();

  auto pkg = pm->packages.Get("Deck");
  const auto h = pkg->Param<Real>("hRK");
  const auto atol = pkg->Param<Real>("atol");
  const auto rtol = pkg->Param<Real>("rtol");
  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");

  const auto gamma_min = pkg->Param<Real>("gamma_min");
  const auto p_BC = pkg->Param<Real>("p_BC");
  const auto p_RE = pkg->Param<Real>("p_RE");

  std::cout << "Device=" << DevExecSpace().name() << std::endl;
  const auto ms = pkg->Param<MollerSource>("MollerSource");
  const auto gce =
      pkg->Param<GuidingCenterEquations<EM_Field, true, false>>("GCE");
  const auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  const auto sa = pkg->Param<SmallAngleCollision<PartialScreening, EnergyScattering, ModifiedCouLog>>("SmallAngleCollision");
  const Real dtSA_min = sa.getSmallAngleCollisionTimestep(momentum_(1.000020));
  const Real dtSA_max = dtLA;
  auto field_interpolation = pkg->Param<EM_Field>("Field");

  Kokkos::Timer timer;

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index,
                                         Kinetic::status>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  auto jre = field_interpolation.getJreDataSubview();

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_r.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        // block and particle indices
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm_r.GetContext(b);
        const auto markers_d = pack_swarm_i.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)&&
            (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE) ) {
          Dim5 X;
          Real t = 0.0;
          X[0] = pack_swarm_r(b, Kinetic::p(), n);
          X[1] = pack_swarm_r(b, Kinetic::xi(), n);
          X[2] = pack_swarm_r(b, Kinetic::R(), n);
          X[3] = pack_swarm_r(b, Kinetic::phi(), n);
          X[4] = pack_swarm_r(b, Kinetic::Z(), n);
          Real w = pack_swarm_r(b, Kinetic::weight(), n);

          bool last_step = false;

          while (last_step == false) {
            Real dtSA =
                sa.getSmallAngleCollisionTimestep(X[0], dtSA_min, dtSA_max);

            if (t + dtSA > dtLA) {
              dtSA = dtLA - t;
              if (dtSA < 1e-16) {
                break;
              }
              last_step = true;
            }

            Kokkos::Array<Dim5, 10> work_d;
            auto ret = solve_dopri5(gce, X, t, t + dtSA, rtol, atol, h, 1e-9,
                         std::numeric_limits<int>::max(), work_d);
            if (ret != SUCCESS) {
              pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
              if ((pack_swarm_i(b, Kinetic::status(), n) & PROTECTED) == 0)
                swarm_d.MarkParticleForRemoval(n);
              break;
            }
            int ii,jj;
            int level = cdg.indicator(X, ii,jj);

            if (level != 2 || X[0] < p_BC) {
              pack_swarm_i(b, Kinetic::status(), n) &= ~Kinetic::ALIVE;
              if ((pack_swarm_i(b, Kinetic::status(), n) & PROTECTED) == 0)
                swarm_d.MarkParticleForRemoval(n);
              break;
            }

            if (X[0] > p_RE) {
              DepositCurrent(X, t, w, jre, dtSA, field_interpolation, cdg);
            }

            sa(X[0], X[1], dtSA, rng_pool);
            t += dtSA;
            if (t > dtLA)
              break;
          }

        	pack_swarm_r(b, Kinetic::p(), n)   = X[0];
        	pack_swarm_r(b, Kinetic::xi(), n)  = X[1];
        	pack_swarm_r(b, Kinetic::R(), n)   = X[2];
        	pack_swarm_r(b, Kinetic::phi(), n) = X[3];
        	pack_swarm_r(b, Kinetic::Z(), n)   = X[4];

          pack_swarm_i(b, Kinetic::will_scatter(), n) =
                 ms(X[0], w, dtLA, gamma_min, rng_pool);
        }
      });
  Kokkos::fence();



  return TaskStatus::complete;
}

TaskStatus CheckScatter(MeshBlock* pmb) {
  auto data = pmb->meshblock_data.Get();
  auto swarm = data->GetSwarmData()->Get("particles");
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index,
                                         Kinetic::status>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(data.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(data.get());

  auto swarm_d = swarm->GetDeviceContext();

  Kokkos::parallel_scan(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      KOKKOS_LAMBDA(const int n, int &running_total, const bool final_pass) {
        const int b = 0;
        if (swarm_d.IsActive(n)&& !swarm_d.IsMarkedForRemoval(n) && (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
          if (pack_swarm_i(b, Kinetic::will_scatter(), n) == 1) {
            running_total += 1;
            if (final_pass) {
              pack_swarm_i(b, Kinetic::secondary_index(), n) = running_total;
            }
          } else {
            if (final_pass) {
              pack_swarm_i(b, Kinetic::secondary_index(), n) = 0;
            }
          }
        }
      });
  Kokkos::fence();

	return TaskStatus::complete;
}

TaskStatus CleanupParticles(MeshBlock* pmb) {
  pmb->meshblock_data.Get()
  ->GetSwarmData()->Get("particles")
  ->RemoveMarkedParticles();
	return TaskStatus::complete;
}

TaskStatus AddSecondaries(MeshBlock* pmb, const Real dtLA,  const Real gamma_min) {
  auto pkg = pmb->packages.Get("Deck");
  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  auto data = pmb->meshblock_data.Get();
  auto swarm = data->GetSwarmData()->Get("particles");
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index,
                                         Kinetic::status>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(data.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(data.get());

  auto swarm_d = swarm->GetDeviceContext();
  int ntot = 0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      KOKKOS_LAMBDA(const int n, int &nnew) {
        const int b = 0;
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)&& (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
          if (pack_swarm_i(b, Kinetic::will_scatter(), n) == 1)
            nnew += 1;
        }
      },
      ntot);
  Kokkos::fence();
  if (ntot > 0) {
    const int oldMaxIndex = pack_swarm_r.GetMaxFlatIndex();
    auto newParticlesContext = swarm->AddEmptyParticles(ntot);
    auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
        swarm_position::x, swarm_position::y, swarm_position::z,
        Kinetic::p, Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z,
        Kinetic::weight>("particles");
    auto desc_swarm_i = parthenon::MakeSwarmPackDescriptor<
        Kinetic::will_scatter, Kinetic::secondary_index, Kinetic::status>("particles");
    pack_swarm_r = desc_swarm_r.GetPack(data.get());
    pack_swarm_i = desc_swarm_i.GetPack(data.get());

    swarm_d = swarm->GetDeviceContext();

    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL, DevExecSpace(), 0,
        newParticlesContext.GetNewParticlesMaxIndex(),
        // new_n ranges from 0 to N_new_particles
        KOKKOS_LAMBDA(const int new_n) {
          // this is the particle index inside the swarm
          const int n = newParticlesContext.GetNewParticleIndex(new_n);
          const int b = 0;
          pack_swarm_i(b, Kinetic::will_scatter(), n) = 0;
          pack_swarm_i(b, Kinetic::secondary_index(), n) = 0;
        });

    parthenon::par_for(
        DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL, DevExecSpace(), 0,
        oldMaxIndex,
        // new_n ranges from 0 to N_new_particles
        KOKKOS_LAMBDA(const int n_primary) {
          const int b = 0;
          // this is the particle index inside the swarm
          if (swarm_d.IsActive(n_primary)&& !swarm_d.IsMarkedForRemoval(n_primary) &&(pack_swarm_i(b, Kinetic::status(), n_primary) & Kinetic::ALIVE))
            if (pack_swarm_i(b, Kinetic::will_scatter(), n_primary) == 1) {
              int new_n =
                  pack_swarm_i(b, Kinetic::secondary_index(), n_primary) - 1;
              const int n = newParticlesContext.GetNewParticleIndex(new_n);
              pack_swarm_r(b, swarm_position::x(), n) =
                  pack_swarm_r(b, swarm_position::x(), n_primary);
              pack_swarm_r(b, swarm_position::y(), n) =
                  pack_swarm_r(b, swarm_position::y(), n_primary);
              pack_swarm_r(b, swarm_position::z(), n) =
                  pack_swarm_r(b, swarm_position::z(), n_primary);
              pack_swarm_r(b, Kinetic::R(), n) =
                  pack_swarm_r(b, Kinetic::R(), n_primary);
              pack_swarm_r(b, Kinetic::phi(), n) =
                  pack_swarm_r(b, Kinetic::phi(), n_primary);
              pack_swarm_r(b, Kinetic::Z(), n) =
                  pack_swarm_r(b, Kinetic::Z(), n_primary);
              Real p = pack_swarm_r(b, Kinetic::p(), n_primary);
              Real xi = pack_swarm_r(b, Kinetic::xi(), n_primary);
              Real w = pack_swarm_r(b, Kinetic::weight(), n_primary);

              LargeAngleCollision(p, xi, w, dtLA, gamma_min, rng_pool);
              if (p > 0.0) {
                pack_swarm_i(b, Kinetic::status(), n) = Kinetic::ALIVE;
                pack_swarm_r(b, Kinetic::p(), n) = p;
                pack_swarm_r(b, Kinetic::xi(), n) = xi;
                pack_swarm_r(b, Kinetic::weight(), n) = w;
              } else {
                pack_swarm_i(b, Kinetic::status(), n) = 0;
                swarm_d.MarkParticleForRemoval(n);
              }
            }
        });

    Kokkos::fence();
  }

	return TaskStatus::complete;
}

void RunawayDriver::PostExecute(parthenon::DriverStatus st) {
  auto pkg = pmesh->packages.Get("Deck");
  auto field_interpolation = pkg->Param<EM_Field>("Field");
  auto jre = field_interpolation.getJreDataSubview();
  auto ts = pkg->Param<int*>("ts");
  const auto filePath = pkg->Param<std::string>("filePath");
  const auto cdg = pkg->Param<ConfigurationDomainGeometry>("CDG");
  const auto dt_cd = pkg->Param<Real>("dt_cd");
  const auto p_RE = pkg->Param<Real>("p_RE");
  dumpToHDF5(field_interpolation, *ts);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  Kokkos::parallel_for(
    PARTHENON_AUTO_LABEL,
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {jre.extent(0), jre.extent(1)}),
    // loop over all particles
    KOKKOS_LAMBDA(int i, int j) {
      jre(i,j,0) /= dt_cd;
      jre(i,j,1) /= dt_cd;
      jre(i,j,2) /= dt_cd;
    });
  Kokkos::fence();
  field_interpolation.interpolate(); // sets jre to electric field

  auto md = pmesh->mesh_data.Get();
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");

  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  Real I_re = 0.0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      // loop over all particles
      KOKKOS_LAMBDA(const int idx, Real &weight) {
        // block and particle indices
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm_r.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n) && (pack_swarm_i(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
          Dim5 X;
          Real t = 0.0;
          X[0] = pack_swarm_r(b, Kinetic::p(), n);
          X[1] = pack_swarm_r(b, Kinetic::xi(), n);
          X[2] = pack_swarm_r(b, Kinetic::R(), n);
          X[3] = pack_swarm_r(b, Kinetic::phi(), n);
          X[4] = pack_swarm_r(b, Kinetic::Z(), n);
          Real w = pack_swarm_r(b, Kinetic::weight(), n);
          if (X[0] > p_RE) {
            weight += getParticleCurrent(X, t, w, field_interpolation);
          }
        }
      },
      I_re);

  Real I_re_integral = 0.0;
  Real I_ohmic = 0.0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL,
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {jre.extent(0), jre.extent(1)}),
      // loop over all particles
      KOKKOS_LAMBDA(int i, int j, Real& integral, Real& integral_ohmic) {
        integral += cdg.dR * cdg.dZ * jre(i,j,1);
        jre(i,j,0) = 0.0;
        jre(i,j,1) = 0.0;
        jre(i,j,2) = 0.0;
        Real t = 0;

        Real R = cdg.R0 + i * cdg.dR;
        Real Z = cdg.Z0 + j * cdg.dZ;
        Dim3 B = {}, curlB = {}, dBdR = {}, dBdZ = {}, E = {};
        Dim5 X = {0.0,0.0,R,0.0,Z};

        auto ret = field_interpolation(X, t, B, curlB, dBdR, dBdZ, E);
        if (ret == SUCCESS) integral_ohmic += cdg.dR * cdg.dZ * curlB[1];
      },
      I_re_integral, I_ohmic);

  Kokkos::fence();

  if (rank == 0) {
    FILE *fout = fopen(filePath.c_str(), "a");
    fprintf(fout, "%20.14le %20.14le %20.14le %20.14le\n", dt_cd * (*ts+1), I_re * pc::qe * pc::c * .5,
        I_re_integral * pc::qe * pc::c * .5, I_ohmic * 5.3  * 2.0 / pc::mu0);
    // Charge * c / minor radius * 1e-6 (to MAmps)"
    fclose(fout);
  }
  *ts += 1;
  Kokkos::fence();
}

TaskCollection RunawayDriver::MakeTaskCollection(BlockList_t &blocks, SimTime tm) {
  TaskCollection tc;
  TaskID none(0);

  auto partitions = pmesh->GetDefaultBlockPartitions();
  int num_partitions = partitions.size();

  // note that task within this region that contains one tasklist per pack
  // could still be executed in parallel
  TaskRegion &single_tasklist_per_pack_region = tc.AddRegion(num_partitions);
  for (int i = 0; i < num_partitions; i++) {
    auto &tl = single_tasklist_per_pack_region[i];
    // Initialize the base MeshData for this partition
    // (this automatically initializes the MeshBlockData objects
    // required by this MeshData object)
    auto &mbase = pmesh->mesh_data.Add("base", partitions[i]);

    // add tasks that are per mesh here
    auto push = tl.AddTask(none, PushParticles, pmesh, tm.dt);
  }

  // these are per block tasklists
  TaskRegion &async_region = tc.AddRegion(blocks.size());
  for (int i = 0; i < blocks.size(); ++i) {
    // required by this MeshData object)
	  auto &pmb = blocks[i];
    auto &tl = async_region[i];
    auto check_scatter = tl.AddTask(none, CheckScatter, pmb.get());
    auto add_secondaries = tl.AddTask(check_scatter, AddSecondaries, pmb.get(), tm.dt, gamma_min);
    auto cleanup = tl.AddTask(add_secondaries, CleanupParticles, pmb.get());
  }

  return tc;
}
}
