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


#include <iostream>
#include <parthenon/driver.hpp>
#include <parthenon/package.hpp>
#include <parthenon_manager.hpp>

#include <Kokkos_Core.hpp>

#include <hFlux/dopri.hpp>
#include <hFlux/FieldInterpolation.hpp>

using namespace parthenon;
using namespace parthenon::driver::prelude;

constexpr bool PartialScreening = true;
constexpr bool EnergyScattering = false;
constexpr bool ModifiedCouLog = true;

#include "kinetic/AnalyticField.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/EM_Field.hpp"
#include "kinetic/CurrentDensity.hpp"

#include "pgen.hpp"

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {

  auto pkg = std::make_shared<StateDescriptor>("Deck");
  const std::string filePath =
      pin->GetOrAddString("Simulation", "file_path", "avalanche.out");
  pkg->AddParam("filePath", filePath);

  const Real final_time = pin->GetOrAddReal("Simulation", "final_time", 5.0);
  const Real timeStep = pin->GetOrAddReal("Simulation", "hRK", 1.e-6);
  const Real atol = pin->GetOrAddReal("Simulation", "atol", 1.e-10);
  const Real rtol = pin->GetOrAddReal("Simulation", "rtol", 1.e-7);
  pkg->AddParam("final_time", final_time);
  pkg->AddParam("hRK", timeStep);
  pkg->AddParam("atol", atol);
  pkg->AddParam("rtol", rtol);

  const Real gamma_min =
      pin->GetOrAddReal("BoundaryConditions", "gamma_min", 1.02);
  const Real p_BC =
      momentum_(pin->GetOrAddReal("BoundaryConditions", "gamma_BC", 1.02));
  const Real p_RE =
      momentum_(pin->GetOrAddReal("BoundaryConditions", "gamma_RE", 1.02));
  pkg->AddParam("gamma_min", gamma_min);
  pkg->AddParam("p_BC", p_BC);
  pkg->AddParam("p_RE", p_RE);

  const Real c_vTe =
      pin->GetOrAddReal("Collisions", "c_vTe", 357.42095304128196);
  const Real Zeff = pin->GetOrAddReal("Collisions", "Zeff", 1.0);
  const int NSA = pin->GetOrAddInteger("Collisions", "NSA", 150);
  const Real Coulog0 = pin->GetOrAddInteger("Collisions", "Coulog0", 10.0);
  const Real k = pin->GetOrAddReal("Collisions", "k", 5.0);
  const Real aI            = pin->GetOrAddReal("Collisions", "aI",              0.3285296762792767);
  const Real FineStructure = pin->GetOrAddReal("Collisions", "FineStructure",   1. / 137. );
  const Real Z0            = pin->GetOrAddReal("Collisions", "Z0",              18.);
  const Real ZI            = pin->GetOrAddReal("Collisions", "ZI",              1.);
  const Real NeI           = Z0 - ZI;
  const Real II            = pin->GetOrAddReal("Collisions", "II",              219.5 / 511.e+03);
  const Real fI            = pin->GetOrAddReal("Collisions", "fI",              1.0);
  const Real PSCoefDnRA    = 1.0 + NeI * fI / (1.0 + ZI * fI);



	SmallAngleCollision<PartialScreening, EnergyScattering, ModifiedCouLog> sa(c_vTe, Zeff, NSA, Coulog0, k,
     aI,
     FineStructure,
     Z0,
     ZI,
     NeI,
     II,
     fI
  );


  pkg->AddParam("SmallAngleCollision", sa);
  MollerSource ms(Coulog0, PSCoefDnRA);
  pkg->AddParam("MollerSource", ms);

  const Real q0 = pin->GetOrAddReal("AnalyticField", "q0", 2.1);
  const Real q2 = pin->GetOrAddReal("AnalyticField", "q2", 2.0);
  const Real R_a = pin->GetOrAddReal("AnalyticField", "AspectRatio", 3.0);
  const Real E_0 = pin->GetOrAddReal("AnalyticField", "E_0", 70.0);
  AnalyticField af(q0, q2, R_a, E_0);
  pkg->AddParam("Field", af);

  int nR_data = 100;
  int nZ_data = 200;

  int nphi_data = 1;
  int nt = 1;
  Real Rmin = 1.525;
  Real Zmin = -2.975;
  Real dR = 0.0345;
  Real dZ = 0.02975;

  EM_Field field_interpolation(nR_data, nZ_data, nphi_data, nt, Rmin, Zmin, dR, dZ, 0.0, 0.0);
  auto field_data = field_interpolation.getDataRef();

  using policy2D = Kokkos::MDRangePolicy<Kokkos::DefaultExecutionSpace, Kokkos::Rank<2>>;
  Kokkos::parallel_for("setfields",
  policy2D({0,0}, {nR_data,nZ_data}),
  KOKKOS_LAMBDA(int i, int j){
    // linearize: row-major numbering
    auto sbv = Kokkos::subview(field_data,
             i, j, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);

    Real R = Rmin + dR * i, Z = Zmin + dZ * j;
    Dim5 X = {};
    X[2] = R; X[4] = Z;
    Dim3 vB = {}, dBdR = {}, dBdZ = {}, curlB = {}, E = {};
    Real t = 0.0;
    af(X, t, vB, curlB, dBdR, dBdZ, E);
    for (int di = 0; di < sbv.extent(1); ++di) {
      for (int k = 0; k < sbv.extent(2); ++k) {
        for (int ti = 0; ti < sbv.extent(3); ++ti) {
          sbv(0,di,k,ti) = vB[di] * R;
          sbv(1,di,k,ti) = E[di];
          sbv(2,di,k,ti) = 0.0;
        }
      }
    }
  });

  field_interpolation.interpolate();

  Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::HostSpace> indicator_h("indicator", nR_data,nZ_data);
  Kokkos::View<Real***, Kokkos::LayoutLeft, Kokkos::HostSpace> jre_h("jre", nR_data, nZ_data, 3);
  for (int j = 0; j < nZ_data; ++j) {
    for (int i = 0; i < nR_data; ++i) {
      if (i > 0 && i < nR_data-1 && j > 0 && j < nZ_data-1)
        indicator_h(i, j) = 2;
      else
        indicator_h(i, j) = -2;
      for (int k = 0; k < 3; ++k)
        jre_h(i, j, k) = 0.0;
    }
  }


  using DevMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
  auto indicator = Kokkos::create_mirror_view_and_copy(DevMemSpace{}, indicator_h);
  Kokkos::View<Real***> jre = Kokkos::create_mirror_view_and_copy(DevMemSpace{}, jre_h);

  const ConfigurationDomainGeometry<Kokkos::View<int**>>
    cdg(Rmin, Zmin, dR, dZ, -3, indicator);

  pkg->AddParam("Field", field_interpolation);
  pkg->AddParam("CDG", cdg);
  pkg->AddParam("Jre", jre);

  const Real c_aw0 = pin->GetOrAddReal("GuidingCenterEquations", "c_aw0",
                                       1.608027383346224e-5);
  const Real ct_a = pin->GetOrAddReal("GuidingCenterEquations", "ct_a", 5.e5);
  const Real alpha0 =
      pin->GetOrAddReal("GuidingCenterEquations", "alpha0", 0.1);
//  GuidingCenterEquations<EM_Field, true, false> gce(field_interpolation, c_aw0, ct_a,
//                                                        alpha0);
  GuidingCenterEquations<EM_Field, true, false> gce(field_interpolation, c_aw0, ct_a,
                                                        alpha0);
  pkg->AddParam("GCE", gce);


  const int npart =
      pin->GetOrAddInteger("ParticleSeed", "num_particles_per_block", 16);
  pkg->AddParam("num_particles_per_block", npart);
  // Initialize random number generator pool
  int rng_seed = pin->GetOrAddInteger("ParticleSeed", "rng_seed", 1234);
  RNGPool rng_pool(rng_seed);
  pkg->AddParam("rng_pool", rng_pool);

  const Real r_0 = pin->GetOrAddReal("ParticleSeed", "r_0", 0.0);
  const Real Rc = pin->GetOrAddReal("ParticleSeed", "Rcenter", 3.0);
  const Real Zc = pin->GetOrAddReal("ParticleSeed", "Zcenter", 0.0);
  pkg->AddParam("r_0", r_0);
  pkg->AddParam("Rcenter", Rc);
  pkg->AddParam("Zcenter", Zc);

  const Real gammamin = pin->GetOrAddReal("ParticleSeed", "gammamin", 10.0);
  pkg->AddParam("pmin", momentum_(gammamin));
  const Real gammamax = pin->GetOrAddReal("ParticleSeed", "gammamax", 20.0);
  pkg->AddParam("pmax", momentum_(gammamax));

  const Real ximin = pin->GetOrAddReal("ParticleSeed", "ximin", -1.0);
  pkg->AddParam("ximin", ximin);
  const Real ximax = pin->GetOrAddReal("ParticleSeed", "ximax", -0.8);
  pkg->AddParam("ximax", ximax);

  Metadata swarm_metadata({Metadata::Provides, Metadata::None});
  pkg->AddSwarm("particles", swarm_metadata);

  Metadata real_swarmvalue_metadata({Metadata::Real});
  pkg->AddSwarmValue(Kinetic::p::name(), "particles",
                     real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::xi::name(), "particles",
                     real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::R::name(), "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::phi::name(), "particles",
                     real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::Z::name(), "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::weight::name(), "particles",
                     real_swarmvalue_metadata);

  Metadata int_swarmvalue_metadata({Metadata::Integer});
  pkg->AddSwarmValue(Kinetic::will_scatter::name(), "particles",
                     int_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::secondary_index::name(), "particles",
                     int_swarmvalue_metadata);

  std::cout << "E_0 " << E_0 << std::endl;
  std::cout << "r_0 " << r_0 << std::endl;
  std::cout << "Coulog0 " << Coulog0 << std::endl;
  std::cout << "c_vTe " << c_vTe << std::endl;

  return pkg;
}

class RunawayDriver : public EvolutionDriver {
public:
	RNGPool rng_pool;
	Real gamma_min;

  RunawayDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh)
      : EvolutionDriver(pin, app_in, pmesh) {
    tm.dt = pin->GetOrAddReal("parthenon/time", "dt", 1e-5);

  	auto pkg = pmesh->packages.Get("Deck");
  	rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  	gamma_min = pkg->Param<Real>("gamma_min");
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      FILE *fout = fopen(pkg->Param<std::string>("filePath").c_str(), "w");
      fclose(fout);
    }
  }

  TaskListStatus Step() {
    PARTHENON_INSTRUMENT
    using DriverUtils::ConstructAndExecuteTaskLists;
    TaskListStatus status = ConstructAndExecuteTaskLists<>(this, tm);
    return status;
  }

  TaskCollection MakeTaskCollection(BlockList_t &blocks, SimTime tm);

};

KOKKOS_INLINE_FUNCTION Real S2(Real x) {
    if (x < 0.5) return 0.75 - x * x;
    else return (3.0 - 2.0 * x) * (3.0 - 2.0 * x) / 8.0;
}

template <class CurrentDensityView, class CDG, class Field>
KOKKOS_INLINE_FUNCTION
void DepositCurrent(const Dim5& X, const Real t, CurrentDensityView jre, const Real time_interval, Field& field, CDG cdg) {

  const Dim5::value_type p = X[0];
  const Dim5::value_type xi = X[1];
  const Dim5::value_type R = X[2];

  Real contribution = -p*xi/gamma_(p) / R / cdg.dR / cdg.dZ / 2.0 / M_PI * time_interval;

  Dim3 B;
  field.evalB(B, X);

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

  const auto filePath = pkg->Param<std::string>("filePath");

  const auto gamma_min = pkg->Param<Real>("gamma_min");
  const auto p_BC = pkg->Param<Real>("p_BC");
  const auto p_RE = pkg->Param<Real>("p_RE");

  std::cout << "Device=" << DevExecSpace().name() << std::endl;

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const auto ms = pkg->Param<MollerSource>("MollerSource");
  const auto af =
      pkg->Param<AnalyticField>("Field");
  const auto gce =
      pkg->Param<GuidingCenterEquations<EM_Field, true, false>>("GCE");
  const auto cdg =
      pkg->Param<ConfigurationDomainGeometry<Kokkos::View<int**>>>("CDG");
  const auto sa =
      pkg->Param<SmallAngleCollision<PartialScreening, EnergyScattering, ModifiedCouLog>>("SmallAngleCollision");
  const Real dtSA_min =
      sa.getSmallAngleCollisionTimestep(momentum_(1.000020));
  const Real dtSA_max = dtLA;

  Kokkos::Timer timer;

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  auto jre = pkg->Param<Kokkos::View<Real***>>("Jre");
  auto field_interpolation = pkg->Param<EM_Field>("Field");

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_r.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        // block and particle indices
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm_r.GetContext(b);
        const auto markers_d = pack_swarm_i.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)) {
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
            solve_dopri5(gce, X, t, t + dtSA, rtol, atol, h, 1e-9,
                         std::numeric_limits<int>::max(), work_d);

            if (X[0] > p_RE) {
              DepositCurrent(X, t, jre, dtSA, af, cdg);
            }



            if (X[0] < p_BC) {
              swarm_d.MarkParticleForRemoval(n);
              markers_d.MarkParticleForRemoval(n);
              break;
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

          pack_swarm_i(b, Kinetic::will_scatter(), n) = 0;

          if (swarm_d.IsMarkedForRemoval(n))
            return;

          pack_swarm_i(b, Kinetic::will_scatter(), n) =
              ms(X[0], w, dtLA, gamma_min, rng_pool);
        }

      });


  Real I_re = 0.0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      // loop over all particles
      KOKKOS_LAMBDA(const int idx, Real &weight) {
        // block and particle indices
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm_r.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)) {
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
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL,
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {100, 200}),
      // loop over all particles
      KOKKOS_LAMBDA(int i, int j, Real& integral) {
        integral += cdg.dR * cdg.dZ * jre(i,j,1) / dtLA;
        jre(i,j,0) = 0.0;
        jre(i,j,1) = 0.0;
        jre(i,j,2) = 0.0;
      },
      I_re_integral);

  Kokkos::fence();

  if (rank == 0) {
    FILE *fout = fopen(filePath.c_str(), "a");
    fprintf(fout, "%20.14le %20.14le \n", I_re, I_re_integral);
    fclose(fout);
  }

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
                                         Kinetic::secondary_index>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(data.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(data.get());

  auto swarm_d = swarm->GetDeviceContext();

  Kokkos::parallel_scan(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      KOKKOS_LAMBDA(const int n, int &running_total, const bool final_pass) {
        const int b = 0;
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)) {
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

TaskStatus AddSecondaries(MeshBlock* pmb, const Real dtLA,  const Real gamma_min, RNGPool rng_pool) {
  auto data = pmb->meshblock_data.Get();
  auto swarm = data->GetSwarmData()->Get("particles");
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(data.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(data.get());

  auto swarm_d = swarm->GetDeviceContext();
  int ntot = 0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      KOKKOS_LAMBDA(const int n, int &nnew) {
        const int b = 0;
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n)) {
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
        Kinetic::will_scatter, Kinetic::secondary_index>("particles");
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
          if (swarm_d.IsActive(n_primary) &&
              !swarm_d.IsMarkedForRemoval(n_primary))
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
              if (p < 0.0) {
                swarm_d.MarkParticleForRemoval(n);
              }
              pack_swarm_r(b, Kinetic::p(), n) = p;
              pack_swarm_r(b, Kinetic::xi(), n) = xi;
              pack_swarm_r(b, Kinetic::weight(), n) = w;
            }
        });
    Kokkos::fence();
  }

	return TaskStatus::complete;
}

TaskStatus CleanupParticles(MeshBlock* pmb) {
  pmb->meshblock_data.Get()
  ->GetSwarmData()->Get("particles")
  ->RemoveMarkedParticles();
	return TaskStatus::complete;
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
    auto add_secondaries = tl.AddTask(check_scatter, AddSecondaries, pmb.get(), tm.dt, gamma_min, rng_pool);
    auto cleanup = tl.AddTask(add_secondaries, CleanupParticles, pmb.get());
  }

  return tc;
}

int main(int argc, char *argv[]) {

  parthenon::ParthenonManager pman;

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

  // Redefine parthenon defaults
  pman.app_input->ProcessPackages = [](std::unique_ptr<ParameterInput> &pin) {
    Packages_t packages;
    packages.Add(Initialize(pin.get()));
    return packages;
  };
  pman.app_input->ProblemGenerator = GenerateParticleRing;

  // call ParthenonInit to set up the mesh
  // scope so that the mesh object, kokkos views, etc, all get cleaned
  // up before kokkos::finalize

  pman.ParthenonInitPackagesAndMesh();
  {
    RunawayDriver driver(pman.pinput.get(), pman.app_input.get(), pman.pmesh.get());
		auto driver_status = driver.Execute();
  }
  pman.ParthenonFinalize();
  return (0);
}
