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
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/EM_Field.hpp"
#include "pgen.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin) {


  const Real B0 = pin->GetOrAddReal("MHD", "B0", 5.3);     ///< On-axis magnetic field [m/s]
  const Real V_A = pin->GetOrAddReal("MHD", "V_A", 1.1e7); ///< Alfven speed [m/s]

  ///< Impurity Parameters
  const Real a = pin->GetOrAddReal("Geometry", "a", 2.0);  ///< Minor radius [m]
  const Real Z0 = pin->GetOrAddReal("ImpurityParameters", "Z0", 10.0); ///< Atomic number of impurity (Z)
  const Real ZI = pin->GetOrAddReal("ImpurityParameters", "ZI", 1.0); ///< Charge of impurity
  const Real NeI = Z0 - ZI; ///< Number of bound electrons
  const Real T_e0 = pin->GetOrAddReal("Distrubution", "T_e0", 10.48); ///< Initial deuterium density [m^-3]
  const Real nD0  = pin->GetOrAddReal("Distrubution", "nD0", 6.8004357342523336e18); ///< Initial deuterium density [m^-3]
  const Real fI = pin->GetOrAddReal("ImpurityParameters", "fI", 133.33);  ///<Fraction of impurity density, normalized to deuterium denstiy (nD0)
  const Real nI = fI*nD0; ///< Impurity density [m^-3]
  const Real n_e0 = nD0 + ZI*nI; ///< Free electron density [m^-3]
  const Real Zeff = pin->GetOrAddReal("Collisions", "Zeff", (ZI*ZI*nI + nD0)/n_e0);
  const Real Coulog0 = pin->GetOrAddInteger("Collisions", "Coulog0", 14.9 - 0.5*log(n_e0/1.0e20) + log(T_e0/1.0e3));
  const Real Rc = pin->GetOrAddReal("ParticleSeed", "Rcenter", 3.1158966549999998e+00);
  const Real Zc = pin->GetOrAddReal("ParticleSeed", "Zcenter", 3.7114360000000002e-01);

  const Real tau_a = 6*M_PI*pc::eps0*pc::me*pc::me * pc::qe*pc::c*pc::c*pc::c /
    (pc::qe*pc::qe*pc::qe*pc::qe*pc::qe*B0*B0); ///<
  const Real tau_c = 4*M_PI*pc::eps0*pc::eps0*pc::me*pc::me*pc::c*pc::c*pc::c/(pc::qe*pc::qe*pc::qe*pc::qe*n_e0*Coulog0); ///< relativistic collision time

  const Real eta = pin->GetOrAddReal("Plasma parameters", "eta", 0.0006101792698927684);
  const Real E_c = pc::me * pc::c / pc::qe / tau_c;
  const Real E_n = B0 * V_A / E_c;
  const Real eta_mu0aVa = eta / pc::mu0 / a / V_A; // converts eta * \curl B to V_A B_0
  const Real etaec_a3VaB0 = eta * pc::qe * pc::c / a / a / a / V_A / B0; // converts eta J to V_A B_0

  auto pkg = std::make_shared<StateDescriptor>("Deck");

  const std::string filePath = pin->GetOrAddString("Simulation", "file_path", "avalanche.out");
  pkg->AddParam("filePath", filePath);

  int * ts = new int;
  *ts = 0;
  pkg->AddParam("ts", ts);

  const Real final_time = pin->GetOrAddReal("Simulation", "final_time", 5.0);
  const Real timeStep = pin->GetOrAddReal("Simulation", "hRK", 1.e-6);
  const Real atol = pin->GetOrAddReal("Simulation", "atol", 1.e-10);
  const Real rtol = pin->GetOrAddReal("Simulation", "rtol", 1.e-7);
  pkg->AddParam("final_time", final_time);
  pkg->AddParam("hRK", timeStep);
  pkg->AddParam("atol", atol);
  pkg->AddParam("rtol", rtol);


  const Real dR = 0.0345;
  const Real dZ = 0.02975;
  const int nR_data = 100;
  const int nZ_data = 200;

  const Real Rmin = 3.05/2.0 + 0.5 * dR;
  const Real Rmax = Rmin + (nR_data-1) * dR;
  const Real Zmin = -5.95/2.0 + 0.5 * dZ;
  const Real Zmax = Zmin + (nZ_data-1) * dZ;

  std::cout << nR_data << " " << nZ_data << std::endl;
  std::cout << Rmin << " " << Rmax << std::endl;
  std::cout << Zmin << " " << Zmax << std::endl;
  std::cout << dR << " " << dZ << std::endl;

  pkg->AddParam("Rmin", Rmin);
  pkg->AddParam("Rmax", Rmax);
  pkg->AddParam("Zmin", Zmin);
  pkg->AddParam("Zmax", Zmax);


  const Real gamma_min =
      pin->GetOrAddReal("BoundaryConditions", "gamma_min", 1.02);
  const Real p_BC =
      momentum_(pin->GetOrAddReal("BoundaryConditions", "gamma_BC", 1.02));
  const Real p_RE =
      momentum_(pin->GetOrAddReal("BoundaryConditions", "gamma_RE", 1.02));
  pkg->AddParam("gamma_min", gamma_min);
  pkg->AddParam("p_BC", p_BC);
  pkg->AddParam("p_RE", p_RE);

  const Real vTe = sqrt(2*T_e0 * pc::qe/pc::me); // Thermal velocity
  const Real c_vTe = pin->GetOrAddReal("Collisions", "c_vTe", pc::c / vTe);
  const int NSA = pin->GetOrAddInteger("Collisions", "NSA", 150);
  const Real k = pin->GetOrAddReal("Collisions", "k", 5.0);
  const Real aI            = pin->GetOrAddReal("Collisions", "aI",              0.3285296762792767);
  const Real FineStructure = 1. / 137.035999;  // Fine Structure constant
  const Real II            = pin->GetOrAddReal("Collisions", "II", 219.5 / pc::me / pc::c / pc::c); // Mean exitation energy

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

  int nphi_data = 1;
  int nt = 2;

  const std::string configurationdomain_file = pin->GetOrAddString("Geometry", "file_path", "../../inputs/AxisSymmetricGeometry.dat");

  ConfigurationDomainGeometry::IndicatorViewType indicator("indicator", nR_data, nZ_data);
  std::ifstream ifs(configurationdomain_file);
  auto indicator_h = Kokkos::create_mirror_view(indicator);
  for (int i = 0; i < nR_data; ++i) {
    for (int j = 0; j < nZ_data; ++j) {
      ifs >> indicator_h(i,j);
    }
  }

  Kokkos::deep_copy(indicator, indicator_h);
  int ic = std::floor((Rc - Rmin) / dR);
  int jc = std::floor((Zc - Zmin) / dZ);

  const ConfigurationDomainGeometry cdg(Rmin, Zmin, dR, dZ, -3, indicator);
  std::cout << "Location of center " << Rc << " " << Zc << std::endl;
  std::cout <<  ic << " " << jc << std::endl;
  std::cout << indicator_h(ic, jc) <<  std::endl;


  const std::string field_file = pin->GetOrAddString("MHD", "file_path", "../../inputs/fields_0.bin");

  EM_Field field_interpolation(nR_data, nZ_data, nphi_data, nt, Rmin, Zmin, dR, dZ, E_n, eta_mu0aVa, etaec_a3VaB0, cdg);
  auto field_data = field_interpolation.getDataRef();
  auto field_data_h = create_mirror_view(Kokkos::HostSpace(), field_data);
  std::ifstream in(field_file, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot open file: " + field_file);

  in.read(reinterpret_cast<char*>(field_data_h.data()), nR_data * nZ_data * nphi_data * nt * 4 * 3 * sizeof(double));
  if (!in) throw std::runtime_error("Failed to read file: " + field_file);

  Kokkos::deep_copy(field_data, field_data_h);

  Kokkos::parallel_for("setfields",
  Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {nR_data,nZ_data}),
  KOKKOS_LAMBDA(int i, int j){
    // linearize: row-major numbering
    auto sbv = Kokkos::subview(field_data,
             i, j, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);

    for (int di = 0; di < sbv.extent(1); ++di) {
      for (int k = 0; k < sbv.extent(2); ++k) {
        for (int ti = 0; ti < sbv.extent(3); ++ti) {
          sbv(1,di,k,ti) = sbv(2,di,k,ti);
          sbv(2,di,k,ti) = 0.0;
          sbv(3,di,k,ti) = 0.0;
        }
      }
    }
  });


  field_interpolation.interpolate();

  pkg->AddParam("CDG", cdg);
  pkg->AddParam("Field", field_interpolation);

  const Real wce0 = pc::qe * B0 / pc::me; // Electron gyrofrequency
  const Real c_aw0 =  pin->GetOrAddReal("GuidingCenterEquations", "c_aw0", pc::c/a/wce0);
  const Real ct_a =   pin->GetOrAddReal("GuidingCenterEquations", "ct_a", pc::c/vTe);
  const Real alpha0 = pin->GetOrAddReal("GuidingCenterEquations", "alpha0", tau_c/tau_a);

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

  pkg->AddParam("Rc", Rc);
  pkg->AddParam("Zc", Zc);

  const Real seed_current = pin->GetOrAddReal("ParticleSeed", "current", 150); // 150 Amps
  pkg->AddParam("seed_current", seed_current * a / pc::qe / pc::c); // Convert from amps

  const Real gammamin = pin->GetOrAddReal("ParticleSeed", "gammamin", 10.0);
  pkg->AddParam("pmin", momentum_(gammamin));
  const Real gammamax = pin->GetOrAddReal("ParticleSeed", "gammamax", 20.0);
  pkg->AddParam("pmax", momentum_(gammamax));

  const Real ximin = pin->GetOrAddReal("ParticleSeed", "ximin", 0.8);
  pkg->AddParam("ximin", ximin);
  const Real ximax = pin->GetOrAddReal("ParticleSeed", "ximax", 1.0);
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

  std::cout << "Coulog0 " << Coulog0 << std::endl;
  std::cout << "c_vTe " << c_vTe << std::endl;

  return pkg;
}

class RunawayDriver : public EvolutionDriver {
public:
	Real gamma_min;

  RunawayDriver(ParameterInput *pin, ApplicationInput *app_in, Mesh *pmesh)
      : EvolutionDriver(pin, app_in, pmesh) {
    tm.dt = pin->GetOrAddReal("parthenon/time", "dt", 1e-5);

  	auto pkg = pmesh->packages.Get("Deck");
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

TaskStatus ComputeParticleWeights(Mesh *pm) {
  auto md = pm->mesh_data.Get();
  auto pkg = pm->packages.Get("Deck");
  const auto field_interpolation = pkg->Param<EM_Field>("Field");

  const Real p_RE = pkg->Param<Real>("p_RE");
  const Real seed_current = pkg->Param<Real>("seed_current");

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");

  Real I_re = 0.0;

  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());

  std::cout << "Calculating current";

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

  Kokkos::fence();
  Real w = seed_current / I_re;
  std::printf("%le %le\n", I_re, w);
  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_r.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        // block and particle indices
        pack_swarm_r(b, Kinetic::weight(), n) = w;
      });

  return TaskStatus::complete;
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
                                         Kinetic::secondary_index>("particles");
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
            auto ret = solve_dopri5(gce, X, t, t + dtSA, rtol, atol, h, 1e-9,
                         std::numeric_limits<int>::max(), work_d);
            if (ret != SUCCESS) {
              swarm_d.MarkParticleForRemoval(n);
              markers_d.MarkParticleForRemoval(n);
              break;
            }
            int ii,jj;
            int level = cdg.indicator(X, ii,jj);

            if (level != 2 || X[0] < p_BC) {
              swarm_d.MarkParticleForRemoval(n);
              markers_d.MarkParticleForRemoval(n);
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

          pack_swarm_i(b, Kinetic::will_scatter(), n) = 0;

          if (swarm_d.IsMarkedForRemoval(n))
            return;

          pack_swarm_i(b, Kinetic::will_scatter(), n) =
              ms(X[0], w, dtLA, gamma_min, rng_pool);
        }

      });

  Kokkos::parallel_for(
      PARTHENON_AUTO_LABEL,
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {jre.extent(0), jre.extent(1)}),
      // loop over all particles
      KOKKOS_LAMBDA(int i, int j) {
        jre(i,j,0) /= dtLA;
        jre(i,j,1) /= dtLA;
        jre(i,j,2) /= dtLA;
      });

  Kokkos::fence();
  field_interpolation.interpolate(); // sets jre to electric field
  auto ts = pkg->Param<int*>("ts");
  if (*ts % 2000 == 0) dumpToHDF5(field_interpolation, (*ts) / 2000);
  *ts += 1;
  Kokkos::fence();


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
    fprintf(fout, "%20.14le %20.14le %20.14le\n", I_re * pc::qe * pc::c * .5,
        I_re_integral * pc::qe * pc::c * .5, I_ohmic * 5.3  * 2.0 / pc::mu0);
    // Charge * c / minor radius * 1e-6 (to MAmps)"
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
    auto add_secondaries = tl.AddTask(check_scatter, AddSecondaries, pmb.get(), tm.dt, gamma_min);
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
  pman.app_input->ProblemGenerator = GenerateParticleCurrentDensity;

  // call ParthenonInit to set up the mesh
  // scope so that the mesh object, kokkos views, etc, all get cleaned
  // up before kokkos::finalize

  pman.ParthenonInitPackagesAndMesh();
  {
    ComputeParticleWeights(pman.pmesh.get());
    RunawayDriver driver(pman.pinput.get(), pman.app_input.get(), pman.pmesh.get());
		auto driver_status = driver.Execute();
  }
  pman.ParthenonFinalize();
  return (0);
}
