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

#include <memory>
#include <numeric>
#include <typeinfo>  //for 'typeid' to work
#include <parthenon/package.hpp>

using namespace parthenon;

#include "kinetic/kinetic.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/EM_Field.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

namespace Kinetic {

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

  printf("cBn = %le, Jn = %le, En = %le, Ec = %le\n", eta_mu0aVa, etaec_a3VaB0, E_n, E_c);

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

  const Real seed_current = pin->GetOrAddReal("ParticleSeed", "current", 150e3); // 150 kAmps
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
  pkg->AddSwarmValue(Kinetic::saved_p::name(), "particles",
                     real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::saved_xi::name(), "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::saved_R::name(), "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::saved_phi::name(), "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::saved_Z::name(), "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::saved_w::name(), "particles", real_swarmvalue_metadata);

  Metadata int_swarmvalue_metadata({Metadata::Integer});
  pkg->AddSwarmValue(Kinetic::will_scatter::name(), "particles",
                     int_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::secondary_index::name(), "particles",
                     int_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::status::name(), "particles",
                     int_swarmvalue_metadata);

  std::cout << "Coulog0 " << Coulog0 << std::endl;
  std::cout << "c_vTe " << c_vTe << std::endl;

  return pkg;
}

auto &GetCoords(std::shared_ptr<MeshBlock> &pmb) { return pmb->coords; }
auto &GetCoords(MeshBlock *pmb) { return pmb->coords; }
auto &GetCoords(Mesh *pm) { return pm->block_list[0]->coords; }

void ComputeParticleWeights(Mesh* pm) {

  auto md = pm->mesh_data.Get();
  auto pkg = pm->packages.Get("Deck");
  const auto field_interpolation = pkg->Param<EM_Field>("Field");

  const Real p_RE = pkg->Param<Real>("p_RE");
  const Real seed_current = pkg->Param<Real>("seed_current");

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      Kinetic::p, Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
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
}

void SaveState(Mesh* pm) {
  auto md = pm->mesh_data.Get();
  auto pkg = pm->packages.Get("Deck");
  const auto field_interpolation = pkg->Param<EM_Field>("Field");

  const Real p_RE = pkg->Param<Real>("p_RE");
  const Real seed_current = pkg->Param<Real>("seed_current");

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      Kinetic::p, Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight,
      Kinetic::saved_p, Kinetic::saved_xi, Kinetic::saved_R, Kinetic::saved_phi, Kinetic::saved_Z, Kinetic::saved_w>(
      "particles");
  auto desc_swarm_i = parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");

  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_r.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        auto [b_r, n_r] = pack_swarm_r.GetBlockParticleIndices(idx);
        auto [b_i, n_i] = pack_swarm_i.GetBlockParticleIndices(idx);
        // block and particle indices

        if (pack_swarm_i(b_i, Kinetic::status(), n_i) & Kinetic::ALIVE) {
          pack_swarm_i(b_i, Kinetic::status(), n_i) |= Kinetic::PROTECTED;
        } else {
          const auto swarm_r = pack_swarm_r.GetContext(b_r);
          const auto swarm_i = pack_swarm_i.GetContext(b_i);

          swarm_r.MarkParticleForRemoval(n_r);
          swarm_i.MarkParticleForRemoval(n_i);
        }

      });

}

void RestoreState(Mesh* pm) {
  auto md = pm->mesh_data.Get();
  auto pkg = pm->packages.Get("Deck");
  const auto field_interpolation = pkg->Param<EM_Field>("Field");

  const Real p_RE = pkg->Param<Real>("p_RE");
  const Real seed_current = pkg->Param<Real>("seed_current");

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      Kinetic::p, Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight,
      Kinetic::saved_p, Kinetic::saved_xi, Kinetic::saved_R, Kinetic::saved_phi, Kinetic::saved_Z, Kinetic::saved_w>(
      "particles");
  auto desc_swarm_i = parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");

  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_r.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        auto [b_r, n_r] = pack_swarm_r.GetBlockParticleIndices(idx);
        auto [b_i, n_i] = pack_swarm_i.GetBlockParticleIndices(idx);
        // block and particle indices

        if (pack_swarm_i(b_i, Kinetic::status(), n_i) & Kinetic::PROTECTED) {
          pack_swarm_r(b_r, Kinetic::p(), n_r)   = pack_swarm_r(b_r, Kinetic::saved_p(), n_r);
          pack_swarm_r(b_r, Kinetic::xi(), n_r)  = pack_swarm_r(b_r, Kinetic::saved_xi(), n_r);
          pack_swarm_r(b_r, Kinetic::R(), n_r)   = pack_swarm_r(b_r, Kinetic::saved_R(), n_r) ;
          pack_swarm_r(b_r, Kinetic::phi(), n_r) = pack_swarm_r(b_r, Kinetic::saved_phi(), n_r);
          pack_swarm_r(b_r, Kinetic::Z(), n_r)   = pack_swarm_r(b_r, Kinetic::saved_Z(), n_r)  ;
          pack_swarm_r(b_r, Kinetic::weight(), n_r)   = pack_swarm_r(b_r, Kinetic::saved_w(), n_r)  ;

          pack_swarm_i(b_i, Kinetic::status(), n_i) |= Kinetic::ALIVE;
        } else {
          const auto swarm_r = pack_swarm_r.GetContext(b_r);
          const auto swarm_i = pack_swarm_i.GetContext(b_i);

          swarm_r.MarkParticleForRemoval(n_r);
          swarm_i.MarkParticleForRemoval(n_i);
        }
      });
}

} // namespace Kinetic
