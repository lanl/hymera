//========================================================================================(
// (C) (or copyright) 2025. Triad National Security, LLC. All rights reserved.
//
// This program was produced under U.S. Government contract 89233218CNA000001 for Los
// Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC // for the U.S. Department of Energy/National Nuclear Security Administration. All rights
// in the program are reserved by Triad National Security, LLC, and the U.S. Department
// of Energy/National Nuclear Security Administration. The Government is granted for
// itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
// license in this material to reproduce, prepare derivative works, distribute copies to
// the public, perform publicly and display publicly, and to permit others to do so.
//========================================================================================
#include <memory>
#include <numeric>
#include <iostream>
#include <format>
#include <random>
#include <typeinfo>  //for 'typeid' to work
#include <parthenon/package.hpp>
#include <hFlux/dopri.hpp>
#include <hFlux/FieldData.hpp>

using namespace parthenon;

#include "kinetic/RunawayDriver.h"
#include "kinetic/kinetic.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/AnalyticField.hpp"
#include "kinetic/FieldEvaluator.hpp"
#include "mhd/mhd.h"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

namespace Kinetic {
std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin, User* mhd_context) {
  /// Physical constants
  static constexpr Real eps0 = pc::eps0; ///< Vacuum permittivity [F / M]
  static constexpr Real c = pc::c;       ///< Speed of light [m/s]
  static constexpr Real mi = pc::amu;    ///< Ion mass = atomic mass unit [kg]
  static constexpr Real mu0 = pc::mu0;   ///< Vacuum permeabulity [N / A^2]
  static constexpr Real me = pc::me;     ///< electron mass [kg]
  static constexpr Real e = pc::qe;      ///< electron charge [C]

  /// Time discretization parameters
  const Real dt_mhd     = pin->GetReal("parthenon/time","dt_force");
  const Real final_time = pin->GetReal("parthenon/time","tlim");

  const int nCorrectorSteps = pin->GetOrAddInteger("Time","nCorrectorSteps", 0);  ///< current deposit timestep for electric field readjustment [s]
  const int nCD = pin->GetOrAddInteger("Time","nCD", 100);  ///< current deposit timestep for electric field readjustment [s]
  Real dtLA_over_tauC = pin->GetOrAddReal("Time","dtLA_over_tauC", 1e-5);  ///< current deposit timestep for electric field readjustment [s]

  const Real timeStep = pin->GetOrAddReal("Simulation", "hRK", 1.e-6);    /// Runge kutta time in tau_c [-]
  const Real atol = pin->GetOrAddReal("Simulation", "atol", 1.e-6);      /// Absoulte tolerance for RK [-]
  const Real rtol = pin->GetOrAddReal("Simulation", "rtol", 1.e-5);       /// Realative toleratnce for RK[ [-]
  const std::string filePath = pin->GetOrAddString("Simulation", "file_path", "current.out");

  /// Reference parameters
  const Real B0 = pin->GetOrAddReal("Reference", "B0", 5.3);  ///< On-axis magnetic field [T]
  const Real a  = pin->GetOrAddReal("Reference", "a", 2.0);    ///< Minor radius [m] and reference length
  const Real R0 = pin->GetOrAddReal("Reference", "R0", 6.0);    ///< Major radius [m]
  const Real nD0 = pin->GetOrAddReal("Reference", "nD0", 1e20);    ///< Deutirium density [m^-3]
  const Real Te0 = pin->GetOrAddReal("Reference", "Te0", 2.0198);  ///< Electron temperature [eV], and plasma temperature single-temperature model

  /// Derived parameters
  const Real VA   = pin->GetOrAddReal("Derived", "VA",  B0 / sqrt(mi * nD0 * mu0)); ///< Alfven velocity [m/s]
  const Real tauA = pin->GetOrAddReal("Derived", "tauA", a / VA);                  ///< Alfven time     [s]
  const Real E0   = pin->GetOrAddReal("Derived", "E0", B0 * VA);                   ///< Reference electric field in MHD [V/m]
  const Real J0   = pin->GetOrAddReal("Derived", "J0", B0 / (mu0 * a));            ///< Reference current density [A/m^2]
  const Real eta0 = pin->GetOrAddReal("Derived", "eta0", a * VA * mu0);            ///< Reference resitivity [Ohm*m]
  const Real eta  = pin->GetOrAddReal("Derived", "eta",  1.0);                     ///< Resitivity scale     [-]
  const Real Re   = pin->GetOrAddReal("Derived", "Re",  200.0);                    ///< Reinolds Number

  ///< Plasma composition parameters
  const Real vTe = sqrt(2.0*Te0*e / me); ///< Thermal velocity
  const Real Z0 = pin->GetOrAddReal("Plasma", "Z0", 10.0); ///< Atomic number of impurity (Z)
  const Real ZI = pin->GetOrAddReal("Plasma", "ZI", 1.0);  ///< Charge of impurity
  const Real fI = pin->GetOrAddReal("Plasma", "fI", 1.0);  ///<Fraction of impurity density, normalized to deuterium denstiy (nD0)
  const Real nI = fI*nD0; ///< Impurity density [m^-3]
  const Real n_e0 = nD0 + ZI*nI; ///< Free electron density [m^-3]
  const Real Zeff = pin->GetOrAddReal("Plasma", "Zeff", (ZI*ZI*nI + nD0)/n_e0); // TODO: why ZI is a square?
  const Real NeI = Z0 - ZI; ///< Number of bound electrons
  const Real Coulog0 = pin->GetOrAddReal("Plasma", "Coulog0", 14.9 - 0.5*log(n_e0/1.0e20) + log(Te0/1.e3));
  const Real Rc = pin->GetOrAddReal("Plasma", "Rc", 3.1158966549999998e+00); ///< Initial guess for magnetic axis, R [-], length normalized
  const Real Zc = pin->GetOrAddReal("Plasma", "Zc", 3.7114360000000002e-01); ///< Initial guess for magnetic axis, Z [-], length normalized


  const Real L11 = 0.58 * 32.0 / (3.0 * M_PI);
  const Real sigmapar = 12.0 * pow(M_PI, 1.5) / sqrt(2.0) * pow(Te0 * e, 1.5) * pow(eps0, 2) / (Zeff * pow(e, 2) * sqrt(me) * Coulog0) * L11;
  const Real etaplasma = pin->GetOrAddReal("Plasma", "etaplasma", 1.0 / sigmapar);

  ///< Geometry parameters
  const Real etawall               = pin->GetOrAddReal("Geometry", "etawall", 4.4e-2);                     ///< Wall resistivity [Ohm*m]
  const Real etawallperp            = pin->GetOrAddReal("Geometry", "etawallperp", etawall);
  const Real etawallphi             = pin->GetOrAddReal("Geometry", "etawallphi", etawall);
  const Real etawallphi_isol_cell   = pin->GetOrAddReal("Geometry", "etawallphi_isol_cell", etawall);
  const Real etasepwal              = pin->GetOrAddReal("Geometry", "etasepwal",  etaplasma);
  const Real etaVV                  = pin->GetOrAddReal("Geometry", "etaVV",  1.30288e-6);
  const Real etaout                 = pin->GetOrAddReal("Geometry", "etaout",  1.30288e-3);

  const Real Rmin                   = pin->GetOrAddReal("Geometry", "rmin",  1.525); ///< Minimum R [-]
  const Real Rmax                   = pin->GetOrAddReal("Geometry", "rmax",  4.975); ///< Maximum R [-]
  const Real Zmin                   = pin->GetOrAddReal("Geometry", "zmin",  -2.975);///<  Minimum Z [-]
  const Real Zmax                   = pin->GetOrAddReal("Geometry", "zmax",   2.975); ///< Maximum Z [-]

  const Real gamma_BC  = pin->GetOrAddReal("BoundaryConditions", "gamma_BC", 1. + 1e-4);   // Assuming E max = 1000
  const Real gamma_min = pin->GetOrAddReal("BoundaryConditions", "gamma_min", 2. * (gamma_BC - 1.) + 1.0);
  const Real p_RE      = 0.0;  // will be computed from maximum electric field
  const Real p_BC      = momentum_(gamma_BC);


  ///< Runaway parameters
  const Real c_vTe = pin->GetOrAddReal("Collisions", "c_vTe", c / vTe); ///< Guiding center equations coefficient [-]
  const int NSA = pin->GetOrAddInteger("Collisions", "NSA", 150);       ///< Number of small angle collisions     [-]
  const Real k = pin->GetOrAddReal("Collisions", "k", 5.0);
  const Real aI            = pin->GetOrAddReal("Collisions", "aI", 0.3285296762792767);  ///<
  const Real FineStructure = 1. / 137.035999;  // Fine Structure constant
  const Real II            = pin->GetOrAddReal("Collisions", "II", 235.2 / pc::eV / pc::me / pc::c / pc::c); // Mean exitation energy
  Real PSCoefDnRA = 1.0;
  if (PartialScreening)
    PSCoefDnRA    = 1.0 + NeI * fI / (1.0 + ZI * fI);

  ///< Numerical paremters
  const Real dampV                  = pin->GetOrAddReal("Numerical", "dampV", 0.01); ///< Stabilization coefficeint for velocity gradient
  const Real itime                  = pin->GetOrAddReal("Numerical", "itime", 0.0); ///< Initial time for mhd counters [sec]
  const int NR                      = pin->GetOrAddInteger("Numerical", "NR", 100);
  const int Nphi                    = pin->GetOrAddInteger("Numerical", "Nphi", 2);
  const int NZ                      = pin->GetOrAddInteger("Numerical", "NZ", 200);

  const Real dR = (Rmax - Rmin) / (Real) NR;
  const Real dZ = (Zmax - Zmin) / (Real) NZ;

  const Real RminCellCenter = Rmin + .5 * dR;
  const Real RmaxCellCenter = Rmax - .5 * dR;
  const Real ZminCellCenter = Zmin + .5 * dZ;
  const Real ZmaxCellCenter = Zmax - .5 * dZ;

  const Real tau_a = 6*M_PI*eps0*pow(me * c, 3) / pow(e,4) / pow(B0,2);     ///< Syncrotron radiation damping time
  const Real tau_c = 4*M_PI*pow(eps0,2)*me*me*c*c*c/(e*e*e*e*n_e0*Coulog0); ///< Relativistic collision time
  const Real Ec = me * c / e / tau_c;                                       ///< Connor-Hastie Electric field
  const Real En = E0 / Ec;
  const Real eta_norm = etaplasma / eta0; // converts eta * \curl B to V_A B_0
  const Real eta_a3VaB0 = etaplasma / pow(a,3) / E0; // converts eta J to V_A B_0

  // Create plotting mesh for interpolated fields
  const int NR_plot = pin->GetOrAddInteger("Output", "NR_plot", 400);
  const int NZ_plot = pin->GetOrAddInteger("Output", "NZ_plot", 800);

  Kokkos::DualView<Real***> Jre_mhd("Jre_mhd", NR, NZ, 3);
  Kokkos::deep_copy(Jre_mhd.view_device(), 0.0);
  Kokkos::deep_copy(Jre_mhd.view_host(), 0.0);

  if (mhd_context != nullptr) {
    /// Initialize MHD context
    mhd_context->mi                     = mi;
    mhd_context->mu0                    = mu0;

    mhd_context->density                = nD0;
    mhd_context->B0                     = B0;
    mhd_context->L0                     = a;
    mhd_context->V_A                    = VA;
    mhd_context->eta0                   = eta0;
    mhd_context->eta                    = eta;
    mhd_context->etawall                = etawall;
    mhd_context->etaplasma              = etaplasma;
    mhd_context->etawallperp            = etawallperp;
    mhd_context->etawallphi             = etawallphi;
    mhd_context->etawallphi_isol_cell   = etawallphi_isol_cell;
    mhd_context->etasepwal              = etasepwal;
    mhd_context->etaVV                  = etaVV;
    mhd_context->etaout                 = etaout;
    mhd_context->dampV                  = dampV;
    mhd_context->rmin                   = Rmin * a;
    mhd_context->rmax                   = Rmax * a;
    mhd_context->phimin                 = 0.0;
    mhd_context->phimax                 = 2.0 * M_PI;
    mhd_context->zmin                   = Zmin * a;
    mhd_context->zmax                   = Zmax * a;;
    mhd_context->dt                     = dt_mhd / tauA;
    mhd_context->ictype                 = 9;
    mhd_context->Nr                     = NR;
    mhd_context->Nphi                   = Nphi;
    mhd_context->Nz                     = NZ;
    mhd_context->Re                     = Re;
    mhd_context->itime                  = itime / tauA;
    mhd_context->ftime                  = final_time / tauA;
    mhd_context->phibtype               = pin->GetOrAddInteger("MHD_Config", "phibtype",  1);
    mhd_context->dr                     = dR;
    mhd_context->dphi                   = (mhd_context->phimax - mhd_context->phimin) / mhd_context->Nphi;
    mhd_context->dz                     = dZ;
    mhd_context->tstype                 = pin->GetOrAddInteger("MHD_Config", "tstype",  2);
    mhd_context->jtype                  = pin->GetOrAddInteger("MHD_Config", "jtype",  2);
    mhd_context->debug                  = pin->GetOrAddInteger("MHD_Config", "debug",  0);
    mhd_context->dump                   = pin->GetOrAddInteger("MHD_Config", "dump",  0);
    mhd_context->savecoords             = pin->GetOrAddInteger("MHD_Config", "savecoords",  0);
    mhd_context->isB                    = NULL;
    mhd_context->isEP                   = NULL;
    mhd_context->istau                  = NULL;
    mhd_context->isV                    = NULL;
    mhd_context->isni                   = NULL;
    mhd_context->isB_boundary           = NULL;
    mhd_context->isE_boundary           = NULL;
    mhd_context->isni_boundary          = NULL;


    // Set default location for input data.
    strcpy(mhd_context->input_folder, pin->GetOrAddString("MHD_Config", "input_folder", "../../inputs/mhd").c_str());
    strcpy(mhd_context->ic_binary_path, pin->GetOrAddString("MHD_Config", "ic_binary_path", "").c_str());
    mhd_context->ic_binary_mode = pin->GetOrAddInteger("MHD_Config", "ic_binary_load", 1) == 1 ? 'l' : 'c';

    mhd_context->Ebc = 0;
    mhd_context->tempdump = 0;
    mhd_context->dumpfreq = std::ceil(mhd_context->ftime / (10.0 * mhd_context->dt));
    mhd_context->testSpGD = 0;
    mhd_context->testSpGDsamerhs = 0;
    mhd_context->oldstep = 0;
    mhd_context->n_record =0;
    mhd_context->n_record_Steady_jRE = 0;

    mhd_context -> jre = wrap_view(Jre_mhd.view_host());

    mhd_initialize(mhd_context);

    /// Do single step before seeding particles
    mhd_step(mhd_context);

  }

  auto pkg = std::make_shared<StateDescriptor>("Deck");

  pkg->AddParam("NR", NR);
  pkg->AddParam("NZ", NZ);

  pkg->AddParam("nCorrectorSteps",  nCorrectorSteps);
  pkg->AddParam("nCD",  nCD);

  pkg->AddParam("tau_c",  tau_c);
  pkg->AddParam("eta_norm",  eta_norm);
  pkg->AddParam("eta_a3VaB0",  eta_a3VaB0);
  pkg->AddParam("En",  En);

  int nLA = std::ceil((dt_mhd / nCD) / (dtLA_over_tauC * tau_c));
  Real my_dtLA_over_tauC = (dt_mhd / nCD / tau_c) / nLA;

  if(Globals::my_rank == 0) std::cout <<
    std::format("Adjusting dtLA to evenly devide RE current deposition step:\n {:g} -> {:g} x {:g} sec\n", dtLA_over_tauC, my_dtLA_over_tauC, tau_c);

  pkg->AddParam("dtLA_over_tauC", my_dtLA_over_tauC);
  pkg->AddParam("nLA", nLA);


  pkg->AddParam("filePath", filePath);
  if (Globals::my_rank == 0) {
    std::ofstream(pkg->Param<std::string>("filePath"));
  }

  pkg->AddParam("gamma_min", gamma_min);
  pkg->AddParam("p_BC", p_BC, Params::Mutability::Mutable);
  pkg->AddParam("p_RE", p_RE, Params::Mutability::Mutable);

  pkg->AddParam("hRK", timeStep);
  pkg->AddParam("atol", atol);
  pkg->AddParam("rtol", rtol);

  pkg->AddParam("Rmin", Rmin);
  pkg->AddParam("Rmax", Rmax);
  pkg->AddParam("Zmin", Zmin);
  pkg->AddParam("Zmax", Zmax);

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

  // Initialize field data structures and geometry
  FieldData_t data(NR, NZ, RminCellCenter, ZminCellCenter, dR, dZ);
  pkg->AddParam("FieldData", data);

  // Initialize geometry
  const std::string configurationdomain_file = pin->GetOrAddString("Geometry", "input_file", "../../inputs/AxisSymmetricGeometry.dat");

  ConfigurationDomainGeometry::IndicatorViewType indicator("indicator", NR, NZ);
  std::ifstream ifs(configurationdomain_file);
  auto indicator_h = Kokkos::create_mirror_view(indicator);
  for (int i = 0; i < NR; ++i) {
    for (int j = 0; j < NZ; ++j) {
      ifs >> indicator_h(i,j);
    }
  }
  Kokkos::deep_copy(indicator, indicator_h);
  ConfigurationDomainGeometry cdg(data.fd_locator, indicator);
  pkg->AddParam("CDG", cdg);

  pkg->AddParam("Jre_mhd", Jre_mhd);
  pkg->AddParam("Jre_mhd", E_base);

  CommunicateBJV(data, mhd_context, FieldComponents::B);
  ComputeBaseElectricField_in_place(data, En, eta_norm, FieldsComponents::B);
  CommunicateBJV(data, mhd_context, FieldComponents::Bt);
  ComputeBaseElectricField_in_place(data, En, eta_norm, FieldsComponents::Bt);

  InterpolateTime(data, dtField);

  auto data_d = data.data.device_view();
  data.data.sync_device();
  Kokkos::parallel_for("Set intial runaway current",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0}, {NR, NZ}),
      KOKKOS_LAMBDA(const int i, const int j) {
        Real mask = indicator(i,j) > 0 > Real(1.0) : Real(0.0);
        Real scale = mask * current_fraction * eta_norm;
        data_d(i,j,FieldComponents::E + 0) -= scale * data.data(i, j, FieldComponents::J + 0);
        data_d(i,j,FieldComponents::E + 1) -= scale * data.data(i, j, FieldComponents::J + 1);
        data_d(i,j,FieldComponents::E + 2) -= scale * data.data(i, j, FieldComponents::J + 2);
      });
  data.data.modify_device();

  View3 E_base("E_base", NR, NZ);
  Kokkos::deep_copy(E_base, 0.0);
  pkg->AddParam("E_base", Jre);

  View3 Jre("Jre", NR, NZ);
  Kokkos::deep_copy(Jre, 0.0);
  pkg->AddParam("Jre", Jre);

  InterpolateHermiteBJE(data, FieldComponents::B);
  InterpolateHermiteBJE(data, FieldComponents::Bt);

  const Real R0_plot = f.hR0 + 1e-10;
  const Real Z0_plot = f.hZ0 + 1e-10;
  const Real dR_plot = (f.nR_hermite_data * f.hR - 2e-10) / static_cast<Real> (NR_plot);
  const Real dZ_plot = (f.nZ_hermite_data * f.hZ - 2e-10) / static_cast<Real> (NZ_plot);

  ParArray1D<Real> hpd_R("Hermite_Field_Plot_data_R", NR_plot);
  ParArray1D<Real> hpd_Z("Hermite_Field_Plot_data_Z", NZ_plot);
  ParArray3D<Real> gce_data("GCE_data", NR_plot, NZ_plot, 5);
  ParArrayND<Real> hpd_F("Hermite_Field_Plot_data_F", NR_plot, NZ_plot, 3, 6);
  ParArrayND<Real> hpd_eval("Hermite_Field_Plot_data_eval", NR_plot, NZ_plot, 3,
      static_cast<size_t>(fid::Count));

  Kokkos::parallel_for("FillGrids", NR_plot,
      KOKKOS_LAMBDA(const int n) {
        hpd_R(n) = R0_plot + n * dR_plot;
      });
  Kokkos::parallel_for("FillGrids", NZ_plot,
      KOKKOS_LAMBDA(const int n) {
        hpd_Z(n) = Z0_plot + n * dZ_plot;
      });

  pkg->AddParam("Hermite_Field_Plot_data_R", hpd_R);
  pkg->AddParam("Hermite_Field_Plot_data_Z", hpd_Z);
  pkg->AddParam("Hermite_Field_Plot_data_F", hpd_F);
  pkg->AddParam("Hermite_Field_Plot_data_eval", hpd_eval);
  pkg->AddParam("GCE_data", gce_data);

  const Real wce0 = pc::qe * B0 / pc::me; // Electron gyrofrequency
  const Real c_aw0 =  pin->GetOrAddReal("GuidingCenterEquations", "c_aw0", pc::c/a/wce0);
  const Real ct_a =   pin->GetOrAddReal("GuidingCenterEquations", "ct_a", pc::c * tau_c / a);
  const Real alpha0 = pin->GetOrAddReal("GuidingCenterEquations", "alpha0", tau_c/tau_a);

  pkg->AddParam("c_aw0", c_aw0);
  pkg->AddParam("ct_a", ct_a);
  pkg->AddParam("alpha0", alpha0);

  int npart =  pin->GetOrAddInteger("ParticleSeed", "num_particles_per_block", 16);
  pkg->AddParam("num_particles_per_block", npart);

// Paticles are counted on Collect current
  pkg->AddParam("num_particles_total", 0, Params::Mutability::Restart);
  const int num_particles_max = pin->GetOrAddInteger("ParticleSeed", "MaxParticles", 500000);

  pkg->AddParam("num_particles_max", num_particles_max);
  // Initialize random number generator pool
  int rng_seed = pin->GetOrAddInteger("ParticleSeed", "rng_seed", 1234) + Globals::my_rank;
  RNGPool rng_pool(rng_seed);
  pkg->AddParam("rng_pool", rng_pool);

  pkg->AddParam("Rc", Rc);
  pkg->AddParam("Zc", Zc);

  const Real seed_current = pin->GetOrAddReal("ParticleSeed", "current", 15e3); // 15 kAmps
  pkg->AddParam("seed_current", seed_current * a); // Convert from amps

  const Real seed_current_fraction = pin->GetOrAddReal("ParticleSeed", "current_fraction", seed_current / 15e6); // Used to determine the initial runaway current to adjust the Electric field.


  const Real gammamin = pin->GetOrAddReal("ParticleSeed", "gammamin", 10.0);
  pkg->AddParam("pmin", momentum_(gammamin));
  const Real gammamax = pin->GetOrAddReal("ParticleSeed", "gammamax", 20.0);
  pkg->AddParam("pmax", momentum_(gammamax));
  const Real ximin = pin->GetOrAddReal("ParticleSeed", "ximin", 0.8);
  pkg->AddParam("ximin", ximin);
  const Real ximax = pin->GetOrAddReal("ParticleSeed", "ximax", 1.0);
  pkg->AddParam("ximax", ximax);

  if (Globals::my_rank == 0) {
    std::ofstream ofs("collision_profiles.dat");
    ofs << std::format("{:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s}",
        "#     p", "gamma", "dtSA", "psi", "CB", "CF", "CouLogee ratio", "probability");
    Real p = momentum_(1. + 2.e-3);
    while (p < pkg->Param<Real>("pmax") + 20.0) {
      auto cc = sa.getCollisionCoefficients(p);
      ofs << std::format("{:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e}",
          p, gamma_(p),
          sa.getSmallAngleCollisionTimestep(p),
          cc.psi, cc.CB, cc.CF, cc.CouLogee_ratio,
          ms.computeProbability(p, 1.0, dtLA_over_tauC, 1.002)
      ) << std::endl;

      p += 1e-2;
    }
  }

  Metadata swarm_metadata({Metadata::Provides, Metadata::None, Metadata::Restart});
  pkg->AddSwarm("particles", swarm_metadata);

  Metadata real_swarmvalue_metadata({Metadata::Real, Metadata::Restart});
  pkg->AddSwarmValue(Kinetic::p::name(), "particles", real_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::xi::name(), "particles", real_swarmvalue_metadata);
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


  Metadata int_swarmvalue_metadata({Metadata::Integer, Metadata::Restart});
  pkg->AddSwarmValue(Kinetic::will_scatter::name(), "particles",
                     int_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::secondary_index::name(), "particles",
                     int_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::status::name(), "particles",
                     int_swarmvalue_metadata);

  int EnableLargeAngleCollisions = pin->GetOrAddInteger("Simulation", "EnableLargeAngleCollisions", 1);
  int EnableSmallAngleCollisions = pin->GetOrAddInteger("Simulation", "EnableSmallAngleCollisions", 1);
  int EnableComputeConservedQuantities = pin->GetOrAddInteger("Simulation", "EnableComputeConservedQuantities", 0);

  pkg->AddParam("EnableLargeAngleCollisions",      EnableLargeAngleCollisions);
  pkg->AddParam("EnableSmallAngleCollisions",      EnableSmallAngleCollisions);
  pkg->AddParam("EnableComputeConservedQuantities",EnableComputeConservedQuantities);

	pkg->AddSwarmValue(p_phi::name(), "particles",real_swarmvalue_metadata);
  pkg->AddSwarmValue(mu::name(), "particles",real_swarmvalue_metadata);

  // This filename is updated with OutputParameters.file_number and file_basename before each restart.
  // It is required to remember the Petsc file name associated with Parthenon restart.
  // Petsc file is where the mhd state is stored, it is separate from Parthenon restart file
  std::string mhd_restart_filename = "";
  pkg->AddParam("mhd_restart_filename", mhd_restart_filename, Params::Mutability::Restart);
  pkg->AddParam("ComputeInitialWeights", 1, Params::Mutability::Restart);
  pkg->AddParam("MyTime", 0.0, Params::Mutability::Restart);
  pkg->AddParam("PredictorIterationNumber", 0, Params::Mutability::Restart);

  if (Globals::my_rank == 0) std::cout << "Init finished\n";

  return pkg;
}

std::shared_ptr<StateDescriptor> InitializeAnalytic(ParameterInput *pin) {
  /// Physical constants
  static constexpr Real eps0 = pc::eps0; ///< Vacuum permittivity [F / M]
  static constexpr Real c = pc::c;       ///< Speed of light [m/s]
  static constexpr Real mi = pc::amu;    ///< Ion mass = atomic mass unit [kg]
  static constexpr Real mu0 = pc::mu0;   ///< Vacuum permeabulity [N / A^2]
  static constexpr Real me = pc::me;     ///< electron mass [kg]
  static constexpr Real e = pc::qe;      ///< electron charge [C]

  /// Time discretization parameters
  const Real dt_mhd     = pin->GetReal("parthenon/time","dt_force");
  const Real final_time = pin->GetReal("parthenon/time","tlim");

  const int nCorrectorSteps = pin->GetOrAddInteger("Time","nCorrectorSteps", 0);  ///< current deposit timestep for electric field readjustment [s]
  const int nCD = pin->GetOrAddInteger("Time","nCD", 100);  ///< current deposit timestep for electric field readjustment [s]
  Real dtLA_over_tauC = pin->GetOrAddReal("Time","dtLA_over_tauC", 1e-5);  ///< current deposit timestep for electric field readjustment [s]

  const Real timeStep = pin->GetOrAddReal("Simulation", "hRK", 1.e-6);    /// Runge kutta time in tau_c [-]
  const Real atol = pin->GetOrAddReal("Simulation", "atol", 1.e-6);      /// Absoulte tolerance for RK [-]
  const Real rtol = pin->GetOrAddReal("Simulation", "rtol", 1.e-5);       /// Realative toleratnce for RK[ [-]
  const std::string filePath = pin->GetOrAddString("Simulation", "file_path", "current.out");

  /// Reference parameters
  const Real B0 = pin->GetOrAddReal("Reference", "B0", 5.3);  ///< On-axis magnetic field [T]
  const Real a  = pin->GetOrAddReal("Reference", "a", 2.0);    ///< Minor radius [m] and reference length
  const Real R0 = pin->GetOrAddReal("Reference", "R0", 6.0);    ///< Major radius [m]
  const Real nD0 = pin->GetOrAddReal("Reference", "nD0", 1e20);    ///< Deutirium density [m^-3]
  const Real Te0 = pin->GetOrAddReal("Reference", "Te0", 2.0198);  ///< Electron temperature [eV], and plasma temperature single-temperature model

  /// Derived parameters
  const Real VA   = pin->GetOrAddReal("Derived", "VA",  B0 / sqrt(mi * nD0 * mu0)); ///< Alfven velocity [m/s]
  const Real tauA = pin->GetOrAddReal("Derived", "tauA", a / VA);                  ///< Alfven time     [s]
  const Real E0   = pin->GetOrAddReal("Derived", "E0", B0 * VA);                   ///< Reference electric field in MHD [V/m]
  const Real J0   = pin->GetOrAddReal("Derived", "J0", B0 / (mu0 * a));            ///< Reference current density [A/m^2]
  const Real eta0 = pin->GetOrAddReal("Derived", "eta0", a * VA * mu0);            ///< Reference resitivity [Ohm*m]
  const Real eta  = pin->GetOrAddReal("Derived", "eta",  1.0);                     ///< Resitivity scale     [-]
  const Real Re   = pin->GetOrAddReal("Derived", "Re",  200.0);                    ///< Reinolds Number

  ///< Plasma composition parameters
  const Real vTe = sqrt(2.0*Te0*e / me); ///< Thermal velocity
  const Real Z0 = pin->GetOrAddReal("Plasma", "Z0", 10.0); ///< Atomic number of impurity (Z)
  const Real ZI = pin->GetOrAddReal("Plasma", "ZI", 1.0);  ///< Charge of impurity
  const Real fI = pin->GetOrAddReal("Plasma", "fI", 100.0);  ///<Fraction of impurity density, normalized to deuterium denstiy (nD0)
  const Real nI = fI*nD0; ///< Impurity density [m^-3]
  const Real n_e0 = nD0 + ZI*nI; ///< Free electron density [m^-3]
  const Real Zeff = pin->GetOrAddReal("Plasma", "Zeff", (ZI*ZI*nI + nD0)/n_e0); // TODO: why ZI is a square?
  const Real NeI = Z0 - ZI; ///< Number of bound electrons
  const Real Coulog0 = pin->GetOrAddReal("Plasma", "Coulog0", 14.9 - 0.5*log(n_e0/1.0e20) + log(Te0/1.e3));
  const Real Rc = pin->GetOrAddReal("Plasma", "Rc", 3.1158966549999998e+00); ///< Initial guess for magnetic axis, R [-], length normalized
  const Real Zc = pin->GetOrAddReal("Plasma", "Zc", 3.7114360000000002e-01); ///< Initial guess for magnetic axis, Z [-], length normalized


  const Real L11 = 0.58 * 32.0 / (3.0 * M_PI);
  const Real sigmapar = 12.0 * pow(M_PI, 1.5) / sqrt(2.0) * pow(Te0 * e, 1.5) * pow(eps0, 2) / (Zeff * pow(e, 2) * sqrt(me) * Coulog0) * L11;
  const Real etaplasma = pin->GetOrAddReal("Plasma", "etaplasma", 1.0 / sigmapar);

  ///< Geometry parameters
  const Real etawall               = pin->GetOrAddReal("Geometry", "etawall", 4.4e-2);                     ///< Wall resistivity [Ohm*m]
  const Real etawallperp            = pin->GetOrAddReal("Geometry", "etawallperp", etawall);
  const Real etawallphi             = pin->GetOrAddReal("Geometry", "etawallphi", etawall);
  const Real etawallphi_isol_cell   = pin->GetOrAddReal("Geometry", "etawallphi_isol_cell", etawall);
  const Real etasepwal              = pin->GetOrAddReal("Geometry", "etasepwal",  etaplasma);
  const Real etaVV                  = pin->GetOrAddReal("Geometry", "etaVV",  1.30288e-6);
  const Real etaout                 = pin->GetOrAddReal("Geometry", "etaout",  1.30288e-3);

  const Real Rmin                   = pin->GetOrAddReal("Geometry", "rmin",  1.525); ///< Minimum R [-]
  const Real Rmax                   = pin->GetOrAddReal("Geometry", "rmax",  4.975); ///< Maximum R [-]
  const Real Zmin                   = pin->GetOrAddReal("Geometry", "zmin",  -2.975);///<  Minimum Z [-]
  const Real Zmax                   = pin->GetOrAddReal("Geometry", "zmax",   2.975); ///< Maximum Z [-]

  const Real gamma_min      = pin->GetOrAddReal("BoundaryConditions", "gamma_min",1.02);
  const Real p_BC = momentum_(pin->GetOrAddReal("BoundaryConditions", "gamma_BC", 1.02));
  const Real p_RE = momentum_(pin->GetOrAddReal("BoundaryConditions", "gamma_RE", 1.02));

  ///< Runaway parameters
  const Real c_vTe = pin->GetOrAddReal("Collisions", "c_vTe", c / vTe); ///< Guiding center equations coefficient [-]
  const int NSA = pin->GetOrAddInteger("Collisions", "NSA", 150);       ///< Number of small angle collisions     [-]
  const Real k = pin->GetOrAddReal("Collisions", "k", 5.0);
  const Real aI            = pin->GetOrAddReal("Collisions", "aI", 0.3285296762792767);  ///<
  const Real FineStructure = 1. / 137.035999;  // Fine Structure constant
  const Real II            = pin->GetOrAddReal("Collisions", "II", 219.5 / pc::eV / pc::me / pc::c / pc::c); // Mean exitation energy

  Real PSCoefDnRA = 1.0;
  if (PartialScreening)
    PSCoefDnRA    = 1.0 + NeI * fI / (1.0 + ZI * fI);

  ///< Numerical paremters
  const Real dampV                  = pin->GetOrAddReal("Numerical", "dampV", 0.01); ///< Stabilization coefficeint for velocity gradient
  const Real itime                  = pin->GetOrAddReal("Numerical", "itime", 0.0); ///< Initial time for mhd counters [sec]
  const int NR                      = pin->GetOrAddInteger("Numerical", "NR", 100);
  const int Nphi                    = pin->GetOrAddInteger("Numerical", "Nphi", 2);
  const int NZ                      = pin->GetOrAddInteger("Numerical", "NZ", 200);

  const Real dR = (Rmax - Rmin) / (Real) NR;
  const Real dZ = (Zmax - Zmin) / (Real) NZ;

  const Real RminCellCenter = Rmin + .5 * dR;
  const Real RmaxCellCenter = Rmax - .5 * dR;
  const Real ZminCellCenter = Zmin + .5 * dZ;
  const Real ZmaxCellCenter = Zmax - .5 * dZ;

  const Real tau_a = 6*M_PI*eps0*pow(me * c, 3) / pow(e,4) / pow(B0,2);     ///< Syncrotron radiation damping time
  const Real tau_c = 4*M_PI*pow(eps0,2)*me*me*c*c*c/(e*e*e*e*n_e0*Coulog0); ///< Relativistic collision time
  const Real Ec = me * c / e / tau_c;                                       ///< Connor-Hastie Electric field
  const Real En = E0 / Ec;
  const Real eta_norm = etaplasma / eta0; // converts eta * \curl B to V_A B_0
  const Real eta_a3VaB0 = etaplasma / pow(a,3) / E0; // converts eta J to V_A B_0

  // Create plotting mesh for interpolated fields
  const int NR_plot = pin->GetOrAddInteger("Output", "NR_plot", 400);
  const int NZ_plot = pin->GetOrAddInteger("Output", "NZ_plot", 800);

  auto pkg = std::make_shared<StateDescriptor>("Deck");


  pkg->AddParam("nCorrectorSteps",  nCorrectorSteps);
  pkg->AddParam("nCD",  nCD);

  pkg->AddParam("tau_c",  tau_c);
  pkg->AddParam("eta_norm",  eta_norm);
  pkg->AddParam("eta_a3VaB0",  eta_a3VaB0);
  pkg->AddParam("En",  En);

  int nLA = std::ceil((dt_mhd / nCD) / (dtLA_over_tauC * tau_c));
  Real my_dtLA_over_tauC = (dt_mhd / nCD / tau_c) / nLA;

  if(Globals::my_rank == 0) std::cout <<
    std::format("Adjusting dtLA to evenly devide RE current deposition step:\n {:g} -> {:g} x {:g} sec\n", dtLA_over_tauC, my_dtLA_over_tauC, tau_c);

  pkg->AddParam("dtLA_over_tauC", my_dtLA_over_tauC);
  pkg->AddParam("nLA", nLA);


  pkg->AddParam("filePath", filePath);
  if (Globals::my_rank == 0) {
    std::ofstream(pkg->Param<std::string>("filePath"));
  }

  pkg->AddParam("gamma_min", gamma_min);
  pkg->AddParam("p_BC", p_BC);
  pkg->AddParam("p_RE", p_RE);

  pkg->AddParam("hRK", timeStep);
  pkg->AddParam("atol", atol);
  pkg->AddParam("rtol", rtol);

  pkg->AddParam("Rmin", Rmin);
  pkg->AddParam("Rmax", Rmax);
  pkg->AddParam("Zmin", Zmin);
  pkg->AddParam("Zmax", Zmax);

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

  FreeGeometry fg;
  pkg->AddParam("CDG", fg);

  const Real q0 = pin->GetOrAddReal("AnalyticField", "q0", 2.1);
  const Real q2 = pin->GetOrAddReal("AnalyticField", "q2", 2.0);
  const Real R_a = pin->GetOrAddReal("AnalyticField", "R_a", 3.0);
  const Real E_0 = pin->GetOrAddReal("AnalyticField", "E_0", 5.0);

  AnalyticField f(q0, q2, R_a, E_0);
  pkg->AddParam("Field", f);

  const Real wce0 = pc::qe * B0 / pc::me; // Electron gyrofrequency
  const Real c_aw0 =  pin->GetOrAddReal("GuidingCenterEquations", "c_aw0", pc::c/a/wce0);
  const Real ct_a =   pin->GetOrAddReal("GuidingCenterEquations", "ct_a", pc::c * tau_c / a);
  const Real alpha0 = pin->GetOrAddReal("GuidingCenterEquations", "alpha0", tau_c/tau_a);

  pkg->AddParam("c_aw0", c_aw0);
  pkg->AddParam("ct_a", ct_a);
  pkg->AddParam("alpha0", alpha0);

  int npart =  pin->GetOrAddInteger("ParticleSeed", "num_particles_per_block", 16);

  pkg->AddParam("num_particles_per_block", npart);
// Paticles are counted on Collect current
  pkg->AddParam("num_particles_total", 0, Params::Mutability::Restart);
  const int num_particles_max = pin->GetOrAddInteger("ParticleSeed", "MaxParticles", 4096 * 5);

  pkg->AddParam("num_particles_max", num_particles_max);
  // Initialize random number generator pool
  int rng_seed = pin->GetOrAddInteger("ParticleSeed", "rng_seed", 1234) + Globals::my_rank;
  RNGPool rng_pool(rng_seed);
  pkg->AddParam("rng_pool", rng_pool);

  pkg->AddParam("Rc", Rc);
  pkg->AddParam("Zc", Zc);

  const Real seed_current = pin->GetOrAddReal("ParticleSeed", "current", 15e3); // 15 kAmps
  pkg->AddParam("seed_current", seed_current * a); // Convert from amps
  const Real gammamin = pin->GetOrAddReal("ParticleSeed", "gammamin", 10.0);
  pkg->AddParam("pmin", momentum_(gammamin));
  const Real gammamax = pin->GetOrAddReal("ParticleSeed", "gammamax", 20.0);
  pkg->AddParam("pmax", momentum_(gammamax));
  const Real ximin = pin->GetOrAddReal("ParticleSeed", "ximin", 0.8);
  pkg->AddParam("ximin", ximin);
  const Real ximax = pin->GetOrAddReal("ParticleSeed", "ximax", 1.0);
  pkg->AddParam("ximax", ximax);

  if (Globals::my_rank == 0) {
    std::ofstream ofs("collision_profiles.dat");
    ofs << std::format("{:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s}",
        "#     p", "gamma", "dtSA", "psi", "CB", "CF", "CouLogee ratio", "probability");
    Real p = momentum_(1. + 2.e-3);
    while (p < pkg->Param<Real>("pmax") + 20.0) {
      auto cc = sa.getCollisionCoefficients(p);
      ofs << std::format("{:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e}",
          p, gamma_(p),
          sa.getSmallAngleCollisionTimestep(p),
          cc.psi, cc.CB, cc.CF, cc.CouLogee_ratio,
          ms.computeProbability(p, 1.0, dtLA_over_tauC, 1.002)
      ) << std::endl;

      p += 1e-2;
    }
  }

  Metadata swarm_metadata({Metadata::Provides, Metadata::None, Metadata::Restart});
  pkg->AddSwarm("particles", swarm_metadata);

  Metadata real_swarmvalue_metadata({Metadata::Real, Metadata::Restart});
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


  Metadata int_swarmvalue_metadata({Metadata::Integer, Metadata::Restart});
  pkg->AddSwarmValue(Kinetic::will_scatter::name(), "particles",
                     int_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::secondary_index::name(), "particles",
                     int_swarmvalue_metadata);
  pkg->AddSwarmValue(Kinetic::status::name(), "particles",
                     int_swarmvalue_metadata);

  int EnableLargeAngleCollisions = pin->GetOrAddInteger("Simulation", "EnableLargeAngleCollisions", 1);
  int EnableSmallAngleCollisions = pin->GetOrAddInteger("Simulation", "EnableSmallAngleCollisions", 1);
  int EnableComputeConservedQuantities = pin->GetOrAddInteger("Simulation", "EnableComputeConservedQuantities", 0);

  pkg->AddParam("EnableLargeAngleCollisions",      EnableLargeAngleCollisions);
  pkg->AddParam("EnableSmallAngleCollisions",      EnableSmallAngleCollisions);
  pkg->AddParam("EnableComputeConservedQuantities",EnableComputeConservedQuantities);

	pkg->AddSwarmValue(p_phi::name(), "particles",real_swarmvalue_metadata);
  pkg->AddSwarmValue(mu::name(), "particles",real_swarmvalue_metadata);

  if (Globals::my_rank == 0) std::cout << "Init finished\n";
  return pkg;
}


void WorkBeforeOutput(Mesh * pm, ParameterInput * pin, SimTime const & tm, User* mhd_context) {
  auto pkg = pm->packages.Get("Deck");

  auto hpd_R = pkg->Param<ParArray1D<Real>>("Hermite_Field_Plot_data_R");
  auto hpd_Z = pkg->Param<ParArray1D<Real>>("Hermite_Field_Plot_data_Z");
  auto hpd_F = pkg->Param<ParArrayND<Real>>("Hermite_Field_Plot_data_F");
  auto hpd_eval = pkg->Param<ParArrayND<Real>>("Hermite_Field_Plot_data_eval").KokkosView();
  auto gce_data = pkg->Param<ParArray3D<Real>>("GCE_data");
  const auto pmin  = pkg->Param<Real>("pmin");
  const auto pmax  = pkg->Param<Real>("pmax");
  const auto ximin = pkg->Param<Real>("ximin");
  const auto ximax = pkg->Param<Real>("ximax");

  const Real p0 = .5 * (pmax + pmin);
  const Real xi0 = .5 * (ximax + ximin);


  const int NR_plot = hpd_R.size();
  const int NZ_plot = hpd_Z.size();

  const auto c_aw0  = pkg->Param<Real>("c_aw0");
  const auto ct_a   = pkg->Param<Real>("ct_a");
  const auto alpha0 = pkg->Param<Real>("alpha0");
  auto data = pkg->Param<FieldData_t>("FieldData");
  auto cdg  = pkg->Param<ConfigurationDomainGeometry>("ConfigurationDomainGeometry");
  FieldEvaluator f{cdg.locator, data};
  GuidingCenterEquations<FieldEvaluator, true, false> gce(f, c_aw0, ct_a, alpha0);

  // Now plot all Hermite fields
  Kokkos::parallel_for("FillInterpolatedData_plot",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {NR_plot,NZ_plot}),
      KOKKOS_LAMBDA(const int i, const int j) {
        Real R = hpd_R(i);
        Real Z = hpd_Z(j);
        Dim5 X = {p0, xi0, R, 0.0, Z};
        Real t = 0.0;
        EvalGCE ev;
        f.eval(ev, R, Z, t);

        for (int k = 0; k < 3; ++k) {
          hpd_F(i,j,k,0) = ev.B[k];
          hpd_F(i,j,k,1) = ev.curlB[k];
          hpd_F(i,j,k,2) = ev.dBdR[k];
          hpd_F(i,j,k,3) = ev.dBdZ[k];
          hpd_F(i,j,k,4) = ev.E[k];
          hpd_F(i,j,k,5) = ev.dbdt[k];
        }

        Dim5 dX = {};
        gce(0.0, X, dX);
        for (int k = 0; k < 5; ++k) {
          gce_data(i,j,k) = dX[k];
        }

      });

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      swarm_position::x, swarm_position::y, swarm_position::z, Kinetic::p,
      Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i =
      parthenon::MakeSwarmPackDescriptor<Kinetic::will_scatter,
                                         Kinetic::secondary_index,
                                         Kinetic::status>("particles");
  auto md = pm->mesh_data.Get();
  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  Real I_re = 0.0;

  const Real p_RE = pkg->Param<Real>("p_RE");
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
            weight += getParticleCurrent(X, t, w, f);
          }
        }
      },
      I_re);


  auto Jre_mhd = pkg->Param<DualView3>("Jre_mhd");

  Jre_mhd.sync_device();
  Jre_d = Jre_mhd.device_view();

  Reak R0 = cdg.indicator_locator.R0;
  Reak Z0 = cdg.indicator_locator.Z0;
  Reak dR = cdg.indicator_locator.dR;
  Reak dZ = cdg.indicator_locator.dZ;
  Reak nR = cdg.indicator_locator.nR;
  Reak nZ = cdg.indicator_locator.nZ;

  data.data.sync_device();
  data_d = data.data.device_view();

  Real I_ohmic_fd = 0.0;
  Real I_ohmic_hermite = 0.0;
  Real I_re_integral = 0.0;
  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL,
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {nR, nZ}),
      // loop over all particles
      KOKKOS_LAMBDA(int i, int j, Real& integral, Real& integral_ohmic_hermite, Real& integral_ohmic_fd) {
        Real R = R0 + i * dR;
        Real Z = Z0 + j * dZ;

        int ii,jj, region;
        cdg.locate_region(R, Z, ii,jj, region);

        if (region > 0) {

          EvalGCE ev;
          Real t = 0.0;

          auto ret = (ev, R, Z, t);

          Real area = dR * dZ;
          integral_ohmic_fd += area * data_d(i, j, FieldComponents::J + 1);
          integral_ohmic_hermite += area * ev.curlB[1];
          integral += area * Jre_d(i,j,1);
        }

      },
      I_re_integral, I_ohmic_hermite, I_ohmic_fd);

  MPI_Allreduce(MPI_IN_PLACE,&I_re,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if (Globals::my_rank == 0) {
    std::cout << std::format("{:20.14e} {:20.14e} {:20.14e}",
        I_re * .5, I_re_integral * 5.3 * 2.0,
        I_ohmic * 5.3  * 2.0 / pc::mu0) << std::endl;
  }
}

void WorkBeforeRestartOutput(Mesh * pm, ParameterInput * pin, OutputParameters * op, User* mhd_context) {
  auto signal = SignalHandler::CheckSignalFlags();
  std::string ext = pin->GetOrAddString("MHD_Config", "file_extention", "dat");

  std::string filename = std::format("{}.{}.",
      op -> file_basename,
      op -> file_id);

  if (signal == SignalHandler::OutputSignal::now) {
    filename.append("now");
  } else if (signal == SignalHandler::OutputSignal::final && op -> file_label_final) {
    filename.append("final");
    // default time based data dump
  } else {
    filename.append(std::format("{:0>{}}", op->file_number, op -> file_number_width));
  }

  filename.append(".r");
  filename.append(ext);

  mhd_savesolution(mhd_context, filename.c_str());

  auto pkg = pm->packages.Get("Deck");
  pkg->UpdateParam("mhd_restart_filename", filename);
}

void WorkBeforeLoop(Mesh * pm, User* mhd_context) {
  if (Globals::is_restart) {
    if (Globals::my_rank == 0) std::cout << "WorkBeforeLoop: reading mhd restart\n";
    auto pkg = pm->packages.Get("Deck");
    auto filename = pkg -> Param<std::string>("mhd_restart_filename");

    mhd_loadsolution(mhd_context, filename.c_str());
  }
}

auto &GetCoords(std::shared_ptr<MeshBlock> &pmb) { return pmb->coords; }
auto &GetCoords(MeshBlock *pmb) { return pmb->coords; }
auto &GetCoords(Mesh *pm) { return pm->block_list[0]->coords; }

TaskStatus SaveState(Mesh* pm) {
  std::cout << "SaveState start" << std::endl;
  auto md = pm->mesh_data.Get();
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
          pack_swarm_r(b_r, Kinetic::saved_p(), n_r)   = pack_swarm_r(b_r, Kinetic::p(), n_r);
          pack_swarm_r(b_r, Kinetic::saved_xi(), n_r)  = pack_swarm_r(b_r, Kinetic::xi(), n_r);
          pack_swarm_r(b_r, Kinetic::saved_R(), n_r)   = pack_swarm_r(b_r, Kinetic::R(), n_r) ;
          pack_swarm_r(b_r, Kinetic::saved_phi(), n_r) = pack_swarm_r(b_r, Kinetic::phi(), n_r);
          pack_swarm_r(b_r, Kinetic::saved_Z(), n_r)   = pack_swarm_r(b_r, Kinetic::Z(), n_r)  ;
          pack_swarm_r(b_r, Kinetic::saved_w(), n_r)   = pack_swarm_r(b_r, Kinetic::weight(), n_r)  ;
          pack_swarm_i(b_i, Kinetic::status(), n_i) |= Kinetic::PROTECTED;
        } else {
          const auto swarm = pack_swarm_i.GetContext(b_i);
          pack_swarm_i(b_i, Kinetic::status(), n_i) &= ~Kinetic::PROTECTED;
          swarm.MarkParticleForRemoval(n_i);
        }
      });
	return TaskStatus::complete;
}


TaskStatus RandomRemove(Mesh* pm) {
  std::cout << "Random remove start" << std::endl;
  auto pkg = pm->packages.Get("Deck");
  int num_particles = 0;

  auto desc_swarm = parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");
  auto md = pm->mesh_data.Get();
  auto pack_swarm = desc_swarm.GetPack(md.get());

  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm.GetMaxFlatIndex() + 1,
      // loop over all particles
      KOKKOS_LAMBDA(const int idx, int& num) {
        // block and particle indices
        auto [b, n] = pack_swarm.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm.GetContext(b);
        if (swarm_d.IsActive(n) && !swarm_d.IsMarkedForRemoval(n) && (pack_swarm(b, Kinetic::status(), n) & Kinetic::ALIVE)) {
          num++;
        }
      },
      num_particles);
  MPI_Allreduce(MPI_IN_PLACE,&num_particles,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  std::cout << "Number of particles " << num_particles << std::endl;
  auto num_particles_max = pkg->Param<int>("num_particles_max");
	if (num_particles_max > num_particles) return TaskStatus::complete;

  auto rng_pool = pkg->Param<Kinetic::RNGPool>("rng_pool");
  auto desc_swarm_i = parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());
  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<Kinetic::weight>("particles");
  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_i.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        auto [b_i, n_i] = pack_swarm_i.GetBlockParticleIndices(idx);
        auto [b_r, n_r] = pack_swarm_r.GetBlockParticleIndices(idx);
        if (pack_swarm_i(b_i, Kinetic::status(), n_i) & Kinetic::ALIVE) {
          auto rng_gen = rng_pool.get_state();
          auto isKilled = (rng_gen.urand(2));
          rng_pool.free_state(rng_gen);
          // block and particle indices

          if (isKilled == 1) {
              pack_swarm_i(b_i, Kinetic::status(), n_i) &= ~Kinetic::ALIVE;
              pack_swarm_i(b_i, Kinetic::status(), n_i) &= ~Kinetic::PROTECTED;
              const auto swarm = pack_swarm_i.GetContext(b_i);
              swarm.MarkParticleForRemoval(n_i);
          } else {
							pack_swarm_r(b_r, Kinetic::weight(), n_r) *= 2.0;
					}
        }
      });

  return TaskStatus::complete;
}

TaskStatus RestoreState(Mesh* pm) {
  auto md = pm->mesh_data.Get();
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
          const auto swarm = pack_swarm_i.GetContext(b_i);
          swarm.MarkParticleForRemoval(n_i);
        }
      });
	return TaskStatus::complete;
}

TaskStatus RestoreJre(Mesh* pm) {
  auto pkg = pm->packages.Get("Deck");

  auto jre_backup = pkg->Param<ParArray3D<Real>>("jre_backup");
  auto f = pkg->Param<EM_Field>("Field");
  auto jre_data = Kokkos::subview(f.data, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, static_cast<size_t>(fid::Jre), 0);

  Kokkos::deep_copy(jre_data, jre_backup);
	return TaskStatus::complete;
}

void ComputeParticleWeights(Mesh* pm) {

  auto md = pm->mesh_data.Get();
  auto pkg = pm->packages.Get("Deck");
  if (pkg->Param<int>("ComputeInitialWeights") == 0) return;
  auto data = pkg->Param<FieldData_t>("FieldData");
  auto cdg  = pkg->Param<ConfigurationDomainGeometry>("ConfigurationDomainGeometry");
  FieldEvaluator f{cdg.hermite_locator, data};

  const Real p_RE = pkg->Param<Real>("p_RE");
  const Real seed_current = pkg->Param<Real>("seed_current");

  auto desc_swarm_r = parthenon::MakeSwarmPackDescriptor<
      Kinetic::p, Kinetic::xi, Kinetic::R, Kinetic::phi, Kinetic::Z, Kinetic::weight>(
      "particles");
  auto desc_swarm_i = parthenon::MakeSwarmPackDescriptor<
      Kinetic::status>(
      "particles");

  Real I_re = 0.0;

  auto pack_swarm_r = desc_swarm_r.GetPack(md.get());
  auto pack_swarm_i = desc_swarm_i.GetPack(md.get());

  if (Globals::my_rank == 0)
    std::cout << "Calculating current: \n";

  Kokkos::parallel_reduce(
      PARTHENON_AUTO_LABEL, pack_swarm_r.GetMaxFlatIndex() + 1,
      // loop over all particles
      KOKKOS_LAMBDA(const int idx, Real &weight) {
        // block and particle indices
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        const auto swarm_d = pack_swarm_r.GetContext(b);
        if (swarm_d.IsActive(n) && (pack_swarm_i(b, Kinetic::status(),n) & Kinetic::ALIVE)) {
          Dim5 X;
          Real t = 0.0;
          X[0] = pack_swarm_r(b, Kinetic::p(), n);
          X[1] = pack_swarm_r(b, Kinetic::xi(), n);
          X[2] = pack_swarm_r(b, Kinetic::R(), n);
          X[3] = pack_swarm_r(b, Kinetic::phi(), n);
          X[4] = pack_swarm_r(b, Kinetic::Z(), n);
          Real w = pack_swarm_r(b, Kinetic::weight(), n);
          if (X[0] > p_RE) {
            weight += getParticleCurrent(X, t, w, f);
          }

        }
      },
      I_re);

  MPI_Allreduce(MPI_IN_PLACE,&I_re,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  Real w = seed_current / I_re;

  if (Globals::my_rank == 0)
    std::cout << std::format("I_re = {:.8E}, w = {:.8E}\n", I_re, w) << std::endl;
  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_r.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        // block and particle indices
        pack_swarm_r(b, Kinetic::weight(), n) = w;
      });

  pkg->UpdateParam("ComputeInitialWeights", 0);
}

TaskStatus MakeOutputs(Outputs* pouts, Mesh* pmesh, ParameterInput* pinput, Real time, int iPR) {
  auto md = pmesh->mesh_data.Get();
  auto pkg = pmesh->packages.Get("Deck");
  pkg->UpdateParam("MyTime", time);
  pkg->UpdateParam("PredictorIterationNumber", iPR);

  pouts->MakeOutputs(pmesh, pinput);
  return TaskStatus::complete;
}

} // namespace Kinetic
