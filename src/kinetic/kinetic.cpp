//========================================================================================
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
#include <typeinfo>  //for 'typeid' to work
#include <parthenon/package.hpp>

using namespace parthenon;

#include "kinetic/RunawayDriver.h"
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
void InitializeMHDConfig(ParameterInput *pin, User* mhd_context) {
  /// Physical constants
  static constexpr Real eps0 = pc::eps0; ///< Vacuum permittivity [F / M]
  static constexpr Real c = pc::c;       ///< Speed of light [m/s]
  static constexpr Real mi = pc::amu;    ///< Ion mass = atomic mass unit [kg]
  static constexpr Real mu0 = pc::mu0;   ///< Vacuum permeabulity [N / A^2]
  static constexpr Real me = pc::me;     ///< electron mass [kg]
  static constexpr Real e = pc::qe;      ///< electron charge [C]

  /// Time discretization parameters
  const Real dt_mhd = pin->GetOrAddReal("Time","dt_mhd", 86.19e-6);       ///< mhd timestep [s]
  const Real dt_cd =  pin->GetOrAddReal("Time","dt_cd",  dt_mhd * 1e-1);  ///< current deposit timestep for electric field readjustment [s]
  const Real dt_LA =  pin->GetOrAddReal("Time","dt_LA",  dt_cd  * 1e-3 ); ///< large-angle collision step [s]
  const Real final_time = pin->GetOrAddReal("Time", "final_time", 1.0);   /// Final time [s]
  const Real timeStep = pin->GetOrAddReal("Simulation", "hRK", 1.e-6);    /// Runge kutta time in tau_c [-]
  const Real atol = pin->GetOrAddReal("Simulation", "atol", 1.e-10);      /// Absoulte tolerance for RK [-]
  const Real rtol = pin->GetOrAddReal("Simulation", "rtol", 1.e-7);       /// Realative toleratnce for RK[ [-]
  /// Derived parameters
  pin->GetOrAddReal("parthenon/time","tlim",0.0);         ///<WARNING: Setting this different in the input file would cause UB!
  pin->GetOrAddReal("parthenon/time","dt_force",dt_LA);   ///<WARNINS: Setting this different in the input file would cause UB!

  /// Reference parameters
  const Real B0 = pin->GetOrAddReal("Reference", "B0", 5.3);  ///< On-axis magnetic field [T]
  const Real a  = pin->GetOrAddReal("Reference", "a", 2.0);    ///< Minor radius [m] and reference length
  const Real R0 = pin->GetOrAddReal("Reference", "R0", 6.0);    ///< Major radius [m]
  const Real nD0 = pin->GetOrAddReal("Reference", "nD0", 1e20);    ///< Deutirium density [m^-3]
  const Real Te0 = pin->GetOrAddReal("Reference", "Te0", 10e3);  ///< Electron temperature [eV], and plasma temperature single-temperature model

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
  const Real Zeff = pin->GetOrAddReal("Plasma", "Zeff", (ZI*ZI*nI + nD0)/n_e0);
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


  ///< Numerical paremters
  const Real dampV                  = pin->GetOrAddReal("Numerical", "dampV", 0.01); ///< Stabilization coefficeint for velocity gradient
  const Real itime                  = pin->GetOrAddReal("Numerical", "itime", 0.0); ///< Initial time for mhd counters
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
  const Real eta_mu0aVa = etaplasma / eta0; // converts eta * \curl B to V_A B_0
  const Real etaec_a3VaB0 = etaplasma * e * c / pow(a,3) / E0; // converts eta J to V_A B_0

  if (mhd_context == NULL) return;
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
  mhd_context->itime                  = itime * tau_c / tauA;
  mhd_context->ftime                  = final_time / tauA;
  mhd_context->phibtype               = pin->GetOrAddInteger("MHD_Config", "phibtype",  1);
  mhd_context->dr                     = dR;
  mhd_context->dphi                   = (mhd_context->phimax - mhd_context->phimin) / mhd_context->Nphi;
  mhd_context->dz                     = dZ;
  mhd_context->pred_loop              = pin->GetOrAddInteger("MHD_Config", "pred_loop",  0);
  mhd_context->tstype                 = pin->GetOrAddInteger("MHD_Config", "tstype",  2);
  mhd_context->jtype                  = pin->GetOrAddInteger("MHD_Config", "jtype",  2);
  mhd_context->adaptdt                = pin->GetOrAddInteger("MHD_Config", "adaptdt",  0);
  mhd_context->debug                  = pin->GetOrAddInteger("MHD_Config", "debug",  0);
  mhd_context->dump                   = pin->GetOrAddInteger("MHD_Config", "dump",  0);
  mhd_context->EnableRelaxation       = pin->GetOrAddInteger("MHD_Config", "EnableRelaxation",  0);
  mhd_context->EnableReadICFromBinary = pin->GetOrAddInteger("MHD_Config", "EnableReadICFromBinary",  1);
  mhd_context->prestep                = pin->GetOrAddInteger("MHD_Config", "prestep",  1);
  mhd_context->savecoords             = pin->GetOrAddInteger("MHD_Config", "savecoords",  0);
  mhd_context->savesol                = pin->GetOrAddInteger("MHD_Config", "savesol",  0);
  mhd_context->delay_kinetic          = pin->GetOrAddReal("MHD_Config", "delay_kinetic",  3);
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


  mhd_context->axis[0] = Rc;
  mhd_context->axis[1] = Zc;

  int nt = 2, ndims = 3;
  mhd_context->jre_data = new double[NR * NZ * ndims * nt];
  mhd_context->jre    = mhd_context->jre_data;
  for (int i = 0; i < NR * NZ * 3 * 2; ++i)
    mhd_context->jre_data[i] = 0.0;
  mhd_context->jreR   = mhd_context->jre;
  mhd_context->jrephi = mhd_context->jre +     NR * NZ;
  mhd_context->jreZ   = mhd_context->jre + 2 * NR * NZ;

  mhd_context->field_data = new double[NR * NZ * 4 * ndims * nt];

  mhd_context->Ebc = 0;
  mhd_context->tempdump = 0;
  mhd_context->dumpfreq = std::ceil(mhd_context->ftime / (10.0 * mhd_context->dt));
  mhd_context->testSpGD = 0;
  mhd_context->testSpGDsamerhs = 0;
  mhd_context->oldstep = 0;
  mhd_context->n_record =0;
  mhd_context->n_record_Steady_jRE = 0;
  mhd_context->CorrectorIdentifier = 1;

  mhd_context->ParticlesCreated = 0;
}

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin, User* mhd_context) {
  /// Physical constants
  static constexpr Real eps0 = pc::eps0; ///< Vacuum permittivity [F / M]
  static constexpr Real c = pc::c;       ///< Speed of light [m/s]
  static constexpr Real me = pc::me;     ///< electron mass [kg]
  static constexpr Real e = pc::qe;      ///< electron charge [C]

  /// Time discretization parameters
  const Real dt_mhd = pin->GetReal("Time","dt_mhd");       ///< mhd timestep [s]
  const Real dt_cd =  pin->GetReal("Time","dt_cd");  ///< current deposit timestep for electric field readjustment [s]
  const Real dt_LA =  pin->GetReal("Time","dt_LA"); ///< large-angle collision step [s]
  const Real final_time = pin->GetReal("Time", "final_time");   /// Final time [s]
  const Real timeStep = pin->GetReal("Simulation", "hRK");    /// Runge kutta time in tau_c [-]
  const Real atol = pin->GetReal("Simulation", "atol");      /// Absoulte tolerance for RK [-]
  const Real rtol = pin->GetReal("Simulation", "rtol");       /// Realative toleratnce for RK[ [-]
  /// Derived parameters
  pin->GetReal("parthenon/time","tlim");         ///<WARNING: Setting this different in the input file would cause UB!
  pin->GetReal("parthenon/time","dt_force");   ///<WARNINS: Setting this different in the input file would cause UB!

  /// Reference parameters
  const Real B0 = pin->GetReal("Reference", "B0");  ///< On-axis magnetic field [T]
  const Real a  = pin->GetReal("Reference", "a");    ///< Minor radius [m] and reference length
  const Real R0 = pin->GetReal("Reference", "R0");    ///< Major radius [m]
  const Real nD0 = pin->GetReal("Reference", "nD0");    ///< Deutirium density [m^-3]
  const Real Te0 = pin->GetReal("Reference", "Te0");  ///< Electron temperature [eV], and plasma temperature single-temperature model

  /// Derived parameters
  const Real VA   = pin->GetReal("Derived", "VA"); ///< Alfven velocity [m/s]
  const Real tauA = pin->GetReal("Derived", "tauA");                  ///< Alfven time     [s]
  const Real E0   = pin->GetReal("Derived", "E0");                   ///< Reference electric field in MHD [V/m]
  const Real J0   = pin->GetReal("Derived", "J0");            ///< Reference current density [A/m^2]
  const Real eta0 = pin->GetReal("Derived", "eta0");            ///< Reference resitivity [Ohm*m]
  const Real eta  = pin->GetReal("Derived", "eta");                     ///< Resitivity scale     [-]
  const Real Re   = pin->GetReal("Derived", "Re");                    ///< Reinolds Number

  ///< Plasma composition parameters
  const Real vTe = sqrt(2.0*Te0 * e / me); ///< Thermal velocity
  const Real Z0 = pin->GetReal("Plasma", "Z0"); ///< Atomic number of impurity (Z)
  const Real ZI = pin->GetReal("Plasma", "ZI");  ///< Charge of impurity
  const Real fI = pin->GetReal("Plasma", "fI");  ///<Fraction of impurity density, normalized to deuterium denstiy (nD0)
  const Real nI = fI*nD0; ///< Impurity density [m^-3]
  const Real n_e0 = nD0 + ZI*nI; ///< Free electron density [m^-3]
  const Real Zeff = pin->GetReal("Plasma", "Zeff");
  const Real NeI = Z0 - ZI; ///< Number of bound electrons
  const Real Coulog0 = pin->GetReal("Plasma", "Coulog0");
  const Real Rc = pin->GetReal("Plasma", "Rc"); ///< Initial guess for magnetic axis, R [-], length normalized
  const Real Zc = pin->GetReal("Plasma", "Zc"); ///< Initial guess for magnetic axis, Z [-], length normalized

  const Real L11 = 0.58 * 32.0 / (3.0 * M_PI);
  const Real sigmapar = 12.0 * pow(M_PI, 1.5) / sqrt(2.0) * pow(Te0 * e, 1.5) * pow(eps0, 2) / (Zeff * pow(e, 2) * sqrt(me) * Coulog0) * L11;
  const Real etaplasma = pin->GetReal("Plasma", "etaplasma");

  const Real Rmin = pin->GetReal("Geometry", "rmin"); ///< Minimum R [-]
  const Real Rmax = pin->GetReal("Geometry", "rmax"); ///< Maximum R [-]
  const Real Zmin = pin->GetReal("Geometry", "zmin");///<  Minimum Z [-]
  const Real Zmax = pin->GetReal("Geometry", "zmax"); ///< Maximum Z [-]

  ///< Runaway parameters
  const Real c_vTe = pin->GetOrAddReal("Collisions", "c_vTe", c / vTe); ///< Guiding center equations coefficient [-]
  const int NSA = pin->GetOrAddInteger("Collisions", "NSA", 150);       ///< Number of small angle collisions     [-]
  const Real k = pin->GetOrAddReal("Collisions", "k", 5.0);
  const Real aI            = pin->GetOrAddReal("Collisions", "aI", 0.3285296762792767);  ///<
  const Real FineStructure = 1. / 137.035999;  // Fine Structure constant
  const Real II            = pin->GetOrAddReal("Collisions", "II", 219.5 / pc::me / pc::c / pc::c); // Mean exitation energy
  const Real PSCoefDnRA    = 1.0 + NeI * fI / (1.0 + ZI * fI);

  ///< Numerical paremters
  const int NR                      = pin->GetInteger("Numerical", "NR");
  const int Nphi                    = pin->GetInteger("Numerical", "Nphi");
  const int NZ                      = pin->GetInteger("Numerical", "NZ");

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
  const Real eta_mu0aVa = etaplasma / eta0; // converts eta * \curl B to V_A B_0
  const Real etaec_a3VaB0 = etaplasma * e * c / pow(a,3) / E0; // converts eta J to V_A B_0

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {

    std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "cBn = " << std::setw(20) << eta_mu0aVa
              << "Jn  = " << std::setw(20) << etaec_a3VaB0
              << "En = "  << std::setw(20) << En
              << "Ec = "  << std::setw(20) << Ec
              << std::endl;
  }

  auto pkg = std::make_shared<StateDescriptor>("Deck");

  pkg->AddParam("dt_LA",  dt_LA / tau_c);
  pkg->AddParam("dt_cd",  dt_cd / tau_c);
  pkg->AddParam("dt_mhd", dt_mhd / tau_c);

  const std::string filePath = pin->GetOrAddString("Simulation", "file_path", "current.out");
  pkg->AddParam("filePath", filePath);
  if (rank == 0) {
    std::ofstream(pkg->Param<std::string>("filePath"));
  }

  const Real gamma_min      = pin->GetOrAddReal("BoundaryConditions", "gamma_min",1.02);
  const Real p_BC = momentum_(pin->GetOrAddReal("BoundaryConditions", "gamma_BC", 1.02));
  const Real p_RE = momentum_(pin->GetOrAddReal("BoundaryConditions", "gamma_RE", 1.02));

  pkg->AddParam("gamma_min", gamma_min);
  pkg->AddParam("p_BC", p_BC);
  pkg->AddParam("p_RE", p_RE);

  int * ts = new int;
  *ts = 0;
  pkg->AddParam("ts", ts);

  pkg->AddParam("final_time", final_time);
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

  int nphi_data = 1;
  int nt = 2;

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

  ConfigurationDomainGeometry cdg(RminCellCenter, ZminCellCenter, dR, dZ, -3, indicator);
  pkg->AddParam("CDG", cdg);
  EM_Field field_interpolation(NR, NZ, nphi_data, nt, RminCellCenter, ZminCellCenter, dR, dZ, En, eta_mu0aVa, etaec_a3VaB0, cdg);
  auto field_data = field_interpolation.getDataRef();
  using Host = Kokkos::HostSpace;
  using Unmanaged = Kokkos::MemoryTraits<Kokkos::Unmanaged>;
  Kokkos::View<Real******, Kokkos::LayoutLeft, Host, Unmanaged> field_data_h(mhd_context->field_data, NR, NZ, 4, 3, 1, 2);
  Kokkos::deep_copy(field_data, field_data_h);
  field_interpolation.interpolate();
  pkg->AddParam("Field", field_interpolation);
  pkg->AddParam("FieldData", field_data_h);

  auto jre = field_interpolation.getJreDataSubview();
  Kokkos::View<double***> jre_backup("jre_backup", jre.extent(0), jre.extent(1), jre.extent(2));
  pkg->AddParam("JreBackup", jre_backup);

  Kokkos::View<Real***, Kokkos::LayoutRight, Host, Unmanaged> jre_h(mhd_context->jre_data, NR, NZ, 3);
  pkg->AddParam("JreData", jre_h);

  const Real wce0 = pc::qe * B0 / pc::me; // Electron gyrofrequency
  const Real c_aw0 =  pin->GetOrAddReal("GuidingCenterEquations", "c_aw0", pc::c/a/wce0);
  const Real ct_a =   pin->GetOrAddReal("GuidingCenterEquations", "ct_a", pc::c/vTe);
  const Real alpha0 = pin->GetOrAddReal("GuidingCenterEquations", "alpha0", tau_c/tau_a);

  GuidingCenterEquations<EM_Field, true, false> gce(field_interpolation, c_aw0, ct_a,
                                                        alpha0);
  pkg->AddParam("GCE", gce);


  const int npart =  pin->GetOrAddInteger("ParticleSeed", "num_particles_per_block", 16);
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


  return pkg;
}

void InitializeDriver(ParthenonManager* man) {
  auto pkg = man->pmesh.get()->packages.Get("Deck");

  auto driver = std::make_shared<RunawayDriver>(man->pinput.get(), man->app_input.get(), man->pmesh.get());
  driver->tm.tlim = 0.0;
  pkg->AddParam("Driver", driver);
}

void Push(ParthenonManager * man) {
  auto pkg = man->pmesh.get()->packages.Get("Deck");
  auto driver = pkg->Param<std::shared_ptr<RunawayDriver>>("Driver");
  const auto dt_mhd = pkg->Param<Real>("dt_mhd");
  const auto dt_cd = pkg->Param<Real>("dt_cd");
  while (driver.get()->tm.tlim < dt_mhd) {
    driver.get()->tm.tlim += dt_cd;
	  auto driver_status = driver.get()->Execute();
  }
}

auto &GetCoords(std::shared_ptr<MeshBlock> &pmb) { return pmb->coords; }
auto &GetCoords(MeshBlock *pmb) { return pmb->coords; }
auto &GetCoords(Mesh *pm) { return pm->block_list[0]->coords; }



void SaveState(Mesh* pm) {
  auto md = pm->mesh_data.Get();
  auto desc_swarm = parthenon::MakeSwarmPackDescriptor<Kinetic::status>("particles");
  auto pack_swarm = desc_swarm.GetPack(md.get());

  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        auto [b, n] = pack_swarm.GetBlockParticleIndices(idx);
        // block and particle indices

        if (pack_swarm(b, Kinetic::status(), n) & Kinetic::ALIVE) {
          pack_swarm(b, Kinetic::status(), n) |= Kinetic::PROTECTED;
        } else {
          const auto swarm = pack_swarm.GetContext(b);
          swarm.MarkParticleForRemoval(n);
        }
      });

  auto pkg = pm->packages.Get("Deck");
  auto jre = pkg->Param<EM_Field>("Field").getJreDataSubview();
  auto jre_backup = pkg->Param<Kokkos::View<Real***>>("JreBackup");
  Kokkos::deep_copy(jre_backup, jre);
}

void RestoreState(Mesh* pm) {
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

  auto pkg = pm->packages.Get("Deck");
  auto jre = pkg->Param<EM_Field>("Field").getJreDataSubview();
  auto jre_backup = pkg->Param<Kokkos::View<Real***>>("JreBackup");
  Kokkos::deep_copy(jre, jre_backup);
}

void ComputeParticleWeights(Mesh* pm) {

  auto md = pm->mesh_data.Get();
  auto pkg = pm->packages.Get("Deck");
  const auto field_interpolation = pkg->Param<EM_Field>("Field");

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

  std::cout << "Calculating current";

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
            weight += getParticleCurrent(X, t, w, field_interpolation);
          }

        }
      },
      I_re);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Allreduce(MPI_IN_PLACE,&I_re,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);


  Kokkos::fence();
  Real w = seed_current / I_re;
  std::cout << "I_re = " << std::setw(20) << I_re
            << "   w = " << std::setw(20) << w;
  parthenon::par_for(DEFAULT_LOOP_PATTERN, PARTHENON_AUTO_LABEL,
                     DevExecSpace(), 0, pack_swarm_r.GetMaxFlatIndex(),
                     // new_n ranges from 0 to N_new_particles
                     KOKKOS_LAMBDA(const int idx) {
        auto [b, n] = pack_swarm_r.GetBlockParticleIndices(idx);
        // block and particle indices
        pack_swarm_r(b, Kinetic::weight(), n) = w;
      });
}

} // namespace Kinetic
