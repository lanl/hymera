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
#include <hFlux/Evaluator.hpp>

using namespace parthenon;

#include <parthenon/driver.hpp>
#include "kinetic/kinetic.hpp"
#include "kinetic/KineticParams.hpp"
#include "kinetic/Log.hpp"
#include "kinetic/DiagnosticsFile.hpp"
#include "kinetic/PushTiming.hpp"
#include "kinetic/GuidingCenterEquations.hpp"
#include "kinetic/LargeAngleCollision.hpp"
#include "kinetic/SmallAngleCollision.hpp"
#include "kinetic/kinetic.hpp"
#include "kinetic/ConfigurationDomainGeometry.hpp"
#include "kinetic/CurrentDensity.hpp"
#include "kinetic/AnalyticField.hpp"
#include "kinetic/FieldEvaluator.hpp"
#include "tasks/Tasks.h"
#include "mhd/mhd.h"
#include "util/H5io.hpp"

using parthenon::constants::SI;
using parthenon::constants::PhysicalConstants;
using pc = PhysicalConstants<SI>;

#include "kinetic/Setup.hpp"

namespace Kinetic {

void SetupMHD(ParameterInput *pin, const KineticParams &prm, User *mhd_context,
              DualView3 &Jre_mhd) {
  /// Physical constants the MHD context needs.
  static constexpr Real mi = pc::amu;    ///< Ion mass = atomic mass unit [kg]
  static constexpr Real mu0 = pc::mu0;   ///< Vacuum permeabulity [N / A^2]

  // Names bound to the run parameters, as Initialize has always used them.
  const Real &dt_mhd = prm.dt_mhd, &final_time = prm.final_time;
  const Real &B0 = prm.B0, &a = prm.a, &nD0 = prm.nD0;
  const Real &VA = prm.VA, &tauA = prm.tauA, &eta0 = prm.eta0, &eta = prm.eta, &Re = prm.Re;
  const Real &etaplasma = prm.etaplasma, &etawall = prm.etawall, &etawallperp = prm.etawallperp;
  const Real &etawallphi = prm.etawallphi, &etawallphi_isol_cell = prm.etawallphi_isol_cell;
  const Real &etasepwal = prm.etasepwal, &etaVV = prm.etaVV, &etaout = prm.etaout;
  const Real &Rmin = prm.Rmin, &Rmax = prm.Rmax, &Zmin = prm.Zmin, &Zmax = prm.Zmax;
  const Real &dampV = prm.dampV, &itime = prm.itime;
  const int &NR = prm.NR, &Nphi = prm.Nphi, &NZ = prm.NZ;
  const Real &dR = prm.dR, &dZ = prm.dZ;

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
  mhd_context->zmax                   = Zmax * a;
  mhd_context->dt                     = dt_mhd / tauA;
  // Initial-condition selector, consumed in src/mhd/ts_functions.c.
  //   9  = EFIT + Grad-Shafranov equilibrium from inputs/mhd/vec{psi,g}*.txt
  //        (production; the only value used by physics runs)
  //   1..8, 10..13 = manufactured analytic solutions with closed forms, used
  //        by the regression harness to measure spatial convergence order.
  //        These read none of the inputs/mhd/*.txt equilibrium data.
  // Default stays 9 so production behaviour is unchanged.
  mhd_context->ictype                 = pin->GetOrAddInteger("MHD_Config", "ictype", 9);
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
  // Per-step diagnostics from Monitor: step norms, max|div B|, toroidal
  // currents. Costs a nested linear solve per step, hence off by default.
  // The regression harness turns this on and uses the output as its
  // per-step fingerprint.
  mhd_context->monitor                = pin->GetOrAddInteger("MHD_Config", "monitor",  0);
  // Advective inertia n_i (V.grad)V in the production momentum residual.
  // Off by default (it used to be multiplied by 0.0).
  mhd_context->inertia                = pin->GetOrAddInteger("MHD_Config", "inertia",  0);
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
  // Number of relaxation steps (100 tau_A each) when the equilibrium is relaxed
  // (ic_binary_load = 0). Read only when present, so that the parameter dump of
  // other runs is unchanged.
  mhd_context->relax_steps = pin->DoesParameterExist("MHD_Config", "relax_steps")
                                 ? pin->GetInteger("MHD_Config", "relax_steps") : 30;

  mhd_context->Ebc = 0;
  mhd_context->tempdump = 0;
  mhd_context->dumpfreq = std::ceil(mhd_context->ftime / (10.0 * mhd_context->dt));
  mhd_context->testSpGD = 0;
  mhd_context->testSpGDsamerhs = 0;
  mhd_context->oldstep = 0;
  mhd_context->n_record =0;
  mhd_context->n_record_Steady_jRE = 0;

  mhd_context -> jre = wrap_view(Jre_mhd.view_host());

  MHD_CHECK(mhd_initialize(mhd_context));

  /// Do single step before seeding particles
  MHD_CHECK(mhd_step(mhd_context));
}

FieldGridResult SetupFieldGrid(ParameterInput *pin, const KineticParams &prm,
                               User *mhd_context, DualView3 &Jre_mhd,
                               std::shared_ptr<StateDescriptor> &pkg) {
  // Names bound to the run parameters, as Initialize has always used them.
  const Real &dt_mhd = prm.dt_mhd;
  const int &nCorrectorSteps = prm.nCorrectorSteps, &nCD = prm.nCD;
  const Real &dtLA_over_tauC = prm.dtLA_over_tauC;
  const Real &timeStep = prm.timeStep, &atol = prm.atol, &rtol = prm.rtol;
  const int &integrator = prm.integrator;
  const std::string &filePath = prm.filePath;
  const Real &Z0 = prm.Z0, &ZI = prm.ZI, &fI = prm.fI, &Zeff = prm.Zeff, &NeI = prm.NeI;
  const Real &Coulog0 = prm.Coulog0;
  Real Rc = prm.Rc, Zc = prm.Zc;  // updated by the magnetic-axis search below
  const Real &Rmin = prm.Rmin, &Rmax = prm.Rmax, &Zmin = prm.Zmin, &Zmax = prm.Zmax;
  const Real &gamma_min = prm.gamma_min, &p_BC = prm.p_BC, &p_RE = prm.p_RE;
  const Real &c_vTe = prm.c_vTe;
  const int &NSA = prm.NSA;
  const Real &k = prm.k, &aI = prm.aI, &FineStructure = prm.FineStructure, &II = prm.II;
  const Real &PSCoefDnRA = prm.PSCoefDnRA;
  const int &NR = prm.NR, &NZ = prm.NZ;
  const Real &dR = prm.dR, &dZ = prm.dZ;
  const Real &RminCellCenter = prm.RminCellCenter, &ZminCellCenter = prm.ZminCellCenter;
  const Real &tau_c = prm.tau_c, &En = prm.En;
  const Real &eta_norm = prm.eta_norm, &eta_a3VaB0 = prm.eta_a3VaB0;
  const Real &E_0 = prm.E_0;

  pkg->AddParam("E_0", E_0);

  pkg->AddParam("NR", NR);
  pkg->AddParam("NZ", NZ);
  pkg->AddParam("NR_maxE", prm.NR_maxE);
  pkg->AddParam("NZ_maxE", prm.NZ_maxE);

  pkg->AddParam("nCorrectorSteps",  nCorrectorSteps);
  pkg->AddParam("nCD",  nCD);

  pkg->AddParam("tau_c",  tau_c);
  pkg->AddParam("eta_norm",  eta_norm);
  pkg->AddParam("eta_a3VaB0",  eta_a3VaB0);
  pkg->AddParam("En",  En);

  int nLA = std::ceil((dt_mhd / nCD) / (dtLA_over_tauC * tau_c));
  Real my_dtLA_over_tauC = (dt_mhd / nCD / tau_c) / nLA;

  Log::Info("Large-angle step adjusted to divide the current-deposition step: {:g} -> {:g} tau_c ({} steps, tau_c = {:g} s)",
            dtLA_over_tauC, my_dtLA_over_tauC, nLA, tau_c);

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
  pkg->AddParam("integrator", integrator);

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

  View3 E_base("E_base", NR, NZ);
  Kokkos::deep_copy(E_base, 0.0);
  pkg->AddParam("E_base", E_base);

  View3 Jre("Jre", NR, NZ);
  Kokkos::deep_copy(Jre, 0.0);
  pkg->AddParam("Jre", Jre);

  pkg->AddParam("Jre_mhd", Jre_mhd);

  // The runaway current of the last MHD step, as passed to MHD; it corrects E in
  // the first current-deposition interval of the next step. Kept across a restart
  // in a side file of the MHD restart file (WorkBeforeRestartOutput).
  View3 Jre_seed("Jre_seed", NR, NZ);
  Kokkos::deep_copy(Jre_seed, 0.0);
  pkg->AddParam("Jre_seed", Jre_seed);

  Real Psi_min = 0.0;
  // Fill the interpolated (Hermite) field grid from the MHD solver. When there
  // is no MHD context (mhd_context == nullptr, e.g. the `profile` executable),
  // this block is skipped and the Hermite grid + magnetic axis are instead
  // loaded from a raw field dump via LoadRawFieldData().
  if (mhd_context != nullptr) {
    CommunicateBJV(data, mhd_context, FieldComponents::B);
    ComputeBaseElectricField_in_place(data, En, eta_norm, FieldComponents::B);
    CommunicateBJV(data, mhd_context, FieldComponents::Bt);
    ComputeBaseElectricField_in_place(data, En, eta_norm, FieldComponents::Bt);

    InterpolateTime(data, dt_mhd / tau_c);

    const Real seed_current_fraction = pin->GetOrAddReal("ParticleSeed", "current_fraction", 1.0e-3); // Used to determine the initial runaway current to adjust the Electric field.

    auto data_d = data.data.view_device();
    data.data.sync_device();
    Kokkos::parallel_for("Set intial runaway current",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {NR, NZ}),
        KOKKOS_LAMBDA(const int i, const int j) {
          Real mask = (indicator(i,j) > 0) ? Real(1.0) : Real(0.0);
          Real scale = mask * seed_current_fraction * eta_norm;
          data_d(i,j,FieldComponents::E + 0) -= scale * data_d(i, j, FieldComponents::J + 0);
          data_d(i,j,FieldComponents::E + 1) -= scale * data_d(i, j, FieldComponents::J + 1);
          data_d(i,j,FieldComponents::E + 2) -= scale * data_d(i, j, FieldComponents::J + 2);
        });
    data.data.modify_device();

    InterpolateHermiteBJE(data, FieldComponents::B);
    InterpolateHermiteBJE(data, FieldComponents::Bt);

    Interpolator<FIELD_SMOOTHNESS, FIELD_FD_STENSIL> itrp;
    itrp.computeFlux(data.hermite_locator,
        data.hermite_data.view_device(),
        data.psi_data.view_device());
    data.psi_data.modify_device();  // Mark psi data as modified

    int accent = 1;
    data.hermite_data.sync_host();
    data.psi_data.sync_host();

    findMagneticAxis(Rc, Zc,
                data.hermite_data.view_host(),
                data.psi_data.view_host(),
                Evaluator{data.hermite_locator},
                accent,
                Psi_min);

    Log::Info("Magnetic axis (R, Z) = ({:g}, {:g}), Psi_min = {:g}", Rc, Zc, Psi_min);
  }

  return {Rc, Zc, Psi_min};
}

void SetupParticles(ParameterInput *pin, const KineticParams &prm, User *mhd_context,
                    const FieldGridResult &axis, std::shared_ptr<StateDescriptor> &pkg) {
  // Names bound to the run parameters, as Initialize has always used them.
  const Real &dt_mhd = prm.dt_mhd;
  const Real &a = prm.a;
  const Real &tau_c = prm.tau_c;
  const int &NR_plot = prm.NR_plot, &NZ_plot = prm.NZ_plot, &Nt_plot = prm.Nt_plot;
  const Real Rc = axis.Rc, Zc = axis.Zc;  // magnetic axis from SetupFieldGrid
  const Real &Psi_min = axis.Psi_min;
  auto data = pkg->Param<FieldData_t>("FieldData");

  pkg->AddParam("Rc", Rc, Params::Mutability::Restart);
  pkg->AddParam("Zc", Zc, Params::Mutability::Restart);

  KineticLateParams late;
  late.ParseSeedPoint(pin, Rc, Zc);
  const Real &Rseed = late.Rseed, &Zseed = late.Zseed;
  pkg->AddParam("Rseed", Rseed, Params::Mutability::Restart);
  pkg->AddParam("Zseed", Zseed, Params::Mutability::Restart);

  if (mhd_context != nullptr) {
    data.psi_data.sync_device();
    Kokkos::parallel_for("Normalize psi",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {data.hermite_locator.nR,data.hermite_locator.nZ}),
    KOKKOS_LAMBDA(int i, int j){
      data.psi_data.view_device()(0,0,i,j) -= Psi_min;
    });
    data.psi_data.modify_device();
  }


  pkg->AddParam("plot_fields", prm.plot_fields);
  pkg->AddParam("substep_snapshots", prm.substep_snapshots);
  if (prm.plot_fields != 0) {
    const Real R0_plot = data.hermite_locator.R0 + 1e-10;
    const Real Z0_plot = data.hermite_locator.Z0 + 1e-10;
    const Real dR_plot = (data.hermite_locator.nR * data.hermite_locator.dR - 2e-10) / static_cast<Real> (NR_plot);
    const Real dZ_plot = (data.hermite_locator.nZ * data.hermite_locator.dZ - 2e-10) / static_cast<Real> (NZ_plot);
    const Real dt_plot = dt_mhd / tau_c / static_cast<Real> (Nt_plot-1);

    ParArray1D<Real> hpd_R("Hermite_Field_Plot_data_R", NR_plot);
    ParArray1D<Real> hpd_t("Hermite_Field_Plot_data_t", Nt_plot);
    ParArray1D<Real> hpd_Z("Hermite_Field_Plot_data_Z", NZ_plot);
    //  ParArrayND<Real> gce_data("GCE_data", NR_plot, NZ_plot, Nt_plot, 5);
    ParArrayND<Real> hpd_F("Hermite_Field_Plot_data_F", NR_plot, NZ_plot, Nt_plot, 3, 7);

    Kokkos::parallel_for("FillGrids_R", NR_plot,
        KOKKOS_LAMBDA(const int n) {
          hpd_R(n) = R0_plot + n * dR_plot;
        });
    Kokkos::parallel_for("FillGrids_Z", NZ_plot,
        KOKKOS_LAMBDA(const int n) {
          hpd_Z(n) = Z0_plot + n * dZ_plot;
        });
    Kokkos::parallel_for("FillGrids_t", Nt_plot,
        KOKKOS_LAMBDA(const int n) {
          hpd_t(n) = n * dt_plot;
        });

    pkg->AddParam("Hermite_Field_Plot_data_R", hpd_R);
    pkg->AddParam("Hermite_Field_Plot_data_Z", hpd_Z);
    pkg->AddParam("Hermite_Field_Plot_data_t", hpd_t);
    pkg->AddParam("Hermite_Field_Plot_data_F", hpd_F);
  }
 // pkg->AddParam("GCE_data", gce_data);

  late.ParseParticles(pin, prm);
  const Real &c_aw0 = late.c_aw0, &ct_a = late.ct_a, &alpha0 = late.alpha0;

  pkg->AddParam("c_aw0", c_aw0);
  pkg->AddParam("ct_a", ct_a);
  pkg->AddParam("alpha0", alpha0);

  const int npart = late.npart;
  pkg->AddParam("num_particles_per_block", npart);

// Paticles are counted on Collect current
  pkg->AddParam("num_particles_total", 0, Params::Mutability::Restart);
  const int num_particles_max = late.num_particles_max;

  pkg->AddParam("num_particles_max", num_particles_max);
  // Initialize random number generator pool
  const int rng_seed = late.rng_seed;
  RNGPool rng_pool(rng_seed);
  pkg->AddParam("rng_pool", rng_pool);

  const Real seed_current = late.seed_current; // 15 kAmps by default
  pkg->AddParam("seed_current", seed_current * a); // Convert from amps

  const Real gammamin = late.gammamin;
  pkg->AddParam("pmin", momentum_(gammamin));
  const Real gammamax = late.gammamax;
  pkg->AddParam("pmax", momentum_(gammamax));
  const Real ximin = late.ximin;
  pkg->AddParam("ximin", ximin);
  const Real ximax = late.ximax;
  pkg->AddParam("ximax", ximax);



  if (npart > 0) {
    Metadata swarm_metadata({Metadata::Provides, Metadata::None, Metadata::Restart});
    pkg->AddSwarm("particles", swarm_metadata);

    Metadata real_swarmvalue_metadata({Metadata::Real, Metadata::Restart});
    pkg->AddSwarmValue(Kinetic::p::name(), "particles",      real_swarmvalue_metadata);
    pkg->AddSwarmValue(Kinetic::xi::name(), "particles",     real_swarmvalue_metadata);
    pkg->AddSwarmValue(Kinetic::R::name(), "particles",      real_swarmvalue_metadata);
    pkg->AddSwarmValue(Kinetic::phi::name(), "particles",    real_swarmvalue_metadata);
    pkg->AddSwarmValue(Kinetic::Z::name(), "particles",      real_swarmvalue_metadata);
    pkg->AddSwarmValue(Kinetic::weight::name(),  "particles", real_swarmvalue_metadata);

    Metadata int_swarmvalue_metadata({Metadata::Integer, Metadata::Restart});
    pkg->AddSwarmValue(Kinetic::will_scatter::name(), "particles",
                       int_swarmvalue_metadata);
    pkg->AddSwarmValue(Kinetic::secondary_index::name(), "particles",
                       int_swarmvalue_metadata);
    pkg->AddSwarmValue(Kinetic::status::name(), "particles",
                       int_swarmvalue_metadata);
	  pkg->AddSwarmValue(p_phi::name(), "particles",real_swarmvalue_metadata);
    pkg->AddSwarmValue(mu::name(), "particles",real_swarmvalue_metadata);
  }
}

void SetupDiagnostics(ParameterInput *pin, const KineticParams &prm,
                      std::shared_ptr<StateDescriptor> &pkg) {
  // Names bound to the run parameters, as Initialize has always used them.
  const Real &dtLA_over_tauC = prm.dtLA_over_tauC;
  using SA = SmallAngleCollision<PartialScreening, EnergyScattering, ModifiedCouLog>;
  auto sa = pkg->Param<SA>("SmallAngleCollision");
  auto ms = pkg->Param<MollerSource>("MollerSource");

  const KineticSwitches sw = KineticSwitches::FromInput(pin);
  const int &EnableLargeAngleCollisions = sw.EnableLargeAngleCollisions;
  const int &EnableSmallAngleCollisions = sw.EnableSmallAngleCollisions;
  const int &EnableComputeConservedQuantities = sw.EnableComputeConservedQuantities;

  pkg->AddParam("EnableLargeAngleCollisions",      EnableLargeAngleCollisions);
  pkg->AddParam("EnableSmallAngleCollisions",      EnableSmallAngleCollisions);
  pkg->AddParam("EnableComputeConservedQuantities",EnableComputeConservedQuantities);

  // Log file for the weight-averaged conservation diagnostic (p_phi, mu).
  const std::string &conservation_log = sw.conservation_log;
  pkg->AddParam("conservation_log", conservation_log);
  if (EnableComputeConservedQuantities == 1 && Globals::my_rank == 0 && sw.legacy_dat == 1) {
    std::ofstream(conservation_log);  // truncate at start of run
  }
  pkg->AddParam("legacy_dat", sw.legacy_dat);
  pkg->AddParam("parallel_push", sw.parallel_push);

  // On restart the file is reopened as is; rows past the restart time are cut in
  // WorkBeforeLoop, once Parthenon has restored the simulation time.
  DiagnosticsFile::Open(pin->GetString("parthenon/job", "problem_id"), Globals::is_restart);
  // Push timing only when the deck asks for it: reading an existing parameter
  // leaves the parameter dump of other runs unchanged.
  if (pin->DoesParameterExist("Simulation", "push_timing") &&
      pin->GetInteger("Simulation", "push_timing") == 1) {
    pkg->AddParam("push_timing", 1);
    PushTiming::enabled = true;
  }
  if (auto *diag = DiagnosticsFile::Get())
    diag->DeclareSeries("conservation", {"t", "p_phi_mean", "mu_mean", "weight_total"},
                        {"s", "", "", ""});

  // Collision coefficients against momentum, written once per run.
  if (Globals::my_rank == 0) {
    std::vector<std::vector<double>> cols(8);
    Real p = momentum_(1. + 2.e-3);
    while (p < pkg->Param<Real>("pmax") + 20.0) {
      auto cc = sa.getCollisionCoefficients(p);
      const Real row[8] = {p, gamma_(p), sa.getSmallAngleCollisionTimestep(p),
                           cc.psi, cc.CB, cc.CF, cc.CouLogee_ratio,
                           ms.computeProbability(p, 1.0, dtLA_over_tauC, 1.002)};
      for (int c = 0; c < 8; ++c) cols[c].push_back(row[c]);
      p += 1e-2;
    }
    if (auto *diag = DiagnosticsFile::Get())
      diag->WriteProfile("collision",
          {"p", "gamma", "dtSA", "psi", "CB", "CF", "CouLogee_ratio", "probability"},
          {"m_e c", "", "tau_c", "", "", "", "", ""}, cols);
    if (sw.legacy_dat == 1) {
      std::ofstream ofs("collision_profiles.dat");
      ofs << std::format("{:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s} {:20s}",
          "#     p", "gamma", "dtSA", "psi", "CB", "CF", "CouLogee ratio", "probability");
      for (size_t i = 0; i < cols[0].size(); ++i)
        ofs << std::format("{:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e} {:20.14e}",
            cols[0][i], cols[1][i], cols[2][i], cols[3][i], cols[4][i], cols[5][i], cols[6][i], cols[7][i]
        ) << std::endl;
    }
  }


  // This filename is updated with OutputParameters.file_number and file_basename before each restart.
  // It is required to remember the Petsc file name associated with Parthenon restart.
  // Petsc file is where the mhd state is stored, it is separate from Parthenon restart file
  std::string mhd_restart_filename = "";
  pkg->AddParam("mhd_restart_filename", mhd_restart_filename, Params::Mutability::Restart);
  pkg->AddParam("ComputeInitialWeights", 1, Params::Mutability::Restart);
  pkg->AddParam("MyTime", 0.0, Params::Mutability::Restart);
  pkg->AddParam("PredictorIterationNumber", 0, Params::Mutability::Restart);
}

std::shared_ptr<StateDescriptor> Initialize(ParameterInput *pin, User* mhd_context) {
  // Read the deck and derive the normalisation (see KineticParams.hpp).
  const KineticParams prm = KineticParams::FromInput(pin);
  Log::Configure(pin);

  DualView3 Jre_mhd("Jre_mhd", prm.NR, prm.NZ, 3);
  Kokkos::deep_copy(Jre_mhd.view_device(), 0.0);
  Kokkos::deep_copy(Jre_mhd.view_host(), 0.0);

  if (mhd_context != nullptr) SetupMHD(pin, prm, mhd_context, Jre_mhd);

  auto pkg = std::make_shared<StateDescriptor>("Deck");
  const FieldGridResult axis = SetupFieldGrid(pin, prm, mhd_context, Jre_mhd, pkg);
  SetupParticles(pin, prm, mhd_context, axis, pkg);
  SetupDiagnostics(pin, prm, pkg);

  Log::Info("Kinetic package initialised");

  return pkg;
}

} // namespace Kinetic
