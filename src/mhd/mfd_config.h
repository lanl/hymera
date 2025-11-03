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

#if !defined(MFD_CONFIG_H)
#define MFD_CONFIG_H

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdraw.h>
#include <petscksp.h>
#include <petscdmstag.h>
#include <petscsys.h>
#include <petscvec.h>
#include <mpi.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <petsc/private/dmstagimpl.h>
#include "KineticContext.h"


/* Shorter, more convenient names for DMStagStencilLocation entries */
#define BACK_DOWN_LEFT   DMSTAG_BACK_DOWN_LEFT
#define BACK_DOWN        DMSTAG_BACK_DOWN
#define BACK_DOWN_RIGHT  DMSTAG_BACK_DOWN_RIGHT
#define BACK_LEFT        DMSTAG_BACK_LEFT
#define BACK             DMSTAG_BACK
#define BACK_RIGHT       DMSTAG_BACK_RIGHT
#define BACK_UP_LEFT     DMSTAG_BACK_UP_LEFT
#define BACK_UP          DMSTAG_BACK_UP
#define BACK_UP_RIGHT    DMSTAG_BACK_UP_RIGHT
#define DOWN_LEFT        DMSTAG_DOWN_LEFT
#define DOWN             DMSTAG_DOWN
#define DOWN_RIGHT       DMSTAG_DOWN_RIGHT
#define LEFT             DMSTAG_LEFT
#define ELEMENT          DMSTAG_ELEMENT
#define RIGHT            DMSTAG_RIGHT
#define UP_LEFT          DMSTAG_UP_LEFT
#define UP               DMSTAG_UP
#define UP_RIGHT         DMSTAG_UP_RIGHT
#define FRONT_DOWN_LEFT  DMSTAG_FRONT_DOWN_LEFT
#define FRONT_DOWN       DMSTAG_FRONT_DOWN
#define FRONT_DOWN_RIGHT DMSTAG_FRONT_DOWN_RIGHT
#define FRONT_LEFT       DMSTAG_FRONT_LEFT
#define FRONT            DMSTAG_FRONT
#define FRONT_RIGHT      DMSTAG_FRONT_RIGHT
#define FRONT_UP_LEFT    DMSTAG_FRONT_UP_LEFT
#define FRONT_UP         DMSTAG_FRONT_UP
#define FRONT_UP_RIGHT   DMSTAG_FRONT_UP_RIGHT

/* Define a structure for user context */
typedef struct{
  PetscReal   density;          /* Density */
  PetscReal   L0;               /* Characteristic length */
  PetscReal   B0;               /* Initial magnetic field magnitude */
  PetscReal   V_A;              /* Alfven speed */
  PetscReal   mu0;              /* Permeability */
  PetscReal   mi;               /* Ion mass */
  PetscReal   eta0;             /* Resistivity constant used in normalization eta0 := (L0 * V_A * mu0) */
  PetscReal   eta;              /* Resistivity normalized wrt eta0 := (L0 * V_A * mu0) */
  PetscReal   etawall;          /* Resistivity between inner Vacuum Vessel wall and plasma wall */
  PetscReal   etawallperp;      /* Poloidal resistivity in between inner Vacuum Vessel wall and plasma wall */
  PetscReal   etawallphi;       /* Toroidal resistivity between inner Vacuum Vessel wall and plasma wall */
  PetscReal   etawallphi_isol_cell;/* Lower toroidal resistivity in isolated cells between inner Vacuum Vessel wall and plasma wall */
  PetscReal   etaplasma;        /* Resistivity inside separatrix */
  PetscReal   etasepwal;        /* Resistivity inside plasma but outside separatrix */
  PetscReal   etaVV;            /* Resistivity inside Vacuum Vessel wall */
  PetscReal   etaout;           /* Resistivity outside outer Vacuum Vessel wall */
  PetscReal   rmin;             /* Minimum radius */
  PetscReal   rmax;             /* Maximum radius */
  PetscReal   phimin;           /* Minimum azimuth */
  PetscReal   phimax;           /* Maximum azimuth */
  PetscReal   zmin;             /* Minimum height */
  PetscReal   zmax;             /* Maximum height */
  PetscInt    Nr;               /* Global number of grid points in r direction */
  PetscInt    Nphi;             /* Global number of grid points in phi direction */
  PetscInt    Nz;               /* Global number of grid points in z direction */
  PetscReal   dt;               /* Length of timestep */
  PetscReal   itime;            /* Initial time */
  PetscReal   ftime;            /* Final time */
  DM          coorda;           /* DM used only to get cell/face/edge/vertex center coordinates */
  PetscInt    ictype;           /* Type of initial conditions */
  PetscInt    phibtype;         /* Boundary Type for phi */
  PetscReal   dr;               /* Radius step size */
  PetscReal   dphi;             /* Azimuth step size */
  PetscReal   dz;               /* Height step size */
  PetscInt    pred_loop;	/* Number of times the predictor step is being executed */
  PetscInt    tstype;	        /* Timestepping method: 1 for Forward Euler, 2 for Backward Euler, 3 for Crank-Nicholson */
  PetscInt    jtype;            /* Jacobian type (0: user provide Jacobian, 1: slow finite difference, 2: fd with coloring) */
  PetscInt    adaptdt;          /* Flag for using adaptive step size */
  PetscInt    n_record;        /* Counter for successful TS steps */
  PetscInt    n_record_Steady_jRE;        /* Counter for TS step where j_RE change is less than 10% */

  PetscInt    debug;            /* Flag for displaying debug information */
  PetscInt    dump;             /* Flag for saving output in vtk files */
  PetscInt    prestep;          /* Flag for activating the prestep to approximate the runaway current contribution */
  Vec 	      oldX;
  Vec         X_star;
  IS          isV;              /* Indexing PETSc object for velocity */
  IS          isni;              /* Indexing PETSc object for ion number density */
  IS 	      isB;              /* Indexing PETSc object for B Field */
  IS          isEP;              /* Indexing PETSc object for Electrostatic Potential */
  IS          istau;              /* Indexing PETSc object for divergence-free tau field */
  IS          isE_boundary;     /* Indexing PETSc object for boundary edges */
  IS          isB_boundary;     /* Indexing PETSc object for boundary faces */
  IS          isni_boundary;     /* Indexing PETSc object for boundary cells */
  PetscReal   *dataC;           /* Array containing the level set function */
  PetscReal   *dataz,*dataphi,*datar;/* Arrays containing the initial B field components */
  PetscReal   *datag,*datapsi;/* Arrays containing the initial psi and G(psi) values */
  PetscInt    numz,numphi,numr;/* Lengths of arrays containing the initial B field components */
  PetscInt    numg,numpsi;/* Lengths of arrays containing the initial psi and G(psi) values */
  PetscInt    Ebc;             /* Type of boundary conditions for E field */
  PetscInt    savecoords;      /* Flag for saving coordinates of cell centers, face centers and edge centers in .m files: if 0 save option disabled, otherwise save enabled */
  PetscInt    savesol;         /* Flag for saving values of magnetic fields at face centers in .m files: if 0 the save option is disabled, otherwise the save option is enabled */
  PetscInt    tempdump;        /* Flag for saving intermediate output in .dat files */
  PetscInt    dumpfreq;        /* Frequency for saving intermediate output in .dat files */
  PetscInt    testSpGD;        /* Flag for testing S the Schur complement matrix vs S - dt*GD the Schur complement matrix augmented by the gradient of derived divergence operator */
  PetscInt    testSpGDsamerhs; /* Flag for setting the right hand side vector while testing S the Schur complement matrix vs S - dt*GD the Schur complement matrix augmented by the gradient of derived divergence operator: 0 for a different right-hand side (S - dt*GD)*E, 1 for the same right-hand side S*E */
  PetscInt    oldstep;           /* Last step in previous simulation when a restart is used
  //PetscInt    dtchange;        /* Flag to indicate if a restart is used with a new time step : olditime...olddt...newitime...newdt...newftime */
  //PetscReal   olddt;           /* Length of old timestep when a restart is used with a new time step */
  PetscScalar Iphi1;             /* Current intensity inside plasma */
  PetscScalar Iphi2;             /* Current intensity outside plasma */
  PetscScalar Iphi3;             /* Current intensity inside the vaccum vessel */
  PetscScalar dampV;             /* Stabilization coefficient in the constraint : ((\nabla x B) x B - dampV V) . e_{R/Z} = 0 */
  PetscScalar Re;                /* Reynolds parameter in the constraint : - (\nabla x B) x B + (V.\nabla) V - 1/(Re) (\nabla^2 V) = 0 */
  Vec         DiagMe;             /* Vector containing the diagonal of mass matrix M_e with material properties in extended DM */
  Vec	      MeVec;		 /* Vector containing the diagonal of mass matrix M_e with material properties in regular DM */
  Vec         DiagMe1;             /* Vector containing the diagonal of mass matrix M_e without material properties in extended DM */
  Vec	      MeVec1;		 /* Vector containing the diagonal of mass matrix M_e without material properties in regular DM */
  TS          ts;             /* time integrator */
  Mat         DiagBlock_V;    /* Diagonal block for Velocity in the user-provided preconditioner */
  Mat         DiagBlock_EP;    /* Diagonal block for Electrostatic Potential in the user-provided preconditioner */
  Mat         DiagBlock_B;    /* Diagonal block for Magnetic Field in the user-provided preconditioner */
  KSP         KSP_V;          /* KSP solver for Velocity in the user-provided preconditioner */
  KSP         KSP_B;          /* KSP solver for Magnetic Field in the user-provided preconditioner */
  KSP         KSP_EP;         /* KSP solver for Electrostatic Potential in the user-provided preconditioner */
  IS          isALL_V;        /* Indexing PETSc object for {Electrostatic Potential,Tau field,Magnetic Field,Ion number density} */
  Mat         OffDiagBlock_U;    /* Upper Off-diagonal block of Jacobian matrix according to the {ETBN,V} partitioning */
  Mat         OffDiagBlock_L;    /* Lower Off-diagonal block of Jacobian matrix according to the {ETBN,V} partitioning */
  PetscScalar NT[6][6];       /* Interpolation array */

  PetscInt EnableRelaxation;
  PetscInt EnableReadICFromBinary;

  double * jre_data;
  double * field_data;
  double * jre;
  double * jreR;
  double * jrephi;
  double * jreZ;
  void* runaway_solver;
  double prev_current;        /* Stores the runaway current of the previous time step */
  double present_current;
  int CorrectorIdentifier;

  char input_folder [PETSC_MAX_PATH_LEN];

  int delay_kinetic;

  int poincare_counter;
  int field_counter;

  KineticContext kctx;
  void* field_interpolation;
  double axis[2];
} User;

typedef struct LocalCoordinate
{
  PetscScalar R;
  PetscScalar Z;
  int iR;
  int iZ;
} tLocalCoordinate;

typedef struct HermiteDivFreeFields
{
  PetscScalar *m_hd_psi_iz;
  PetscScalar *m_hd_chi_iz;
  PetscScalar *m_hd_psi_ir;

  int m_nRR;
  int m_nZR;
  int m_nRZ;
  int m_nZZ;

  PetscScalar *m_hp_RR;
  PetscScalar *m_hp_RZ;
  PetscScalar *m_hp_ZR;
  PetscScalar *m_hp_ZZ;
} tHermiteDivFreeFields;
#endif /* defined(MFD_CONFIG_H) */
