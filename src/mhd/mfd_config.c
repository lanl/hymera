#include <petscsys.h>
#include <petscdm.h>
#include <petscts.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>

#include "mfd_config.h"

PetscErrorCode AppCtxView(MPI_Comm comm, const User *ctx) {
  PetscFunctionBeginUser;
  if (!ctx) PetscFunctionReturn(PETSC_ERR_ARG_NULL);

  PetscCall(PetscPrintf(comm,"================ AppCtx dump ================\n"));

  PetscCall(PetscPrintf(comm,"-- Physical / normalization parameters --\n"));
  PetscCall(PetscPrintf(comm,"density      = %g\n",(double)ctx->density));
  PetscCall(PetscPrintf(comm,"L0           = %g\n",(double)ctx->L0));
  PetscCall(PetscPrintf(comm,"B0           = %g\n",(double)ctx->B0));
  PetscCall(PetscPrintf(comm,"V_A         = %g\n",(double)ctx->V_A));
  PetscCall(PetscPrintf(comm,"mu0          = %g\n",(double)ctx->mu0));
  PetscCall(PetscPrintf(comm,"mi           = %g\n",(double)ctx->mi));
  PetscCall(PetscPrintf(comm,"eta0         = %g\n",(double)ctx->eta0));
  PetscCall(PetscPrintf(comm,"eta          = %g\n",(double)ctx->eta));

  PetscCall(PetscPrintf(comm,"etawall      = %g\n",(double)ctx->etawall));
  PetscCall(PetscPrintf(comm,"etawallperp  = %g\n",(double)ctx->etawallperp));
  PetscCall(PetscPrintf(comm,"etawallphi   = %g\n",(double)ctx->etawallphi));
  PetscCall(PetscPrintf(comm,"etawallphi_isol_cell = %g\n",(double)ctx->etawallphi_isol_cell));
  PetscCall(PetscPrintf(comm,"etaplasma    = %g\n",(double)ctx->etaplasma));
  PetscCall(PetscPrintf(comm,"etasepwal    = %g\n",(double)ctx->etasepwal));
  PetscCall(PetscPrintf(comm,"etaVV        = %g\n",(double)ctx->etaVV));
  PetscCall(PetscPrintf(comm,"etaout       = %g\n",(double)ctx->etaout));

  PetscCall(PetscPrintf(comm,"\n-- Geometry / mesh --\n"));
  PetscCall(PetscPrintf(comm,"rmin   = %g\n",(double)ctx->rmin));
  PetscCall(PetscPrintf(comm,"rmax   = %g\n",(double)ctx->rmax));
  PetscCall(PetscPrintf(comm,"phimin = %g\n",(double)ctx->phimin));
  PetscCall(PetscPrintf(comm,"phimax = %g\n",(double)ctx->phimax));
  PetscCall(PetscPrintf(comm,"zmin   = %g\n",(double)ctx->zmin));
  PetscCall(PetscPrintf(comm,"zmax   = %g\n",(double)ctx->zmax));

  PetscCall(PetscPrintf(comm,"Nr     = %d\n",(int)ctx->Nr));
  PetscCall(PetscPrintf(comm,"Nphi   = %d\n",(int)ctx->Nphi));
  PetscCall(PetscPrintf(comm,"Nz     = %d\n",(int)ctx->Nz));

  PetscCall(PetscPrintf(comm,"dr     = %g\n",(double)ctx->dr));
  PetscCall(PetscPrintf(comm,"dphi   = %g\n",(double)ctx->dphi));
  PetscCall(PetscPrintf(comm,"dz     = %g\n",(double)ctx->dz));

  PetscCall(PetscPrintf(comm,"\n-- Time stepping --\n"));
  PetscCall(PetscPrintf(comm,"dt      = %g\n",(double)ctx->dt));
  PetscCall(PetscPrintf(comm,"itime   = %g\n",(double)ctx->itime));
  PetscCall(PetscPrintf(comm,"ftime   = %g\n",(double)ctx->ftime));
  PetscCall(PetscPrintf(comm,"tstype  = %d\n",(int)ctx->tstype));
  PetscCall(PetscPrintf(comm,"n_record = %d\n",(int)ctx->n_record));
  PetscCall(PetscPrintf(comm,"n_record_Steady_jRE = %d\n",(int)ctx->n_record_Steady_jRE));
  PetscCall(PetscPrintf(comm,"oldstep = %d\n",(int)ctx->oldstep));

  PetscCall(PetscPrintf(comm,"\n-- Initial/boundary conditions --\n"));
  PetscCall(PetscPrintf(comm,"ictype   = %d\n",(int)ctx->ictype));
  PetscCall(PetscPrintf(comm,"phibtype = %d\n",(int)ctx->phibtype));
  PetscCall(PetscPrintf(comm,"Ebc      = %d\n",(int)ctx->Ebc));

  PetscCall(PetscPrintf(comm,"\n-- Flags --\n"));
  PetscCall(PetscPrintf(comm,"debug       = %d\n",(int)ctx->debug));
  PetscCall(PetscPrintf(comm,"dump        = %d\n",(int)ctx->dump));
  PetscCall(PetscPrintf(comm,"monitor     = %d\n",(int)ctx->monitor));
  PetscCall(PetscPrintf(comm,"prestep     = %d\n",(int)ctx->prestep));
  PetscCall(PetscPrintf(comm,"savecoords  = %d\n",(int)ctx->savecoords));
  PetscCall(PetscPrintf(comm,"savesol     = %d\n",(int)ctx->savesol));
  PetscCall(PetscPrintf(comm,"tempdump    = %d\n",(int)ctx->tempdump));
  PetscCall(PetscPrintf(comm,"dumpfreq    = %d\n",(int)ctx->dumpfreq));
  PetscCall(PetscPrintf(comm,"testSpGD    = %d\n",(int)ctx->testSpGD));
  PetscCall(PetscPrintf(comm,"testSpGDsamerhs = %d\n",(int)ctx->testSpGDsamerhs));
  PetscCall(PetscPrintf(comm,"ic_binary_mode       = %c\n",ctx->ic_binary_mode));
  PetscCall(PetscPrintf(comm,"ic_binary_path       = %s\n",ctx->ic_binary_path));

  PetscCall(PetscPrintf(comm,"\n-- Counters / currents --\n"));
  PetscCall(PetscPrintf(comm,"Iphi1 = %g\n",(double)PetscRealPart(ctx->Iphi1)));
  PetscCall(PetscPrintf(comm,"Iphi2 = %g\n",(double)PetscRealPart(ctx->Iphi2)));
  PetscCall(PetscPrintf(comm,"Iphi3 = %g\n",(double)PetscRealPart(ctx->Iphi3)));
  PetscCall(PetscPrintf(comm,"dampV = %g\n",(double)PetscRealPart(ctx->dampV)));
  PetscCall(PetscPrintf(comm,"Re    = %g\n",(double)PetscRealPart(ctx->Re)));

  PetscCall(PetscPrintf(comm,"\n-- Grid / coord DM --\n"));
  PetscCall(PetscPrintf(comm,"coorda        = %s\n", ctx->coorda ? "set" : "NULL"));

  PetscCall(PetscPrintf(comm,"\n-- PETSc Vec/IS/Mat/TS objects (pointer presence only) --\n"));
  PetscCall(PetscPrintf(comm,"X0          = %s\n", ctx->X0 ? "set" : "NULL"));

  PetscCall(PetscPrintf(comm,"isV           = %s\n", ctx->isV ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"isni          = %s\n", ctx->isni ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"isB           = %s\n", ctx->isB ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"isEP          = %s\n", ctx->isEP ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"istau         = %s\n", ctx->istau ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"isE_boundary  = %s\n", ctx->isE_boundary ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"isB_boundary  = %s\n", ctx->isB_boundary ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"isni_boundary = %s\n", ctx->isni_boundary ? "set" : "NULL"));

  PetscCall(PetscPrintf(comm,"DiagMe        = %s\n", ctx->DiagMe ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"MeVec         = %s\n", ctx->MeVec ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"DiagMe1       = %s\n", ctx->DiagMe1 ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"MeVec1        = %s\n", ctx->MeVec1 ? "set" : "NULL"));

  PetscCall(PetscPrintf(comm,"ts            = %s\n", ctx->ts ? "set" : "NULL"));

  PetscCall(PetscPrintf(comm,"DiagBlock_V   = %s\n", ctx->DiagBlock_V ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"DiagBlock_EP  = %s\n", ctx->DiagBlock_EP ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"DiagBlock_B   = %s\n", ctx->DiagBlock_B ? "set" : "NULL"));

  PetscCall(PetscPrintf(comm,"KSP_V         = %s\n", ctx->KSP_V ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"KSP_B         = %s\n", ctx->KSP_B ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"KSP_EP        = %s\n", ctx->KSP_EP ? "set" : "NULL"));

  PetscCall(PetscPrintf(comm,"isALL_V       = %s\n", ctx->isALL_V ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"OffDiagBlock_U= %s\n", ctx->OffDiagBlock_U ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"OffDiagBlock_L= %s\n", ctx->OffDiagBlock_L ? "set" : "NULL"));

  PetscCall(PetscPrintf(comm,"\n-- Raw data arrays (pointer presence only) --\n"));
  PetscCall(PetscPrintf(comm,"dataC      = %s\n", ctx->dataC ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"dataz      = %s\n", ctx->dataz ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"dataphi    = %s\n", ctx->dataphi ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"datar      = %s\n", ctx->datar ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"datag      = %s\n", ctx->datag ? "set" : "NULL"));
  PetscCall(PetscPrintf(comm,"datapsi    = %s\n", ctx->datapsi ? "set" : "NULL"));

  PetscCall(PetscPrintf(comm,"numz       = %d\n",(int)ctx->numz));
  PetscCall(PetscPrintf(comm,"numphi     = %d\n",(int)ctx->numphi));
  PetscCall(PetscPrintf(comm,"numr       = %d\n",(int)ctx->numr));
  PetscCall(PetscPrintf(comm,"numg       = %d\n",(int)ctx->numg));
  PetscCall(PetscPrintf(comm,"numpsi     = %d\n",(int)ctx->numpsi));

  PetscCall(PetscPrintf(comm,"\n-- Input folder --\n"));
  PetscCall(PetscPrintf(comm,"input_folder = %s\n", ctx->input_folder));

  PetscCall(PetscPrintf(comm,"================ End AppCtx dump ============\n"));

  PetscFunctionReturn(PETSC_SUCCESS);
}

