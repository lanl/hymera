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

  PetscPrintf(comm,"================ AppCtx dump ================\n");

  PetscPrintf(comm,"-- Physical / normalization parameters --\n");
  PetscPrintf(comm,"density      = %g\n",(double)ctx->density);
  PetscPrintf(comm,"L0           = %g\n",(double)ctx->L0);
  PetscPrintf(comm,"B0           = %g\n",(double)ctx->B0);
  PetscPrintf(comm,"V_A         = %g\n",(double)ctx->V_A);
  PetscPrintf(comm,"mu0          = %g\n",(double)ctx->mu0);
  PetscPrintf(comm,"mi           = %g\n",(double)ctx->mi);
  PetscPrintf(comm,"eta0         = %g\n",(double)ctx->eta0);
  PetscPrintf(comm,"eta          = %g\n",(double)ctx->eta);

  PetscPrintf(comm,"etawall      = %g\n",(double)ctx->etawall);
  PetscPrintf(comm,"etawallperp  = %g\n",(double)ctx->etawallperp);
  PetscPrintf(comm,"etawallphi   = %g\n",(double)ctx->etawallphi);
  PetscPrintf(comm,"etawallphi_isol_cell = %g\n",(double)ctx->etawallphi_isol_cell);
  PetscPrintf(comm,"etaplasma    = %g\n",(double)ctx->etaplasma);
  PetscPrintf(comm,"etasepwal    = %g\n",(double)ctx->etasepwal);
  PetscPrintf(comm,"etaVV        = %g\n",(double)ctx->etaVV);
  PetscPrintf(comm,"etaout       = %g\n",(double)ctx->etaout);

  PetscPrintf(comm,"\n-- Geometry / mesh --\n");
  PetscPrintf(comm,"rmin   = %g\n",(double)ctx->rmin);
  PetscPrintf(comm,"rmax   = %g\n",(double)ctx->rmax);
  PetscPrintf(comm,"phimin = %g\n",(double)ctx->phimin);
  PetscPrintf(comm,"phimax = %g\n",(double)ctx->phimax);
  PetscPrintf(comm,"zmin   = %g\n",(double)ctx->zmin);
  PetscPrintf(comm,"zmax   = %g\n",(double)ctx->zmax);

  PetscPrintf(comm,"Nr     = %d\n",(int)ctx->Nr);
  PetscPrintf(comm,"Nphi   = %d\n",(int)ctx->Nphi);
  PetscPrintf(comm,"Nz     = %d\n",(int)ctx->Nz);

  PetscPrintf(comm,"dr     = %g\n",(double)ctx->dr);
  PetscPrintf(comm,"dphi   = %g\n",(double)ctx->dphi);
  PetscPrintf(comm,"dz     = %g\n",(double)ctx->dz);

  PetscPrintf(comm,"\n-- Time stepping --\n");
  PetscPrintf(comm,"dt      = %g\n",(double)ctx->dt);
  PetscPrintf(comm,"itime   = %g\n",(double)ctx->itime);
  PetscPrintf(comm,"ftime   = %g\n",(double)ctx->ftime);
  PetscPrintf(comm,"tstype  = %d\n",(int)ctx->tstype);
  PetscPrintf(comm,"pred_loop = %d\n",(int)ctx->pred_loop);
  PetscPrintf(comm,"adaptdt = %d\n",(int)ctx->adaptdt);
  PetscPrintf(comm,"n_record = %d\n",(int)ctx->n_record);
  PetscPrintf(comm,"n_record_Steady_jRE = %d\n",(int)ctx->n_record_Steady_jRE);
  PetscPrintf(comm,"oldstep = %d\n",(int)ctx->oldstep);

  PetscPrintf(comm,"\n-- Initial/boundary conditions --\n");
  PetscPrintf(comm,"ictype   = %d\n",(int)ctx->ictype);
  PetscPrintf(comm,"phibtype = %d\n",(int)ctx->phibtype);
  PetscPrintf(comm,"Ebc      = %d\n",(int)ctx->Ebc);

  PetscPrintf(comm,"\n-- Flags --\n");
  PetscPrintf(comm,"debug       = %d\n",(int)ctx->debug);
  PetscPrintf(comm,"dump        = %d\n",(int)ctx->dump);
  PetscPrintf(comm,"prestep     = %d\n",(int)ctx->prestep);
  PetscPrintf(comm,"savecoords  = %d\n",(int)ctx->savecoords);
  PetscPrintf(comm,"savesol     = %d\n",(int)ctx->savesol);
  PetscPrintf(comm,"tempdump    = %d\n",(int)ctx->tempdump);
  PetscPrintf(comm,"dumpfreq    = %d\n",(int)ctx->dumpfreq);
  PetscPrintf(comm,"testSpGD    = %d\n",(int)ctx->testSpGD);
  PetscPrintf(comm,"testSpGDsamerhs = %d\n",(int)ctx->testSpGDsamerhs);
  PetscPrintf(comm,"EnableRelaxation       = %d\n",(int)ctx->EnableRelaxation);
  PetscPrintf(comm,"EnableReadICFromBinary = %d\n",(int)ctx->EnableReadICFromBinary);
  PetscPrintf(comm,"delay_kinetic          = %d\n",(int)ctx->delay_kinetic);
  PetscPrintf(comm,"poincare_counter       = %d\n",(int)ctx->poincare_counter);
  PetscPrintf(comm,"field_counter          = %d\n",(int)ctx->field_counter);
  PetscPrintf(comm,"CorrectorIdentifier    = %d\n",(int)ctx->CorrectorIdentifier);

  PetscPrintf(comm,"\n-- Counters / currents --\n");
  PetscPrintf(comm,"Iphi1 = %g\n",(double)PetscRealPart(ctx->Iphi1));
  PetscPrintf(comm,"Iphi2 = %g\n",(double)PetscRealPart(ctx->Iphi2));
  PetscPrintf(comm,"Iphi3 = %g\n",(double)PetscRealPart(ctx->Iphi3));
  PetscPrintf(comm,"dampV = %g\n",(double)PetscRealPart(ctx->dampV));
  PetscPrintf(comm,"Re    = %g\n",(double)PetscRealPart(ctx->Re));
  PetscPrintf(comm,"prev_current    = %g\n",ctx->prev_current);
  PetscPrintf(comm,"present_current = %g\n",ctx->present_current);

  PetscPrintf(comm,"\n-- Grid / coord DM --\n");
  PetscPrintf(comm,"coorda        = %s\n", ctx->coorda ? "set" : "NULL");

  PetscPrintf(comm,"\n-- PETSc Vec/IS/Mat/TS objects (pointer presence only) --\n");
  PetscPrintf(comm,"oldX          = %s\n", ctx->oldX ? "set" : "NULL");
  PetscPrintf(comm,"X_star        = %s\n", ctx->X_star ? "set" : "NULL");

  PetscPrintf(comm,"isV           = %s\n", ctx->isV ? "set" : "NULL");
  PetscPrintf(comm,"isni          = %s\n", ctx->isni ? "set" : "NULL");
  PetscPrintf(comm,"isB           = %s\n", ctx->isB ? "set" : "NULL");
  PetscPrintf(comm,"isEP          = %s\n", ctx->isEP ? "set" : "NULL");
  PetscPrintf(comm,"istau         = %s\n", ctx->istau ? "set" : "NULL");
  PetscPrintf(comm,"isE_boundary  = %s\n", ctx->isE_boundary ? "set" : "NULL");
  PetscPrintf(comm,"isB_boundary  = %s\n", ctx->isB_boundary ? "set" : "NULL");
  PetscPrintf(comm,"isni_boundary = %s\n", ctx->isni_boundary ? "set" : "NULL");

  PetscPrintf(comm,"DiagMe        = %s\n", ctx->DiagMe ? "set" : "NULL");
  PetscPrintf(comm,"MeVec         = %s\n", ctx->MeVec ? "set" : "NULL");
  PetscPrintf(comm,"DiagMe1       = %s\n", ctx->DiagMe1 ? "set" : "NULL");
  PetscPrintf(comm,"MeVec1        = %s\n", ctx->MeVec1 ? "set" : "NULL");

  PetscPrintf(comm,"ts            = %s\n", ctx->ts ? "set" : "NULL");

  PetscPrintf(comm,"DiagBlock_V   = %s\n", ctx->DiagBlock_V ? "set" : "NULL");
  PetscPrintf(comm,"DiagBlock_EP  = %s\n", ctx->DiagBlock_EP ? "set" : "NULL");
  PetscPrintf(comm,"DiagBlock_B   = %s\n", ctx->DiagBlock_B ? "set" : "NULL");

  PetscPrintf(comm,"KSP_V         = %s\n", ctx->KSP_V ? "set" : "NULL");
  PetscPrintf(comm,"KSP_B         = %s\n", ctx->KSP_B ? "set" : "NULL");
  PetscPrintf(comm,"KSP_EP        = %s\n", ctx->KSP_EP ? "set" : "NULL");

  PetscPrintf(comm,"isALL_V       = %s\n", ctx->isALL_V ? "set" : "NULL");
  PetscPrintf(comm,"OffDiagBlock_U= %s\n", ctx->OffDiagBlock_U ? "set" : "NULL");
  PetscPrintf(comm,"OffDiagBlock_L= %s\n", ctx->OffDiagBlock_L ? "set" : "NULL");

  PetscPrintf(comm,"\n-- Raw data arrays (pointer presence only) --\n");
  PetscPrintf(comm,"dataC      = %s\n", ctx->dataC ? "set" : "NULL");
  PetscPrintf(comm,"dataz      = %s\n", ctx->dataz ? "set" : "NULL");
  PetscPrintf(comm,"dataphi    = %s\n", ctx->dataphi ? "set" : "NULL");
  PetscPrintf(comm,"datar      = %s\n", ctx->datar ? "set" : "NULL");
  PetscPrintf(comm,"datag      = %s\n", ctx->datag ? "set" : "NULL");
  PetscPrintf(comm,"datapsi    = %s\n", ctx->datapsi ? "set" : "NULL");

  PetscPrintf(comm,"numz       = %d\n",(int)ctx->numz);
  PetscPrintf(comm,"numphi     = %d\n",(int)ctx->numphi);
  PetscPrintf(comm,"numr       = %d\n",(int)ctx->numr);
  PetscPrintf(comm,"numg       = %d\n",(int)ctx->numg);
  PetscPrintf(comm,"numpsi     = %d\n",(int)ctx->numpsi);

  PetscPrintf(comm,"jre_data   = %s\n", ctx->jre_data ? "set" : "NULL");
  PetscPrintf(comm,"field_data = %s\n", ctx->field_data ? "set" : "NULL");
  PetscPrintf(comm,"jre        = %s\n", ctx->jre ? "set" : "NULL");
  PetscPrintf(comm,"jreR       = %s\n", ctx->jreR ? "set" : "NULL");
  PetscPrintf(comm,"jrephi     = %s\n", ctx->jrephi ? "set" : "NULL");
  PetscPrintf(comm,"jreZ       = %s\n", ctx->jreZ ? "set" : "NULL");

  PetscPrintf(comm,"\n-- Misc / other pointers --\n");
  PetscPrintf(comm,"field_interpolation = %s\n", ctx->field_interpolation ? "set" : "NULL");
  PetscPrintf(comm,"manager             = %s\n", ctx->manager ? "set" : "NULL");
  PetscPrintf(comm,"axis[0], axis[1]    = %g, %g\n",
              ctx->axis[0], ctx->axis[1]);
  PetscPrintf(comm,"ParticlesCreated    = %d\n", ctx->ParticlesCreated);

  PetscPrintf(comm,"\n-- NT interpolation array --\n");
  {
    int i,j;
    for (i=0; i<6; i++) {
      PetscPrintf(comm,"NT[%d][:] =", i);
      for (j=0; j<6; j++) {
        PetscPrintf(comm," %g",(double)PetscRealPart(ctx->NT[i][j]));
      }
      PetscPrintf(comm,"\n");
    }
  }

  PetscPrintf(comm,"\n-- Input folder --\n");
  PetscPrintf(comm,"input_folder = %s\n", ctx->input_folder);

  PetscPrintf(comm,"================ End AppCtx dump ============\n");

  PetscFunctionReturn(PETSC_SUCCESS);
}

