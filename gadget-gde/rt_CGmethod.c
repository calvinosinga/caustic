#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef RADTRANSFER

#define MAX_ITER 10000
#define ACCURACY 1.0e-2
#define EPSILON 1.0e-5
#define tiny 1e-8

/*structures for radtransfer*/
struct radtransferdata_in
{
  int NodeList[NODELISTLENGTH];
  MyDouble Pos[3];
  MyFloat Hsml;
  MyFloat ET[6];
  double Kappa, Lambda;
  MyFloat Mass, Density;
}
 *RadTransferDataIn, *RadTransferDataGet;

struct radtransferdata_out
{
  double Out, Sum;
}
 *RadTransferDataResult, *RadTransferDataOut;

static double *XVec;
static double *QVec, *DVec, *Residue, *Zvec;
static double *Kappa, *Lambda, *Diag, *Diag2;
static double c_light, dt, a3inv, hubble_a;

void radtransfer(void)
{
  int i, j, iter;
  double alpha_cg, beta, delta_old, delta_new, sum, min_diag, glob_min_diag, max_diag, glob_max_diag;
  double nH;
  double rel, res, maxrel, glob_maxrel;
  double DQ;

  c_light = C / All.UnitVelocity_in_cm_per_s;

  /*  the actual time-step we need to do */
  dt = (All.Radiation_Ti_endstep - All.Radiation_Ti_begstep) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      a3inv = 1 / (All.Time * All.Time * All.Time);
      hubble_a = hubble_function(All.Time);
      /* in comoving case, timestep is dloga at this point. Convert to dt */
      dt /= hubble_a;
    }
  else
    {
      a3inv = hubble_a = 1.0;
    }

  XVec = (double *) mymalloc("XVec", N_gas * sizeof(double));
  QVec = (double *) mymalloc("QVec", N_gas * sizeof(double));
  DVec = (double *) mymalloc("DVec", N_gas * sizeof(double));
  Residue = (double *) mymalloc("Residue", N_gas * sizeof(double));
  Kappa = (double *) mymalloc("Kappa", N_gas * sizeof(double));
  Lambda = (double *) mymalloc("Lambda", N_gas * sizeof(double));
  Diag = (double *) mymalloc("Diag", N_gas * sizeof(double));
  Zvec = (double *) mymalloc("Zvec", N_gas * sizeof(double));
  Diag2 = (double *) mymalloc("Diag2", N_gas * sizeof(double));

  for(i = 0; i < N_BINS; i++)
    {
      /* initialization for the CG method */
      for(j = 0; j < N_gas; j++)
	if(P[j].Type == 0)
	  {
	    XVec[j] = 0;
	    QVec[j] = 0;
	    DVec[j] = 0;
	    Residue[j] = 0;
	    Kappa[j] = 0;
	    Lambda[j] = 0;
	    Diag[j] = 0;
	    Zvec[j] = 0;
	    Diag2[j] = 0;

	    XVec[j] = SphP[j].n_gamma[i];

#ifndef UM_CHEMISTRY
	    nH = HYDROGEN_MASSFRAC * SphP[j].d.Density / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
	    Kappa[j] = a3inv * (SphP[j].HI + tiny) * nH * rt_sigma_HI[i];

#if defined(RT_INCLUDE_HE) && defined(RT_MULTI_FREQUENCY)
	    Kappa[j] += a3inv * ((SphP[j].HeI + tiny) * nH * rt_sigma_HeI[i] +
				(SphP[j].HeII + tiny) * nH * rt_sigma_HeII[i]);
#endif	    
	    
#else
	    double denom = 4 * SphP[j].HeI + 
	      1 * SphP[j].HII + 
	      4 * SphP[j].HeIII + 
	      1 * SphP[i].HI + 
	      4 * SphP[j].HeII + 
	      1 * SphP[j].HM +
	      2 * SphP[j].H2II + 
	      2 * SphP[j].H2I +	    
	      5. * SphP[j].HeHII + 
	      3. * SphP[j].HD + 
	      2. * SphP[j].DI + 
	      2. * SphP[j].DII;
	    
	    nH = SphP[j].d.Density * All.UnitMass_in_g / All.HubbleParam / PROTONMASS / denom;
	    nH /= All.UnitLength_in_cm * All.UnitLength_in_cm / All.HubbleParam / All.HubbleParam;
	    
	    Kappa[j] = a3inv * (SphP[j].HI * nH * sigma24[i] +
				SphP[j].HeI * nH * sigma26[i] +
				SphP[j].HeII * nH * sigma25[i] +
				SphP[j].HM * nH * sigma27[i] +
                                SphP[j].H2I * nH * (sigma29[i] + sigma31[i]) +
                                SphP[j].H2II * nH * (sigma28[i] + sigma30[i])) + tiny;
#endif
	    
	    if(All.ComovingIntegrationOn)
	      Kappa[j] *= All.Time;

#ifdef RADTRANSFER_FLUXLIMITER
	    /* now calculate flux limiter */

	    if(SphP[j].n_gamma[i] > 0)
	      {
		double R = sqrt(SphP[j].Grad_ngamma[0][i] * SphP[j].Grad_ngamma[0][i] +
				SphP[j].Grad_ngamma[1][i] * SphP[j].Grad_ngamma[1][i] +
				SphP[j].Grad_ngamma[2][i] * SphP[j].Grad_ngamma[2][i]) / (SphP[j].n_gamma[i] *
											  Kappa[j]);

		if(All.ComovingIntegrationOn)
		  R /= All.Time;

		R *= 0.1;

		Lambda[j] = (1 + R) / (1 + R + R * R);
		if(Lambda[j] < 1e-100)
		  Lambda[j] = 0;
	      }
	    else
	      Lambda[j] = 1.0;
#endif

	    /* add the source term */
	    SphP[j].n_gamma[i] += dt * SphP[j].Je[i] * P[j].Mass;
	  }
      
      radtransfer_matrix_multiply(XVec, Residue, Diag);

      /* Let's take the diagonal matrix elements as Jacobi preconditioner */

      for(j = 0, min_diag = MAX_REAL_NUMBER, max_diag = -MAX_REAL_NUMBER; j < N_gas; j++)
	if(P[j].Type == 0)
	  {
	    Residue[j] = SphP[j].n_gamma[i] - Residue[j];

	    /* note: in principle we would have to substract the w_ii term, but this is always zero */
	    if(Diag[j] < min_diag)
	      min_diag = Diag[j];
	    if(Diag[j] > max_diag)
	      max_diag = Diag[j];

	    Zvec[j] = Residue[j] / Diag[j];
	    DVec[j] = Zvec[j];
	  }

      MPI_Allreduce(&min_diag, &glob_min_diag, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&max_diag, &glob_max_diag, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      delta_new = radtransfer_vector_multiply(Zvec, Residue);
      delta_old = delta_new;

#ifndef RT_SILENT
      if(ThisTask == 0)
	{ 
	  printf("Begin N_BIN %d/%d\n", i+1, N_BINS);
	  printf("\nBegin CG iteration\nmin-diagonal=%g, max-diagonal=%g\n",
	       glob_min_diag, glob_max_diag);
	}
#endif


      /* begin the CG method iteration */
      iter = 0;

      do
	{
	  radtransfer_matrix_multiply(DVec, QVec, Diag2);

	  DQ = radtransfer_vector_multiply(DVec, QVec);
	  if(DQ == 0)
	    alpha_cg = 0;
	  else
	    alpha_cg = delta_new / DQ;


	  for(j = 0, maxrel = 0; j < N_gas; j++)
	    {
	      XVec[j] += alpha_cg * DVec[j];
	      Residue[j] -= alpha_cg * QVec[j];

	      Zvec[j] = Residue[j] / Diag[j];

	      rel = fabs(alpha_cg * DVec[j]) / (XVec[j] + 1.0e-10);
	      if(rel > maxrel)
		maxrel = rel;
	    }

	  delta_old = delta_new;
	  delta_new = radtransfer_vector_multiply(Zvec, Residue);

	  sum = radtransfer_vector_sum(XVec);
	  res = radtransfer_vector_sum(Residue);

	  if(delta_old)
	    beta = delta_new / delta_old;
	  else
	    beta = 0;

	  for(j = 0; j < N_gas; j++)
	    DVec[j] = Zvec[j] + beta * DVec[j];

	  MPI_Allreduce(&maxrel, &glob_maxrel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#ifndef RT_SILENT
	  if(ThisTask == 0)
	    {
	      printf("radtransfer: iter=%3d  |res|/|x|=%12.6g  maxrel=%12.6g  |x|=%12.6g | res|=%12.6g\n",
		     iter, res / sum, glob_maxrel, sum, res);
	      fflush(stdout);
	    }
#endif
	  iter++;

	  if(iter >= MAX_ITER)
            terminate("failed to converge radtransfer\n");
	}
      while((res > ACCURACY * sum && iter < MAX_ITER) || iter < 2);

#ifdef RT_SILENT
      if(ThisTask == 0)
	{
	  printf("%d iterations performed\n", iter);
	  fflush(stdout);
	}
#endif

      /* update the intensity */
      for(j = 0; j < N_gas; j++)
	if(P[j].Type == 0)
	  {
	    if(XVec[j] < 0)
	      XVec[j] = 0;

	    SphP[j].n_gamma[i] = XVec[j];
	  }
    }

  myfree(Diag2);
  myfree(Zvec);
  myfree(Diag);
  myfree(Lambda);
  myfree(Kappa);
  myfree(Residue);
  myfree(DVec);
  myfree(QVec);
  myfree(XVec);

}

/* internal product of two vectors */
double radtransfer_vector_multiply(double *a, double *b)
{
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      sum += a[i] * b[i];

  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}


double radtransfer_vector_sum(double *a)
{
  int i;
  double sum, sumall;

  for(i = 0, sum = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      sum += fabs(a[i]);

  MPI_Allreduce(&sum, &sumall, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sumall;
}


/* this function computes the vector b(out) given the vector x(in) such as Ax = b, where A is a matrix */
void radtransfer_matrix_multiply(double *in, double *out, double *sum)
{
  int i, j, k, ngrp, dummy, ndone, ndone_flag;
  int recvTask, nexport, nimport, place;
  double ainv, dt;

  /* allocate buffers to arrange communication */

  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct radtransferdata_in) +
					     sizeof(struct radtransferdata_out) +
					     sizemax(sizeof(struct radtransferdata_in),
						     sizeof(struct radtransferdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  dt = (All.Radiation_Ti_endstep - All.Radiation_Ti_begstep) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      ainv = 1.0 / All.Time;
      /* in comoving case, timestep is dloga at this point. Convert to dt */
      dt /= hubble_function(All.Time);
    }
  else
    {
      ainv = 1.0;
    }

  i = 0;

  do				/* communication loop */
    {

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      for(nexport = 0; i < N_gas; i++)
	{
	  if(P[i].Type == 0)
	    if(radtransfer_evaluate(i, 0, in, out, sum, &nexport, Send_count) < 0)
	      break;
	}

#ifdef MYSORT
      mysort_dataindex(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#else
      qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
#endif

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      RadTransferDataGet =
	(struct radtransferdata_in *) mymalloc("RadTransferDataGet",
					       nimport * sizeof(struct radtransferdata_in));
      RadTransferDataIn =
	(struct radtransferdata_in *) mymalloc("RadTransferDataIn",
					       nexport * sizeof(struct radtransferdata_in));

      /* prepare particle data for export */

      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  for(k = 0; k < 3; k++)
	    {
	      RadTransferDataIn[j].Pos[k] = P[place].Pos[k];
	      RadTransferDataIn[j].ET[k] = SphP[place].ET[k];
	      RadTransferDataIn[j].ET[k + 3] = SphP[place].ET[k + 3];
	    }
	  RadTransferDataIn[j].Hsml = PPP[place].Hsml;
	  RadTransferDataIn[j].Kappa = Kappa[place];
	  RadTransferDataIn[j].Lambda = Lambda[place];
	  RadTransferDataIn[j].Mass = P[place].Mass;
	  RadTransferDataIn[j].Density = SphP[place].d.Density;

	  memcpy(RadTransferDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&RadTransferDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct radtransferdata_in), MPI_BYTE,
			       recvTask, TAG_RT_A,
			       &RadTransferDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct radtransferdata_in), MPI_BYTE,
			       recvTask, TAG_RT_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      myfree(RadTransferDataIn);
      RadTransferDataResult =
	(struct radtransferdata_out *) mymalloc("RadTransferDataResult",
						nimport * sizeof(struct radtransferdata_out));
      RadTransferDataOut =
	(struct radtransferdata_out *) mymalloc("RadTransferDataOut",
						nexport * sizeof(struct radtransferdata_out));

      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
	radtransfer_evaluate(j, 1, in, out, sum, &dummy, &dummy);

      if(i < N_gas)
	ndone_flag = 0;
      else
	ndone_flag = 1;

      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      /* get the result */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&RadTransferDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct radtransferdata_out),
			       MPI_BYTE, recvTask, TAG_RT_B,
			       &RadTransferDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct radtransferdata_out),
			       MPI_BYTE, recvTask, TAG_RT_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  out[place] += RadTransferDataOut[j].Out;
	  sum[place] += RadTransferDataOut[j].Sum;
	}

      myfree(RadTransferDataOut);
      myfree(RadTransferDataResult);
      myfree(RadTransferDataGet);

    }
  while(ndone < NTask);

  /* do final operations on results */
  for(i = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      {
	/* divide c_light by a to get comoving speed of light (because kappa is comoving) */
	if((1 + dt * c_light * ainv * Kappa[i] + sum[i]) < 0)
	  {
	    printf("1 + sum + rate= %g   sum=%g rate=%g i =%d\n",
		   1 + dt * c_light * ainv * Kappa[i] + sum[i], sum[i], dt * c_light * ainv * Kappa[i], i);
	    endrun(11111111);
	  }

	sum[i] += 1.0 + dt * c_light * ainv * Kappa[i];

	out[i] += in[i] * sum[i];
      }

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);
}

/* this function evaluates parts of the matrix A */
int radtransfer_evaluate(int target, int mode, double *in, double *out, double *sum, int *nexport,
			 int *nsend_local)
{
  int startnode, numngb, listindex = 0;
  int j, n, k;
  MyFloat *ET_aux, ET_j[6], ET_i[6], ET_ij[6];
  MyFloat kappa_i, kappa_j, kappa_ij;

#ifdef RADTRANSFER_FLUXLIMITER
  MyFloat lambda_i, lambda_j;
#endif
  MyDouble *pos;
  MyFloat mass, mass_i, rho, rho_i;
  double sum_out = 0, sum_w = 0, fac = 0;

  double dx, dy, dz;
  double h_j, hinv, hinv4, h_i;
  double dwk_i, dwk_j, dwk;
  double r, r2, r3inv, u, ainv, dt;

  dt = (All.Radiation_Ti_endstep - All.Radiation_Ti_begstep) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      ainv = 1.0 / All.Time;
      /* in comoving case, timestep is dloga at this point. Convert to dt */
      dt /= hubble_function(All.Time);
    }
  else
    {
      ainv = 1.0;
    }

  if(mode == 0)
    {
      ET_aux = SphP[target].ET;
      pos = P[target].Pos;
      h_i = PPP[target].Hsml;
      kappa_i = Kappa[target];
#ifdef RADTRANSFER_FLUXLIMITER
      lambda_i = Lambda[target];
#endif
      mass_i = P[target].Mass;
      rho_i = SphP[target].d.Density;
    }
  else
    {
      ET_aux = RadTransferDataGet[target].ET;
      pos = RadTransferDataGet[target].Pos;
      h_i = RadTransferDataGet[target].Hsml;
      kappa_i = RadTransferDataGet[target].Kappa;
#ifdef RADTRANSFER_FLUXLIMITER
      lambda_i = RadTransferDataGet[target].Lambda;
#endif
      mass_i = RadTransferDataGet[target].Mass;
      rho_i = RadTransferDataGet[target].Density;
    }

#ifdef RADTRANSFER_MODIFY_EDDINGTON_TENSOR
  /*modify Eddington tensor */
  ET_i[0] = 2 * ET_aux[0] - 0.5 * ET_aux[1] - 0.5 * ET_aux[2];
  ET_i[1] = 2 * ET_aux[1] - 0.5 * ET_aux[2] - 0.5 * ET_aux[0];
  ET_i[2] = 2 * ET_aux[2] - 0.5 * ET_aux[0] - 0.5 * ET_aux[1];

  for(k = 3; k < 6; k++)
    ET_i[k] = 2.5 * ET_aux[k];
#else
  for(k = 0; k < 6; k++)
    ET_i[k] = ET_aux[k];
#endif

  if(mode == 0)
    {
      startnode = All.MaxPart;
    }
  else
    {
      startnode = RadTransferDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb = ngb_treefind_pairs(pos, h_i, target, &startnode, mode, nexport, nsend_local);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = Ngblist[n];
	      {
		dx = pos[0] - P[j].Pos[0];
		dy = pos[1] - P[j].Pos[1];
		dz = pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
		if(dx > boxHalf_X)
		  dx -= boxSize_X;
		if(dx < -boxHalf_X)
		  dx += boxSize_X;
		if(dy > boxHalf_Y)
		  dy -= boxSize_Y;
		if(dy < -boxHalf_Y)
		  dy += boxSize_Y;
		if(dz > boxHalf_Z)
		  dz -= boxSize_Z;
		if(dz < -boxHalf_Z)
		  dz += boxSize_Z;
#endif
		r2 = dx * dx + dy * dy + dz * dz;
		r = sqrt(r2);
		r3inv = 1.0 / (r2 * r);
		h_j = PPP[j].Hsml;

		if(r > 0 && (r < h_i || r < h_j))
		  {
		    mass = P[j].Mass;
		    rho = SphP[j].d.Density;
		    kappa_j = Kappa[j];
#ifdef RADTRANSFER_FLUXLIMITER
		    lambda_j = Lambda[j];
#endif

#ifdef RADTRANSFER_MODIFY_EDDINGTON_TENSOR
		    ET_aux = SphP[j].ET;

		    /*modify Eddington tensor */
		    ET_j[0] = 2 * ET_aux[0] - 0.5 * ET_aux[1] - 0.5 * ET_aux[2];
		    ET_j[1] = 2 * ET_aux[1] - 0.5 * ET_aux[2] - 0.5 * ET_aux[0];
		    ET_j[2] = 2 * ET_aux[2] - 0.5 * ET_aux[0] - 0.5 * ET_aux[1];

		    for(k = 3; k < 6; k++)
		      ET_j[k] = 2.5 * ET_aux[k];
#else
		    for(k = 0; k < 6; k++)
		      ET_j[k] = SphP[j].ET[k];
#endif

		    for(k = 0; k < 6; k++)
		      ET_ij[k] = 0.5 * (ET_i[k] + ET_j[k]);

		    if(r < h_i)
		      {
			hinv = 1.0 / h_i;
			hinv4 = hinv * hinv * hinv * hinv;
			u = r * hinv;

			if(u < 0.5)
			  dwk_i = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			else
			  dwk_i = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		      }
		    else
		      dwk_i = 0;

		    if(r < h_j)
		      {
			hinv = 1.0 / h_j;
			hinv4 = hinv * hinv * hinv * hinv;
			u = r * hinv;

			if(u < 0.5)
			  dwk_j = hinv4 * u * (KERNEL_COEFF_3 * u - KERNEL_COEFF_4);
			else
			  dwk_j = hinv4 * KERNEL_COEFF_6 * (1.0 - u) * (1.0 - u);
		      }
		    else
		      dwk_j = 0;

		    kappa_ij = 0.5 * (1 / kappa_i + 1 / kappa_j);
		    dwk = 0.5 * (dwk_i + dwk_j);
		    mass = 0.5 * (mass + mass_i);
		    rho = 0.5 * (rho + rho_i);

		    double tensor = (ET_ij[0] * dx * dx + ET_ij[1] * dy * dy + ET_ij[2] * dz * dz
				     + 2.0 * ET_ij[3] * dx * dy + 2.0 * ET_ij[4] * dy * dz +
				     2.0 * ET_ij[5] * dz * dx);

		    if(tensor > 0)
		      {
			fac = -2.0 * dt * c_light * ainv * (mass / rho) * kappa_ij * dwk * r3inv * tensor;

#ifdef RADTRANSFER_FLUXLIMITER
			fac *= 0.5 * (lambda_i + lambda_j);
#endif

			sum_out -= fac * in[j];

			sum_w += fac;
		      }
		  }
	      }
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = RadTransferDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;
	    }
	}
    }

  if(mode == 0)
    {
      out[target] = sum_out;
      sum[target] = sum_w;
    }
  else
    {
      RadTransferDataResult[target].Out = sum_out;
      RadTransferDataResult[target].Sum = sum_w;
    }

  return 0;
}

/* this function sets up simple initial conditions for a single source in a uniform field of gas with constant density*/
void radtransfer_set_simple_inits(void)
{
  int i, j;

  for(i = 0; i < N_gas; i++)
    if(P[i].Type == 0)
      {
	for(j = 0; j < N_BINS; j++)
	  SphP[i].n_gamma[j] = 0.0;
	
#ifndef UM_CHEMISTRY
	/* in code units */
	SphP[i].HII = tiny;
	SphP[i].HI = 1.0 - SphP[i].HII;
	SphP[i].elec = SphP[i].HII;

#ifdef RT_INCLUDE_HE
	double fac = (1-HYDROGEN_MASSFRAC)/4.0/HYDROGEN_MASSFRAC;

	SphP[i].HeIII = tiny * fac;
	SphP[i].HeII = tiny * fac;
	SphP[i].HeI = (1.0 - SphP[i].HeII - SphP[i].HeIII) * fac;

	SphP[i].elec += SphP[i].HeII + 2.0 * SphP[i].HeIII;
#endif
#endif
      }
}

void rt_get_sigma(void)
{
  double fac = 1.0 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;

#ifndef RT_MULTI_FREQUENCY
  rt_sigma_HI[0] = 6.3e-18 * fac;
  
#else 
  int i, j, integral;
  double e, d_nu, e_start, e_end;
  double sum_HI_sigma, sum_HI_G;
  double hc, T_eff, I_nu;
  double sig, f, fac_two;
#ifdef RT_INCLUDE_HE
  double sum_HeI_sigma, sum_HeII_sigma;
  double sum_HeI_G, sum_HeII_G;
#endif

  T_eff = All.star_Teff;
  hc = C * PLANCK;

  integral = 10000;

  fac_two = ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;

  nu[0] = 13.6;
  nu[1] = 24.6;
  nu[2] = 54.4;
  nu[3] = 70.0;

  sum_HI_sigma = 0.0;
  sum_HI_G = 0.0;
#ifdef RT_INCLUDE_HE
  sum_HeI_G = sum_HeII_G = 0.0;
  sum_HeI_sigma = 0.0;
  sum_HeII_sigma = 0.0;
#endif
  
  for(i = 0; i < N_BINS; i++)
    {
      e_start = nu[i];
      
      if(i == N_BINS - 1)
	e_end = 500.0;
      else
	e_end = nu[i+1];
      
      d_nu = (e_end - e_start) / (float)(integral - 1);
      
      rt_sigma_HI[i] = 0.0;
      G_HI[i] = 0.0;

#ifdef RT_INCLUDE_HE
      rt_sigma_HeI[i] = 0.0;
      rt_sigma_HeII[i] = 0.0;
      G_HeI[i] = G_HeII[i] = 0.0;
#endif	  
	
      for(j = 0; j < integral; j++)
	{
	  e = e_start + j * d_nu;
	  
	  I_nu = 2.0 * pow(e * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc)
	    / (exp(e * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0);
	  
	  if(nu[i] >= 13.6)
	    {
	      f = sqrt((e / 13.6) - 1.0);
	  
	      if(j == 0)
		sig = 6.3e-18;
	      else
		sig = 6.3e-18 * pow(13.6 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));
	      
	      rt_sigma_HI[i] += d_nu * sig * I_nu / e;

	      sum_HI_sigma += d_nu * I_nu / e; 

	      G_HI[i] += d_nu * sig * (e - 13.6) * I_nu / e;

	      sum_HI_G += d_nu * sig * I_nu / e;
	    }
	  
#ifdef RT_INCLUDE_HE
	  if(nu[i] >= 24.6)
	    {
	      f = sqrt((e / 24.6) - 1.0);
	  
	      if(j == 0)
		sig = 7.83e-18;
	      else
		sig = 7.83e-18 * pow(24.6 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));
	      
	      rt_sigma_HeI[i] += d_nu * sig * I_nu / e;
	      
              sum_HeI_sigma += d_nu * I_nu / e;

              G_HeI[i] += d_nu * sig * (e - 24.6) * I_nu / e;

              sum_HeI_G += d_nu * sig * I_nu / e;
	    }
	  
	  if(nu[i] >= 54.4)
	    {
	      f = sqrt((e / 54.4) - 1.0);
	  
	      if(j == 0)
		sig = 1.58e-18;
	      else
		sig = 1.58e-18 * pow(54.4 / e, 4) * exp(4 - (4 * atan(f) / f)) / (1.0 - exp(-2 * M_PI / f));
	      
	      rt_sigma_HeII[i] += d_nu * sig * I_nu / e;

              sum_HeII_sigma += d_nu * I_nu / e;

              G_HeII[i] += d_nu * sig * (e - 54.4) * I_nu / e;

              sum_HeII_G += d_nu * sig * I_nu / e;
	    }
#endif
	}
    }
  for(i = 0; i < N_BINS; i++)
    {      
      if(nu[i] >= 13.6)
	{
	  rt_sigma_HI[i] *= fac / sum_HI_sigma;
	  G_HI[i] *= fac_two / sum_HI_G;
	}
      
#ifdef RT_INCLUDE_HE
      if(nu[i] >= 24.6)
	{
	  rt_sigma_HeI[i] *= fac / sum_HeI_sigma;
	  G_HeI[i] *= fac_two / sum_HeI_G;
	}
      
      if(nu[i] >= 54.4)
	{
	  rt_sigma_HeII[i] *= fac / sum_HeII_sigma;         
	  G_HeII[i] *= fac_two / sum_HeII_G;
	}
#endif
    }

  if(ThisTask == 0)
    for(i = 0; i < N_BINS; i++)
      printf("%g %g | %g %g | %g %g\n", 
	     rt_sigma_HI[i]/fac, G_HI[i]/fac_two, 
	     rt_sigma_HeI[i]/fac, G_HeI[i]/fac_two, 
	     rt_sigma_HeII[i]/fac, G_HeII[i]/fac_two);
  
#endif
}

#ifdef RT_MULTI_FREQUENCY
#if defined(EDDINGTON_TENSOR_STARS) || defined(EDDINGTON_TENSOR_SFR)
void rt_get_lum_stars(void)
{
  int i;
  double T_eff, R_eff, I_nu;
  double hc, sum, d_nu;
#ifndef UM_CHEMISTRY
  int j, integral;
  double e, e_start, e_end;
  integral = 10000;
#endif

  T_eff = All.star_Teff;
  R_eff = 7e11;
  hc = C * PLANCK;

  for(i = 0, sum = 0; i < N_BINS; i++)
    {
#ifndef UM_CHEMISTRY
      e_start = nu[i];

      if(i == N_BINS - 1)
        e_end = 500.0;
      else
        e_end = nu[i+1];
      
      d_nu = (e_end - e_start) / (float)(integral - 1);

      lum[i] = 0.0;

      for(j = 0; j < integral; j++)
        {
          e = e_start + j*d_nu;
	  
	  I_nu = 2.0 * pow(e * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc)
	    / (exp(e * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0);

	  lum[i] += 4.0 * M_PI * R_eff * R_eff * M_PI * I_nu / e * d_nu / PLANCK; // number/s
	}
#else
      I_nu = 2.0 * pow(nu[i] * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc)
	/ (exp(nu[i] * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0);

      lum[i] = 4.0 * M_PI * R_eff * R_eff * M_PI * I_nu / PLANCK; // number/s

      if(i == N_BINS - 1)
	d_nu = 500.0 - nu[i];
      else
	d_nu = nu[i+1] - nu[i];
#endif

#ifndef UM_CHEMISTRY
      sum += lum[i];
#else
      sum += lum[i] * d_nu / nu[i];
#endif
      
      lum[i] *= All.UnitTime_in_s / All.HubbleParam; //number/time
    }
  
#ifdef RT_TEST_SST
#ifdef EDDINGTON_TENSOR_STARS
  for(i = 0; i < N_BINS; i++)
    lum[i] *= All.IonizingLumPerSolarMass / sum;
#endif
#ifdef EDDINGTON_TENSOR_SFR
  for(i = 0; i < N_BINS; i++)
    lum[i] *= All.IonizingLumPerSFR / sum;
#endif
#endif
  
#ifdef RT_POPIII
  T_eff = All.star_Teff_popIII;
  R_eff = 7e11;
  hc = C * PLANCK;

  for(i = 0, sum = 0; i < N_BINS; i++)
    {
#ifndef UM_CHEMISTRY
      e_start = nu[i];

      if(i == N_BINS - 1)
        e_end = 500.0;
      else
        e_end = nu[i+1];
      
      d_nu = (e_end - e_start) / (float)(integral - 1);

      lum_popIII[i] = 0.0;

      for(j = 0; j < integral; j++)
        {
          e = e_start + j*d_nu;
	  
	  I_nu = 2.0 * pow(e * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc)
	    / (exp(e * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0);

	  lum_popIII[i] += 4.0 * M_PI * R_eff * R_eff * M_PI * I_nu / e * d_nu / PLANCK; // number/s
	}
#else
      I_nu = 2.0 * pow(nu[i] * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc)
	/ (exp(nu[i] * ELECTRONVOLT_IN_ERGS / (BOLTZMANN * T_eff)) - 1.0);

      lum_popIII[i] = 4.0 * M_PI * R_eff * R_eff * M_PI * I_nu / PLANCK; // number/s

      if(i == N_BINS - 1)
	d_nu = 500.0 - nu[i];
      else
	d_nu = nu[i+1] - nu[i];
#endif

#ifndef UM_CHEMISTRY
      sum += lum_popIII[i];
#else
      sum += lum_popIII[i] * d_nu / nu[i];
#endif
      
      lum_popIII[i] *= All.UnitTime_in_s / All.HubbleParam; //number/time
    }
  
#ifdef RT_TEST_SST
#ifdef EDDINGTON_TENSOR_STARS
  for(i = 0; i < N_BINS; i++)
    lum_popIII[i] *= All.IonizingLumPerSolarMassPopIII / sum;
#endif
#ifdef EDDINGTON_TENSOR_SFR
  for(i = 0; i < N_BINS; i++)
    lum_popIII[i] *= All.IonizingLumPerSFR / sum;
#endif
#endif
 
#endif // RT_POPIII
 
  if(ThisTask == 0)
    {
      fprintf(FdStar, "T_eff %g\n", All.star_Teff);   
      fflush(FdStar);
      for(i = 0; i < N_BINS; i++)
	{
	  fprintf(FdStar, "%d %g %g\n", i, nu[i], lum[i] / All.UnitTime_in_s * All.HubbleParam);   
	  fflush(FdStar);
	}
#ifdef RT_POPIII
      fprintf(FdStar, "T_eff_popIII %g\n", All.star_Teff_popIII);
      fflush(FdStar);
      for(i = 0; i < N_BINS; i++)
        {
          fprintf(FdStar, "%d %g %g\n", i, nu[i], lum_popIII[i] / All.UnitTime_in_s * All.HubbleParam);
          fflush(FdStar);
        }
#endif
    }

}
#endif
#endif

#if defined(EDDINGTON_TENSOR_GAS) && defined(RT_MULTI_FREQUENCY)
void rt_get_lum_gas(int target, double *je)
{
  int j;
  double temp, molecular_weight;
  double kT, hc, BB_l, BB_r, next;
  double sigma_SB, R_eff;
  double u_cooling, u_BB;
  double dt, dtime, fac;
  double d_nu;

  sigma_SB = 5.6704e-5;

  dt = (All.Radiation_Ti_endstep - All.Radiation_Ti_begstep) * All.Timebase_interval;

  if(All.ComovingIntegrationOn)
    {
      dtime = dt / hubble_function(All.Time);
      a3inv = 1.0 / All.Time / All.Time / All.Time;
    }
  else
    {
      dtime = dt;
      a3inv = 1.0;
    }

  R_eff = 3.0 / 4.0 / M_PI * pow(P[target].Mass / SphP[target].d.Density * a3inv, 1. / 3.);
  R_eff *= All.UnitLength_in_cm / All.HubbleParam; //cm

  molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[target].elec);

  temp =  SphP[target].Entropy *
    pow(SphP[target].d.Density * a3inv, GAMMA_MINUS1) *
    (molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam) /
    (BOLTZMANN / All.UnitEnergy_in_cgs * All.HubbleParam);

  kT = temp * BOLTZMANN;
  hc = C * PLANCK;

  u_cooling = rt_get_cooling_rate(target, SphP[target].Entropy) *
    dtime / (SphP[target].d.Density * a3inv);

  if(!(dtime > 0))
    u_cooling = 0.0;

  u_BB = sigma_SB * pow(temp, 4) * dtime * All.UnitTime_in_s / All.HubbleParam; // erg/cm^2
  u_BB /= All.UnitEnergy_in_cgs / All.HubbleParam; //energy/cm^2
  u_BB *= 4.0 * M_PI * R_eff * R_eff; //energy
  u_BB /= P[target].Mass; //energy/mass

  if(u_BB > 0)
    fac = - u_cooling / u_BB;
  else
    fac = 0.0;

  for(j = 0; j < N_BINS; j++)
    {
      if(j == N_BINS - 1)
        next = 500;
      else
        next = nu[j+1];

      d_nu = next - nu[j];

      BB_r = 2.0 * pow(next * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc)
        / (exp(next * ELECTRONVOLT_IN_ERGS / (kT)) - 1.0);

      BB_l = 2.0 * pow(nu[j] * ELECTRONVOLT_IN_ERGS, 3) / (hc * hc)
        / (exp(nu[j] * ELECTRONVOLT_IN_ERGS / (kT)) - 1.0);

      je[j] += 4.0 * M_PI * R_eff * R_eff * M_PI * 0.5 * (BB_l / nu[j] + BB_r / next) * d_nu / PLANCK; // number/s

      je[j] *= All.UnitTime_in_s / All.HubbleParam; //number/time

      je[j] *= fac;
    }

}
#endif

#endif
