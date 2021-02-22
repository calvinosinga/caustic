#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "allvars.h"
#include "proto.h"
#include "kernel.h"
#ifdef NUM_THREADS
#include <pthread.h>
#endif

#ifdef INTER_SPH_LOOP

#ifdef NUM_THREADS
extern pthread_mutex_t mutex_nexport;
extern pthread_mutex_t mutex_partnodedrift;
#define LOCK_NEXPORT     pthread_mutex_lock(&mutex_nexport);
#define UNLOCK_NEXPORT   pthread_mutex_unlock(&mutex_nexport);
#else
#define LOCK_NEXPORT
#define UNLOCK_NEXPORT
#endif

static double a2inv;
#define NV_MYSIGN(x) (( x > 0 ) - ( x < 0 ))

struct kernel_addSPH
{
  double dx, dy, dz;
  double r;
  double wk_i, wk_j, dwk_i, dwk_j;
  double h_i;
};

struct addSPHdata_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
#ifndef DONOTUSENODELIST
  int NodeList[NODELISTLENGTH];
#endif
}
 *AddSPHDataIn, *AddSPHDataGet;

struct addSPHdata_out
{
#ifdef VISCOSITY_SUPPRESSION
  MyFloat NV_R;
#endif
}
 *AddSPHDataResult, *AddSPHDataOut;


static inline void particle2in_addSPH(struct addSPHdata_in *in, int i);
static inline void out2particle_addSPH(struct addSPHdata_out *out, int i, int mode);

static inline void particle2in_addSPH(struct addSPHdata_in *in, int i)
{
  int k;

  for(k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];
  in->Hsml = PPP[i].Hsml;
}

static inline void out2particle_addSPH(struct addSPHdata_out *out, int i, int mode)
{
#ifdef VISCOSITY_SUPPRESSION
  ASSIGN_ADD(SphP[i].NV_R, out->NV_R, mode);
#endif
}

void additionalSPH_calc(void)
{
  int i, j, k, ngrp, ndone, ndone_flag;
  int recvTask, place;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0, timenetwork = 0;
  double timecomp, timecomm, timewait, tstart, tend, t0, t1;
  int save_NextParticle;
  long long n_exported = 0;

  if(All.ComovingIntegrationOn)
    {
      a2inv= All.Time * All.Time;
    }
  else
    {
      a2inv = 1.0;
    }

  /* allocate buffers to arrange communication */

  int NTaskTimesNumPart;

  NTaskTimesNumPart = maxThreads * NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct addSPHdata_in) +
					     sizeof(struct addSPHdata_out) +
					     sizemax(sizeof(struct addSPHdata_in),
						     sizeof(struct addSPHdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  CPU_Step[CPU_HYDMISC] += measure_time();
  t0 = second();

  NextParticle = FirstActiveParticle;	/* beginn with this index */

  do
    {

      BufferFullFlag = 0;
      Nexport = 0;
      save_NextParticle = NextParticle;

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      tstart = second();

#ifdef NUM_THREADS
      pthread_t mythreads[NUM_THREADS - 1];
      int threadid[NUM_THREADS - 1];
      pthread_attr_t attr;

      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      pthread_mutex_init(&mutex_nexport, NULL);
      pthread_mutex_init(&mutex_partnodedrift, NULL);

      TimerFlag = 0;

      for(j = 0; j < NUM_THREADS - 1; j++)
	{
	  threadid[j] = j + 1;
	  pthread_create(&mythreads[j], &attr, addSPH_evaluate_primary, &threadid[j]);
	}
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
#ifdef _OPENMP
	int mainthreadid = omp_get_thread_num();
#else
	int mainthreadid = 0;
#endif
	addSPH_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
      }

#ifdef NUM_THREADS
      for(j = 0; j < NUM_THREADS - 1; j++)
	pthread_join(mythreads[j], NULL);
#endif


      tend = second();
      timecomp1 += timediff(tstart, tend);

      if(BufferFullFlag)
	{
	  int last_nextparticle = NextParticle;

	  NextParticle = save_NextParticle;

	  while(NextParticle >= 0)
	    {
	      if(NextParticle == last_nextparticle)
		break;

	      if(ProcessedFlag[NextParticle] != 1)
		break;

	      ProcessedFlag[NextParticle] = 2;

	      NextParticle = NextActiveParticle[NextParticle];
	    }

	  if(NextParticle == save_NextParticle)
	    {
	      /* in this case, the buffer is too small to process even a single particle */
	      endrun(12998);
	    }

	  int new_export = 0;

	  for(j = 0, k = 0; j < Nexport; j++)
	    if(ProcessedFlag[DataIndexTable[j].Index] != 2)
	      {
		if(k < j + 1)
		  k = j + 1;

		for(; k < Nexport; k++)
		  if(ProcessedFlag[DataIndexTable[k].Index] == 2)
		    {
		      int old_index = DataIndexTable[j].Index;

		      DataIndexTable[j] = DataIndexTable[k];
		      DataNodeList[j] = DataNodeList[k];
		      DataIndexTable[j].IndexGet = j;
		      new_export++;

		      DataIndexTable[k].Index = old_index;
		      k++;
		      break;
		    }
	      }
	    else
	      new_export++;

	  Nexport = new_export;

	}

      n_exported += Nexport;

      for(j = 0; j < NTask; j++)
	Send_count[j] = 0;
      for(j = 0; j < Nexport; j++)
	Send_count[DataIndexTable[j].Task]++;

      MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);

      tstart = second();

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timewait1 += timediff(tstart, tend);

      for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	{
	  Nimport += Recv_count[j];

	  if(j > 0)
	    {
	      Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
	      Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
	    }
	}

      AddSPHDataGet = (struct addSPHdata_in *) mymalloc("AddSPHDataGet", Nimport * sizeof(struct addSPHdata_in));
      AddSPHDataIn = (struct addSPHdata_in *) mymalloc("AddSPHDataIn", Nexport * sizeof(struct addSPHdata_in));

      /* prepare particle data for export */

      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  particle2in_addSPH(&AddSPHDataIn[j], place);
#ifndef DONOTUSENODELIST
	  memcpy(AddSPHDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
#endif

	}

      /* exchange particle data */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&AddSPHDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct addSPHdata_in), MPI_BYTE,
			       recvTask, TAG_INTERLOOP_A,
			       &AddSPHDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct addSPHdata_in), MPI_BYTE,
			       recvTask, TAG_INTERLOOP_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm1 += timediff(tstart, tend);

      myfree(AddSPHDataIn);
      AddSPHDataResult =
	(struct addSPHdata_out *) mymalloc("AddSPHDataResult", Nimport * sizeof(struct addSPHdata_out));
      AddSPHDataOut =
	(struct addSPHdata_out *) mymalloc("AddSPHDataOut", Nexport * sizeof(struct addSPHdata_out));


      report_memory_usage(&HighMark_addSPH, "SPH_INTERMEDIATE_LOOP");

      /* now do the particles that were sent to us */

      tstart = second();

      NextJ = 0;

#ifdef NUM_THREADS
      for(j = 0; j < NUM_THREADS - 1; j++)
	pthread_create(&mythreads[j], &attr, addSPH_evaluate_secondary, &threadid[j]);
#endif
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
#ifdef _OPENMP
	int mainthreadid = omp_get_thread_num();
#else
	int mainthreadid = 0;
#endif
	addSPH_evaluate_secondary(&mainthreadid);
      }

#ifdef NUM_THREADS
      for(j = 0; j < NUM_THREADS - 1; j++)
	pthread_join(mythreads[j], NULL);

      pthread_mutex_destroy(&mutex_partnodedrift);
      pthread_mutex_destroy(&mutex_nexport);
      pthread_attr_destroy(&attr);
#endif

      tend = second();
      timecomp2 += timediff(tstart, tend);

      if(NextParticle < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

      tstart = second();
      MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      tend = second();
      timewait2 += timediff(tstart, tend);

      /* get the result */
      tstart = second();
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  recvTask = ThisTask ^ ngrp;
	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* send the results */
		  MPI_Sendrecv(&AddSPHDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct addSPHdata_out),
			       MPI_BYTE, recvTask, TAG_INTERLOOP_B,
			       &AddSPHDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct addSPHdata_out),
			       MPI_BYTE, recvTask, TAG_INTERLOOP_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = second();
      timecommsumm2 += timediff(tstart, tend);

      /* add the result to the local particles */
      tstart = second();
      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  out2particle_addSPH(&AddSPHDataOut[j], place, 1);
	}
      tend = second();
      timecomp1 += timediff(tstart, tend);

      myfree(AddSPHDataOut);
      myfree(AddSPHDataResult);
      myfree(AddSPHDataGet);
    }
  while(ndone < NTask);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

  /* do final operations on results */

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 0)
      {
#ifdef VISCOSITY_SUPPRESSION
	SphP[i].NV_R /= SphP[i].d.Density;

	double NV_dummy = fabs(2.0 * pow(1.0 - SphP[i].NV_R,4.0) * SphP[i].NV_DivVel);
	double NV_limiter = NV_dummy*NV_dummy / (NV_dummy*NV_dummy + SphP[i].NV_trSSt); 
  
        double NV_dt =  (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval;
	if(All.ComovingIntegrationOn)
	  NV_dt /= hubble_function(All.Time);

	double NV_A = NV_limiter * DMAX(-SphP[i].NV_dt_DivVel, 0.0);

	SphP[i].alphaloc = All.ViscosityAMax * PPP[i].Hsml*PPP[i].Hsml * NV_A / 
                          (0.25 * SphP[i].MaxSignalVel*SphP[i].MaxSignalVel + PPP[i].Hsml*PPP[i].Hsml * NV_A);

        if(SphP[i].alpha < SphP[i].alphaloc)
	  SphP[i].alpha = SphP[i].alphaloc;
	else if (SphP[i].alpha > SphP[i].alphaloc)
	  SphP[i].alpha = SphP[i].alphaloc + (SphP[i].alpha - SphP[i].alphaloc) * exp (-NV_dt / PPP[i].Hsml * 2.0 * VISCOSITY_SUPPRESSION * 0.5 * SphP[i].MaxSignalVel);

	if(SphP[i].alpha < All.ViscosityAMin)
	  SphP[i].alpha = All.ViscosityAMin;
#endif
     }
  /* collect some timing information */

  t1 = WallclockTime = second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_HYDCOMPUTE] += timecomp;
  CPU_Step[CPU_HYDWAIT] += timewait;
  CPU_Step[CPU_HYDCOMM] += timecomm;
  CPU_Step[CPU_HYDNETWORK] += timenetwork;
  CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm + timenetwork);
}


int addSPH_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		   int *ngblist)
{
  int startnode, numngb, listindex = 0;
  int j, k, n, l;

  /* avoid unused variable compiler warnings */
  (void) l;
  (void) k;

  double hinv, hinv3, hinv4, r2, u;

  struct kernel_addSPH kernel;
  struct addSPHdata_in local;
  struct addSPHdata_out out;
  memset(&out, 0, sizeof(struct addSPHdata_out));

  if(mode == 0)
    particle2in_addSPH(&local, target);
  else
    local = AddSPHDataGet[target];

  kernel.h_i = local.Hsml;

  /* Now start the actual SPH computation for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = AddSPHDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  numngb =
	    ngb_treefind_pairs_threads(local.Pos, kernel.h_i, target, &startnode, mode, exportflag,
				       exportnodecount, exportindex, ngblist);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = ngblist[n];

	      kernel.dx = local.Pos[0] - P[j].Pos[0];
	      kernel.dy = local.Pos[1] - P[j].Pos[1];
	      kernel.dz = local.Pos[2] - P[j].Pos[2];
#ifdef PERIODIC			/*  now find the closest image in the given box size  */
	      kernel.dx = NEAREST_X(kernel.dx);
	      kernel.dy = NEAREST_Y(kernel.dy);
	      kernel.dz = NEAREST_Z(kernel.dz);
#endif
	      
	      r2 = kernel.dx * kernel.dx + kernel.dy * kernel.dy + kernel.dz * kernel.dz;

	      if(r2 < kernel.h_i * kernel.h_i)
		{
		  kernel.r = sqrt(r2);
		  if(kernel.r > 0)
		    {
		      if(kernel.r < kernel.h_i)
			{
			  kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
			  u = kernel.r * hinv;
			  kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, -1);
			}
		      else
			{
			  kernel.wk_i = 0;
			}

#ifdef VISCOSITY_SUPPRESSION
		      out.NV_R += NV_MYSIGN(SphP[j].NV_DivVel) * P[j].Mass * kernel.wk_i;
#endif

		    }
		}
	    }
	}

#ifndef DONOTUSENODELIST
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = AddSPHDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
#endif
    }

  /* Now collect the result at the right place */
  if(mode == 0)
    out2particle_addSPH(&out, target, 0);
  else
    AddSPHDataResult[target] = out;

  return 0;
}


void *addSPH_evaluate_primary(void *p)
{
  int thread_id = *(int *) p;
  int i, j;
  int *exportflag, *exportnodecount, *exportindex, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;
  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

  while(1)
    {
      int exitFlag = 0;
      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	if(BufferFullFlag != 0 || NextParticle < 0)
	  {
	    exitFlag = 1;
	  }
	else
	  {
	    i = NextParticle;
	    ProcessedFlag[i] = 0;
	    NextParticle = NextActiveParticle[NextParticle];
	  }
      }
      UNLOCK_NEXPORT;
      if(exitFlag)
	break;

      if(P[i].Type == 0)
	{
	  if(addSPH_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;
}

void *addSPH_evaluate_secondary(void *p)
{
  int thread_id = *(int *) p;
  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;

  while(1)
    {
      LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
      {
	j = NextJ;
	NextJ++;
      }
      UNLOCK_NEXPORT;

      if(j >= Nimport)
	break;

      addSPH_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;
}

#endif /* INTER_SPH_LOOP */
