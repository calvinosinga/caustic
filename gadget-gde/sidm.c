#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <sys/stat.h>
#include <sys/types.h>


#include "allvars.h"
#include "proto.h"
#include "kernel.h"

#ifdef SIDM

#define SCATTER_STATS	       	     //do some statistics and checks
//#define ENERGY_STATS               //check kinetic energy before and after scatter
#define MAX_NGBS 44
#define SIDM_INELASTIC_MODEL_3
//cross section interpolation (if vel-dep.)
#define CROSS_VBINS 1000000          //bins
#define CROSS_VMIN  1e-2             //v_min [internal units]
#define CROSS_VMAX  1e5              //v_max [internal units]

/* speed of light in cm/s */
#define SPEED_OF_LIGHT 2.99792458e10

/* find smoothing length */
void sidm_Hsml(void);
int sidm_Hsml_evaluate(int target, int mode, int *nexport, int *nsend_local);

/* construct ngb list for scatter particles (<<NumPart) */
void sidm_NgbList(void);
int sidm_NgbList_evaluate(int target, int mode, int *nexport, int *nsend_local);

/* do actual scatter for scatter particles (<<NumPart) */
void sidm_Scatter(void);
int sidm_Scatter_evaluate(int target, int mode, int *nexport, int *nsend_local);

/* modified tree find for DM particles */
int sidm_ngb_treefind_variable(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode, int *nexport, int *nsend_local);

/*init/reinit data structures */
void sidm_ReInit(void);

/* check for scatter */
void sidm_check_particle_scatter(void);

/* select scatter partner */
void sidm_AssignScatterPartner(void);

/* stat functions */
double get_kinetic_energy(void);
void scatter_stats(void);

/* set up cross section interpolation table */
void init_cross_table_elastic(void);
void init_cross_table_inelastic(void);

/* scatter lists for particles that scatter */
int *ScatterParticleList_1, *ScatterParticleList_2;

/* numper of local particles that scatter */
int NumScatterParticles;

/* derived constants */
/* conversion factors internal units*/
static double CrossUnitFac;
#ifdef SIDM_INELASTIC
/* speed of light in internal units */
static double SpeedOfLight;
/* velocity kick due to excited state */
static double SplittingVelocity;
#endif
/* array for cross section interpolation */
#ifndef CONST_CROSS
static double Dvlog;
static double CrossTable[CROSS_VBINS];
#endif

typedef struct
{
  MyIDType NgbIDs;
  MyDouble P0j_half;
  float Distance;
} ngb_entry;

struct ngb_data
{
  int offset;
  ngb_entry E[MAX_NGBS];
}
 *NgbData;

struct densdata_in
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  float Hsml;
  int NodeList[NODELISTLENGTH];
  MyIDType ScatterID;
}
 *DensDataIn, *DensDataGet;

struct hsmldata_out
{
  float Density;
  float VelDisp, Vx, Vy, Vz;
  MyDouble PSum;
  int Ngb;
}
 *HsmlDataResult, *HsmlDataOut;

struct ngbdata_out
{
  int Ngb;
  ngb_entry E[MAX_NGBS];
}
 *NgbDataResult, *NgbDataOut;

struct scatterdata_out
{
  float Vx, Vy, Vz;
  int DidScatterInStep;
}
 *ScatterDataResult, *ScatterDataOut;

static MyFloat *DM_Vx, *DM_Vy, *DM_Vz;

/* ID compare function */
int ID_cmp(const void *a, const void *b)
{
  const MyIDType *ia = (const MyIDType *) a;
  const MyIDType *ib = (const MyIDType *) b;
  if(*ia < *ib)
    return -1;
  if(*ia > *ib)
    return +1;
  return 0;
}

/* ngb struct compare, increasing distance */
int ngb_entry_distance_cmp(const void *a, const void *b)
{
  const ngb_entry *ia = (const ngb_entry *) a;
  const ngb_entry *ib = (const ngb_entry *) b;
  if(ia->Distance < ib->Distance)
    return -1;
  if(ia->Distance > ib->Distance)
    return +1;
  return 0;
}


/* main loop for scattering, called in run.c */
void sidm_DoScatter()
{
  int n;

  sidm_ReInit();

  sidm_Hsml();

  ScatterParticleList_1 = (int *) mymalloc("ScatterParticleList_1", NumPart * sizeof(int));
  ScatterParticleList_2 = (int *) mymalloc("ScatterParticleList_2", NumPart * sizeof(int));

  sidm_check_particle_scatter();

  NgbData = (struct ngb_data *) mymalloc("NgbData", NumScatterParticles * sizeof(struct ngb_data));

  for(n = 0; n < NumScatterParticles; n++)
    memset(&NgbData[n], 0, MAX_NGBS * sizeof(struct ngb_data));

  sidm_NgbList();

  sidm_AssignScatterPartner();

  myfree(NgbData);

#ifndef SIDM_DO_NOT_SCATTER
  sidm_Scatter();
#endif

  myfree(ScatterParticleList_2);
  myfree(ScatterParticleList_1);
}


/* assign scatter partner to particles supposed to scatter */
void sidm_AssignScatterPartner(void)
{
  int n, p, done_flag, i, numdup;
  MyDouble PSum_partial = 0;
  MyIDType *scatterlist_local, *scatterlist_global, *duplicate;
  int *count, *offset, tot_count;
  int TotNumScatterParticles=2*NumScatterParticles;
  
  scatterlist_local = mymalloc("scatterlist_local", TotNumScatterParticles*sizeof(MyIDType));

  /* select scatter partner */
  for(n = 0; n < NumScatterParticles; n++)
    {
      //sort by increasing distance
      qsort(&NgbData[n].E[0], NgbData[n].offset, sizeof(ngb_entry), ngb_entry_distance_cmp); 

      done_flag = 0;
      for(p = 0; p < NgbData[n].offset; p++)
	{
	  PSum_partial += NgbData[n].E[p].P0j_half;
	  if(PSum_partial > P[ScatterParticleList_1[n]].RandX)
	    {
	      P[ScatterParticleList_1[n]].ScatterID = NgbData[n].E[p].NgbIDs;
	      scatterlist_local[n]=NgbData[n].E[p].NgbIDs;
	      scatterlist_local[n+NumScatterParticles]=P[ScatterParticleList_1[n]].ID;
	      done_flag = 1;
	      break;
	    }
	}
      if(done_flag == 0)
	{
	  printf("no scatter partner found id=%d offset=%d numbgb=%g\n", P[ScatterParticleList_1[n]].ID,
		 NgbData[n].offset, P[ScatterParticleList_1[n]].sidm_NumNgb);
	  endrun(1012);
	}
    }
  
  
  count = mymalloc("count", sizeof(int) * NTask);
  offset = mymalloc("offset", sizeof(int) * NTask);
  
  MPI_Allgather(&TotNumScatterParticles, 1, MPI_INT, count, 1, MPI_INT, MPI_COMM_WORLD);
  
  for(i = 0, tot_count = 0, offset[0] = 0; i < NTask; i++)
    {
      tot_count += count[i];
      if(i > 0)
        offset[i] = offset[i - 1] + count[i - 1];
    }
  
  for(i = 0; i < NTask; i++)
    {
      count[i]*=sizeof(MyIDType);
      offset[i]*=sizeof(MyIDType);
    }
  
  scatterlist_global = mymalloc("scatterlist_global", tot_count * sizeof(MyIDType));
  MPI_Allgatherv(scatterlist_local, TotNumScatterParticles * sizeof(MyIDType), MPI_BYTE, scatterlist_global, count, offset, MPI_BYTE, MPI_COMM_WORLD);
  
  
  duplicate = mymalloc("duplicate", tot_count * sizeof(MyIDType));
  qsort(scatterlist_global, tot_count, sizeof(MyIDType), ID_cmp);
  
  numdup=0;
  for (i=1; i < tot_count; i++)
    {
      if (scatterlist_global[i-1]==scatterlist_global[i])
	{
	  duplicate[numdup++]=scatterlist_global[i];
	}
    }

  for (n=0; n<NumScatterParticles; n++)   
    {
      for (i=0; i<numdup; i++) 
	{
	  if ((P[ScatterParticleList_1[n]].ScatterID==duplicate[i]) || (P[ScatterParticleList_1[n]].ID==duplicate[i]))
	    P[ScatterParticleList_1[n]].ShouldScatterInStep=-1;
	}
    }
  
  myfree(duplicate);
  myfree(scatterlist_global);
  myfree(offset);
  myfree(count);
  myfree(scatterlist_local);
}

/* calculate smoothing length */
void sidm_Hsml(void)
{
  long long ntot;
  int i, j, ndone, ndone_flag, npleft, dummy, iter = 0;
  MyFloat *Left, *Right;
  char *Todo;
  int ngrp, recvTask, place, nexport, nimport;
  double t0, t1;


  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc("Ngblist", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct densdata_in) + sizeof(struct hsmldata_out) +
					     sizemax(sizeof(struct densdata_in),
						     sizeof(struct hsmldata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  Left = mymalloc("Left", sizeof(MyFloat) * NumPart);
  Right = mymalloc("Right", sizeof(MyFloat) * NumPart);
  Todo = mymalloc("Todo", sizeof(char) * NumPart);

  DM_Vx = mymalloc("DM_Vx", sizeof(MyFloat) * NumPart);
  DM_Vy = mymalloc("DM_Vy", sizeof(MyFloat) * NumPart);
  DM_Vz = mymalloc("DM_Vz", sizeof(MyFloat) * NumPart);

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      Left[i] = Right[i] = 0;
      Todo[i] = 1;
    }

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {
      t0 = second();

      i = FirstActiveParticle;


      do
	{
	  for(j = 0; j < NTask; j++)
	    {
	      Send_count[j] = 0;
	      Exportflag[j] = -1;
	    }

	  /* do local particles and prepare export list */

	  for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	    {
	      if(Todo[i])
	        if((1 << P[i].Type) & (SIDM))
		  {
		    if(sidm_Hsml_evaluate(i, 0, &nexport, Send_count) < 0)
		      break;
		  }
	    }

#ifdef OMP_SORT
	  omp_qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
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

	  DensDataGet =
	    (struct densdata_in *) mymalloc("	  DensDataGet", nimport * sizeof(struct densdata_in));
	  DensDataIn =
	    (struct densdata_in *) mymalloc("	  DensDataIn", nexport * sizeof(struct densdata_in));

	  /* prepare particle data for export */
	  for(j = 0; j < nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

              DensDataIn[j].Pos[0] = P[place].Pos[0];
              DensDataIn[j].Pos[1] = P[place].Pos[1];
              DensDataIn[j].Pos[2] = P[place].Pos[2];
              DensDataIn[j].Vel[0] = P[place].Vel[0];
              DensDataIn[j].Vel[1] = P[place].Vel[1];
              DensDataIn[j].Vel[2] = P[place].Vel[2];
              DensDataIn[j].Hsml = P[place].sidm_Hsml;

	      memcpy(DensDataIn[j].NodeList,
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
		      MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &DensDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

	  myfree(DensDataIn);
	  HsmlDataResult =
	    (struct hsmldata_out *) mymalloc("	  HsmlDataResult", nimport * sizeof(struct hsmldata_out));
	  HsmlDataOut =
	    (struct hsmldata_out *) mymalloc("	  HsmlDataOut", nexport * sizeof(struct hsmldata_out));


	  /* now do the particles that were sent to us */
	  for(j = 0; j < nimport; j++)
	    sidm_Hsml_evaluate(j, 1, &dummy, &dummy);

	  if(i < 0)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

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
		      MPI_Sendrecv(&HsmlDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct hsmldata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &HsmlDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct hsmldata_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }

      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
        {
          place = DataIndexTable[j].Index;

          P[place].sidm_Density += HsmlDataOut[j].Density;
          P[place].sidm_VelDisp += HsmlDataOut[j].VelDisp;
          P[place].sidm_PSum += HsmlDataOut[j].PSum;
          P[place].sidm_NumNgb += HsmlDataOut[j].Ngb;
          DM_Vx[place] += HsmlDataOut[j].Vx;
          DM_Vy[place] += HsmlDataOut[j].Vy;
          DM_Vz[place] += HsmlDataOut[j].Vz;
        }


	  myfree(HsmlDataOut);
	  myfree(HsmlDataResult);
	  myfree(DensDataGet);
	}
      while(ndone < NTask);


      /* do final operations on results */
      for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
	{
	  /* now check whether we had enough neighbours */
          if(!((1 << P[i].Type) & (SIDM)))
	    continue;

	  if(Todo[i])
	    {
              if((P[i].sidm_NumNgb < (All.DesNumNgb-All.MaxNumNgbDeviation) || P[i].sidm_NumNgb > (All.DesNumNgb+All.MaxNumNgbDeviation)) &&
                 ((Right[i] - Left[i]) > 1.0e-4 * Left[i] || Left[i] == 0 || Right[i] == 0))
		{
		  /* need to redo this particle */
		  npleft++;

		  if(P[i].sidm_NumNgb < All.DesNumNgb)
		    Left[i] = DMAX(P[i].sidm_Hsml, Left[i]);
		  else
		    {
		      if(Right[i] != 0)
			{
			  if(P[i].sidm_Hsml < Right[i])
			    Right[i] = P[i].sidm_Hsml;
			}
		      else
			Right[i] = P[i].sidm_Hsml;
		    }

		  if(iter >= MAXITER - 10)
		    {
		      printf
			("i=%d task=%d ID=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (int) P[i].ID, P[i].sidm_Hsml, Left[i], Right[i],
			 (double) P[i].sidm_NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1],
			 P[i].Pos[2]);
		      fflush(stdout);
		    }

		  if(Right[i] > 0 && Left[i] > 0)
		    P[i].sidm_Hsml = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
		  else
		    {
		      if(Right[i] == 0 && Left[i] == 0)
			endrun(8187);	/* can't occur */

		      if(Right[i] == 0 && Left[i] > 0)
			P[i].sidm_Hsml *= 1.26;

		      if(Right[i] > 0 && Left[i] == 0)
			P[i].sidm_Hsml /= 1.26;
		    }
		}
	      else
		Todo[i] = 0;
	    }
	}

      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
	{
	  iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("ngb iteration %d: need to repeat for %d%09d particles. (took %g sec)\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000), timediff(t0, t1));
	      fflush(stdout);
	    }

	  if(iter > MAXITER)
	    {
	      printf("failed to converge in neighbour iteration in density()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}
    }
  while(ntot > 0);


  /* final operations on velocity dispersion */
  for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
    {
      if(!((1 << P[i].Type) & (SIDM)))
        continue;

      DM_Vx[i] /= P[i].sidm_NumNgb;
      DM_Vy[i] /= P[i].sidm_NumNgb;
      DM_Vz[i] /= P[i].sidm_NumNgb;
      P[i].sidm_VelDisp /= P[i].sidm_NumNgb;

      P[i].sidm_VelDisp = sqrt(P[i].sidm_VelDisp - DM_Vx[i] * DM_Vx[i] - DM_Vy[i] * DM_Vy[i] - DM_Vz[i] * DM_Vz[i]);
    }


  myfree(DM_Vz);
  myfree(DM_Vy);
  myfree(DM_Vx);

  myfree(Todo);
  myfree(Right);
  myfree(Left);

  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);
}


int sidm_Hsml_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;
  int startnode, numngb, ngb, listindex = 0;
  double wk, rho, h, h2, hinv, hinv3;
  double r, r2, u, mass_j, v2, vx, vy, vz;
  MyDouble *pos, *vel;
  MyDouble phys_rho, phys_rel_vel, phys_rel_velx, phys_rel_vely, phys_rel_velz, PSum;
  MyDouble arho, avel, hubble_a;
  MyDouble dx, dy, dz;


  /* to physical */
  if(All.ComovingIntegrationOn) 
    {
      hubble_a = hubble_function(All.Time);
      arho = 1. / (All.Time * All.Time * All.Time);
      avel = 1. / All.Time;
    }
  else
    {
      hubble_a=0.0;
      arho = 1;
      avel = 1;
    }

  arho *= All.HubbleParam * All.HubbleParam;

  rho = 0;
  numngb = 0;
  v2 = vx = vy = vz = 0;
  PSum = 0.0;

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = P[target].Vel;
      h = P[target].sidm_Hsml;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
    }

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DensDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  ngb = sidm_ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nsend_local);

	  if(ngb < 0)
	    return -1;

	  for(n = 0; n < ngb; n++)
	    {
	      j = Ngblist[n];

	      h2 = h * h;

	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC
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


              if((r2 < h2) && (r2 > 0))
                {
                  hinv = 1. / h;
                  hinv3 = hinv * hinv * hinv;

                  r = sqrt(r2);

                  u = r * hinv;

                  if(u < 0.5)
                    {
                      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
                    }
                  else
                    {
                      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
                    }

                  mass_j = P[j].Mass;

                  rho += (mass_j * wk);
 
                  phys_rel_velx = avel * (P[j].Vel[0] - vel[0]) + hubble_a*All.Time*dx;
                  phys_rel_vely = avel * (P[j].Vel[1] - vel[1]) + hubble_a*All.Time*dy;
                  phys_rel_velz = avel * (P[j].Vel[2] - vel[2]) + hubble_a*All.Time*dz;

                  phys_rel_vel = sqrt(phys_rel_velx*phys_rel_velx + phys_rel_vely*phys_rel_vely + phys_rel_velz*phys_rel_velz);

                  phys_rho = mass_j * wk;

                  phys_rho *= arho;

                  //add to scatter sum; divide by two because of double counting
                  PSum += phys_rho * sidm_cross_sigma(phys_rel_vel) * phys_rel_vel / 2.0;

                  vx += P[j].Vel[0];
                  vy += P[j].Vel[1];
                  vz += P[j].Vel[2];

                  v2 += P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2];

                  numngb++;
                }
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = DensDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    {
      P[target].sidm_Density = rho;
      P[target].sidm_VelDisp = v2;
      P[target].sidm_PSum = PSum;
      DM_Vx[target] = vx;
      DM_Vy[target] = vy;
      DM_Vz[target] = vz;
      P[target].sidm_NumNgb = numngb;
    }
  else
    {
      HsmlDataResult[target].Density = rho;
      HsmlDataResult[target].VelDisp = v2;
      HsmlDataResult[target].PSum = PSum;
      HsmlDataResult[target].Vx = vx;
      HsmlDataResult[target].Vy = vy;
      HsmlDataResult[target].Vz = vz;
      HsmlDataResult[target].Ngb = numngb;
    }

  return 0;
}


/* construct neighbour list */
void sidm_NgbList(void)
{
  int i, j, dummy;
  int ndone_flag, ndone;
  int ngrp, sendTask, recvTask, place, nexport, nimport;


  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc("NgbList", NumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct densdata_in) + sizeof(struct ngbdata_out) +
                                             sizemax(sizeof(struct densdata_in),
                                                     sizeof(struct ngbdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));


  i = FirstActiveParticle;
  do
    {

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
	{
	  if(P[i].ShouldScatterInStep > 0)
	    if(sidm_NgbList_evaluate(i, 0, &nexport, Send_count) < 0)
	      break;
	}


#ifdef OMP_SORT
      omp_qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
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

      DensDataGet = (struct densdata_in *) mymalloc("DensDataGet", nimport * sizeof(struct densdata_in));
      DensDataIn = (struct densdata_in *) mymalloc("DensDataIn", nexport * sizeof(struct densdata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  DensDataIn[j].Pos[0] = P[place].Pos[0];
	  DensDataIn[j].Pos[1] = P[place].Pos[1];
	  DensDataIn[j].Pos[2] = P[place].Pos[2];
	  DensDataIn[j].Vel[0] = P[place].Vel[0];
	  DensDataIn[j].Vel[1] = P[place].Vel[1];
	  DensDataIn[j].Vel[2] = P[place].Vel[2];
	  DensDataIn[j].Hsml = P[place].sidm_Hsml;

	  memcpy(DensDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &DensDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      myfree(DensDataIn);

      NgbDataResult =
	(struct ngbdata_out *) mymalloc("	  NgbDataResult", nimport * sizeof(struct ngbdata_out));
      NgbDataOut =
	(struct ngbdata_out *) mymalloc("	  NgbDataOut", nexport * sizeof(struct ngbdata_out));


      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
	sidm_NgbList_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

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
		  MPI_Sendrecv(&NgbDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct ngbdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &NgbDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct ngbdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = ScatterParticleList_2[DataIndexTable[j].Index];
	  memcpy(&NgbData[place].E[NgbData[place].offset], &NgbDataOut[j].E[0],
		 NgbDataOut[j].Ngb * sizeof(ngb_entry));
	  NgbData[place].offset += NgbDataOut[j].Ngb;

	}


      myfree(NgbDataOut);
      myfree(NgbDataResult);
      myfree(DensDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);

}


int sidm_NgbList_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;
  int startnode, numngb, ngb, listindex = 0;
  double h, h2;
  double r, r2;
  MyDouble *pos, *vel;
  MyDouble dx, dy, dz;
  ngb_entry e[MAX_NGBS];
  double wk, hinv, hinv3, u;
  MyDouble phys_rho, phys_rel_vel, phys_rel_velx, phys_rel_vely, phys_rel_velz;
  MyDouble arho, avel, hubble_a;

  /* to physical */
  if(All.ComovingIntegrationOn) 
    {
      hubble_a = hubble_function(All.Time);
      arho = 1. / (All.Time * All.Time * All.Time);
      avel = 1. / All.Time;
    }
  else
    {
      hubble_a=0.0;
      arho = 1;
      avel = 1;
    }

  arho *= All.HubbleParam * All.HubbleParam;


  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = P[target].Vel;
      h = P[target].sidm_Hsml;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
    }

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DensDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  ngb = sidm_ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nsend_local);

	  if(ngb < 0)
	    return -1;

	  for(n = 0; n < ngb; n++)
	    {
	      j = Ngblist[n];
	      h2 = h * h;

	      dx = pos[0] - P[j].Pos[0];
	      dy = pos[1] - P[j].Pos[1];
	      dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC
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
	      if((r2 < h2) && (r2 > 0))
		{
		  hinv = 1. / h;
		  hinv3 = hinv * hinv * hinv;

		  r = sqrt(r2);

		  u = r * hinv;

		  if(u < 0.5)
		    {
		      wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
		    }
		  else
		    {
		      wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);
		    }

                  phys_rel_velx = avel * (P[j].Vel[0] - vel[0]) + hubble_a*All.Time*dx;
                  phys_rel_vely = avel * (P[j].Vel[1] - vel[1]) + hubble_a*All.Time*dy;
                  phys_rel_velz = avel * (P[j].Vel[2] - vel[2]) + hubble_a*All.Time*dz;

                  phys_rel_vel = sqrt(phys_rel_velx*phys_rel_velx + phys_rel_vely*phys_rel_vely + phys_rel_velz*phys_rel_velz);

		  phys_rho = P[j].Mass * wk;

		  phys_rho *= arho;

		  e[numngb].P0j_half = phys_rho * sidm_cross_sigma(phys_rel_vel) * phys_rel_vel / 2.0;
		  e[numngb].NgbIDs = P[j].ID;
		  e[numngb].Distance = r;

		  numngb++;

		}
	    }

	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = DensDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }
  if(mode == 0)
    {
      memcpy(&NgbData[ScatterParticleList_2[target]].E[NgbData[ScatterParticleList_2[target]].offset], &e[0],
	     numngb * sizeof(ngb_entry));
      NgbData[ScatterParticleList_2[target]].offset = numngb;
    }
  else
    {
      memcpy(&NgbDataResult[target].E[0], &e[0], numngb * sizeof(ngb_entry));
      NgbDataResult[target].Ngb = numngb;
    }

  return 0;
}


/* do scatter */
void sidm_Scatter(void)
{
  int i, j, dummy;
  int ndone_flag, ndone;
  int ngrp, sendTask, recvTask, place, nexport, nimport;
#ifdef ENERGY_STATS
  double kin_energy_before, kin_energy_after;
  kin_energy_before = get_kinetic_energy();
#endif

  /* allocate buffers to arrange communication */
  Ngblist = (int *) mymalloc("NgbList", NumPart * sizeof(int));


  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
                                             sizeof(struct densdata_in) + sizeof(struct scatterdata_out) +
                                             sizemax(sizeof(struct densdata_in),
                                                     sizeof(struct scatterdata_out))));

  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  i = FirstActiveParticle;	/* first particle for this task */

  do
    {
      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */

      for(nexport = 0; i >= 0; i = NextActiveParticle[i])
        if((1 << P[i].Type) & (SIDM))
	  {
	    if(P[i].ShouldScatterInStep > 0)
	      {
		if(sidm_Scatter_evaluate(i, 0, &nexport, Send_count) < 0)
		  break;
	      }
	  }


#ifdef OMP_SORT
      omp_qsort(DataIndexTable, nexport, sizeof(struct data_index), data_index_compare);
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

      DensDataGet = (struct densdata_in *) mymalloc("DensDataGet", nimport * sizeof(struct densdata_in));
      DensDataIn = (struct densdata_in *) mymalloc("DensDataIn", nexport * sizeof(struct densdata_in));

      /* prepare particle data for export */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;

	  DensDataIn[j].Pos[0] = P[place].Pos[0];
	  DensDataIn[j].Pos[1] = P[place].Pos[1];
	  DensDataIn[j].Pos[2] = P[place].Pos[2];
	  DensDataIn[j].Vel[0] = P[place].Vel[0];
	  DensDataIn[j].Vel[1] = P[place].Vel[1];
	  DensDataIn[j].Vel[2] = P[place].Vel[2];
	  DensDataIn[j].Hsml = P[place].sidm_Hsml;
	  DensDataIn[j].ScatterID = P[place].ScatterID;

	  memcpy(DensDataIn[j].NodeList,
		 DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	}

      /* exchange particle data */
      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ ngrp;

	  if(recvTask < NTask)
	    {
	      if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		{
		  /* get the particles */
		  MPI_Sendrecv(&DensDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A,
			       &DensDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct densdata_in), MPI_BYTE,
			       recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      myfree(DensDataIn);

      ScatterDataResult =
	(struct scatterdata_out *) mymalloc("	  ScatterDataResult", nimport * sizeof(struct scatterdata_out));
      ScatterDataOut =
	(struct scatterdata_out *) mymalloc("	  ScatterDataOut", nexport * sizeof(struct scatterdata_out));


      /* now do the particles that were sent to us */
      for(j = 0; j < nimport; j++)
	sidm_Scatter_evaluate(j, 1, &dummy, &dummy);

      if(i < 0)
	ndone_flag = 1;
      else
	ndone_flag = 0;

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
		  MPI_Sendrecv(&ScatterDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct scatterdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B,
			       &ScatterDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct scatterdata_out),
			       MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}

      /* add the result to the local particles */
      for(j = 0; j < nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  P[place].DidScatterInStep += ScatterDataOut[j].DidScatterInStep;
	  if(ScatterDataOut[j].DidScatterInStep > 0)
	    {
	      /* assign new velocity */
	      P[place].Vel[0] = ScatterDataOut[j].Vx;
	      P[place].Vel[1] = ScatterDataOut[j].Vy;
	      P[place].Vel[2] = ScatterDataOut[j].Vz;
	    }
	}
      myfree(ScatterDataOut);
      myfree(ScatterDataResult);
      myfree(DensDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Ngblist);


#ifdef ENERGY_STATS
  kin_energy_after = get_kinetic_energy();

  if(ThisTask == 0)
      printf("kinetic energy: before=%g after=%g before-after=%g\n", kin_energy_before, kin_energy_after,
	     kin_energy_before - kin_energy_after);
#endif

#ifdef SCATTER_STATS
  scatter_stats();
#endif
}


int sidm_Scatter_evaluate(int target, int mode, int *nexport, int *nsend_local)
{
  int j, n;
  int startnode, numngb, ngb, listindex = 0;
  double h;
  MyDouble *pos, *vel;
  MyDouble theta, phi;
  MyDouble xe = 0, ye = 0, ze = 0;
  int didScatterInStep = 0;
  double vcm_x = 0, vcm_y = 0, vcm_z = 0;
  double dvx, dvy, dvz, dv = 0;
  MyIDType scatterID;
  MyDouble arho, avel, hubble_a;
  MyDouble dx, dy, dz;
  MyDouble new_velx_internal, new_vely_internal, new_velz_internal;
#ifdef SIDM_INELASTIC
  float randnr;
#endif 

  /* to physical */
  if(All.ComovingIntegrationOn) 
    {
      hubble_a = hubble_function(All.Time);
      arho = 1 / (All.Time * All.Time * All.Time);
      avel = 1 / All.Time;
    }
  else
    {
      hubble_a=0.0;
      arho = 1;
      avel = 1;
    }  

  if(mode == 0)
    {
      pos = P[target].Pos;
      vel = P[target].Vel;
      h = P[target].sidm_Hsml;
      scatterID = P[target].ScatterID;
    }
  else
    {
      pos = DensDataGet[target].Pos;
      vel = DensDataGet[target].Vel;
      h = DensDataGet[target].Hsml;
      scatterID = DensDataGet[target].ScatterID;
    }

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = DensDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  numngb = 0;

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  ngb = sidm_ngb_treefind_variable(pos, h, target, &startnode, mode, nexport, nsend_local);

	  if(ngb < 0)
	    return -1;

	  for(n = 0; n < ngb; n++)
	    {
	      j = Ngblist[n];

	      if(P[j].ID == scatterID)
		{
		  didScatterInStep++;
		  P[j].DidScatterInStep++;

		  dx = pos[0] - P[j].Pos[0];
		  dy = pos[1] - P[j].Pos[1];
		  dz = pos[2] - P[j].Pos[2];

#ifdef PERIODIC
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
		  /* physical center-of-mass velocity */
		  vcm_x = ( avel*vel[0] + hubble_a*All.Time*(pos[0])    +    avel*P[j].Vel[0] + hubble_a*All.Time*(P[j].Pos[0]) ) / 2.0;
		  vcm_y = ( avel*vel[1] + hubble_a*All.Time*(pos[1])    +    avel*P[j].Vel[1] + hubble_a*All.Time*(P[j].Pos[1]) ) / 2.0;
		  vcm_z = ( avel*vel[2] + hubble_a*All.Time*(pos[2])    +    avel*P[j].Vel[2] + hubble_a*All.Time*(P[j].Pos[2]) ) / 2.0;

		  /* physical relative velocity */
		  dvx = avel*(vel[0] - P[j].Vel[0]) + hubble_a*All.Time*dx;
		  dvy = avel*(vel[1] - P[j].Vel[1]) + hubble_a*All.Time*dy;
		  dvz = avel*(vel[2] - P[j].Vel[2]) + hubble_a*All.Time*dz;

		  dv = sqrt(dvx * dvx + dvy * dvy + dvz * dvz);


#ifdef SIDM_INELASTIC
#ifdef SIDM_INELASTIC_MODEL_1
		  dv+=+SplittingVelocity; //exo-thermic only
#endif 
#ifdef SIDM_INELASTIC_MODEL_2
                  dv+=-SplittingVelocity; //endo-thermic only
#endif
#ifdef SIDM_INELASTIC_MODEL_3
		  randnr=get_random_number(P[j].ID)-0.5;
                  dv+=randnr/fabs(randnr)*SplittingVelocity; //exo/endo-thermic randomly
#endif
#endif
                  dv /= 2.0;

		  /* unit sphere vector for isotropic scattering */
		  phi = 2 * M_PI * get_random_number(P[j].ID+1);
		  theta = acos(get_random_number(P[j].ID+2) * 2 - 1);

		  xe = sin(theta) * cos(phi);
		  ye = sin(theta) * sin(phi);
		  ze = cos(theta);

#ifdef SIDM_DEBUG
		  printf("EVT: %g   %g %g %g   %g %g %g   %g %g %g   %g %g %g   %g %g %g   %g %g %g\n",
			 All.Time,
			 vcm_x + dv * xe, vcm_y + dv * ye, vcm_z + dv * ze,
			 vcm_x - dv * xe, vcm_y - dv * ye, vcm_z - dv * ze,
			 P[j].Vel[0], P[j].Vel[1], P[j].Vel[2],
			 vel[0], vel[1], vel[2],
			 P[j].Pos[0], P[j].Pos[1], P[j].Pos[2], pos[0], pos[1], pos[2]);
		  fflush(stdout);
#endif
		  /* physical velocity to internal velocity */
		  new_velx_internal = 1./avel*( (vcm_x + dv * xe) - hubble_a*All.Time*(P[j].Pos[0]) );
		  new_vely_internal = 1./avel*( (vcm_y + dv * ye) - hubble_a*All.Time*(P[j].Pos[1]) );
		  new_velz_internal = 1./avel*( (vcm_z + dv * ze) - hubble_a*All.Time*(P[j].Pos[2]) );
		  /* assign new velocity */
		  P[j].Vel[0] = new_velx_internal;
		  P[j].Vel[1] = new_vely_internal;
		  P[j].Vel[2] = new_velz_internal;
		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = DensDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }
  if(mode == 0)
    {
      P[target].DidScatterInStep = didScatterInStep;
      if(didScatterInStep > 0)
	{
          /* physical velocity to internal velocity */
	  new_velx_internal = 1./avel*( (vcm_x - dv * xe) - hubble_a*All.Time*(pos[0]) );
	  new_vely_internal = 1./avel*( (vcm_y - dv * ye) - hubble_a*All.Time*(pos[1]) );
          new_velz_internal = 1./avel*( (vcm_z - dv * ze) - hubble_a*All.Time*(pos[2]) ); 
	  /* assign new velocity */
	  P[target].Vel[0] = new_velx_internal;
	  P[target].Vel[1] = new_vely_internal;
	  P[target].Vel[2] = new_velz_internal;
	}
    }
  else
    {
      /* physical velocity to internal velocity */
      new_velx_internal = 1./avel*( (vcm_x - dv * xe) - hubble_a*All.Time*(pos[0]) );
      new_vely_internal = 1./avel*( (vcm_y - dv * ye) - hubble_a*All.Time*(pos[1]) );
      new_velz_internal = 1./avel*( (vcm_z - dv * ze) - hubble_a*All.Time*(pos[2]) );
      /* assign new velocity */
      ScatterDataResult[target].DidScatterInStep = didScatterInStep;
      ScatterDataResult[target].Vx = new_velx_internal;
      ScatterDataResult[target].Vy = new_vely_internal;
      ScatterDataResult[target].Vz = new_velz_internal;

    }
  return 0;
}


/* treefind for dark matter particles, type in SIDM encoded */
int sidm_ngb_treefind_variable(MyDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			       int *nexport, int *nsend_local)
{
  int numngb, no, p, task, nexport_save;
  struct NODE *current;
  // cache some global vars locally for improved compiler alias analysis
  int maxPart = All.MaxPart;
  int maxNodes = MaxNodes;
  int bunchSize = All.BunchSize;
  integertime ti_Current = All.Ti_Current;
  MyDouble dist, dx, dy, dz;
#ifdef PERIODIC
  MyDouble xtmp;
#endif

  nexport_save = *nexport;

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < maxPart)		/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

          if(!((1 << P[p].Type) & (SIDM)))
	    continue;

	  if(P[p].Ti_current != ti_Current)
	    drift_particle(p, ti_Current);

	  dist = hsml;

	  dx = NGB_PERIODIC_LONG_X(P[p].Pos[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(P[p].Pos[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(P[p].Pos[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  Ngblist[numngb++] = p;
	}
      else
	{
	  if(no >= maxPart + maxNodes)	/* pseudo particle */
	    {
	      if(mode == 1)
		endrun(12312);

	      if(target >= 0)	/* if no target is given, export will not occur */
		{
		  if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
		    {
		      Exportflag[task] = target;
		      Exportnodecount[task] = NODELISTLENGTH;
		    }

		  if(Exportnodecount[task] == NODELISTLENGTH)
		    {
		      if(*nexport >= bunchSize)
			{
			  *nexport = nexport_save;
			  if(nexport_save == 0)
			    endrun(13004);	/* in this case, the buffer is too small to process even a single particle */
			  for(task = 0; task < NTask; task++)
			    nsend_local[task] = 0;
			  for(no = 0; no < nexport_save; no++)
			    nsend_local[DataIndexTable[no].Task]++;
			  return -1;
			}
		      Exportnodecount[task] = 0;
		      Exportindex[task] = *nexport;
		      DataIndexTable[*nexport].Task = task;
		      DataIndexTable[*nexport].Index = target;
		      DataIndexTable[*nexport].IndexGet = *nexport;
		      *nexport = *nexport + 1;
		      nsend_local[task]++;
		    }

		  DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] =
		    DomainNodeIndex[no - (maxPart + maxNodes)];

		  if(Exportnodecount[task] < NODELISTLENGTH)
		    DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;
		}

	      no = Nextnode[no - maxNodes];
	      continue;
	    }

	  current = &Nodes[no];

	  if(mode == 1)
	    {
	      if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
		{
		  *startnode = -1;
		  return numngb;
		}
	    }

	  if(current->Ti_current != ti_Current)
	    force_drift_node(no, ti_Current);

	  if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
	    {
	      if(current->u.d.mass)	/* open cell */
		{
		  no = current->u.d.nextnode;
		  continue;
		}
	    }

	  no = current->u.d.sibling;	/* in case the node can be discarded */

	  dist = hsml + 0.5 * current->len;
	  dx = NGB_PERIODIC_LONG_X(current->center[0] - searchcenter[0]);
	  if(dx > dist)
	    continue;
	  dy = NGB_PERIODIC_LONG_Y(current->center[1] - searchcenter[1]);
	  if(dy > dist)
	    continue;
	  dz = NGB_PERIODIC_LONG_Z(current->center[2] - searchcenter[2]);
	  if(dz > dist)
	    continue;
	  /* now test against the minimal sphere enclosing everything */
	  dist += FACT1 * current->len;
	  if(dx * dx + dy * dy + dz * dz > dist * dist)
	    continue;

	  no = current->u.d.nextnode;	/* ok, we need to open the node */
	}
    }
  *startnode = -1;
  return numngb;
}

/* check which particle to scatter */
void sidm_check_particle_scatter(void)
{
  int i;
  MyDouble dtime, xran, PSum, dt, time_hubble_a, hubble_a;

  NumScatterParticles = 0;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {

      dt = (P[i].TimeBin ? (((integertime) 1) << P[i].TimeBin) : 0) * All.Timebase_interval;

      if(All.ComovingIntegrationOn) 
	{
	  hubble_a = hubble_function(All.Time);
	  time_hubble_a = All.Time * hubble_a;
	  dtime = All.Time * dt / time_hubble_a;
	}
      else
        {
	 dtime = dt;
        }

      dtime *= 1.0 / All.HubbleParam;

      xran = get_random_number(P[i].ID);
      PSum = P[i].sidm_PSum * dtime;

      if(xran < PSum)
	{
	  P[i].sidm_NumTotalScatter++;
	  P[i].ShouldScatterInStep = 1;
          P[i].RandX = xran;
	  ScatterParticleList_1[NumScatterParticles] = i;
	  ScatterParticleList_2[i] = NumScatterParticles;
	  NumScatterParticles++;
	}
    }
}


/* first init from init.c */
void sidm_Init_Particles(void)
{
  int i;
  
  mpi_printf("SIDM: Init Particles...\n");
  if (All.DesNumNgb+All.MaxNumNgbDeviation>=MAX_NGBS)
   { 
    terminate("All.DesNumNgb+All.MaxNumNgbDeviation>=MAX_NGBS\n");
   }

  for(i = 0; i < NumPart; i++)
    {
      P[i].ShouldScatterInStep = 0;
      P[i].DidScatterInStep = 0;
      P[i].ScatterID = 0;
      P[i].RandX = 0;
      P[i].sidm_PSum = 0.0;
      P[i].sidm_NumTotalScatter = 0;
      P[i].sidm_Hsml = All.SofteningTable[P[i].Type];
      P[i].sidm_Density = 0.0;
      P[i].sidm_VelDisp = 0.0;
      P[i].sidm_NumNgb = 0;
    }
  mpi_printf("done.\n");
}

void sidm_Init_CrossSection(void)
{
  mpi_printf("SIDM: Init CrossSection...\n");
  CrossUnitFac=All.UnitMass_in_g / (All.UnitLength_in_cm * All.UnitLength_in_cm);
#ifdef SIDM_INELASTIC
  SpeedOfLight=SPEED_OF_LIGHT / All.UnitVelocity_in_cm_per_s;
  SplittingVelocity=sqrt(All.Splitting)*SpeedOfLight;
#endif
#ifndef CONST_CROSS
#ifdef SIDM_INELASTIC
  init_cross_table_inelastic();
#else
  init_cross_table_elastic();
#endif
#endif
 mpi_printf("done.\n");
}

/* reinit before each scatter loop */
void sidm_ReInit(void)
{
  int i;

  for(i = 0; i < NumPart; i++)
    {
      P[i].ShouldScatterInStep = 0;
      P[i].DidScatterInStep = 0;
      P[i].ScatterID = 0;
    }
}

/* kinetic energy of all particles */
double get_kinetic_energy(void)
{
  double global_kin_energy = 0, local_kin_energy = 0;
  int i;

  for(i = 0; i < NumPart; i++)
    local_kin_energy +=
      0.5 * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

  MPI_Allreduce(&local_kin_energy, &global_kin_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return global_kin_energy;
}

/* some scatter statistics */
void scatter_stats(void)
{
  int global_ShouldScatterInStep = 0, local_ShouldScatterInStep = 0;
  int global_DidScatterInStep = 0, local_DidScatterInStep = 0;
  int global_rejected = 0, local_rejected = 0;

  int i;

  for(i = 0; i < NumPart; i++)
    {
      if (P[i].ShouldScatterInStep>0)
	local_ShouldScatterInStep += 1;

      if (P[i].ShouldScatterInStep<0)
	local_rejected += 1;

      local_DidScatterInStep += (int) P[i].DidScatterInStep;
      
      if((P[i].ShouldScatterInStep > 0) && (P[i].DidScatterInStep == 0))
	{
	  printf("failed scatter: %d %d %d %d\n", P[i].DidScatterInStep, P[i].ID, P[i].ScatterID,
		 P[i].ShouldScatterInStep);
	  endrun(1012);
	}
      
      if(P[i].DidScatterInStep > 1)
	{
	  printf("multiple scatter: %d %d %d %d\n", P[i].DidScatterInStep, P[i].ID, P[i].ScatterID,
		 P[i].ShouldScatterInStep);
	  endrun(1012);
	}
    }
  
  MPI_Allreduce(&local_ShouldScatterInStep, &global_ShouldScatterInStep, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_DidScatterInStep, &global_DidScatterInStep, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_rejected, &global_rejected, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
  if(ThisTask == 0)
    printf("active scatters: ShouldScatterInStep = %08d  DidScatterInStep = %08d  rejected=%08d  (%06.2f percent)\n",
	   global_ShouldScatterInStep, global_DidScatterInStep, global_rejected, 100.0*global_rejected/(2.0*global_ShouldScatterInStep+1e-5));
}

#ifdef CONST_CROSS
/* cross section in code units, constant */
MyDouble sidm_cross_sigma(MyDouble rel_vel)
{
  return CrossUnitFac * All.CrossSectionPerMass_in_cgs;
}

#else
/* cross section in code units, from table*/
MyDouble sidm_cross_sigma(MyDouble rel_vel)
{
  int bin=(int)(log(rel_vel/CROSS_VMIN)/Dvlog-log(CROSS_VMIN));

  if (bin<0)
   bin=0;
  if (bin>CROSS_VBINS-1)
   bin=CROSS_VBINS-1;
  return CrossTable[bin];
}

/* set up cross section interpolation table for inelastic case */
void init_cross_table_inelastic(void)
{
}

/* set up cross section interpolation table for elastic case */
void init_cross_table_elastic(void)
{
 int i;
 Dvlog=log(CROSS_VMAX/CROSS_VMIN)/CROSS_VBINS;
 double rel_vel, beta, sigma_norm;
 FILE *fp=fopen("sidm_cross.txt", "w");

 for (i=0; i < CROSS_VBINS; i++)
   {
    rel_vel=exp(Dvlog*(i+0.5) + log(CROSS_VMIN));
    beta=M_PI*(All.PeakSigma/rel_vel)*(All.PeakSigma/rel_vel);
    sigma_norm=0.0;
    if (beta<0.1)
     sigma_norm=(4.0*M_PI/22.7) * beta*beta * log(1.0+1.0/beta);
    if ((beta<1.0e3) && (beta>0.1))
     sigma_norm=(8.0*M_PI/22.7) * beta*beta / (1.0 + 1.5*pow(beta, 1.65));
    if (beta>1.0e3)
     sigma_norm=(M_PI/22.7) * pow(log(beta) + 1.0 - 0.5/log(beta), 2.0);

    CrossTable[i]=CrossUnitFac * All.CrossSectionPerMass_in_cgs*sigma_norm;
    fprintf(fp, "%d %g %g %g %g %g\n", i, rel_vel, beta, sigma_norm, CrossTable[i], CrossTable[i]/CrossUnitFac);
   }
 fclose(fp);
}

#endif

#endif

