#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef VORONOI
#include "voronoi.h"


void voronoi_hydro_force(void)
{
  int i, j, q1, q2, li, ri, timebin1, timebin2;
  double length, pressure1, pressure2, cx, cy, cz, c, fac1, fac2, ex, ey, ez, forcex, forcey, forcez, vdotr2;
  double mass1, mass2, dEdt, fviscx, fviscy, fviscz, w, dens1, dens2, volume1, volume2;

  double hubble_a, hubble_a2, fac_mu, atime;
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)
  MyFloat f,h,h1,h2;
  double DivVel1, DivVel2, CurlVel1, CurlVel2;
#endif
#ifdef VORONOI_SHAPESCHEME
  double w1, w2;
#endif
#ifdef VORONOI_TIME_DEP_ART_VISC
  double damp1, damp2, delta_damp1, delta_damp2;
  double min_damp, v_sig, v_perp, partofcell1, partofcell2, dim_fac;
  min_damp = All.ArtBulkViscConst / 16.0;
#endif
#ifdef KH_MIXING
  double v_chg[3];
  double mix_fac, entr_chg, tkh, tkh_inv, deltav2, al_fac, energy_change, entropy_corr;
  double  entr_chg1, entr_chg2;
#endif
  MyFloat *vel1, *vel2, *center1, *center2;
  point *p1, *p2;


  if (All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration */
      hubble_a = hubble_function(All.Time);
      hubble_a2 = All.Time * All.Time * hubble_a;

      fac_mu = pow(All.Time, 3 * (GAMMA - 1) / 2) / All.Time;
      atime = All.Time;
    }
  else
    hubble_a = hubble_a2 = atime = fac_mu = 1.0;


  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	{
	  SphP[i].e.DtEntropy = 0;

	  for(j = 0; j < 3; j++)
	    SphP[i].a.dHydroAccel[j] = 0;
	}
    }


  voronoi_exchange_ghost_variables();

  for(i = 0; i < Nvf; i++)
    {
      p1 = &DP[VF[i].p1];
      p2 = &DP[VF[i].p2];

      q1 = li = p1->index;
      q2 = ri = p2->index;

      if(!
	 ((li >= 0 && li < N_gas && p1->task == ThisTask) || (ri >= 0 && ri < N_gas && p2->task == ThisTask)))
	continue;

      if(li < 0 || ri < 0)
	continue;

      if(li >= N_gas && p1->task == ThisTask)
	li -= N_gas;

      if(ri >= N_gas && p2->task == ThisTask)
	ri -= N_gas;

      if(p1->task == ThisTask)
	{
	  pressure1 = SphP[li].Pressure;
	  dens1 = SphP[li].d.Density;
	  vel1 = SphP[li].VelPred;
	  mass1 = P[li].Mass;
	  center1 = SphP[li].Center;
	  volume1 = SphP[li].Volume;
	  timebin1 = P[li].TimeBin;
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)
	  DivVel1 = SphP[li].v.DivVel;
	  CurlVel1 = SphP[li].r.CurlVel;
#endif
#ifdef VORONOI_TIME_DEP_ART_VISC
	  damp1 = SphP[li].alpha;
#endif
#ifdef VORONOI_SHAPESCHEME
	  w1 = SphP[li].W;
#endif
	}
      else
	{
	  pressure1 = PrimExch[q1].Pressure;
	  dens1 = PrimExch[q1].Density;
	  vel1 = PrimExch[q1].VelPred;
	  mass1 = PrimExch[q1].Mass;
	  center1 = PrimExch[q1].Center;
	  volume1 = PrimExch[q1].Volume;
	  timebin1 = PrimExch[q1].TimeBin;
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)
	  DivVel1 = PrimExch[q1].DivVel;
	  CurlVel1 = PrimExch[q1].CurlVel;
#endif
#ifdef VORONOI_TIME_DEP_ART_VISC
	  damp1 = PrimExch[q1].Alpha;
#endif
#ifdef VORONOI_SHAPESCHEME
	  w1 = PrimExch[q1].W;
#endif
	}

      if(p2->task == ThisTask)
	{
	  pressure2 = SphP[ri].Pressure;
	  dens2 = SphP[ri].d.Density;
	  vel2 = SphP[ri].VelPred;
	  mass2 = P[ri].Mass;
	  center2 = SphP[ri].Center;
	  volume2 = SphP[ri].Volume;
	  timebin2 = P[ri].TimeBin;
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)
	  DivVel2 = SphP[ri].v.DivVel;
	  CurlVel2 = SphP[ri].r.CurlVel;
#endif
#ifdef VORONOI_TIME_DEP_ART_VISC
	  damp2 = SphP[ri].alpha;
#endif
#ifdef VORONOI_SHAPESCHEME
	  w2 = SphP[ri].W;
#endif
	}
      else
	{
	  pressure2 = PrimExch[q2].Pressure;
	  dens2 = PrimExch[q2].Density;
	  vel2 = PrimExch[q2].VelPred;
	  mass2 = PrimExch[q2].Mass;
	  center2 = PrimExch[q2].Center;
	  volume2 = PrimExch[q2].Volume;
	  timebin2 = PrimExch[q2].TimeBin;
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)
	  DivVel2 = PrimExch[q2].DivVel;
	  CurlVel2 = PrimExch[q2].CurlVel;
#endif
#ifdef VORONOI_TIME_DEP_ART_VISC
	  damp2 = PrimExch[q2].Alpha;
#endif
#ifdef VORONOI_SHAPESCHEME
	  w2 = PrimExch[q2].W;
#endif
	}

#ifdef TWODIMS
      h1 = sqrt(volume1 / M_PI);
      h2 = sqrt(volume2 / M_PI);
#else
      h1 = pow(volume1 / (4.0 / 3 * M_PI), 1.0 / 3);
      h2 = pow(volume2 / (4.0 / 3 * M_PI), 1.0 / 3);
#endif
      h = 0.5 * (h1 + h2);

      cx = p2->x - p1->x;
      cy = p2->y - p1->y;
      cz = p2->z - p1->z;

      c = sqrt(cx * cx + cy * cy + cz * cz);	/* distance of the two points */
      length = VF[i].area;	/* length/area of common face */


      fac1 = 0.5 * (pressure2 + pressure1) * length / c;
      fac2 = (pressure2 - pressure1) * length / c;


      ex = VF[i].cx - 0.5 * (p1->x + p2->x);
      ey = VF[i].cy - 0.5 * (p1->y + p2->y);
      ez = VF[i].cz - 0.5 * (p1->z + p2->z);

#ifdef VORONOI_SHAPESCHEME
			double Stffnss = All.VoronoiStiffNess;
			double Rndnss  = All.VoronoiRoundNess;
#endif

      double ee = sqrt(ex * ex  + ey * ey + ez * ez);

      /* put in a limiter for highly distorted cells */
      if(ee > 0.5 * c)
	fac2 *= 0.5 / (ee/c);

///////////////////// HERE VACUUM /////////////
#ifdef VORONOI_VACUUM
#ifdef TWODIMS
			double VACUUM_AREA = 4.0 * All.VoronoiVacuumLength;
			double VACUUM_VOL  = 0.2 * pow(0.1*All.VoronoiVacuumLength,2); //minimum 0.16*...//
#else
			double VACUUM_AREA = pow(4.0*All.VoronoiVacuumLength,2);
			double VACUUM_VOL  = 0.2 * pow(0.1*All.VoronoiVacuumLength,3);
#endif

			if (c > 2.0* All.VoronoiVacuumLength)
					continue;
// 			if (length >= VACUUM_AREA )
			if ((length >= VACUUM_AREA ) || (sqrt(ex * ex + ey * ey + ez * ez) > 0.25*All.VoronoiVacuumLength)|| (volume1 > VACUUM_VOL || volume2 > VACUUM_VOL))
// 					|| vol1 > VACUUM_VOL || vol2 > VACUUM_VOL)
				{
// 					fprintf(stderr, "\n really happens\n");
					ex = 0.0;
					ey = 0.0;
					ez = 0.0;
// 					length = VACUUM_AREA;
      fac1 = 0.5 * (pressure2 - pressure1) * length / c;
      fac2 = (pressure2 - pressure1) * length / c;
#ifdef VORONOI_SHAPESCHEME
					Stffnss = 0.0;
					Rndnss  = 0.0;
#endif
				}
#endif
    
      forcex = -fac1 * cx - fac2 * ex;
      forcey = -fac1 * cy - fac2 * ey;
      forcez = -fac1 * cz - fac2 * ez;


      /* calculate viscous force */

      vdotr2 = cx * (vel2[0] - vel1[0]) 
	     + cy * (vel2[1] - vel1[1]) 
	     + cz * (vel2[2] - vel1[2]);

      if(All.ComovingIntegrationOn)
	vdotr2 += hubble_a2 * (c * c);	/* corresponds to physical r * dr/dt  */

      if(vdotr2 < 0)
	{
	  w = vdotr2 / c;   /* velocity of approach (physical times a), projected onto separation vector */ 
 
	  double csound, pvisc, pvisc_max, dens_avg;

	  /* calculate viscous force */

	  csound = 0.5 * (sqrt(GAMMA * pressure1 / dens1) 
			+ sqrt(GAMMA * pressure2 / dens2));
	  
	  //	  dens_avg = 0.5 * (dens1 + dens2)
	  dens_avg = 2.0 * (dens1 * dens2 / (dens1 + dens2));
#ifndef VORONOI_TIME_DEP_ART_VISC
	  pvisc = All.ArtBulkViscConst * dens_avg * (-w * csound + 2.0 * w * w);
#else
	  pvisc =0.5 * (damp1 + damp2) * dens_avg * (-w * csound + 2.0 * w * w);

	  /* now update damp  */
	  if (All.ComovingIntegrationOn)
	    DivVel1 = DivVel1 / (atime * atime) + hubble_a;
	  if (All.ComovingIntegrationOn)
	    DivVel2 = DivVel2 / (atime * atime) + hubble_a;

	  /* source term */
	  if ((DivVel1 > 0.0) && (damp1 < 2.0 * All.ArtBulkViscConst))
	    delta_damp1 += dim_fac * partofcell1 * 0.50 * 16.0 * DivVel1;
	  if ((DivVel2 > 0.0) && (damp2 < 2.0 * All.ArtBulkViscConst))
	    delta_damp2 += dim_fac * partofcell2 * 0.50 * 16.0 * DivVel2;
	  /* dim_fac*32.0 (limit 2.0) works for sedov2D und sod2D3D slightly overdamped */

	  /* decay term */
	  delta_damp1 -=
	    0.1 * v_sig * (damp1 - min_damp) * 4.0 * 0.5 * partofcell1;
	  delta_damp2 -=
	    0.1 * v_sig * (damp2 - min_damp) * 4.0 * 0.5 * partofcell2;

	  /* spatial smoothing */
	  delta_damp1 +=
	    2.0 * 0.25 * (damp2 - damp1) * v_sig * 0.5 * partofcell1;
	  delta_damp2 -=
	    2.0 * 0.25 * (damp2 - damp1) * v_sig * 0.5 * partofcell2;

	  if (DivVel1 > 0.0)
	    delta_damp1 +=
	      (damp2 - damp1) * v_sig * v_perp * 0.5 * partofcell1;
	  if (DivVel2 > 0.0)
	    delta_damp2 -=
	      (damp2 - damp1) * v_sig * v_perp * 0.5 * partofcell2;
#endif
#ifdef VORONOI_BALSARA
	  f = 0.5 * fabs(DivVel1) /
	    (fabs(DivVel1) + fabs(CurlVel1) +
	     0.00000001 * csound / fac_mu / h);
	  f += 0.5 * fabs(DivVel2) /
	    (fabs(DivVel2) + fabs(CurlVel2) +
	     0.00000001 * csound / fac_mu / h);

	  pvisc *= f;
# else
	  f = 1.0;
#endif
  

	  /* now limit the viscosity such that it can not lead to a reversal of the relative velocity of the two particles */

	  int dtbin = IMAX(timebin1, timebin2);
	  double dt = (dtbin ? (1 << dtbin) : 0) * All.Timebase_interval;
	  
	  if(All.ComovingIntegrationOn)
	    dt /= (hubble_a * pow(All.Time, 3 * (GAMMA - 1)));

	  if(dt > 0)
	    pvisc_max = -w / (length * dt  * (1/mass1 + 1/mass2));
	  else
	    pvisc_max = 0;

	  if(pvisc > pvisc_max)
	    pvisc = pvisc_max;
 
	  fviscx = -pvisc * length * cx / c;
	  fviscy = -pvisc * length * cy / c;
	  fviscz = -pvisc * length * cz / c;

	  /* rate at which energy is dissipated */
	  dEdt = -w * pvisc * length;


#ifdef KH_MIXING
  double dvx, dvy, dvz;
      entr_chg = 0.0;
      mix_fac = 1.0;
      al_fac = 0.0;
      dvx = vel2[0] - vel1[0];
      dvy = vel2[1] - vel1[1];
      dvz = vel2[2] - vel1[2];
      deltav2 = dvx * dvx + dvy * dvy + dvz * dvz;
      deltav2 -= w * w / (c * c);
#endif
 	}
      else
	{
	  fviscx = 0;
	  fviscy = 0;
	  fviscz = 0;
	  dEdt = 0;
	}


#ifdef VORONOI_SHAPESCHEME
      double s1x = VF[i].cx - p1->x;
      double s1y = VF[i].cy - p1->y;
      double s1z = VF[i].cz - p1->z;

      double s2x = VF[i].cx - p2->x;
      double s2y = VF[i].cy - p2->y;
      double s2z = VF[i].cz - p2->z;

      double d1x = NEAREST_X(p1->x - center1[0]);
      double d1y = NEAREST_Y(p1->y - center1[1]);
      double d1z = NEAREST_Z(p1->z - center1[2]);

      double d2x = NEAREST_X(p2->x - center2[0]);
      double d2y = NEAREST_Y(p2->y - center2[1]);
      double d2z = NEAREST_Z(p2->z - center2[2]);

      double d1 = d1x * d1x + d1y * d1y + d1z * d1z;
      double d2 = d2x * d2x + d2y * d2y + d2z * d2z;

      double square1 = 1 + Stffnss * d1 * pow(volume1, -2.0 / DIMS);
      double square2 = 1 + Stffnss * d2 * pow(volume2, -2.0 / DIMS);

      double curly1 = 1 + Rndnss * (w1 * pow(volume1, -2.0 / DIMS) - SHAPE_FAC);
      double curly2 = 1 + Rndnss * (w2 * pow(volume2, -2.0 / DIMS) - SHAPE_FAC);


      double Q1 = Stffnss / GAMMA_MINUS1 * pressure1 * pow(volume1, 1.0 - 2.0 / DIMS) * curly1;
      double Q2 = Stffnss / GAMMA_MINUS1 * pressure2 * pow(volume2, 1.0 - 2.0 / DIMS) * curly2;

      double L1 = Rndnss / GAMMA_MINUS1 * pressure1 * pow(volume1, 1.0 - 2.0 / DIMS) * square1;
      double L2 = Rndnss / GAMMA_MINUS1 * pressure2 * pow(volume2, 1.0 - 2.0 / DIMS) * square2;

      pressure1 *= (square1 * curly1 +
		    Stffnss * (2.0 / DIMS) / GAMMA_MINUS1 * d1 / pow(volume1,
										  2.0 / DIMS) * curly1 +
		    Rndnss * (2.0 / DIMS) / GAMMA_MINUS1 * w1 / pow(volume1,
										  2.0 / DIMS) * square1);

      pressure2 *= (square2 * curly2 +
		    Stffnss * (2.0 / DIMS) / GAMMA_MINUS1 * d2 / pow(volume2,
										  2.0 / DIMS) * curly2 +
		    Rndnss * (2.0 / DIMS) / GAMMA_MINUS1 * w2 / pow(volume2,
										  2.0 / DIMS) * square2);

      pressure1 += 2 * d1 * Q1 / volume1 + w1 * L1 / volume1;
      pressure2 += 2 * d2 * Q2 / volume2 + w2 * L2 / volume2;

      fac1 = 0.5 * (pressure2 + pressure1) * length / c;
      fac2 = (pressure2 - pressure1) * length / c;

      /* put in a limiter for highly distorted cells */
      if(ee > 0.5 * c)
	fac2 *= 0.5 / (ee / c);

      forcex = -fac1 * cx - fac2 * ex;
      forcey = -fac1 * cy - fac2 * ey;
      forcez = -fac1 * cz - fac2 * ez;

      double g1x = NEAREST_X(VF[i].cx - center1[0]);
      double g1y = NEAREST_Y(VF[i].cy - center1[1]);
      double g1z = NEAREST_Z(VF[i].cz - center1[2]);

      double g2x = NEAREST_X(VF[i].cx - center2[0]);
      double g2y = NEAREST_Y(VF[i].cy - center2[1]);
      double g2z = NEAREST_Z(VF[i].cz - center2[2]);

      double f1x = L1 / volume1 * g1x;
      double f1y = L1 / volume1 * g1y;
      double f1z = L1 / volume1 * g1z;

      double f2x = L2 / volume2 * g2x;
      double f2y = L2 / volume2 * g2y;
      double f2z = L2 / volume2 * g2z;


      double e1x = Q1 / volume1 * d1x;
      double e1y = Q1 / volume1 * d1y;
      double e1z = Q1 / volume1 * d1z;

      double e2x = Q2 / volume2 * d2x;
      double e2y = Q2 / volume2 * d2y;
      double e2z = Q2 / volume2 * d2z;

      double prod2 = 2 * ((s1x * e1x + s1y * e1y + s1z * e1z) - (s2x * e2x + s2y * e2y + s2z * e2z));

      prod2 += ((g2x * f2x + g2y * f2y + g2z * f2z) - (g1x * f1x + g1y * f1y + g1z * f1z));
      prod2 += (VF[i].T_xx + VF[i].T_yy + VF[i].T_zz) * (L2 / volume2 - L1 / volume1);


      double edx = 2 * (e1x - e2x + f2x - f1x);
      double edy = 2 * (e1y - e2y + f2y - f1y);
      double edz = 2 * (e1z - e2z + f2z - f1z);

      double tfx = VF[i].T_xx * edx + VF[i].T_xy * edy + VF[i].T_xz * edz;
      double tfy = VF[i].T_xy * edx + VF[i].T_yy * edy + VF[i].T_yz * edz;
      double tfz = VF[i].T_xz * edx + VF[i].T_yz * edy + VF[i].T_zz * edz;

#ifndef TWODIMS
      double ggx = (L2 / volume2 - L1 / volume1) * VF[i].g_x;
      double ggy = (L2 / volume2 - L1 / volume1) * VF[i].g_y;
      double ggz = (L2 / volume2 - L1 / volume1) * VF[i].g_z;
#else
      double ggx = 0, ggy = 0, ggz = 0;
#endif

      double fac3 = 1.0;

      /* put in a limiter for highly distorted cells */
      if(ee > 0.5 * c)
	fac3 = 0.5 / (ee / c);

      double fshape1_x = fac3 * length / c * (tfx + s1x * prod2 + ggx);
      double fshape1_y = fac3 * length / c * (tfy + s1y * prod2 + ggy);
      double fshape1_z = fac3 * length / c * (tfz + s1z * prod2 + ggz);

      double fshape2_x = fac3 * length / c * (-tfx - s2x * prod2 - ggx);
      double fshape2_y = fac3 * length / c * (-tfy - s2y * prod2 - ggy);
      double fshape2_z = fac3 * length / c * (-tfz - s2z * prod2 - ggz);
#endif

#ifdef VORONOI_MESHRELAX
      fviscx = fviscy = fviscz = dEdt = 0;
#endif

      if(p1->task == ThisTask && q1 >= 0 && q1 < N_gas)
	{
	  if(TimeBinActive[P[q1].TimeBin])
	    {
	      SphP[q1].a.dHydroAccel[0] += (forcex + fviscx) / mass1;
	      SphP[q1].a.dHydroAccel[1] += (forcey + fviscy) / mass1;
	      SphP[q1].a.dHydroAccel[2] += (forcez + fviscz) / mass1;

#ifdef VORONOI_SHAPESCHEME
	      SphP[q1].a.dHydroAccel[0] += fshape1_x / mass1;
	      SphP[q1].a.dHydroAccel[1] += fshape1_y / mass1;
	      SphP[q1].a.dHydroAccel[2] += fshape1_z / mass1;
	      SphP[q1].e.DtEntropy +=  dEdt * GAMMA_MINUS1 
		/ pow(dens1, GAMMA_MINUS1) / (mass1+mass2) 
		/ (square1 * curly1)/ hubble_a2;
#else
	      SphP[q1].e.DtEntropy +=  dEdt * GAMMA_MINUS1 
		/ pow(dens1, GAMMA_MINUS1) / (mass1+mass2)/ hubble_a2;
#endif
	    }
	}

      if(p2->task == ThisTask && q2 >= 0 && q2 < N_gas)
	{
	  if(TimeBinActive[P[q2].TimeBin])
	    {
	      SphP[q2].a.dHydroAccel[0] -= (forcex + fviscx) / mass2;
	      SphP[q2].a.dHydroAccel[1] -= (forcey + fviscy) / mass2;
	      SphP[q2].a.dHydroAccel[2] -= (forcez + fviscz) / mass2;

#ifdef VORONOI_SHAPESCHEME
	      SphP[q2].a.dHydroAccel[0] += fshape2_x / mass2;
	      SphP[q2].a.dHydroAccel[1] += fshape2_y / mass2;
	      SphP[q2].a.dHydroAccel[2] += fshape2_z / mass2;
	      SphP[q2].e.DtEntropy +=  dEdt * GAMMA_MINUS1 
		/ pow(dens2, GAMMA_MINUS1) / (mass2+mass1) 
		/ (square2 * curly2)/ hubble_a2;
#else
	      SphP[q2].e.DtEntropy +=  dEdt * GAMMA_MINUS1 
		/ pow(dens2, GAMMA_MINUS1) / (mass2+mass1)/ hubble_a2;
#endif
	    }
	}
    }


#ifdef VORONOI_SHAPESCHEME
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	{
	  double Q =
	    All.VoronoiStiffNess / GAMMA_MINUS1 * SphP[i].Pressure * pow(SphP[i].Volume, 1.0 - 2.0 / DIMS);
#ifdef VORONOI_VACUUM
#ifdef TWODIMS
			double VACUUM_VOL  = pow(0.07*All.VoronoiVacuumLength,2);
#else
			double VACUUM_VOL  = pow(0.07*All.VoronoiVacuumLength,3);
#endif
		if (SphP[i].Volume > VACUUM_VOL)
				Q = 0.0;
#endif
	  double dx = P[i].Pos[0] - SphP[i].Center[0];
	  double dy = P[i].Pos[1] - SphP[i].Center[1];
	  double dz = P[i].Pos[2] - SphP[i].Center[2];

	  SphP[i].a.dHydroAccel[0] += (-2 * Q * dx) / P[i].Mass;
	  SphP[i].a.dHydroAccel[1] += (-2 * Q * dy) / P[i].Mass;
	  SphP[i].a.dHydroAccel[2] += (-2 * Q * dz) / P[i].Mass;
	}
    }
#endif


#ifdef VORONOI_MESHRELAX
  voronoi_meshrelax();

  myfree(Grad);
  myfree(GradExch);
#endif

  myfree(PrimExch);

  myfree(List_P);
  myfree(ListExports);

  myfree(DT);
  myfree(DP - 5);
  myfree(VF);			/* free the list of faces */
}

#endif
