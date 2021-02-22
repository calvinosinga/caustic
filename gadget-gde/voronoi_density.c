#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef VORONOI
#include "voronoi.h"

void voronoi_density(void)
{
  int i;
  int dt_step;
  double dt_entr;


  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
	{
      SphP[i].d.Density = P[i].Mass / SphP[i].Volume;
      SphP[i].v.DivVel = 0;

      dt_step = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0);
      dt_entr = (All.Ti_Current - (P[i].Ti_begstep + dt_step / 2)) * All.Timebase_interval;


#if defined(VORONOI_MESHRELAX)

#ifdef VORONOI_MESHRELAX_KEEPRESSURE
      SphP[i].Entropy = SphP[i].Pressure / (GAMMA_MINUS1 * SphP[i].d.Density);
#else 
      SphP[i].Pressure = GAMMA_MINUS1 * SphP[i].Entropy * SphP[i].d.Density;
#endif

#else
      SphP[i].Pressure = (SphP[i].Entropy + SphP[i].e.DtEntropy * dt_entr) * pow(SphP[i].d.Density, GAMMA);
#endif

      SphP[i].MaxSignalVel = sqrt(GAMMA * SphP[i].Pressure / SphP[i].d.Density);
	}
    }


#if defined(VORONOI_CFL_COND) || defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC) || defined(VORONOI_ENERGYGRAD)  
  voronoi_exchange_ghost_variables();
  //  voronoi_setup_exchange();
  int q1,q2,li,ri;
  double c,cx, cy, cz;
#if defined(VORONOI_CFL_COND)
  double v_dummy, c_sound;
#endif
  double mass1, mass2, dens1, dens2;
#if defined(VORONOI_ENERGYGRAD)
  double entr1, entr2;
#endif
  double h1,h2;
  double invvol1, invvol2;
  point *p1, *p2;


#if defined(VORONOI_ENERGYGRAD)  
			SphP[i].gradu[0] =  0.0;
			SphP[i].gradu[1] =  0.0;
			SphP[i].gradu[2] =  0.0;
#endif	

#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)
  double  ex, ey, ez, dvx,dvy,dvz, length;
  MyFloat *vel1, *vel2;
for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) 
	{
	if(P[i].Type == 0)
			{
			SphP[i].v.DivVel = 0.0;
			SphP[i].r.Rot[2] = 0.0;
			SphP[i].r.Rot[0] = 0.0;
			SphP[i].r.Rot[1] = 0.0;
			SphP[i].r.CurlVel= 0.0;

// #if defined(KH_MIXING)  
// 			SphP[i].dV_dt = 0.0;
// #endif
#if defined(VORONOI_PARTIAL)  
 			SphP[i].surface =  0.0;
#endif	
			}
	}
#endif


for(i = 0; i < Nvf; i++)
{
      p1 = &DP[VF[i].p1];
      p2 = &DP[VF[i].p2];

      
      q1 = li = p1->index;
      q2 = ri = p2->index;

      if(!((li >= 0 && li < N_gas && p1->task == ThisTask) ||
	   (ri >= 0 && ri < N_gas && p2->task == ThisTask)))
	continue;

      //      if(!((li >= 0) &&
      //   (ri >= 0)))
      //continue;

      if(li >= N_gas && p1->task == ThisTask)
	li -= N_gas;

      if(ri >= N_gas && p2->task == ThisTask)
	ri -= N_gas;

      if(p1->task == ThisTask)
	{
          dens1 = SphP[li].d.Density;
          mass1 = P[li].Mass;
#if defined(VORONOI_ENERGYGRAD)
          entr1 = SphP[li].Entropy;
#endif  
          invvol1 = dens1 / mass1;
// fprintf(stdout,"\n !! invvol1 %f dens1  %f\n", invvol1,SphP[li].d.Density);
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)         
          vel1 = SphP[li].VelPred;
#endif   			
	}
      else
	{
	  dens1 = PrimExch[q1].Density;
	  mass1 = PrimExch[q1].Mass;
#if defined(VORONOI_ENERGYGRAD)
	  entr1 = PrimExch[q1].Entropy;
#endif  
	  invvol1 =  dens1 / mass1;
// 	  fprintf(stdout,"\n invvol1 %f \n", invvol1);
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)         
	  vel1 = PrimExch[q1].VelPred;
#endif   			  
	}

      if(p2->task == ThisTask)
	{
          dens2 = SphP[ri].d.Density;
          mass2 = P[ri].Mass;
#if defined(VORONOI_ENERGYGRAD)
          entr2 = SphP[ri].Entropy;
#endif  
          invvol2 =  dens2 / mass2;
// fprintf(stdout,"\n !! invvol2 %f dens2  %f\n", invvol2,SphP[ri].d.Density);
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)         
          vel2 = SphP[ri].VelPred;
#endif  
	}
      else
	{
	  dens2 = PrimExch[q2].Density;
	  mass2 = PrimExch[q2].Mass;
#if defined(VORONOI_ENERGYGRAD)
	  entr2 = PrimExch[q2].Entropy;
#endif  
	  invvol2 = dens2 / mass2;
// 	  fprintf(stdout,"\n invvol2 %f \n", invvol2);
if (mass1 == 0.0) fprintf(stdout,"\n invvol1 %f mass1 %f \n", invvol1,mass1);
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)         
	  vel2 = PrimExch[q2].VelPred;
#endif  
	}

if (invvol1 == 0.0) fprintf(stdout,"\n invvol1 %f mass1 %f \n", invvol1,mass1);
if (invvol2 == 0.0) fprintf(stdout,"\n invvol2 %f mass2 %f \n", invvol2,mass2);
      cx = p2->x - p1->x;
      cy = p2->y - p1->y;
      cz = p2->z - p1->z;

#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)   
			dvx = (vel2[0] - vel1[0]);
			dvy = (vel2[1] - vel1[1]);
			dvz = (vel2[2] - vel1[2]);
  
      ex = VF[i].cx - 0.5 * (p1->x + p2->x);
      ey = VF[i].cy - 0.5 * (p1->y + p2->y);
      ez = VF[i].cz - 0.5 * (p1->z + p2->z);
      length = VF[i].area;										/* length/area of common face */
#endif 
      
	h1 = 0.0;
	h2 = 0.0;
#ifdef TWODIMS
			if (invvol1 > 0.0) 
			h1 = sqrt(1.0/invvol1 / M_PI);
			if (invvol2 > 0.0)
			h2 = sqrt(1.0/invvol2 / M_PI);
#else
			if (invvol1 > 0.0)
			h1 = pow(1.0/invvol1 / (4.0 / 3 * M_PI), 1.0 / 3);
			if (invvol2 > 0.0)
			h2 = pow(1.0/invvol2 / (4.0 / 3 * M_PI), 1.0 / 3);
#endif
      c = sqrt(cx * cx + cy * cy + cz * cz   + 0.000001 * h1 * h2);	/* distance of the two points */


  if(p1->task == ThisTask && q1 >= 0 && q1 < N_gas)
//     if( q1 >= 0 && q1 < N_gas)
	{
		if(TimeBinActive[P[q1].TimeBin]) {
#if defined(VORONOI_CFL_COND)
			c_sound = sqrt(GAMMA * SphP[q1].Pressure / dens1);
			v_dummy = c_sound* h1 / c;
// 			v_dummy *= 2.0; 
#ifdef TWODIMS
                        if (((2.0*length) < c) && ((8.0*length) > c) ) v_dummy = c_sound* h2 / (2.0*length);
                        if (((2.0*length) < c) && ((8.0*length) <=c) ) v_dummy *= 8.0;
#else
                        if ((sqrt(2.0*length) < c) && (sqrt(32.0*length) > c) ) v_dummy = c_sound* h2 / sqrt(2.0*length);
                        if ((sqrt(2.0*length) < c) && (sqrt(32.0*length) <=c) ) v_dummy *= 8.0;
#endif
			if (SphP[q1].MaxSignalVel < v_dummy) 	SphP[q1].MaxSignalVel = v_dummy;
#endif
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)  
 			SphP[q1].v.DivVel +=  -((0.5 *cx+ex) * dvx + (0.5 *cy+ey) * dvy + (0.5 *cz+ez) * dvz)/ c *length*invvol1;
			SphP[q1].r.Rot[2] += 0.5 *(cx * dvy - cy * dvx)/ c *length *invvol1;
			SphP[q1].r.Rot[0] += 0.5 *(cy * dvz - cz * dvy)/ c *length *invvol1;
			SphP[q1].r.Rot[1] += 0.5 *(cz * dvx - cx * dvz)/ c *length *invvol1;
// 			SphP[q1].r.Rot[2] += 0.5 *(ex * dvy - ey * dvx)/ c *length *invvol1;
// 			SphP[q1].r.Rot[0] += 0.5 *(ey * dvz - ez * dvy)/ c *length *invvol1;
// 			SphP[q1].r.Rot[1] += 0.5 *(ez * dvx - ex * dvz)/ c *length *invvol1;
// 			fprintf(stdout,"\n SphP[q1].v.DivVel %f  1.0/ c *length*invvol1 %f  c %e 1 %e  1 %e 1 %e\n",  -((0.5 *cx+ex) * dvx + (0.5 *cy+ey) * dvy + (0.5 *cz+ez) * dvz)/ c *length*invvol1,1.0/ c *length*invvol1,c, h1 * h2,invvol1,1.0/invvol2);
#endif		
// #if defined(KH_MIXING)  
//  			SphP[q1].dV_dt +=  -(1/invvol2 - 1/invvol1)*((0.5 *cx+ex) * vel1[0] + (0.5 *cy+ey) * vel1[1] + (0.5 *cz+ez) * vel1[2])/ c *length*invvol1;
// #endif	
#if defined(VORONOI_ENERGYGRAD)  
			double u_diff = -(entr1*pow(dens1,GAMMA_MINUS1)-entr2*pow(dens2,GAMMA_MINUS1))/GAMMA_MINUS1;
			SphP[q1].gradu[0] +=  (0.5 *cx) * u_diff / c *length*invvol1;
			SphP[q1].gradu[1] +=  (0.5 *cy) * u_diff / c *length*invvol1;
			SphP[q1].gradu[2] +=  (0.5 *cz) * u_diff / c *length*invvol1;
#endif	
#if defined(VORONOI_PARTIAL)  
 			SphP[q1].surface +=  length;
#endif	
		}
	}
	if(p2->task == ThisTask && q2 >= 0 && q2 < N_gas)
// 	if(q2 >= 0 && q2 < N_gas)
	{
		if(TimeBinActive[P[q2].TimeBin]) {
#if defined(VORONOI_CFL_COND)
			c_sound = sqrt(GAMMA * SphP[q2].Pressure / dens2);
			v_dummy = c_sound* h2 / c;
// 			v_dummy *= 2.0;
#ifdef TWODIMS
                        if (((2.0*length) < c) && ((8.0*length) > c) ) v_dummy = c_sound* h2 / (2.0*length);
                        if (((2.0*length) < c) && ((8.0*length) <=c) ) v_dummy *= 8.0;
#else
                        if ((sqrt(2.0*length) < c) && (sqrt(32.0*length) > c) ) v_dummy = c_sound* h2 / sqrt(2.0*length);
                        if ((sqrt(2.0*length) < c) && (sqrt(32.0*length) <=c) ) v_dummy *= 8.0;
#endif

			if (SphP[q2].MaxSignalVel < v_dummy) 	SphP[q2].MaxSignalVel = v_dummy;
#endif
#if defined(VORONOI_BALSARA)  || defined(VORONOI_TIME_DEP_ART_VISC)     		
			SphP[q2].v.DivVel +=  -(( 0.5 *cx-ex) * dvx + ( 0.5 *cy-ey) * dvy + (0.5 *cz-ez) * dvz)/ c *length*invvol2;
			SphP[q2].r.Rot[2] += 0.5 *(cx * dvy - cy * dvx)/ c *length *invvol2;
			SphP[q2].r.Rot[0] += 0.5 *(cy * dvz - cz * dvy)/ c *length *invvol2;
			SphP[q2].r.Rot[1] += 0.5 *(cz * dvx - cx * dvz)/ c *length *invvol2;
// 			SphP[q2].r.Rot[2] -= 0.5 *(ex * dvy - ey * dvx)/ c *length *invvol2;
// 			SphP[q2].r.Rot[0] -= 0.5 *(ey * dvz - ez * dvy)/ c *length *invvol2;
// 			SphP[q2].r.Rot[1] -= 0.5 *(ez * dvx - ex * dvz)/ c *length *invvol2;

// 			fprintf(stdout,"\n SphP[q2].v.DivVel %f  1.0/ c *length*invvol2 %f c %e 1 %e  1 %e 1 %e\n",-(( 0.5 *cx-ex) * dvx + ( 0.5 *cy-ey) * dvy + (0.5 *cz-ez) * dvz)/ c *length*invvol2 ,1.0/ c *length*invvol2,c , h1 * h2,1.0/invvol1,1.0/invvol2);
#endif
// #if defined(KH_MIXING)  
//  			SphP[q2].dV_dt +=  (1/invvol2 - 1/invvol1)*((0.5 *cx-ex) * vel2[0] + (0.5 *cy-ey) * vel2[1] + (0.5 *cz-ez) * vel2[2])/ c *length*invvol1;
// #endif	
#if defined(VORONOI_ENERGYGRAD)  
			double u_diff = (entr1*pow(dens1,GAMMA_MINUS1)-entr2*pow(dens2,GAMMA_MINUS1))/GAMMA_MINUS1;
			SphP[q2].gradu[0] -=  -(0.5 *cx) * u_diff / c *length*invvol2;
			SphP[q2].gradu[1] -=  -(0.5 *cy) * u_diff / c *length*invvol2;
			SphP[q2].gradu[2] -=  -(0.5 *cz) * u_diff / c *length*invvol2;
#endif	
#if defined(VORONOI_PARTIAL)  
 			SphP[q2].surface +=  length;
#endif	
		}
	}




} /* endfor(i = 0; i < Nvf; i++) */

for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) 
	{
	if(P[i].Type == 0)
			{
			SphP[i].r.CurlVel = sqrt(SphP[i].r.Rot[0] * SphP[i].r.Rot[0] +
			SphP[i].r.Rot[1] * SphP[i].r.Rot[1] +
			SphP[i].r.Rot[2] * SphP[i].r.Rot[2]);
			}
	}

#endif
}
#endif
