#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#include "cooling.h"

#if defined(RT_COOLING_PHOTOHEATING)

/* rate1 : photoheating for a blackbody spectrum */
/* rate2 : recombination cooling rate */
/* rate3 : collisional ionization cooling rate */
/* rate4 : collisional excitation cooling rate */
/* rate5 : Bremsstrahlung cooling rate */

/* now do the heating (note: we know how many photons we absorbed) */

#ifndef RT_MULTI_FREQUENCY
double rt_DoHeating(int i, double dt_internal)
{
  double sigma, a3inv, nH;
  double hubble_a, rate, du, de, c_light;
  double n_gamma, nHI;
  double E;
  
  if(All.ComovingIntegrationOn)
    {
      a3inv = 1.0 / All.Time / All.Time / All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    {
      a3inv = hubble_a = 1.0;
    }
  
  c_light = C / All.UnitVelocity_in_cm_per_s;
  
  nH = HYDROGEN_MASSFRAC * SphP[i].d.Density * a3inv / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
  nHI = SphP[i].HI * nH;
  
  sigma = 1.63e-18 / All.UnitLength_in_cm / All.UnitLength_in_cm * All.HubbleParam * All.HubbleParam;
  n_gamma = SphP[i].n_gamma[0] / P[i].Mass * a3inv;
  E = 30.0 * ELECTRONVOLT_IN_ERGS / All.UnitEnergy_in_cgs * All.HubbleParam;
  rate = nHI * c_light * E * sigma * n_gamma;
  
  du = rate * dt_internal / hubble_a / (SphP[i].d.Density * a3inv);
  de = du * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
  
  return de / dt_internal;
}

#else

double rt_DoHeating(int i, double dt_internal)
{
  int j;
  double a3inv, nH;
  double hubble_a, rate, du, de;
  double nHI, n_gamma;
  double c_light;

#ifdef RT_INCLUDE_HE
  double nHeI, nHeII, nHeIII;
#endif

  c_light = C / All.UnitVelocity_in_cm_per_s;
  
  if(All.ComovingIntegrationOn)
    {
      a3inv = All.Time / All.Time / All.Time;
      hubble_a = hubble_function(All.Time);
    }
  else
    {
      a3inv = hubble_a = 1.0;
    }
  
  nH = HYDROGEN_MASSFRAC * SphP[i].d.Density * a3inv / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;
  nHI = SphP[i].HI * nH;
  
#ifdef RT_INCLUDE_HE
  nHeI = SphP[i].HeI * nH;
  nHeII = SphP[i].HeII * nH;
  nHeIII = SphP[i].HeIII * nH;
#endif

  for(j = 0, rate = 0; j < N_BINS; j++)
    {
      n_gamma = SphP[i].n_gamma[j] / P[i].Mass * a3inv;

      if(nu[j] >= 13.6)
	rate += nHI * c_light * rt_sigma_HI[j] * G_HI[j] * n_gamma;
      
#ifdef RT_INCLUDE_HE
      if(nu[j] >= 24.6)
	rate += nHeI * c_light * rt_sigma_HeI[j] * G_HeI[j] * n_gamma;
	
      if(nu[j] >= 54.4)
	rate += nHeII * c_light * rt_sigma_HeII[j] * G_HeII[j] * n_gamma;
#endif
    }
  
  du = rate * dt_internal / hubble_a / (SphP[i].d.Density * a3inv);
  de = du * GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);

  return de / dt_internal;
}

#endif

double rt_DoCooling(int i, double dt_internal)
{
  double dtime, a3inv, iter;
  double lambda, du, de;
  double u_old, u_lower, u_upper, ratefact, u;
  double fac_u_to_entr, entropy;

  if(All.ComovingIntegrationOn)
    {
      dtime = dt_internal / hubble_function(All.Time);
      a3inv = 1.0 / All.Time / All.Time / All.Time;
    }
  else
    {
      dtime = dt_internal;
      a3inv = 1.0;
    }

  fac_u_to_entr = GAMMA_MINUS1 / pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1);
    
  entropy = SphP[i].Entropy;

  /* do the cooling */
  lambda = rt_get_cooling_rate(i, entropy);
  du = lambda * dtime / (SphP[i].d.Density * a3inv);
  de = du * fac_u_to_entr;

  if(fabs(de) < 0.2 * entropy)
    {
      /* cooling is slow, we can do it explicitly */
      return de / dt_internal;
    }
  else
    {
      /* rapid cooling. Better calculate an implicit solution, which is determined by bisection */
      u_old = entropy / fac_u_to_entr;
      u_lower = u_old / sqrt(1.1);
      u_upper = u_old * sqrt(1.1);
      ratefact = dtime / (SphP[i].d.Density * a3inv);
      iter = 0;

      /* bracketing */
      while(u_lower - u_old - ratefact * rt_get_cooling_rate(i, u_lower * fac_u_to_entr) > 0)
        {
          u_upper /= 1.1;
          u_lower /= 1.1;


          if(iter++ >= 1000) //MAXITER)
            terminate("bracketing failure");
        }

      /* bisection */
      iter = 0;
      do
        {
          u = 0.5 * (u_lower + u_upper);

          if(u - u_old - ratefact * rt_get_cooling_rate(i, u * fac_u_to_entr) > 0)
            u_upper = u;
          else
            u_lower = u;

          du = u_upper - u_lower;

          iter++;

          if(iter >= (MAXITER - 10))
            printf("u= %g\n", u);

          if(iter >= MAXITER)
            terminate("convergence failure");
        }
      while(fabs(du / u) > 1.0e-6);

      du = u - u_old;

      return du * fac_u_to_entr / dt_internal;
    }

}

/* returns cooling rate */
double rt_get_cooling_rate(int i, double entropy)
{
  double Lambda;
  double temp, molecular_weight;
  double a3inv;
  double nH;
  double rate2, rate3, rate4, rate5;
  double de2, de3, de4, de5;
#ifdef RT_INCLUDE_HE
  double rateHe2, rateHe3, rateHe4, rateHe5;
  double deHe2, deHe3, deHe4, deHe5;
#endif

  double fac = All.UnitTime_in_s / pow(All.UnitLength_in_cm, 3) / All.UnitEnergy_in_cgs * pow(All.HubbleParam, 3);

  if(All.ComovingIntegrationOn)
    a3inv = 1 / All.Time / All.Time / All.Time;
  else
    a3inv = 1;

  nH = HYDROGEN_MASSFRAC * SphP[i].d.Density * a3inv / PROTONMASS * All.UnitMass_in_g / All.HubbleParam;	//physical
  molecular_weight = 4 / (1 + 3 * HYDROGEN_MASSFRAC + 4 * HYDROGEN_MASSFRAC * SphP[i].elec);

  temp = entropy * pow(SphP[i].d.Density * a3inv, GAMMA_MINUS1) *
    molecular_weight * PROTONMASS / All.UnitMass_in_g * All.HubbleParam /
    BOLTZMANN * All.UnitEnergy_in_cgs / All.HubbleParam;

  /* all rates in erg cm^3 s^-1 in code units */
  /* recombination cooling rate */
  rate2 = 8.7e-27 * sqrt(temp) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7)) * fac;
  de2 = SphP[i].HII * nH * SphP[i].elec * nH * rate2;

  /* collisional ionization cooling rate */
  rate3 = 1.27e-21 * sqrt(temp) * exp(-157809.1 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
  de3 = SphP[i].HI * nH * SphP[i].elec * nH * rate3;

  /* collisional excitation cooling rate */
  rate4 = 7.5e-19 / (1.0 + sqrt(temp / 1e5)) * exp(-118348 / temp) * fac;
  de4 = SphP[i].HI * nH * SphP[i].elec * nH * rate4;

  /* Bremsstrahlung cooling rate */
  rate5 = 1.42e-27 * sqrt(temp) * fac;
  de5 = SphP[i].HII * nH * SphP[i].elec * nH * rate5;

  Lambda = de2 + de3 + de4 + de5;

  /* inverse Compton cooling rate */
  if(All.ComovingIntegrationOn)
    Lambda += 5.406e-36 * SphP[i].elec * (temp - (2.73 / All.Time)) / pow(All.Time, 4) * fac;
  
#ifdef RT_INCLUDE_HE
  /* recombination cooling rate */
  rateHe2 = 1.55e-26 * pow(temp, 0.3647) * fac;
  deHe2 = SphP[i].HeII * nH * SphP[i].elec * nH * rateHe2;

  rateHe2 = 3.48e-26 * sqrt(temp) * pow(temp / 1e3, -0.2) / (1.0 + pow(temp / 1e6, 0.7)) * fac;
  deHe2 += SphP[i].HeIII * nH * SphP[i].elec * nH * rateHe2;

  /* collisional ionization cooling rate */
  rateHe3 = 9.38e-22 * sqrt(temp) * exp(-285335.4 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
  deHe3 = SphP[i].HeI * nH * SphP[i].elec * nH * rateHe3;

  rateHe3 = 4.95e-22 * sqrt(temp) * exp(-631515 / temp) / (1.0 + sqrt(temp / 1e5)) * fac;
  deHe3 += SphP[i].HeII * nH * SphP[i].elec * nH * rateHe3;

  rateHe3 = 5.01e-27 * pow(temp, -0.1687) / (1.0 + sqrt(temp / 1e5)) * exp(-55338 / temp) * fac;
  rateHe3 *= pow(All.HubbleParam / All.UnitLength_in_cm, 3);
  deHe3 += SphP[i].HeII * nH * SphP[i].elec * nH * SphP[i].elec * nH * rateHe3;

  /* collisional excitation cooling rate */
  rateHe4 = 5.54e-17 * pow(temp, -0.397) / (1.0 + sqrt(temp / 1e5)) * exp(-473638 / temp) * fac;
  deHe4 = SphP[i].HeII * nH * SphP[i].elec * nH * rateHe4;

  rateHe4 = 9.10e-27 * pow(temp, -0.1687) / (1.0 + sqrt(temp / 1e5)) * exp(-13179 / temp) * fac;
  rateHe4 *= pow(All.HubbleParam / All.UnitLength_in_cm, 3);
  deHe4 += SphP[i].HeII * nH * SphP[i].elec * nH * SphP[i].elec * nH * rateHe4;

  /* Bremsstrahlung cooling rate */
  rateHe5 = 1.42e-27 * sqrt(temp) * fac;
  deHe5 = (SphP[i].HeII + 4.0 * SphP[i].HeIII + SphP[i].elec) * nH * rateHe5;

  Lambda += deHe2 + deHe3 + deHe4 + deHe5;
#endif

  return -Lambda;
}



#endif
