#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


#ifdef RADTRANSFER
#ifdef EDDINGTON_TENSOR_GAS
#ifdef RT_MULTI_FREQUENCY

void gas_lum(void)
{
  int i, j;
  double *je;

#ifndef EDDINGTON_TENSOR_STARS
  for(j = 0; j < N_gas; j++)
    if(P[j].Type == 0)
      for(i = 0; i < N_BINS; i++)
	SphP[j].Je[i] = 0.0;
#endif

  je = (double *) mymalloc("je", N_BINS * sizeof(double));

  for(j = 0; j < N_gas; j++)
    if(P[j].Type == 0)
      {
	rt_get_lum_gas(j, je);

	for(i = 0; i < N_BINS; i++)
	  SphP[j].Je[i] += je[i];
      }

  myfree(je);
}

#endif
#endif
#endif
