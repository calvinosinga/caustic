#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef SUB_TURB_DRIVING


static int Nsubs;

static struct subtable_data
{
  double Msub, Vmax, RadVmax, Conc, Pos[3], Vel[3], Acc[3];
}
*SubTable;


static double rhocrit;


void sub_turb_read_table(void)
{
  FILE *fd;
  int i;

  if(!(fd = fopen("substructures.txt", "r")))
    terminate("can't open file");

  fscanf(fd, " %d ", &Nsubs);

  SubTable = mymalloc("SubTable", Nsubs * sizeof(struct subtable_data));

  for(i = 0; i < Nsubs; i++)
    {
      if(fscanf(fd, " %lg %lg %lg  %lg   %lg %lg %lg  %lg %lg %lg", &SubTable[i].Msub, &SubTable[i].Vmax,
          &SubTable[i].RadVmax, &SubTable[i].Conc, &SubTable[i].Pos[0], &SubTable[i].Pos[1], &SubTable[i].Pos[2],
          &SubTable[i].Vel[0], &SubTable[i].Vel[1], &SubTable[i].Vel[2]) != 10)
        terminate("read problem");

      if(ThisTask == 0)
	printf(" %10g %10g %10g  %10g   %10g %10g %10g  %10g %10g %10g\n", SubTable[i].Msub, SubTable[i].Vmax,
	       SubTable[i].RadVmax, SubTable[i].Conc, SubTable[i].Pos[0], SubTable[i].Pos[1], SubTable[i].Pos[2],
	       SubTable[i].Vel[0], SubTable[i].Vel[1], SubTable[i].Vel[2]);
    }
  fclose(fd);

  if(ThisTask == 0)
    printf("\nHave read %d subhalo perturbers\n\n", Nsubs);
  
  rhocrit = 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
}



void sub_turb_add_forces(void)
{
  int i, n;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
        {
          for(n = 0; n < Nsubs; n++)
            {
              double dx = P[i].Pos[0] - SubTable[n].Pos[0];
              double dy = P[i].Pos[1] - SubTable[n].Pos[1];
              double dz = P[i].Pos[2] - SubTable[n].Pos[2];

              double r = sqrt(dx * dx + dy * dy + dz * dz);

	      double rfid = 2 * SubTable[n].RadVmax;
	      
	      rfid = 20.0;

	      if(r < rfid)
		{
		  double vsub = sqrt(pow(SubTable[n].Vel[0], 2) + pow(SubTable[n].Vel[1], 2) + pow(SubTable[n].Vel[2], 2));

		  double aram = 0.25 * vsub * vsub / (2 * rfid);

                  P[i].g.GravAccel[0] += aram * SubTable[n].Vel[0] / vsub;
		  P[i].g.GravAccel[1] += aram * SubTable[n].Vel[1] / vsub;
		  P[i].g.GravAccel[2] += aram * SubTable[n].Vel[2] / vsub;
                }
            }
        }
    }
}








void sub_turb_add_forces_dm_gravity(void)
{
  int i, n;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type == 0)
        {
          for(n = 0; n < Nsubs; n++)
            {
              double dx = P[i].Pos[0] - SubTable[n].Pos[0];
              double dy = P[i].Pos[1] - SubTable[n].Pos[1];
              double dz = P[i].Pos[2] - SubTable[n].Pos[2];

              double r = sqrt(dx * dx + dy * dy + dz * dz);

              double m = sub_turb_enclosed_mass(r, SubTable[n].Msub, SubTable[n].Vmax, SubTable[n].RadVmax,
                  SubTable[n].Conc);

	      //	      printf("SubTable[n].RadVmax = %g         m=%g\n", SubTable[n].RadVmax, m);

              if(r > 0)
                {
                  r += 0.01 * SubTable[n].RadVmax / 2.163;

                  P[i].g.GravAccel[0] += -All.G * m * dx / (r * r * r);
                  P[i].g.GravAccel[1] += -All.G * m * dy / (r * r * r);
                  P[i].g.GravAccel[2] += -All.G * m * dz / (r * r * r);
                }
            }
        }
    }
}


double sub_turb_enclosed_mass(double r, double msub, double vmax, double radvmax, double c)
{
  double rs = radvmax / 2.163;
  double delta_c = 200.0 / 3 * c * c * c / (log(1 + c) - c / (1 + c));
  double x = r / rs;

  double m = 4 * M_PI * delta_c * rhocrit * rs * rs * rs * (log(1 + x) - x / (1 + x));

  if(m > msub)
    m = msub;

  return m;
}


void sub_turb_move_perturbers(double t0, double t1)
{
  double dt = t1 - t0;
  int i, j;

  for(i=0; i< Nsubs; i++)
    sub_turb_parent_halo_accel(SubTable[i].Pos[0], SubTable[i].Pos[1], SubTable[i].Pos[2], SubTable[i].Acc);

  for(i=0; i< Nsubs; i++)
    for(j=0;j<3;j++)
      SubTable[i].Vel[j] += 0.5 * dt * SubTable[i].Acc[j];

  for(i=0; i< Nsubs; i++)
    for(j=0;j<3;j++)
      SubTable[i].Pos[j] += dt * SubTable[i].Vel[j];

  for(i=0; i< Nsubs; i++)
    sub_turb_parent_halo_accel(SubTable[i].Pos[0], SubTable[i].Pos[1], SubTable[i].Pos[2], SubTable[i].Acc);

  for(i=0; i< Nsubs; i++)
    for(j=0;j<3;j++)
      SubTable[i].Vel[j] += 0.5 * dt * SubTable[i].Acc[j];


  if(ThisTask == 0)
    {
      for(i=0; i< Nsubs; i++)
	{
	  double r = sqrt(SubTable[i].Pos[0] * SubTable[i].Pos[0] + 
			  SubTable[i].Pos[1] * SubTable[i].Pos[1] + 
			  SubTable[i].Pos[2] * SubTable[i].Pos[2]);
	  double v = sqrt(SubTable[i].Vel[0] * SubTable[i].Vel[0] + 
			  SubTable[i].Vel[1] * SubTable[i].Vel[1] + 
			  SubTable[i].Vel[2] * SubTable[i].Vel[2]);
	  printf("i=%d r=%g  v=%g\n", i, r, v);
	}
    }

}



void sub_turb_parent_halo_accel(double dx, double dy, double dz, double *acc)
{
  double r, m;

  r = sqrt(dx * dx + dy * dy + dz * dz);

  /*
  if(r > ISO_R200)
    m = ISO_M200;
  else
  */
  m = ISO_M200 * r / ISO_R200;

  if(r > 0)
    {
      acc[0] = -All.G * m * dx / r / (r * r + ISO_Eps * ISO_Eps);
      acc[1] = -All.G * m * dy / r / (r * r + ISO_Eps * ISO_Eps);
      acc[2] = -All.G * m * dz / r / (r * r + ISO_Eps * ISO_Eps);
    }
  else
    {
      acc[0] = 0;
      acc[1] = 0;
      acc[2] = 0;
    }
}


#endif
