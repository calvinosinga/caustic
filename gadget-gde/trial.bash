#!/bin/bash
mpicc -Wall -g -I $FFTW_INCDIR -c -o -pm_periodic.o pm_periodic.c
