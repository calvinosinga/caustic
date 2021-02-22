#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        PERIODIC\n"
"        PMGRID=512\n"
"        PEANOHILBERT\n"
"        WALLCLOCK\n"
"        MYSORT\n"
"        NO_ISEND_IRECV_IN_DOMAIN\n"
"        FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG\n"
"        NOTYPEPREFIX_FFTW\n"
"        DEBUG\n"
"        DISTORTIONTENSORPS\n"
"        OUTPUT_DISTORTIONTENSORPS\n"
"        OUTPUT_TIDALTENSORPS\n"
"        CAUSTIC_FINDER=2\n"
"        OUTPUT_LAST_CAUSTIC\n"
"        COMOVING_DISTORTION\n"
"\n");
}
