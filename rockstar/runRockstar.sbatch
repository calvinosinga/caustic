#!/bin/bash

#SBATCH --job-name=RS_25
#SBATCH --output=log_script_run%j.dat
#SBATCH --time=16:00:00
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=2000

# ---------------------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------------------

box="cbol_l0025_n0256_gde"
mdef="200m"

# ---------------------------------------------------------------------------------------
# CONSTANTS
# ---------------------------------------------------------------------------------------

out_dir="/n/home05/diemer/Simulations/Box_${box}/Rockstar${mdef}/"
exe="/n/home05/diemer/Software/rockstar/rockstar"

# ---------------------------------------------------------------------------------------
# SCRIPT
# ---------------------------------------------------------------------------------------

time_begin=`date +%s`

echo "SBatch file: Executing server..."
${exe} -c config.cfg >& log_server_run${SLURM_JOBID}.dat &

# The server will create a copy of the config file which will be used by the readers
# and writers. So we wait until the server has successfully started, and move the 
# config file to the run directory.
echo "SBatch file: Waiting for automatic config file to appear"
while [ ! -e "${out_dir}/auto-rockstar.cfg" ]
do
    sleep 1
done
echo "SBatch file: Found config file."

# Check if a run was previously started. If so, re-start.
if [ -f "_RUN_STARTED" ]
then
    echo "SBatch file: Re-starting readers."
    srun ${exe} -c ${out_dir}/restart.cfg >& log_client_run${SLURM_JOBID}.dat
else
    echo "SBatch file: Starting readers."
    touch _RUN_STARTED
    srun ${exe} -c ${out_dir}/auto-rockstar.cfg >& log_client_run${SLURM_JOBID}.dat
fi

# Finish script, write time it took.
rm _RUN_STARTED
time_end=`date +%s`
time_spent=$[$time_end-$time_begin]
time_min=$(echo "scale=2; {r = $time_spent/60}; r" | bc -l)
echo "SBatch file: Done, script ran $time_min minutes."
