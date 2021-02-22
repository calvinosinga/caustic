#!/bin/bash

#SBATCH --job-name=test-gadget
#SBATCH --output=out_script_L0063_N0256_CBol_gde_run%j.dat
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=sandyb
#SBATCH --account=pi-kravtsov
#SBATCH --constraint=ib
#SBATCH --shared

# Constants for this run

run="L0063_N0256_CBol_gde"
nproc=16

# General constants

folder="Box_${run}"

#exeGadget="/home/diemer/Software/uberl-gadget2/L-Gadget2"
exeGadget="/home/diemer/Software/gadget-gde/P-Gadget3"
pathSim="/home/diemer/Simulations/3_Sim/"

# Now run Gadget. uberLGadget outputs a flag file called restartme if a
# run aborts.

time_begin=`date +%s`

if [ -f "_RUN_STARTED" ]
then
    echo "Re-starting LGadget"
    mpirun -np $nproc $exeGadget ${pathSim}${folder}/lgadget_${run}.param 1 >> ${pathSim}${folder}/out_gadget_${run}_run${SLURM_JOBID}.dat 2>&1
else
    #echo "Running 2LPTic"
    #mpirun -np $nproc /home/diemer/Software/2LPTic/2LPTic /home/diemer/Simulations/2_ICs/${folder}/2lptic_${run}.param >> /home/diemer/Simulations/2_ICs/${folder}/output_2lptic.dat 2>&1

    echo "Running LGadget"
    touch _RUN_STARTED
    mpirun -np $nproc $exeGadget ${pathSim}${folder}/lgadget_${run}.param >> ${pathSim}${folder}/out_gadget_${run}_run${SLURM_JOBID}.dat 2>&1
fi

# Finish, print time

time_end=`date +%s`
time_spent=$[$time_end-$time_begin]
time_min=$(echo "scale=2; {r = $time_spent/60}; r" | bc -l)
echo "Done, script ran $time_min minutes."
