#!/bin/bash

#SBATCH --job-name=CT_25
#SBATCH --output=log_trees_run%j.dat
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10GB

# ---------------------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------------------

box="L0025_N0256_GDE"
mdef="200m"

# ---------------------------------------------------------------------------------------
# CONSTANTS
# ---------------------------------------------------------------------------------------

out_dir="/n/home05/diemer/Simulations/Box_${box}/Rockstar${mdef}/"

# ---------------------------------------------------------------------------------------
# SCRIPT
# ---------------------------------------------------------------------------------------

# These are the commands which need to be executed to generate both the 
# merger trees, and the hlist files from the merger trees.

# This first commands generates the merger_tree.cfg automatic config file.

echo "SBatch script: Generating config file..."
perl /n/home05/diemer/Software/rockstar/scripts/gen_merger_cfg.pl config_${mdef}.cfg

# The rest needs to be executed in the consistent trees directory
echo "SBatch script: Going to consistent trees directory..."
cd /n/home05/diemer/Software/consistent-trees

echo "SBatch script: Running merger trees..."
perl do_merger_tree.pl ${out_dir}outputs/merger_tree.cfg

echo "SBatch script: Creating the combined halo lists..."
perl halo_trees_to_catalog.pl ${out_dir}outputs/merger_tree.cfg

echo "SBatch script: Creating directory structure..."

cd ${out_dir}

mkdir hlist_raw
mkdir hlist_combined

mv *.ascii hlist_raw
mv out*.list hlist_combined

echo "SBatch script: Done."
