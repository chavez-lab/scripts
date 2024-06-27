#!/bin/bash

#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # NOT MULTITHREADED
#SBATCH -t 0-12:00 # Runtime in D-HH:MM
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=hicpro2juicebox
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent

# Usage:
# sbatch hicpro2juicebox.sh -i <allValidPairs> -g <genome.sizes>


############################
# Begin script
############################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate hic-pro

script=/expanse/lustre/projects/csd677/collab/bin/HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh
j=/expanse/lustre/projects/csd677/collab/bin/juicer_tools_1.22.01.jar

if [ $# -eq 0 ]; then
        $script -h
else
        time ( $script -j $j $@ )
fi

echo "Done!"
