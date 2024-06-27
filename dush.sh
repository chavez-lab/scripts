#!/bin/bash

#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 0-04:00 # Runtime in D-HH:MM
#SBATCH --mem=32G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-du-sh.o%j # File to which STDOUT will be written
#SBATCH --job-name=du-sh
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
### Set this to the working directory

# Submit this sbatch job with
# sbatch this_script.sh

#################################
# Run your code below this line #
#################################

cd /expanse/lustre/projects/csd677/
du -sh

