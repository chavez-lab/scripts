#!/bin/bash

#SBATCH --nodes=1 # All cores on one machine
#SBATCH --account=csd677 # Lukas lab allocation
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-24:00 # Runtime in D-HH:MM
#SBATCH --mem=12G # Memory pool for all cores
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=test
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
##SBATCH --array=1-1%1
### Set this to the working directory
#SBATCH -D .

# Submit this sbatch job with
# sbatch this_script.sh

#################################
# Run your code below this line #
#################################

# Load a module (see all with module avail
module load cpu/0.15.4  gcc/9.2.0  bwa/0.7.17

# Load your conda env

bwa mem -t $SLURM_CPUS_PER_TASK $@
