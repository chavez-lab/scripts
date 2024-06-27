#!/bin/bash

#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-00:10 # Runtime in D-HH:MM
##SBATCH --mem=124000 # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o exc # File to which STDOUT will be written
#SBATCH -e excerr # File to which STDERR will be written
#SBATCH --job-name=EXC
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
### Set this to the working directory

# Submit this sbatch job with
# sbatch this_script.sh

#################################
# Run your code below this line #
#################################

DIR=$1
rm -r $DIR
