#!/bin/bash

#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-4:00 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o sbatch-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=AA_amplified_intervals.py
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH --array=1-1%1
### Set this to the working directory

source $HOME/.bash_profile
conda activate py27 # contains python2.7 mosek scipy numpy pysam matplotlib
AA_DIR=/oasis/projects/nsf/ddp360/collab/bin/AmpliconArchitect
AA_DATA_REPO=$AA_DIR/data_repo
export AA_DATA_REPO

python2 $AA_DIR/AmpliconArchitect-jluebeck/src/amplified_intervals.py $@
