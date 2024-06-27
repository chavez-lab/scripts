#!/bin/bash

#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=compute
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH -t 0-2:00 # Runtime in D-HH:MM
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-coverage-.o%j # File to which STDOUT will be written
#SBATCH --job-name=coverage
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
### Set this to the working directory

# Submit this sbatch job with
# sbatch this_script.sh

#######################

#source ~/.bashrc
source ~/miniconda3/etc/profile.d/conda.sh
conda activate deeptools

BAMFILE=$1
OUTFILE=$(basename $BAMFILE .bam)
OUTFILE=${OUTFILE}.bw
bamCoverage -b $BAMFILE -o $OUTFILE -p max --normalizeUsing RPKM
#bamCoverage -b $BAMFILE -o $OUTFILE -r chr8:128732756:129122669 -p max --binSize 2000 --normalizeUsing RPKM
