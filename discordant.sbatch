#!/bin/bash

#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-2:00 # Runtime in D-HH:MM
#SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=discordant
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
### Set this to the working directory
#SBATCH -D .

# Submit this sbatch job with
# sbatch this_script.sh

#######################

module load cpu/0.15.4  gcc/9.2.0 samtools/1.10

BAMFILE=$1
OUTFILE=$(basename $BAMFILE .bam)
OUTFILE=${OUTFILE}_discordant.bam
samtools view -F 1806 -o ${OUTFILE} ${BAMFILE}
samtools index ${OUTFILE}
echo "Done"
