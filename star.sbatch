#!/bin/bash

#SBATCH --nodes=1 # All cores on one machine
#SBATCH --account=csd677 # Lukas lab allocation
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 0-4:00 # Runtime in D-HH:MM
#SBATCH --mem=64G # Memory pool for all cores
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=star
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH -D .

# See a nice tutorial at 
# https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
#################################
# Run your code below this line #
#################################

module load cpu/0.15.4  gcc/10.2.0 star/2.7.9a

STAR $@

## Usage
## Example using ENCODE standard args
## https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

#STAR --runMode genomeGenerate \
#	--runThreadN 8 \
#	--genomeDir $genomeDir \
#	--outSAMtype BAM SortedByCoordinate \
#	--readFilesIn \
#	--outFilterType BySJout \
#	--outFilterMultimapNmax 20 \
#	--alignSJoverhangMin 8 \
#	--alignSJDBoverhangMin 1 \
#	--outFilterMismatchNmax 999 \
#	--outFilterMismatchNoverReadLmax 0.04 \
#	--alignIntronMin 20 \
#	--alignIntronMax 1000000 \
#	--alignMatesGapMax 1000000
