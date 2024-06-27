#!/bin/bash

#SBATCH --account=csd677
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=4
#SBATCH -t 0-01:00 # Runtime in D-HH:MM
##SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=fingerprint
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
##SBATCH --array=1-1%1

# Submit this sbatch job with
# sbatch this_script.sh

# Usage: sbatch fingerprint.sh [bam_file]
# output to current directory.
#################################
# Run your code below this line #
#################################

#source ~/.bashrc
source ~/miniconda3/etc/profile.d/conda.sh
conda activate fingerprint
LUKAS_HOME=/expanse/lustre/projects/csd677/collab

input_bam=$1
snp_bed=$LUKAS_HOME/bin/fingerprint/snp138Common.n1000.vh20140318.bed

echo "Fingerprinting $input_bam."

SAMPLENAME=$(basename $input_bam .bam)
$LUKAS_HOME/bin/fingerprint/bsnp.py $snp_bed $input_bam \
	> ${SAMPLENAME}.fp
#	> ${LUKAS_HOME}/projects/fingerprints/medulloblastoma/${SAMPLENAME}.fp

conda deactivate
