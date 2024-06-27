#!/bin/bash

#SBATCH --account=csd677
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH -t 0-10:00 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=icgc-download
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
##SBATCH --array=1-1%1

# Submit this sbatch job with
# sbatch this_script.sh

# Usage: sbatch icgc-download.sh [uid] [sample_id]
# output to $lukas_home/data/icgc.
#################################
# Run your code below this line #
#################################

source ~/.bashrc
SCORE_DOWNLOAD_CLIENT=/expanse/lustre/projects/csd677/collab/bin/score-client-5.0.0/bin/score-client
OUT_DIR=/expanse/lustre/projects/csd677/collab/data/icgc
ICGC_UID=${1}
SAMPLE_ID=${2}
if [ ! -z ${SAMPLE_ID} ]; then
	OUT_DIR=${OUT_DIR}/${SAMPLE_ID}
	mkdir ${OUT_DIR}
fi

echo "Downloading ${ICGC_UID} as ${SAMPLE_ID}.bam..."

time ( \
${SCORE_DOWNLOAD_CLIENT} download --object-id ${ICGC_UID} \
	--output-dir ${OUT_DIR} \
)

# TODO rename file
if [ ! -z ${SAMPLE_ID} ]; then
	mv ${OUT_DIR}/*.bam ${OUT_DIR}/${SAMPLE_ID}.bam
	mv ${OUT_DIR}/*.bam.bai ${OUT_DIR}/${SAMPLE_ID}.bam.bai
fi

