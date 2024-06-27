#!/bin/bash

#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-2:00 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o sbatch-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=AA_seed_from_coverage
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH --array=1-1%1
### Set this to the working directory

BAM=$1
COV=$2
REF=$3

SAMPLE=$(basename $BAM .bam)
echo "Getting AA seed file ${SAMPLE}_AA_seed.bdg from input ${BAM} at coverage depth ${COV}"

module load bedtools

# coverage across genome
bedtools genomecov -bg -ibam ${BAM} > ${SAMPLE}.bdg

# coverage greater than threshold
awk -v cov="${COV}" '($4 > cov)' ${SAMPLE}.bdg > temp.bdg

# remove mitochondrial
awk '($1 != "chrM")' temp.bdg > ${SAMPLE}_coverage_over_${COV}.bdg

# merge
#bedtools merge -i ${SAMPLE}_coverage_over_${COV}.bdg -d 1000 > ${SAMPLE}_AA_seed.bdg

# Perform AA preprocessing - merge, min length, etc filters
awk 'BEGIN{OFS="\t"}; {print $1,$2,$3,"*",$4}' ${SAMPLE}_coverage_over_${COV}.bdg \
	> ${SAMPLE}_coverage_over_${COV}.bed
conda activate py27 # contains python2.7 mosek scipy numpy pysam matplotlib
AA_DIR=/oasis/projects/nsf/ddp360/collab/bin/AmpliconArchitect
AA_DATA_REPO=$AA_DIR/data_repo
export AA_DATA_REPO
python2 $AA_DIR/AmpliconArchitect-jluebeck/src/amplified_intervals.py \
	--bed ${SAMPLE}_coverage_over_${COV}.bed \
	--bam ${BAM} \
	--ref ${REF} \
	--gain ${COV} \
	--out ${SAMPLE}_AA_seed.bed

rm temp.bdg

# get discordant reads because we'll use them later
module load samtools
samtools view -F 1806 -b -o ${SAMPLE}_discordant.bam ${BAM}

echo "Done!"
