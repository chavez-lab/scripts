#!/bin/bash

#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH -t 0-01:00 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o logs/out%a # File to which STDOUT will be written
#SBATCH -e logs/err%a # File to which STDERR will be written
#SBATCH --job-name=EXC
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=yil031@ucsd.edu # Email to which notifications will be sent
#SBATCH --array=1-1%1 # 1 task, 1 maximum task.

## See https://slurm.schedmd.com/sbatch.html for SBATCH docs.

# Can use --threads up to 23. 
# Usage: sbatch $lukas_home/collab/scripts/encode_ChIP_se.sh $BAM_FILE

##############################

module load bedtools
module load samtools
module load picard

#set -o pipefail
#function fail {
#	echo "$@" >&2
#	exit 1
#}

# =============================
# Remove unmapped, mate unmapped not primary alignment, reads failing platform
# ==================

RAW_BAM_FILE=$1
OUTPUT_DIR=$(dirname $1)
RAW_BAM_FILE_NAME=$(basename $1 .bam)

cd $OUTPUT_DIR
echo $OUTPUT_DIR
echo $PWD

FILT_BAM_FILE="${RAW_BAM_FILE_NAME}.filt.srt.bam"
MAPQ_THRESH=30

samtools view -F 1804 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} -o ${FILT_BAM_FILE} #\ 
#	|| fail "samtools view failed. Did not generate ${FILT_BAM_FILE}"

# ========================
# Mark duplicates
# ======================

TMP_FILT_BAM_FILE="${RAW_BAM_FILE_NAME}.dupmark.bam"
MARKDUP="${PICARDHOME}/picard.jar MarkDuplicates"
DUP_FILE_QC="${RAW_BAM_FILE_NAME}.dup.qc" # QC file


echo $FILT_BAM_FILE
echo $TMP_FILT_BAM_FILE
echo $MARKDUP
echo $DUP_FILE_QC

java -Xmx4G -jar ${MARKDUP} I=${FILT_BAM_FILE} O=${TMP_FILT_BAM_FILE} M=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false  #\ 
#	|| fail "picard MarkDuplicates failed. I=${FILT_BAM_FILE} O=${TMP_FILT_BAM_FILE}"

mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

# ============================
# Remove duplicates
# Index final position sorted BAM
#=============================
FINAL_BAM_PREFIX="${RAW_BAM_FILE_NAME}.filt.nodup.srt"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

# Final bam contains only reads with flags 1804 unset.
samtools view -F 1804 -b ${FILT_BAM_FILE} -o ${FINAL_BAM_FILE} 
#|| fail "Samtools view failed. Did not generate ${FINAL_BAM_FILE}"

# Index Final BAM file
# Need to sort first, then index
# Possible theyre already sorted, try without
#samtools sort -n --threads 23 ${FINAL_BAM_FILE}  
samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE} 
#|| fail "samtools index the final bam failed."

# =============================
# Compute library complexity
# =============================
# sort by position and strand
# Obtain unique count statistics

PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

# PBC File output
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
bedtools bamtobed -i ${FILT_BAM_FILE} | \
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | \
	grep -v 'chrM\|MT' | sort | uniq -c | \
	awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
#|| fail "PBC calculations failed."

rm ${FILT_BAM_FILE}

