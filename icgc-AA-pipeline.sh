#!/bin/bash

#SBATCH --partition=shared
#SBATCH --account=ddp360
#SBATCH --nodes=1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --time=0-04:00 # Runtime in D-HH:MM
##SBATCH --mem=1G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH --output=slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=icgc-AA-pipeline
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent

# Submit this sbatch job with
# sbatch this_script.sh

# Usage: sbatch icgc_AA_pipeline.sh [uid] [sample_id]
# Bam file to $lukas_home/data/icgc
# Analysis to $lukas_home/projects/ecDNA/[sample_id]
#################################
# Run your code below this line #
#################################

source ~/.bashrc
ICGC_UID=${1}
SAMPLE_ID=${2}

LUKAS_HOME=/oasis/projects/nsf/ddp360/collab
SCORE_DOWNLOAD_CLIENT=${LUKAS_HOME}/bin/score-client-5.0.0/bin/score-client
DATA_DIR=${LUKAS_HOME}/data/icgc
ANALYSIS_DIR=${LUKAS_HOME}/projects/ecDNA/AmpliconArchitect/${SAMPLE_ID}

function get_job_id {
	if [[ "$1" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
	echo "${BASH_REMATCH[1]}"
else
	echo "Sbatch failed"
	exit 1
fi
}

# Download
j1="$(sbatch --output=slurm-${SAMPLE_ID}.o%j -D ${DATA_DIR} \
$LUKAS_HOME/scripts/icgc-download.sh $ICGC_UID $SAMPLE_ID \
| grep -o '[[:digit:]]*')"	
#j1=$(get_job_id "$j1")
echo $j1

BAMFILE=${DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.bam
mkdir ${ANALYSIS_DIR}

# PrepareAA
j2="$(sbatch --dependency=afterok:${j1} \
	--output=slurm-PAA-${SAMPLE_ID}.o%j \
	-D ${ANALYSIS_DIR} \
$LUKAS_HOME/scripts/PrepareAA-cnvkit.sh \
	-s $SAMPLE_ID \
	--sorted_bam ${BAMFILE} \
	--ref GRCh37 \
| grep -o '[[:digit:]]*')"
#j2=$(get_job_id "$j2")
echo $j2

# Fingerprint
j3="$(sbatch --dependency=afterok:${j1} \
	--output=slurm-${SAMPLE_ID}.o%j \
	-D ${LUKAS_HOME}/projects/fingerprints/medulloblastoma \
$LUKAS_HOME/scripts/fingerprint.sh \
	${BAMFILE} \
| grep -o '[[:digit:]]*')"
echo $j3

# AmpliconArchitect
j4="$(sbatch --dependency=afterok:${j2} \
	--output=slurm-AA-${SAMPLE_ID}.o%j \
	-D ${ANALYSIS_DIR} \
${LUKAS_HOME}/scripts/AmpliconArchitect.sh \
	${BAMFILE} \
	${ANALYSIS_DIR}/${SAMPLE_ID}_AA_CNV_SEEDS.bed \
	GRCh37 \
| grep -o '[[:digit:]]*')"
echo $j4

# Delete bam file
#j5="$(sbatch --dependency=afterok:${j4}:${j3} $LUKAS_HOME/scripts/icgc-rm.sh ${DATA_DIR}/${SAMPLE_ID})"

echo "Done!"
