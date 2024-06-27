#!/bin/bash
#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH -t 0-00:10 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o exc%a # File to which STDOUT will be written
#SBATCH -e excerr%a # File to which STDERR will be written
#SBATCH --job-name=EXC
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH --array=1-1%1 # 1 task, 1 maximum task.

## See https://slurm.schedmd.com/sbatch.html for SBATCH docs.

##################
### Run script ###
##################

# Runs deeptools to generate heatmaps of coverage at specific sites.

# Example usage:
# /oasis/projects/nsf/ddp360/collab/scripts/deeptools._heatmap_at_bed.sh 
# 	       -b /oasis/projects/nsf/ddp360/collab/anno/mm10/refTSS_v3.1_mouse_coordinate.mm10.bed 
# 	       -o OPTION_PREFIX -w MB018_filterdup.pileup.bw -O OUTPUT_DIR

source $HOME/.bashrc

BED=""
PREFIX=""
BW=""
OUTPUT_DIR=$(pwd)

# Parse arguments
while [ "$1" != "" ]; do
	case $1 in
		-b | --bed )		shift
					BED=$1
					shift
					;;
		-o | --output_prefix )	shift
					PREFIX=$1
					shift
					;;
		-w | --bw )		shift
					BW=$1
					shift
					;;
		-O | --output_dir )	shift 
					OUTPUT_DIR=$1
					shift
					;;
	esac
done
if [$PREFIX = ""]
then
	PREFIX=$(basename $BW _treat_pileup.bw)
fi
echo "bed file: ${BED}"
echo "prefix: ${PREFIX}" 
echo "bw files: $BW"

START=$SECONDS

conda activate deeptools


OUT=$OUTPUT_DIR/$PREFIX

if [ "$BED" != "" ]; then
# Plot matrix
computeMatrix reference-point\
	--referencePoint TSS\
	-a 2000 -b 2000\
	-R $BED\
	-S $BW \
	-p max\
	-bs 50\
	--missingDataAsZero\
	-o ${OUT}_mat.gz &&

plotHeatmap -m ${OUT}_mat.gz\
	-o ${OUT}_peaks_sorted.svg\
	--heatmapHeight 100\
	--boxAroundHeatmaps no\
	--outFileSortedRegions ${OUT}_peaks_sorted.bed\
	--sortUsing max 

fi


conda deactivate

DIFF=$(($SECONDS - $START))
echo "deeptools heatmap took $DIFF seconds."
