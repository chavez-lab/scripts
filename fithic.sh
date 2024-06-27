#!/bin/bash
#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH -t 0-12:00 # Runtime in D-HH:MM
#SBATCH --mem=32G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=fithic
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent


source ~/miniconda3/etc/profile.d/conda.sh
conda activate fithic

SAMPLE=$1
RES=$2

INDIR=/expanse/lustre/projects/csd677/collab/projects/medulloblastoma_hic/fithic/$SAMPLE

# get rid of contigs

## awk $1 >> look at first column. !~ RegEx match, inverted. /_/ Any field with '_'
zcat $INDIR/inputs_${RES}/fithic.interactionCounts.gz \
	| awk '$1 !~ /_|^M$|^Y$|^EBV$/' \
	| gzip > $INDIR/inputs_${RES}/fithic.interactionCounts.chrom.gz
zcat $INDIR/inputs_${RES}/fithic.fragmentMappability.gz \
	| awk '$1 !~ /_|^M$|^Y$|^EBV$/' \
	| gzip > $INDIR/inputs_${RES}/fithic.fragmentMappability.chrom.gz
zcat $INDIR/inputs_${RES}/fithic.biases.gz \
	| awk '$1 !~ /_|^M$|^Y$|^EBV$|^$/' \
	| gzip > $INDIR/inputs_${RES}/fithic.biases.chrom.gz

time (
fithic \
	-f $INDIR/inputs_${RES}/fithic.fragmentMappability.chrom.gz \
	-i $INDIR/inputs_${RES}/fithic.interactionCounts.chrom.gz \
	-o $INDIR \
	-r $RES \
	-t $INDIR/inputs_${RES}/fithic.biases.chrom.gz \
	-p 2 \
	-l $SAMPLE \
	-L 100000 \
	-x interOnly \
	--biasUpperBound 10000 \
) >> $INDIR/fithic.log 2>&1

conda deactivate
