#!/bin/bash
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH -t 0-12:00 # Runtime in D-HH:MM
#SBATCH --mem=32G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=fithic_input_from_hicpro
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent


source ~/miniconda3/etc/profile.d/conda.sh
conda activate fithic2

SAMPLE=$1
RES=$2

INDIR=/oasis/projects/nsf/ddp360/collab/projects/medulloblastoma_hic/fithic/$SAMPLE

# get rid of nonlocalized genome contigs
zcat $INDIR/inputs/fithic.interactionCounts.gz | awk '$1 !~ /_/' | gzip > $INDIR/inputs/fithic.interactionCounts.chrom.gz
zcat $INDIR/inputs/fithic.fragmentMappability.gz | awk '$1 !~ /_/' | gzip > $INDIR/inputs/fithic.fragmentMappability.chrom.gz

time (
fithic \
	-f $INDIR/inputs/fithic.fragmentMappability.chrom.gz \
	-i $INDIR/inputs/fithic.interactionCounts.chrom.gz \
	-o $INDIR \
	-r $RES \
	-t $INDIR/inputs/fithic.biases.gz \
	-p 2 \
	-l $SAMPLE \
	-L 100000 \
	-x All \
) > $INDIR/fithic.log 2>&1

conda deactivate

