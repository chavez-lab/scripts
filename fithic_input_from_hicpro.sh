#!/bin/bash
#SBATCH --account=csd677
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH -t 0-04:00 # Runtime in D-HH:MM
##SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=fithic_input_from_hicpro
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent


source ~/miniconda3/etc/profile.d/conda.sh
conda activate hicpro
PREP_SCRIPT=/expanse/lustre/projects/csd677/collab/bin/HiC-Pro_3.0.0/bin/utils/hicpro2fithic.py

SAMPLE=$1
RES=$2

#INDIR=/oasis/projects/nsf/ddp360/collab/projects/medulloblastoma_hic/hicpro/$SAMPLE/results/hic_results/matrix/$SAMPLE
INDIR=/expanse/lustre/projects/csd677/collab/projects/medulloblastoma_hic/hicpro/$SAMPLE/hic_results/matrix/$SAMPLE
OUTDIR=/expanse/lustre/projects/csd677/collab/projects/medulloblastoma_hic/fithic/$SAMPLE
mkdir -p $OUTDIR/input_$RES

# hack hack
#SAMPLE=rawdata

time (
$PREP_SCRIPT \
	-i $INDIR/raw/$RES/${SAMPLE}_${RES}.matrix \
	-b $INDIR/raw/$RES/${SAMPLE}_${RES}_abs.bed \
	-s $INDIR/iced/$RES/${SAMPLE}_${RES}_iced.matrix.biases \
	-r $RES \
	-o $OUTDIR/input_$RES
) > $OUTDIR/fithic_input_from_hicpro.log 2>&1

conda deactivate

