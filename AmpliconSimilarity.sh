≈#!/bin/bash

#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-06:00 # Runtime in D-HH:MM
#SBATCH --mem=32G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-AmpliconClassifier-.o%j # File to which STDOUT will be written
##SBATCH -e excerr%a # File to which STDERR will be written
#SBATCH --job-name=AmpliconClassifier
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH -D .
##SBATCH --array=1-1%1
### Set this to the working directory

# Performs environment setup and sbatch args for running amplicon_classifier.
# Submit this sbatch job with
# sbatch AmpliconSimilarity.sh [params passed to amplicon_classifier.py]

## Copypasta from AC parames:
#usage: amplicon_classifier.py [-h] [-i INPUT] [-c CYCLES] [-g GRAPH] --ref {hg19,GRCh37,hg38,GRCh38,mm10,GRCm38}
#                              [--min_cn_flow MIN_CN_FLOW] [--min_size MIN_SIZE] [-o O]
#                              [--plotstyle {grouped,individual,noplot}] [--force] [--add_chr_tag]
#                              [--report_complexity] [--verbose_classification] [--annotate_cycles_file]
#                              [--no_LC_filter] [--decomposition_strictness DECOMPOSITION_STRICTNESS] [-v]

#################################
# Run your code below this line #
#################################

source /home/ssridhar/miniconda3/etc/profile.d/conda.sh
conda activate ampliconclassifier

AA_dir=/expanse/lustre/projects/csd677/collab/bin/AmpliconArchitect
AA_DATA_REPO=${AA_dir}/data_repo
export AA_DATA_REPO

python /expanse/lustre/projects/csd677/collab/bin/AmpliconClassifier/amplicon_similarity.py $@
