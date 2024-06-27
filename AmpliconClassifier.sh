#!/bin/bash

#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH --mem=32G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-AmpliconClassifier-.o%j # File to which STDOUT will be written
##SBATCH -e excerr%a # File to which STDERR will be written
#SBATCH --job-name=AmpliconClassifier
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH -D .
##SBATCH --array=1-1%1
### Set this to the working directory

# Submit this sbatch job with
# sbatch AmpliconClassifier.sh --ref [hg19, GRCh37, or GRCh38] \
#	--cycles [/path/to/amplicon_cycles.txt] --graph [/path/to/amplicon_graph.txt] > \
#	classifier_stdout.log

#################################
# Run your code below this line #
#################################

dir=${1} # should contain _amplicon[X]_edges.txt

source /home/ssridhar/miniconda3/etc/profile.d/conda.sh
conda activate ampliconclassifier

AA_dir=/expanse/lustre/projects/csd677/collab/bin/AmpliconArchitect
AA_DATA_REPO=${AA_dir}/data_repo
export AA_DATA_REPO

python /expanse/lustre/projects/csd677/collab/scripts/AmpliconClassifier.py $@
