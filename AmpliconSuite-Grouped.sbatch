#!/bin/bash

## SBATCH settings for pines
#SBATCH --partition=normal
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH --mem=128G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=AmpliconSuite-Grouped
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH -D . # Set this to the working directory

# Run AmpliconArchitect in grouped samples mode. See https://github.com/AmpliconSuite/AmpliconSuite-pipeline?tab=readme-ov-file#--grouped-analysis-of-related-samples-groupedanalysisampsuitepy
# Use:
# sbatch AmpliconSuite-Grouped.sbatch [input_file.txt] [output_directory]

#################################
# Run your code below this line #
#################################

# Required on expanse, not on pines
#source ~/.bashrc
#module load singularitypro

LUKAS_HOME=/mnt/beegfs/shares/chavez_lab/expanse
AA=$LUKAS_HOME/bin/AmpliconArchitect/AmpliconArchitect/src/AmpliconArchitect.py
AAG=$LUKAS_HOME/bin/AmpliconSuite-pipeline/GroupedAnalysisAmpSuite.py
AA_DATA_REPO=$LUKAS_HOME/bin/AmpliconArchitect/data_repo
export AA_DATA_REPO
AA_SRC=/home/programs/AmpliconArchitect-master/src # path to AA in Docker container
export AA_SRC
AC_SRC=/home/programs/AmpliconClassifier-main/src # path to AC in Docker container
export AC_SRC


input_file=$1
default_out=$(dirname $input_file)
output_dir=${2:-$default_out}
reference=${3:-'hg38'}
threads=${SLURM_CPUS_PER_TASK:-1}

echo "Running AA in grouped analysis mode on $threads threads on samples in $input_file. Output to $output_dir."

singularity exec \
	-B $LUKAS_HOME \
	docker://jluebeck/prepareaa:latest \
$AAG --input ${input_file} --out ${output_dir} --nthreads ${threads} \
	--ref ${reference} \
	--cngain 5

