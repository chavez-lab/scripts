#!/bin/bash

#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH --mem=32G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-AA-%j.out # File to which STDOUT will be written
##SBATCH -e excerr%a # File to which STDERR will be written
#SBATCH --job-name=AmpliconArchitect
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
##SBATCH --array=1-1%1
### Set this to the working directory
#SBATCH -D .
# Submit this sbatch job with
# sbatch this_script.sh

#################################
# Run your code below this line #
#################################

source ~/.bashrc
module load singularitypro

AA_dir=/expanse/lustre/projects/csd677/collab/bin/AmpliconArchitect
AA=$AA_dir/AmpliconArchitect/src/AmpliconArchitect.py
AA_DATA_REPO=${AA_dir}/data_repo
export AA_DATA_REPO

input_bam=$1
input_bed=$2
default_out=$(dirname $input_bed)
output_prefix=${4:-$default_out}
reference=${3:-'hg19'}

echo "Running AA on $input_bam and $input_bed. Output to $output_prefix."
SAMPLE=$(basename $input_bam)

mkdir $output_prefix/AA_results
rm $output_prefix/AA_results/*

singularity exec \
	-B $(dirname $input_bam) \
	-B $(dirname $input_bed) \
	-B $HOME/mosek \
	-B $output_prefix/AA_results \
	-B $AA_dir \
	$AA_dir/AmpliconArchitect.simg \
$AA --bam ${input_bam} --bed ${input_bed} --out ${output_prefix}/AA_results/${SAMPLE} \
	--ref ${reference}
