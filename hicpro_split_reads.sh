#!/bin/bash
#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-08:00 # Runtime in D-HH:MM
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
##SBATCH -e excerr%a # File to which STDERR will be written
#SBATCH --job-name=split_reads.sh
#SBATCH -D .
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
##SBATCH --array=1-1%1 # 1 task, 1 maximum task. 

# USAGE
# hicpro_split_reads.sh [--results_folder dir] [--nreads nreads] file.fastq
## See https://slurm.schedmd.com/sbatch.html for SBATCH docs.

##############################

source ~/miniconda3/etc/profile.d/conda.sh
conda activate hic-pro
script=/oasis/projects/nsf/ddp360/collab/bin/HiC-Pro_2.11.4/bin/utils/split_reads.py

if [ $# -eq 0 ]; then
	$script -h
else
	time ( $script $@ )
fi

