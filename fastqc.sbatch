#!/bin/bash
#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH -t 0-04:00 # Runtime in D-HH:MM
#SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
##SBATCH -e excerr%a # File to which STDERR will be written
#SBATCH --job-name=fastqc.sh
#SBATCH -D .
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
##SBATCH --array=1-1%1 # 1 task, 1 maximum task. 

## See https://slurm.schedmd.com/sbatch.html for SBATCH docs.

# Not tested, but you should be able to run up to 23 files simultaneously with -t 23.

##############################

module load cpu/0.15.4  gcc/10.2.0 fastqc

if [ $# -eq 0 ]; then
	fastqc -h
else
	time ( fastqc $@ )
fi

