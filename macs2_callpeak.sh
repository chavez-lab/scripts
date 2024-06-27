#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-08:00 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=macs2
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH -D .
## See https://slurm.schedmd.com/sbatch.html for SBATCH docs.

######################################

#module load macs2
#source ~/.bashrc
source ~/miniconda3/etc/profile.d/conda.sh
conda activate macs2
time ( macs2 callpeak $@ )
conda deactivate
