#! /bin/bash

#SBATCH --partition=normal
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 0-2:00 # Runtime in D-HH:MM
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
##SBATCH --mem-per-cpu=8G
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=rsem-calculate-expression
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
##SBATCH --array=1-1%1


#################################
# usage: sbatch rsem-prepare-reference.sbatch <BAM>
# assumes rsem is installed via mamba
source ~/mambaforge/etc/profile.d/conda.sh
conda activate rsem

INPUTS=$@
rsem-calculate-expression --num-threads $SLURM_CPUS_PER_TASK $INPUTS
