#!/bin/bash

#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH --mem 128G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=bowtie_pe.sh
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH -D .

# Usage:
# bowtie2.sh <fastq_r1> fastq_r2> <dir/of/bt2/indices> [prefix]


############################
# Begin script
############################

R1=$1
R2=$2
IDX=$3
PREFIX=${4:-$(basename "${R1%.*}" )} # Set prefix to $4 if it exists, otherwise root of $1
echo "Write to ${PREFIX}.bam"
BAM=${PREFIX}.bam

module load cpu/0.15.4  gcc/9.2.0 bowtie2/2.4.1 samtools

#source $HOME/.bashrc

mkdir -p tmp
bowtie2 --local --time --mm --threads $SLURM_CPUS_PER_TASK -x $IDX --phred33 \
	-1 ${R1} -2 ${R2} | \
	samtools view -u | samtools sort -T tmp -o ${BAM} -
echo "Done!"
