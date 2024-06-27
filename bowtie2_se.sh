#!/bin/bash

#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH -t 0-00:45 # Runtime in D-HH:MM
#SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o exc%a # File to which STDOUT will be written
#SBATCH -e excerr%a # File to which STDERR will be written
#SBATCH --job-name=EXC
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH --array=1-1%1 # 1 task, 1 maximum task.

## See https://slurm.schedmd.com/sbatch.html for SBATCH docs.

# Should be able to run parallel on 23 threads.

####################################################

# Usage:
# bowtie2.sh <fastq.fastq> <dir/of/bt2/indices> [prefix]
# sbatch $lukas_home/scripts/bowtie2_se.sh \
#	../trimmomatic/CL100145633_L2_AE19071898-571_1.fastq \
#	$lukas_home/anno/mm10/mm10


FASTQ=$1
IDX=$2
PREFIX=${3:-$(basename "${FASTQ%.*}" )} # Set prefix to $3 if it exists, otherwise root of $1
echo $PREFIX

bam=$PREFIX.bam

module load bowtie2
module load samtools


time ( bowtie2 --mm --threads 23 -x $IDX -U $FASTQ | \
	samtools view -Su - | samtools sort -o $bam - \
)
