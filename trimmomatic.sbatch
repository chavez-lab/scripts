#!/bin/bash

#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=23
#SBATCH -t 0-24:00 # Runtime in D-HH:MM
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
##SBATCH -e excerr%a # File to which STDERR will be written
#SBATCH --job-name=trimmomatic.sh
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
##SBATCH --array=1-1%1 # 1 task, 1 maximum task.
#SBATCH -D .

## See https://slurm.schedmd.com/sbatch.html for SBATCH docs.
## Trimmomatic settings from ENCODE specifications for ChIP-seq.

# Can use --threads up to 23. 
# Usage: sbatch $lukas_home/collab/scripts/trimmomatic.sh SE 
# 	-threads 23 -phred33 \
#	$lukas_home/collab/data/Bagdachi_p53/ChIPseq/CL100145633_L2_AE19071898-571_1.fq.gz \
#	./CL100145633_L2_AE19071898-571_1.fastq \
#	ILLUMINACLIP:/opt/biotools/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 \
#	LEADING:15 \
#	TRAILING:15 \
#	SLIDINGWINDOW:4:15

##############################

module load trimmomatic
# module display trimmomatic: 
#setenv           TRIMMOMATICHOME /opt/biotools/trimmomatic
#append-path      CLASSPATH /opt/biotools/trimmomatic/trimmomatic-0.36.jar

time ( trimmomatic $@ )
#time ( java -jar $CLASSPATH $@ )
