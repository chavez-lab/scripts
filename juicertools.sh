#!/bin/bash

#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o sbatch-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=juicertools_eigenvalue
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent

export JUICERTOOLS=/oasis/projects/nsf/ddp360/collab/bin/juicer_tools_1.19.02.jar

java -jar $JUICERTOOLS $@
