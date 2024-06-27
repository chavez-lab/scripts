#!/bin/bash
#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-08:00 # Runtime in D-HH:MM
##SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-md5sums-.%j # File to which STDOUT will be written
##SBATCH -e excerr%a # File to which STDERR will be written
#SBATCH --job-name=md5sums
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent

# Takes an input file of the format
# md5sum file
# and checks each md5sum against the record md5sum

# Usage: check_md5sums_against_file.sh <file>

INPUT=$1

cat $INPUT | while read line || [[ -n $line ]];
do
	line=($line)
	md5=${line[0]}
	filename=${line[1]}
	echo "Checking $filename ..."
	md5tocheck=($(md5sum $filename))
	if [ $md5 == $md5tocheck ] 
	then
		echo "md5sums checked."
	else
		echo "md5sums do not match:"
		echo "file md5sum: $md5"
		echo "record md5sum: $md5tocheck"
	fi
done
