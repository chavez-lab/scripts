#! /bin/bash
#SBATCH --account=ddp360 # DO NOT CHANGE
#SBATCH --partition=gpu-shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks-per-node=6
#SBATCH --gres=gpu:k80:1
#SBATCH -t 0-4:00 # Runtime in D-HH:MM
#SBATCH --mem=16G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT B$##SBATCH --mem-per-cpu=8G
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=hiccups
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent

# Environment
source ~/.bash_profile # need $JUICERTOOLS
module load cuda/7.5

# Params
hic_file=$1
default_out=$(pwd)
out_dir=${2:-$default_out}


time ( \
java -Xms512m -Xmx14g \
	-jar $JUICERTOOLS hiccups \
	--threads 0 --ignore-sparsity -k KR -f .2 \
	-t 0.02,1.5,1.75,2 -d 20000,20000,50000 -r 5000,10000,25000 -p 4,2,1 -i 7,5,3 \
	$hic_file $out_dir ) \
> $out_dir/hiccups.log 2>&1	
