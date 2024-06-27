#! /bin/bash
#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 0-2:00 # Runtime in D-HH:MM
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
##SBATCH --mem-per-cpu=8G
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=PrepareAA
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH -D .


#################################
# usage: PrepareAA-cnvkit.sh -s [sample_name] --sorted_bam [bamfile] --ref [hg19,GRCh37,GRCh38]

source $HOME/.bashrc
source ~/miniconda3/etc/profile.d/conda.sh
conda activate py2

sample=${1}

AA_DIR=/expanse/lustre/projects/csd677/collab/bin/AmpliconArchitect

export AA_SRC=$AA_DIR/AmpliconArchitect/src
export AA_DATA_REPO=$AA_DIR/data_repo

python2 /expanse/lustre/projects/csd677/collab/bin/PrepareAA/PrepareAA.py \
	-t 16 \
        --cngain 4 \
        --cnvkit_dir /home/ochapman/miniconda3/envs/cnvkit/bin \
	--rscript_path /home/ochapman/miniconda3/envs/cnvkit/bin/Rscript \
	--python3_path /home/ochapman/miniconda3/envs/cnvkit/bin \
	$@

conda deactivate
