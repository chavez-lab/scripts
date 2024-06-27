#! /bin/bash
#SBATCH --partition=normal
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -t 0-2:00 # Runtime in D-HH:MM
#SBATCH --mem=64G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
##SBATCH --mem-per-cpu=8G
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=rsem-prepare-reference
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
##SBATCH --array=1-1%1


#################################
# usage: sbatch rsem-prepare-reference.sbatch
source ~/mambaforge/etc/profile.d/conda.sh
conda activate rsem

#gtf=/oasis/projects/nsf/ddp360/collab/anno/hg38/Annotation/Genes/genes.gtf
#fasta=/oasis/projects/nsf/ddp360/collab/anno/hg38/Sequence/WholeGenomeFasta/genome.fa
#out=/oasis/projects/nsf/ddp360/collab/anno/hg38/RSEM/hg38

dir=/mnt/beegfs/shares/chavez_lab/expanse/anno/hs37d5
gtf=$dir/genes.gtf
fasta=$dir/hs37d5.fa
out=$dir/RSEM/hs37d5

mkdir -p $dir/RSEM

rsem-prepare-reference --gtf ${gtf} -p 16 --star ${fasta} ${out}
