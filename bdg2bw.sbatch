#!/bin/bash
#SBATCH --account=csd677 # DO NOT CHANGE
#SBATCH --partition=shared
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -t 0-01:00 # Runtime in D-HH:MM
#SBATCH --mem=8G # Memory pool for all cores (see also --mem-per-cpu) SHOULD NOT BE USED
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=bdg2bw
##SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH -D .
## See https://slurm.schedmd.com/sbatch.html for SBATCH docs.

##############################

source $HOME/.bash_profile
module load bedtools2

# check commands: slopBed, bedGraphToBigWig and bedClip
# /expanse/lustre/projects/csd677/collab/bin/UCSC_GenomeBrowser_Tools

which bedtools &>/dev/null || { echo "bedtools not found! Add to PATH: /expanse/lustre/projects/csd677/collab/bin/UCSC_GenomeBrowser_Tools"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }

# end of checking

if [ $# -lt 2 ];then
    echo "Need 2 parameters! <bedgraph> <chrom info>"
    exit
fi

F=$1
G=$2

PREFIX=$(basename "${F%.*}" )

bedtools slop -i ${F} -g ${G} -b 0 | bedClip stdin ${G} ${PREFIX}.clip

LC_COLLATE=C sort -k1,1 -k2,2n ${PREFIX}.clip > ${PREFIX}.sort.clip

bedGraphToBigWig ${PREFIX}.sort.clip ${G} ${PREFIX}.bw

rm -f ${PREFIX}.clip ${PREFIX}.sort.clip
