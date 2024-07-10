#!/bin/bash

## SBATCH settings for pines
#SBATCH --partition=normal
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -t 2-00:00 # Runtime in D-HH:MM
#SBATCH --mem=128G # Memory pool for all cores
#SBATCH -o slurm-%x-%j.out # File to which STDOUT will be written
#SBATCH --job-name=AmpliconSuite-Grouped
#SBATCH --mail-type=FAIL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=ochapman@ucsd.edu # Email to which notifications will be sent
#SBATCH -D . # Set this to the working directory

# Wrapper around AmpliconSuite-pipeline. Sets environment variables, runs the docker container and provides reasonable default values for --nthreads and --output_directory.
# See https://github.com/AmpliconSuite/AmpliconSuite-pipeline

# Usage:
# sbatch AmpliconSuite-pipeline.sbatch --sample_name test --bam test.bam
# Copy pasta params for ASP:
# AmpliconSuite-pipeline.py [-h] [-v]
#                           [--download_repo {hg19,GRCh37,GRCh38,mm10,GRCh38_viral,hg19_indexed,GRCh37_indexed,GRCh38_indexed,mm10_indexed,GRCh38_viral_indexed} [{hg19,GRCh37,GRCh38,mm10,GRCh38_viral,hg19_indexed,GRCh37_indexed,GRCh38_indexed,mm10_indexed,GRCh38_viral_indexed} ...]]
#                           [-o PATH] [-s STR] [-t INT] [--run_AA] [--run_AC] [--ref STR] [--cngain FLOAT] [--cnsize_min INT] [--downsample FLOAT] [--rscript_path PATH] [--python3_path PATH] [--aa_python_interpreter PATH] [--sv_vcf FILE] [--sv_vcf_no_filter]
#                           [--AA_src PATH] [--AA_runmode STR] [--AA_extendmode STR] [--AA_insert_sdevs FLOAT] [--normal_bam FILE] [--ploidy FLOAT] [--purity FLOAT] [--cnvkit_segmentation STR] [--no_filter] [--no_QC] [--sample_metadata FILE]
#                           [--samtools_path SAMTOOLS_PATH] [--bam FILE | --fastqs TWO FILES TWO FILES | --completed_AA_runs PATH] [--cnv_bed FILE | --cnvkit_dir PATH | --completed_run_metadata FILE | --align_only]

# Requirements
# singularity
# NB. A valid mosek license file at ~/mosek/mosek.lic is required to run AmpliconArchitect. See https://github.com/AmpliconSuite/AmpliconSuite-pipeline?tab=readme-ov-file#installation

#################################
# Run your code below this line #
#################################

LUKAS_HOME=/mnt/beegfs/shares/chavez_lab/expanse
ASP=/home/programs/AmpliconSuite-pipeline-master/AmpliconSuite-pipeline.py # path to ASP in Docker
AA_DATA_REPO=$LUKAS_HOME/bin/AmpliconArchitect/data_repo
export AA_DATA_REPO
AA_SRC=/home/programs/AmpliconArchitect-master/src # path to AA in Docker container
export AA_SRC
AC_SRC=/home/programs/AmpliconClassifier-main # path to AC in Docker container
export AC_SRC

# function to parse arguments
parse_arguments() {
  # Loop through all arguments
  while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
      -s|--sample_name)
      sample_name="$2"
      shift # Skip the argument value
      ;;
      *)  # Unknown option
      ;;
    esac
    shift # Move to the next argument
  done
}

parse_arguments "$@"

# function to check if a given flag is present.
# call as `check_flag -t --nthreads ${args[@]}`
check_flag() {
  local flag1="$1"
  local flag2="$2"
  shift 2
  for arg in "${@}"; do
    if [[ "$arg" == "$flag1" || "$arg" == "$flag2" ]]; then
      return 1 # success
    fi
  done
  return 0 # failure
}

# set number of threads if not specified
flag_to_check="--nthreads"
args=("$@") # Capture all arguments in an array

if [[ ! $(check_flag "-t" "--nthreads" "${args[@]}") ]] ; then # if --nthreds is not set
  if [[ ! -v SLURM_CPUS_PER_TASK ]]; then # if SLURM_CPUS_PER_TASK is not set
    SLURM_CPUS_PER_TASK=1 # set SLURM_CPUS_PER_TASK
  fi
  args+=("--nthreads" ${SLURM_CPUS_PER_TASK}) # set --nthreads
fi

if ! [[ $(check_flag "-o" "--output_directory" ${args[@]}) ]] ; then
  args+=("--output_directory" ${sample_name})
fi

echo "Running command:\n AmpliconSuite-pipeline.py ${args[@]}"

singularity exec \
	-B $LUKAS_HOME \
	docker://jluebeck/prepareaa:latest \
$ASP ${args[@]}
