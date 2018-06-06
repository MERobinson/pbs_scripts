#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default args
workdir=$PWD
outdir=''
check='yes'

# help message
help_message="
Wrapper to run FASTQC.

usage:
    $(basename "$0") [-options] -f <FASTQ>
required arguments:
    -f|--fastq : FASTQ files to analyse, comma separated if multiple
optional arguments:
    -o|--outdir : output directory for metric files (default = pwd)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling

"

# parse arguments
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -f|--fastq)
            fastq=$2
            shift
            ;;
		-o|--outdir)
            outdir=$2
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        *)
            printf "ERROR: Undefined argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
	shift
done

# check outdir exists
mkdir -p $workdir/$outdir

# check fastq and split
if [[ -z ${fastq:-} ]]; then
	printf "\nERROR: --fastq argument required\n"
	echo "$help_message"; exit 1
else
    IFS="," read -r -a fastq_array <<< $fastq
fi

# check readable unless flagged off
if [[ ${check:-} = yes ]]; then
    for fastq in ${fastq_array[@]}; do
        if [[ ! -r ${workdir}/${fastq} ]]; then
            printf "\nERROR: FASTQ is not readable: %s/%s\n" $workdir $bam
            echo "$help_message"; exit 1
        fi
    done
fi

# set arguments
for fastq in ${fastq_array[@]}; do
    fastq_base=$(basename $fastq)
    cp_arg="${fastq_arg:-}cp $workdir/$fastq .;"
    fq_arg="${fq_arg:-}$fastq_base "
done 

scr_name=$(basename "$0" .sh)

# run job
jobid=$(cat <<- EOS | qsub -N $scr_name -
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:ncpus=8:mem=8gb
		#PBS -j oe
		#PBS -q med-bio
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# load modules
		module load fastqc/0.11.2
	
		# copy to scratch
		$cp_arg

		# run
		fastqc --noextract -t 8 $fq_arg

		# copy metrics to outdir
		cp *fastqc.zip $workdir/$outdir/
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
	EOS
)
# echo job id and exit
echo "JOBID: $jobid"
exit 0
