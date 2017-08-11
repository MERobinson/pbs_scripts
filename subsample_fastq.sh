#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=""
check="yes"

# help message
help_message="
usage:
    bash $(basename "$0") [-options] -f1 <FASTQ1> -s <NUMBER>
purpose:
    # quick bash script to subsample from a fastq (or pair of)
required arguments:
    -f1|--fastq1 : input FASTQ file
    -s|--sub_n : amount of reads to sub (treated as fraction if < 1)
optional arguments:
    -f2|--fastq2 : mate FASTQ if PE (default = NULL)
    -n|--name : prefix for output files [string] (default = extracted from bam)
    -o|--outdir : output directory [string] (default = '.')
    -l|--logdir : output dir for log files [string] (default = --outdir)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list of PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # example depenency list: 'afterok:123456,afterok:123457'
output:
    # outputs subsampled FASTQ to --outdir & log file to --logdir
"

# parse arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -f1|--fastq1)
            fq1=$2
            shift
            ;;
        -f2|--fastq2)
            fq2=$2
            shift
            ;;
        -s|--sub_n)
            sub_n=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -l|--logdir)
            logdir=$2
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
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check whether mate provided
if [[ -z ${fq2:-} ]]; then
    echo "Running in SE mode"
else
    echo "Running in PE mode"
fi

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arg
if [[ -z ${fq1:-} ]]; then
	printf "\nERROR: no --fastq1 argument provided\n"
	echo "$help_message"; exit 1
elif [[ -z ${sub_n:-} ]]; then
    printf "\nERROR: no --sub_n argument provided\n"
    echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$fq1 ]]; then
        printf "\nERROR: FASTQ file is not readable: %s/%s\n" $workdir $fq1
        echo "$help_message"; exit 1
    fi
    if [[ ! -z ${fq2:-} ]] & [[ ! -r $workdir/$fq2 ]]; then
        printf "\nERROR: FASTQ file is not readable: %s/%s\n" $workdir $fq2
        echo "$help_message"; exit 1
    fi
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
	name=$(basename ${fq1})
    name=${name%%.*}
fi

# get basenames/prefix
fq1_base=$(basename "$fq1")
if [[ ! -z ${fq2:-} ]]; then fq2_base=$(basename "$fq2"); fi

# setup output dir
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# run subsampling
jobid=$(cat <<- EOS | qsub -N $name.subsample -
		#!/bin/bash
		#PBS -l walltime=02:00:00
		#PBS -l select=1:mem=20gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.subsample_fastq.log
		${depend:-}

		# copy files to scratch
		cp $workdir/$fq1 .
		cp $workdir/$fq2 .

		# sub
		paste $fq1_base $fq2_base |\
			awk '{ printf("%s",\$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' |\
			shuf |\
			head -n $sub_n |\
			sed 's/\t\t/\n/g' |\
			awk -F '\t' '{ print \$1 > "${name}.R1.fq"; print \$2 > "${name}.R2.fq" }'

		cp ${name}.R1.fq $workdir/$outdir
		cp ${name}.R2.fq $workdir/$outdir
		
        ls -lhAR
	EOS
)
echo "JOBID: $jobid"
