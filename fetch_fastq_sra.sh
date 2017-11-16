#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
sradir='~/sradata/sra'
date=`date '+%Y-%m-%d %H:%M:%S'`
scr_name=$(basename "$0")

# help message
help_message="
Wrapper to pre-fetch and extract FASTQ data using SRA-tools.

usage:
    bash $scr_name [-options] -s <SRA_accession>
required arguments:
    -s|--sra : SRA accession number
optional arguments:
    -n|--name : name prefix for output files (default = SRA accession)
    -o|--outdir : output directory (default = PWD)
    -ld|--logdir : output directory for log files (default = --outdir)
    -fq|--fqdir : output directory for fastq files (defalt = --outdir)
    -sd|--sradir : directory setup for sra toolkit (default = ~/sradata/sra)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    # all paths (except sradir) should be given relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified
    # make sure sradir points to the directory setup for SRA toolkit

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -s|--sra)
            sra=$2
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
        -ld|--logdir)
            logdir=$2
            shift
            ;;
        -fd|--fqdir)
            fqdir=$2
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

# check/parse required arguments
if [[ -z ${sra:-} ]]; then
    printf "\nERROR: --sra argument required.\n"
    echo "$help_message"; exit 1
fi

# set name if not provided
if [[ -z "${name:-}" ]]; then
    name=${sra}
fi

# set output directories
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
if [[ -z "${fqdir:-}" ]]; then
    fqdir=$outdir
fi

# create outdir if needed
mkdir -p $workdir/$fqdir
mkdir -p $workdir/$logdir

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$workdir/$logdir/$name.$scr_name.out.log

# job
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=02:00:00
		#PBS -l select=1:mem=10gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.sra
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.fetch_fastq_sra.log

		echo START: \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log
	
		# load modules
		module load sra-toolkit/2.8.1-3

		# prefetch
		if [[ ! -e $workdir/$sradir/$sra.sra ]]; then
			prefetch --max-size 100G $sra >> $out_log
		fi

		# dump
		cp $sradir/$sra.sra .
		fastq-dump -v -I -W -B --skip-technical --outdir . \
			--gzip --split-files $sra.sra >> $out_log

		# copy to outdir
		cp $sra*.fastq.gz $workdir/$outdir/

		echo END: \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
