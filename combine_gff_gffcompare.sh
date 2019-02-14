#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='on'
date=`date '+%Y-%m-%d %H:%M:%S'`

# help message
help_message="
Wrapper to merge/compare transcripts with gffcompare

usage:
    bash $(basename $0) [-options] --input <GFF,list>
required arguments:
    -i|--input : comma separated list of input gff files
optional arguments:
    -R|--ref : reference GTF/GFF file (default = NULL)
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -l|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on|off] (default = on)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--input)
            gtf_list=$2
            shift
            ;;
        -R|--ref)
            ref_gtf=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -l|--logdir)
            logdir=$2
            shift
            ;;
        -n|--name)
            name=$2
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
            printf "ERROR: Unrecognised argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z ${gtf_list:-} ]]; then
    printf "\nERROR: --input argument required\n"
    echo "$help_message"; exit 1
else
    IFS=',' read -r -a gtf_array <<< $gtf_list
fi

# check files 
if [[ "${check:-}" = 'on' ]]; then
    for gtf in ${gtf_array[@]}; do
        if [[ ! -r $gtf ]]; then
            printf "\nERROR: GTF cannot be read: %s/%s\n" $workdir $gtf
            echo "$help_message"; exit 1
        fi
    done
    if [[ -n "${ref_gtf:-}" ]] && [[ ! -r "${ref_gtf:-}" ]]; then
        printf "\nERROR: ref GTF cannot be read: %s/%s\n" $workdir ${ref_gtf:-}
        echo "$help_message"; exit 1
    fi
fi

# set output directories
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set basenames
for gtf in ${gtf_array[@]}; do
    gtf_cp="${gtf_cp:-}; cp $workdir/$gtf* ."
    gtf_arg="${gtf_arg:-}$(basename $gtf) "
done
gtf_cp=${gtf_cp#*;}

# set optional arguments
if [[ -n "${ref_gtf:-}" ]]; then
    ref_cp="cp $workdir/$ref_gtf* ."
    ref_base=$(basename "$ref_gtf")
    ref_arg="-r ${ref_base}"
fi

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    name=$(basename "${gtf_array[0]}")
    name="${name%%.*}.gffcmp"
fi

# set command
gffcmp_cmd=("gffcompare"
            "-V -o $name"
            "${ref_arg:-}"
            "${gtf_arg:-}")

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=20gb:ncpus=20
		#PBS -j oe
		#PBS -N gffcmp.$name
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load modules
		module load gffcompare/0.10.1
		module load samtools/1.2

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log
		
		# copy files to scratch
		${ref_cp:-} &>> $out_log
		${gtf_cp:-} &>> $out_log

		# run stringtie
		printf "\nRunning gffcompare:\n" >> $out_log 
		${gffcmp_cmd[@]} &>> $out_log
		cp $name* $workdir/$outdir/ &>> $out_log

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR &>> $out_log
		ls -lhAR 
		cp $out_log $workdir/$logdir/
        
	EOS
) 
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0 
