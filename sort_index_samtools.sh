#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'
date=`date '+%Y-%m-%d %H:%M:%S'`
scr_name=$(basename "$0")
threads=8
name_sort='off'

# help message
help_message="
Simple wrapper to samtools sort 

usage:
    bash $scr_name [-options] -i <BAM>
required arguments:
    -i|--input : input file [BAM|SAM|..]
optional arguments:
    -n|--name : name prefix for output files (default = extracted from first bam)
    --outdir : outut directory for merged BAM (default = PWD)
    --qcdir : output directory for qc metrics (default = --outdir)
    --logdir : output directory for log files (default = --outdir)
    --name_sort : whether to sort by name [on|off] (default = off)
    -t|--threads : number of threads (default = 8)
    --depend : comma-sep list of job dependencies - format afterok:<jobid>,afterok:<jobid>
    --check : whether to check files prior to running job [yes,no] (default = yes)
additional_info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # bam/log/qc output directories inherit from --outdir unless specified

"

# parse command line arg
while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
        -i|--input)
            input=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        --name_sort)
            name_sort=$2
            shift
            ;;
        -t|--threads)
            threads=$2
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
        -q|--qcdir)
            qcdir=$2
            shift
            ;;
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "${help_message}"; exit 1
            ;;
    esac
    shift
done

# check required arg provided
if [[ -z "${input:-}" ]]; then
    printf "\nERROR: --input argument required\n"
    echo "${help_message}"; exit 1
fi
if [[ ${check} = "on" ]]; then
    if [[ ! -r ${workdir}/${input} ]]; then
        printf "\nERROR: input file not readable: %s/%s\n" ${workdir} ${input}
        echo "${help_message}"; exit 1
    fi
fi

# if no name provided extract from bam
if [[ -z ${name:-} ]]; then
    name=$(basename ${input})
    name="${name%%.*}"
fi

# set output directories
if [[ -z ${logdir:-} ]]; then
    logdir=${outdir}
fi
if [[ -z "${qcdir:-}" ]]; then
    qcdir=${outdir}
fi

# create required dirs
mkdir -p ${workdir}/${outdir}
mkdir -p ${workdir}/${logdir}
mkdir -p ${workdir}/${qcdir}

# parse flag args
if [[ ${name_sort} = 'on' ]]; then
    namesort_arg="-n"
fi

# set commands
sort_cmd=("samtools sort ${namesort_arg:-}"
         "-@ $threads"
         "-o ${name}.srt.bam"
         "$(basename ${input})")
index_cmd=("samtools index ${name}.bam")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log

# write job script
job=$(cat <<- EOS | qsub -N sort.$name - 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=8gb:ncpus=8
		#PBS -j oe
		#PBS -o $std_log
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# load modules
		module load anaconda3/personal
		source activate sambedtools

		cp ${workdir}/${input%%.*}* . 
		${sort_cmd[@]}
		mv ${name}.srt.bam ${name}.bam 
		${index_cmd[@]}
		cp ${name}.bam* ${workdir}/$outdir/

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` 
		ls -lhAR 
	EOS
)
echo "JOBID: $job"
exit 0
