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
format='bam'
wheader='off'
header='off'
refinfo='off'
count='off'
valid='off'
sam='off'
threads=8

# help message
help_message="
Simple wrapper to Sambamba view

usage:
    bash $scr_name [-options] -i <BAM>
required arguments:
    -i|--input : input file [BAM|SAM|..]
optional arguments:
    -n|--name : name prefix for output files (default = extracted from first bam)
    --outdir : outut directory for merged BAM (default = PWD)
    --qcdir : output directory for qc metrics (default = --outdir)
    --logdir : output directory for log files (default = --outdir)
    -F|--filter : custom filter - make sure quoted
    -f|--format : output format [sam|bam|cram|json|unpack] (default = bam)
    -h|--with-header : include header if SAM format [on|off] (default = off)
    -H|--header : only output header [on|off] (default = off)
    -I|--reference-info : only output ref names and lengths [on|off] (default = off)
    -L|--regions : output only reads overlapping regions [BED]
    -c|--count : only count records [on|off] (default = off)
    -v|--valid : output only valid anlignments [on|off] (default = off)
    -t|--threads : number of threads (default = 8)
    -s|--subsample : fraction of reads to subsample 
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
        -F|--filter)
            filter="--filter=$2"
            shift
            ;;
        -f|--format)
            format="$2"
            shift
            ;;
        -h|--with-header)
            wheader=$2
            shift
            ;;
        -H|--header)
            header=$2
            shift
            ;;
        -I|--reference-info)
            refinfo=$2
            shift
            ;;
        -L|--regions)
            regions=$2
            shift
            ;;
        -c|--count)
            count=$2
            shift
            ;;
        -v|--valid)
            valid=$2
            shift
            ;;
        -t|--threads)
            threads=$2
            shift
            ;;
        -s|--subsample)
            subsample="--subsample=$2"
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
    elif [[ -n ${regions:-} ]] && [[ ! -r ${workdir}/${regions} ]]; then
        printf "\nERROR: input regions not readable: %s/%s\n" ${workdir} ${regions}
        echo "${help_message}"; exit 1
    fi
fi

# set region arg if given
if [[ -n ${regions:-} ]]; then
    regions_cp="cp ${workdir}/${regions}* ."
    regions_arg="--regions=$(basename ${regions})"
fi

# if no name provided extract from first bam
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
if [[ ${wheader} = 'on' ]]; then
    wheader_arg="--with-header"
fi
if [[ ${header} = 'on' ]]; then
    header_arg="--header"
fi
if [[ ${refinfo} = 'on' ]]; then
    refinfo_arg="--reference-info"
fi
if [[ ${count} = 'on' ]]; then
    count_arg="--count"
fi
if [[ ${valid} = 'on' ]]; then
    wheader_arg="--valid"
fi
if [[ ${sam} = 'on' ]]; then
    sam_arg="--sam-input"
fi

# set commands
view_cmd=("sambamba view ${filter:-} ${subsample:-}"
         "--format=$format ${wheader_arg:-} ${header_arg:-}"
         "${refinfo_arg:-} ${count_arg:-} ${valid_arg:-} ${sam_arg:-}"
         "-o ${name}.${format} ${regions_arg:-}"
         "$(basename ${input})")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=8gb:ncpus=8
		#PBS -j oe
		#PBS -N view.$name
		#PBS -o $std_log
		###PBS -q med-bio
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# load modules
		module load anaconda3/personal
		source activate sambedtools

		# copy accross bam and fasta
		cp ${workdir}/${input%%.*}* . 
		${regions_cp:-}

		# alignment metrics
		${view_cmd[@]} 
		samtools index ${name}.${format}
		cp ${name}.${format} ${workdir}/$outdir/

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` 
		ls -lhAR 
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
