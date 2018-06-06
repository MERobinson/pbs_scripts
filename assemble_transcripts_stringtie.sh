#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='on'
date=`date '+%Y-%m-%d %H:%M:%S'`
exprs='off'
merge='off'
abund='off'

# help message
help_message="
Wrapper to assemble transcripts with Stringtie

usage:
    bash $(basename $0) [-options] --input <BAM>
required arguments:
    -i|--input : input bam/gtf file(s)
optional arguments:
    -m|--merge : whether to run in gtf merging mode [on|off] (default = off)
    -G|--gtf : reference GTF/GFF file (default = NULL)
    -s|--strand : strand orientation of reads [rf|fr] (default = NULL)
    -e|--exprs : whether to limit to input gtf [on|off] (default = off)
    -b|--ballgown : name of ballgown output dir (default = NULL)
    -a|--abund : whether to measure abundances [on|off] (default = off)
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -l|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on|off] (default = on)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    !! if run in merge mode input should be a comma separated list of gtf,
       if not, input should be a single bam file !!
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified
    # ballgown data will only be output if --ballgown output dir argument given

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--input)
            input=$2
            shift
            ;;
        -m|--merge)
            merge=$2
            shift
            ;;
        -G|--gtf)
            gtf=$2
            shift
            ;;
        -s|--strand)
            strand=$2
            shift
            ;;
        -e|--exprs)
            exprs=$2
            shift
            ;;
        -b|--ballgown)
            bgdir=$2
            shift
            ;;
        -a|--abund)
            abund=$2
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

# check required arguments
if [[ -z ${input:-} ]]; then
    printf "\nERROR: --input argument required\n"
    echo "$help_message"; exit 1
fi
if [[ ${merge} = 'on' ]]; then
    printf "\nRunning in transcript merge mode\n"
    IFS=',' read -r -a merge_array <<< $input
else
    printf "\nRunning in transcript assembly mode\n"
    bam=$input
fi

# check files 
if [[ "${check:-}" = 'on' ]]; then
    if [[ -n ${bam:-} ]] && [[ ! -r $bam ]]; then
        printf "\nERROR: BAM cannot be read: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    elif [[ -n ${gtf:-} ]] && [[ ! -r ${gtf:-} ]]; then
        printf "\nERROR: GTF/GFF cannot be read: %s/%s\n" $workdir ${gtf:-}
        echo "$help_message"; exit 1
    elif [[ -n ${merge_array:-} ]]; then
        for f in ${merge_array[@]}; do
            if [[ ! -r $f ]]; then
                printf "\nERROR: GTF cannot be read: %s/%s\n" $workdir ${f:-}
                echo "$help_message"; exit 1
            fi
        done
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
if [[ -n ${merge_array:-} ]]; then
    for f in ${merge_array[@]}; do
        merge_cp="${merge_cp:-}cp ${workdir}/${f} .;"
        merge_input="${merge_input:-}$(basename $f) "
    done
    merge_cp=${merge_cp%;*}
    merge_input="${merge_input}--merge"
else
    bam_base=$(basename "$bam")
    bam_prefix=${bam%%.*}
    bam_cp="cp ${workdir}/${bam_prefix}* ."
fi

# extract filename prefix if not provided
if [[ -z ${name:-} ]] && [[ -n ${bam:-} ]]; then
    name=${bam_base%%.*}
elif [[ -z ${name:-} ]] && [[ -n ${merge_array} ]]; then
    name="merged"
fi

# set optional arguments
if [[ -n "${gtf:-}" ]]; then
    gtf_cp="cp $workdir/$gtf* ."
    gtf_base=$(basename "$gtf")
    gtf_arg="-G ${gtf_base}"
fi
if [[ -n ${strand:-} ]]; then
    strand_arg="--${strand}"
fi
if [[ -n ${bgdir:-} ]]; then
    mkdir -p $bgdir
    bg_arg="-b $name"
    bg_cp="cp -r $name $workdir/$bgdir/"
fi
if [[ ${exprs} = 'on' ]]; then
    expr_arg="-e"
fi
if [[ ${abund} = 'on' ]]; then
    abund_arg="-A ${name}.abund.tab"
    abund_cp="cp ${name}.abund.tab $workdir/$outdir"
fi

# set command
stringtie_cmd=("stringtie"
               "${bam_base:-}"
               "${merge_input:-}"
               "-v"
               "-o ${name}.gtf"
               "${abund_arg:-}"
               "${gtf_arg:-}"
               "${strand_arg:-}"
               "${bg_arg:-}"
               "${expr_arg:-}")

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
		#PBS -N stringtie.$name
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load modules
		module load stringtie/1.3.3
		module load samtools/1.2

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log
		
		# copy files to scratch
		${bam_cp:-} &>> $out_log
		${gtf_cp:-} &>> $out_log
		${merge_cp:-} &>> $out_log

		# run stringtie
		printf "\nRunning stringtie:\n" >> $out_log 
		${stringtie_cmd[@]} &>> $out_log
		cp $name.gtf $workdir/$outdir/ &>> $out_log
		${abund_cp:-} &>> $out_log
		${bg_cp:-}

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
