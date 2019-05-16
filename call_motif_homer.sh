#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default args
workdir=$PWD
outdir=''
check='on'
date=`date '+%Y-%m-%d %H:%M:%S'`
scr_name=$(basename "$0")
size='-size 200'

# help message
help_message="
Wrapper to run HOMERs findMotifsGenome.pl

usage:
    bash $scr_name [-options] -b <BED> -g <GENOME>
required args:
    -b|--bed : BED file of regions to analyse 
    -g|--genome : genome version [STRING]
optional args:
    --size : frag size [INT|'given'] (default = 200)
    -bg|--background : background regions [BED] (default = NULL)
    --mknown : known motif file to use in place of in-built DB (default = NULL)
    -n|--name : output filename prefix (default = FASTA prefix)
    -o|--outdir : output directory for results (default = PWD)
    -l|--logdir : output directory for logs (default = --outdir)
    --check : whther to check input files [on|off] (default = on)
    --depend : list of PBD dependencies (default = NULL)
additional info:
    # all paths should be relative to working dir
    # check and depend options used for job scheduling
    # log output dir inherits from --outdir unless specified
    # see HOMER docs for more info: http://homer.ucsd.edu/homer/ngs/peakMotifs.html

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bed)
            bed=$2
            shift
            ;;
        -g|--genome)
            genome=$2
            shift
            ;;
        -bg|--background)
            background=$2
            shift
            ;;
        --size)
            size="-size ${2}"
            shift
            ;;
        --mknown)
            mknown=$2
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
    esac
    shift
done

# check required args
if [[ -z ${bed:-} ]]; then
    printf "\nERROR: --bed argument is required\n"
    echo "$help_message"; exit 1
elif [[ -z ${genome:-} ]]; then
    printf "\nERROR: --genome argument if required\n"
    echo "$help_message"; exit 1
fi

# check files
if [[ $check = 'on' ]]; then
    if [[ ! -r $workdir/$bed ]]; then
        printf "\nERROR: input BED file cannot be read: %s/%s\n" $workdir $bed
        echo "$help_message"; exit 1
    elif [[ -n ${background:-} ]] && [[ ! -r $workdir/$background ]]; then
        printf "\nERROR: input motif file cannot be read: %s/%s\n" $workdir $background
        echo "$help_message"; exit 1
    fi
fi

# set outdirs
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set name
if [[ -z ${name:-} ]]; then
    name=$(basename "$bed")
    name=${name%%.*}.ame
fi

# handle optional args
if [[ -n ${background:-} ]]; then
    bg_cp="cp $workdir/$background ."
    bg_arg="-bg $(basename $background)"
fi
if [[ -n ${mknown:-} ]]; then
    mknown_cp="cp $workdir/$mknown ."
    mknown_arg="-mknown $(basename $mknown)"
fi

# set command
homer_cmd=("findMotifsGenome.pl $(basename ${bed}) ${genome}"
           "homer_output ${size} ${bg_arg:-} ${mknown_arg:-}")

# set log filenames
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=10gb:ncpus=4
		#PBS -j oe
		#PBS -N homer.${name}
		#PBS -q med-bio
		#PBS -o ${std_log}
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		module load homer/4.9

		cp $workdir/$bed . &>> $out_log
		${bg_cp:-} &>> $out_log
		${mknown_cp:-} &>> $out_log

		${homer_cmd[@]} &>> $out_log

		cp homer_output/* $workdir/$outdir/ &>> $out_log

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
