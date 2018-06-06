#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default args
workdir=$PWD
check='on'
length_correction='off'
date=`date '+%Y-%m-%d %H:%M:%S'`
scr_name=$(basename "$0")

# help message
help_message="
Wrapper to run AME from the MEME suite to look for enriched motifs

usage:
    bash $scr_name [-options] -s <SEQ> -m <MOTIF>
required args:
    -s|--seq : FASTA file of regions to analyse 
    -m|--motif : MEME formatted motif file
optional args:
    -c|--control : FASTA of control regions to compare against (default = NULL)
    --method : test to use [fisher|ranksum|mhg|4dmhg|spearman|linreg] (defualt = fisher)
    --scoring : scoring method [avg|max|sum|totalhits] (default = totalhits)
    --length-correction : whether to correct for length bias [on|off] (default = off)
    -n|--name : output filename prefix (default = FASTA prefix)
    -o|--outdir : output directory for results (default = PWD)
    -l|--logdir : output directory for logs (default = --outdir)
    --check : whther to check input files [on|off] (default = on)
    --depend : list of PBD dependencies (default = NULL)
additional info:
    # all paths should be relative to working dir
    # check and depend options used for job scheduling
    # log output dir inherits from --outdir unless specified
    # see AME docs for further info: http://meme-suite.org/doc/ame.html?man_type=cmd

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -s|--seq)
            fasta=$2
            shift
            ;;
        -m|--motif)
            motif=$2
            shift
            ;;
        -c|--control)
            control=$2
            shift
            ;;
        --method)
            method="--method $2"
            shift
            ;;
        --scoring)
            scoring="--scoring $2"
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -ld|--logdir)
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
if [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --seq argument is required\n"
    echo "$help_message"; exit 1
elif [[ -z ${motif:-} ]]; then
    printf "\nERROR: --motif argument if required\n"
    echo "$help_message"; exit 1
fi

# check files
if [[ $check = 'on' ]]; then
    if [[ ! -r $fasta ]]; then
        printf "\nERROR: input FASTA file cannot be read: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    elif [[ ! -r $motif ]]; then
        printf "\nERROR: input motif file cannot be read: %s/%s\n" $workdir $motif
        echo "$help_message"; exit 1
    elif [[ -n ${control:-} ]] && [[ ! -r ${control:-} ]]; then
        printf "\nERROR: input control FASTA cannot be read: %s/%s\n" $workdir $control
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
    name=$(basename "$fasta")
    name=${name%%.*}
fi

# handle optional args
if [[ -n ${control:-} ]]; then
    control_arg="--control $(basename $control)"
    control_cp="cp $workdir/$control ."
fi
if [[ $length_correction = 'on' ]]; then
    correct="--length-correction"
fi

# set command
ame_cmd=("ame -o tmp ${control_arg:-} ${correct:-} ${method:-} ${scoring:-}"
         "$(basename $fasta) $(basename $motif)")

# set log filenames
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=10gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.ame
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		module load meme/4.12.0 &>> $out_log

		cp $workdir/$fasta . &>> $out_log
		cp $workdir/$motif . &>> $out_log
		${control_cp:-}

		${ame_cmd[@]} &>> $out_log

		cp tmp/ame.html $workdir/$outdir/$name.html &>> $out_log
		cp tmp/ame.txt $workdir/$outdir/$name.txt &>> $out_log

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
