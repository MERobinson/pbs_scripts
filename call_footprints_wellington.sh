#!/usr/bin/env bash
set -o pipefail
set -o nounset 

# default arg
workdir=$PWD
outdir=''
check='on'

# help message
help_message="
Wrapper to call footprints from DNase data with Wellington

usage:
    bash $(basename "$0") [-options] -r <BED> -b <BAM>
required arguments:
    -r|--regions : accessible regions to profile [BED]
    -b|--bam : BAM files to find footprints in [BAM]
optional arguments:
    -n|--name : output filename prefix (default = extracted from input file)
    -o|--outdir : output directory name (default = PWD)
    --fdrlimit : minimum pval to consider sig for FDR (default = -20)
    --logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on,off] (default = on)
    --depend : list of PBS dependencies (default = NULL)
additional info:

"

# parse command line arguments
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-r|--regions)
		    regions=$2
		    shift
		    ;;
		-b|--bam)
		    bam=$2
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
        --fdrlimit)
            fdrlimit="-fdrlimit $2"
            shift
            ;;
        --logdir)
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
		    echo "ERROR: Unrecognised argument: %s %s" $1 $2
		    echo "$help_message"; exit 1
		    ;;
	esac
	shift
done

# check required args
if [[ -z ${regions:-} ]]; then
    printf "\nERROR: --region argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${bam:-} ]]; then
    printf "\nERROR: --bam argument required\n"
    echo "$help_message"; exit 1
fi

# set optional args
if [[ -z ${name:-} ]]; then
    name=$(basename "${bam}")
    name=${name%%.*}
fi
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# check files
if [[ ${check} = 'on' ]]; then
    if [[ ! -r "${workdir}/${regions}" ]]; then
        printf "\nERROR: input file cannot be read: %s/%s\n" $workdir $regions
        echo "$help_message"; exit 1
    elif [[ ! -r "${workdir}/${bam}" ]]; then
        printf "\nERROR: index cannot be read: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    fi
fi

# set commands
pydnase_cmd=("wellington_footprints.py ${fdr_limit:-} -p 12" 
             "preprocessed.bed $(basename ${bam}) ${name}")

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$workdir/$logdir/$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:ncpus=12:mem=20gb
		#PBS -j oe
		#PBS -N fp.$name
		#PBS -q med-bio
		#PBS -o ${std_log}
		${depend:-}

		# load modules
		module load anaconda3/personal
		source activate pydnase

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy input & index to scratch
		cp -L ${workdir}/${bam}* .
		cp -L ${workdir}/${regions}* .

		# pre process bed
		cut -f 1-3 $(basename $regions) > preprocessed.bed 

		# run wellington
		mkdir -p ${name}
		${pydnase_cmd[@]} &>> $out_log

		# copy results to output directory
		cp -r ${name}/* $workdir/$outdir/
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR &>> $out_log
		ls -lhAR 
		cp -r $out_log $workdir/$logdir/
	EOS
)
echo "$script" > $pbs_log

# submit job, echo id and exit
jobid=$(qsub "$pbs_log")
echo "JOBID: $jobid"
exit 0
