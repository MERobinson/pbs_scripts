#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=""
check="yes"
date=`date '+%Y-%m-%d %H:%M:%S'`
scr_name=$(basename "$0")
padding=0
bin_length=1000
merge_rule='OVERLAPPING_ONLY'

# help message
help_message="
Wrapper to GATK4 PreprocessIntervals - to prepare bins for coverage collection.

usage:
    bash $scr_name [-options] -F <FASTA> -i <intervals>
required arguments:
    -i|--intervals : interval list to process
    -F|--fasta : whole genome FASTA file
optional arguments:
    -p|--padding : length (in bp) to add to interval ends post-merge (default = 0)
    -b|--bin_length : length of bins to split intervals into (default = 1000)
    -m|--merge_rule : interval merging rule [see gatk doc] (default = OVERLAPPING_ONLY)
    -o|--outdir : output directory for processed intervals (default = PWD)
    -ld|--logdir : output directory for log files (default = --outdir)
    -n|--name : name of output file [string] (default = extracted from input file)
    --check : whether to check input files [yes,no] (default = 'yes')
    --depend : dependency list to pass to PBS script (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # depend arguments should have PBS format, e.g. 'afterok:123456,afterok:123457'

"

# parse arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--intervals)
            intervals=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -p|--padding)
            padding=$2
            shift
            ;;
        -b|--bin_length)
            bin_length=$2
            shift
            ;;
        -m|--merge_rule)
            merge_rule=$2
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
        *)
            printf "\nERROR: Undefined argument provided\n"
            echo "$help_message"; exit 2
            ;;
    esac
    shift
done

# set outdirs
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arg
if [[ -z ${intervals:-} ]]; then
	printf "\nERROR: --intervals argument is required\n"
	echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument is required\n"
    echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$intervals ]]; then
        printf "\nERROR: Intervals file is not readable: %s/%s\n" $workdir $intervals
        echo "$help_message"; exit 1
    elif [[ ! -r $workdir/$fasta ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    fi
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
    name=$(basename $intervals)
    name="${name%%.*}.pro"
fi

# get basename/prefix
fasta_base=$(basename "$fasta")
fasta_prefix=${fasta%%.*}
intervals_base=$(basename "$intervals")

# setup output directories
mkdir -p $workdir/$logdir
mkdir -p $workdir/$outdir

# set commands
prepro_command=("gatk PreprocessIntervals"
                "-R $fasta_base"
                "-L $intervals_base"
                "-O $name.interval_list"
                "--bin-length $bin_length"
                "--padding $padding"
                "--interval-merging-rule $merge_rule")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=8gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.prepro
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}		

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# load modules
		module load gatk/4.0
		module load java/jdk-8u66
		module load samtools/1.2

		# copy files to scratch
		cp -rL $workdir/$intervals* . &>> $out_log
		cp -rL $workdir/$fasta_prefix* . &>> $out_log

		# preprocess
		${prepro_command[@]} &>> $out_log

		# copy output file to outdir 
		cp $name.interval_list* $workdir/$outdir/  &>> $out_log
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		printf "JOBID: %s\n" \${PBS_JOBID:-} >> $out_log
		ls -lhAR >> $out_log
		cp $out_log $workdir/$logdir/
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
