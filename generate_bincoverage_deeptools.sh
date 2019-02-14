#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
quality=10
dedup='yes'
extend=200
binsize=1000
raw='no'

# help message
help_message="
Wrapper to run DeepTools multiBamSummary to generate genome coverage data.

usage:
    bash $(basename "$0") [-options] -b <BAM>
required arguments:
    -b|--bam : comma seperated list of bam files
optional arguments:
    -o|--outdir : out directory for .npz file (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output file (default = BAM filename)
    -r|--region : region to analyse (e.g. 'chr1' 'chr1:1000-100000') (default = NULL)
    -e|--extend : length to extend reads (bp) (default = 200)
    -q|--quality : min alignment quality to filter reads (default = 10)
                   (to disable simply set to 0)
    -d|--dedup : whether to filter duplicates [yes,no] (default = yes)
    -s|--binsize : binning for track output (default = 1000)
    --raw : whether to output raw counts matrix [yes|no] (default = no)

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam_list=$2
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
        -r|--region)
            region=$2
            shift
            ;;
        -e|--extend)
            extend=$2
            shift
            ;;
        -q|--quality)
            quality=$2
            shift
            ;;
        -d|--dedup)
            dedup=$2
            shift
            ;;
        -s|--binsize)
            binsize=$2
            shift
            ;;
        --raw)
            raw=$2
            shift
            ;;
        *)
            printf "Error: Unrecognised argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check bam argument is provided and split
if [[ -z "${bam_list:-}" ]]; then
    printf "\nERROR: Input BAM file required\n"
    echo "$help_message"; exit 1
else
    IFS=',' read -r -a bam_array <<< $bam_list
fi

# check bam files exist and add to args
for bam in ${bam_array[@]}; do
    if [[ ! -f "$bam" ]]; then
        printf "\nERROR: input BAM file does not exist: %s/%s/%s\n" $bam
        echo "$help_message"; exit 1
    else
        bam_base=$(basename "$bam")
        bam_arg="${bam_arg:-} $bam_base"
        bam_cp="${bam_cp:-}; cp $workdir/$bam* ."
    fi
done
bam_arg=${bam_arg#* }
bam_cp=${bam_cp#*; }

# if no name provided extract from first bam
if [[ -z "$name" ]]; then
    name=$(basename "${bam_array[1]}")
    name=${name%%.*}
fi

# set optional arg
if [[ $dedup = yes ]]; then
    dedup_arg="--samFlagExclude 1024"
fi
if [[ -n ${binsize:-} ]]; then
    bs_arg="--binSize $binsize"
fi
if [[ -n ${region:-} ]]; then
    region_arg="--region $region"
fi
if [[ -n ${extend:-} ]]; then
    extend_arg="--extendReads $extend"
fi
if [[ -n ${quality:-} ]]; then
    qual_arg="--minMappingQuality $quality"
fi
if [[ $raw = 'yes' ]]; then
    raw_arg="--outRawCounts $name.bin_counts.tab"
fi

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set log files
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# set commands
multibam_cmd=("multiBamSummary bins --bamfiles $bam_arg --outFileName $name.npz"
              "--numberOfProcessors 18 ${region_arg:-} ${extend_arg:-}"
              "${qual_arg:-} ${dedup_arg:-} ${bs_arg:-} ${raw_arg:-}")

# run job
script=$(cat <<- EOS
	#!/bin/bash
	#PBS -l walltime=24:00:00
	#PBS -l select=1:mem=12gb:ncpus=18
	#PBS -j oe
	#PBS -N $name.multibam
	#PBS -q med-bio
	#PBS -o $std_log
	
	# load modules
	module load samtools/1.2
	module load anaconda
	source activate my_deeptools-2.4.2    

	printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

	# copy to scratch
	$bam_cp &>> $out_log

	# generate coverage
	${multibam_cmd[@]} &>> $out_log
	
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
