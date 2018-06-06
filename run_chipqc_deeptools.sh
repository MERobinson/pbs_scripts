#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
quality=10
dedup='yes'
extend=150
binsize=500
log_file=$(basename $0 .sh).$(date "+%Y-%m-%d").log

# help message
help_message="
Pipeline to run various ChIP QC tools from DeepTools library.

usage:
    bash $(basename "$0") [-options] -b <list,of,bam>
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
    -s|--binsize : binning for track output (default = 500)

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
        *)
            printf "\nError: Unrecognised argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# set logdir
if [[ -n ${logdir:-} ]]; then
    logdir=$outdir
fi

# check bam argument is provided and split
if [[ -z "${bam_list:-}" ]]; then
    printf "\nERROR: Input BAM file required\n"
    echo "$help_message"; exit 1
fi

# if no name provided extract from first bam
if [[ -z "$name" ]]; then
    name=$(basename "${bam_array[1]}")
    name=${name%%.*}.deeptools
fi

# set optional arg
if [[ -n $dedup = yes ]]; then
    dedup_arg="--dedup $dedup"
fi
if [[ -n ${binsize:-} ]]; then
    bs_arg="--binsize $binsize"
fi
if [[ -n ${region:-} ]]; then
    region_arg="--region $region"
fi
if [[ -n ${extend:-} ]]; then
    extend_arg="--extend $extend"
fi
if [[ -n ${quality:-} ]]; then
    qual_arg="--qual $quality"
fi

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# run multiBamSummary
multibam_cmd=("bash $scriptdir/generate_bincoverage_deeptools.sh --bam $bam_arg"
              "--outdir $outdir --logdir $logdir --name $name ${region_arg:-}"
              "--binsize $binsize ${dedup_arg:-} ${qual_arg:-} ${extend_arg:-}")
printf "\nRunning multiBamSummary with command:\n" >> $log_file
printf "\t%s\n" ${multibam_cmd[@]} >> $log_file
${multibam_cmd[@]}

# plot heatmap & pca

