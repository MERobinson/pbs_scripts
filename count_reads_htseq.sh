#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'
date=`date '+%Y-%m-%d %H:%M:%S'`
format='bam'
order='name'
strand='yes'
minqual='10'
feature='exon'
mode='union'

# help message
help_message="
Wrapper to align reads with BWA

usage:
    bash $(basename "$0" .sh) [-options] -fq1 <FASTQ>
required arguments:
    -i|--input : input bam file
    -G|--gtf : reference gtf file
optional arguments:
    -f|--format : input file format [sam|bam] (default = bam)
    -r|--order : sort order of inut file [name|pos] (default = name)
    -s|--strand : strandedness of reads [no|yes|reverse] (default = yes)
    -a|--minqual : minimum alignment quality [numeric] (default = 10)
    -t|--feature : feature type to count [exon|gene] (default = exon)
    -d|--add_attr : additional attributes to add e.g. gene_name (default = NULL)
    -m|--mode : count mode [union|intersection-strict|intersection-nonempty] (defualt = union)
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -l|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log output directory inherits from --outdir unless specified

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--input)
            bam=$2
            shift
            ;;
        -G|--gtf)
            gtf=$2
            shift
            ;;
        -f|--format)
            format=$2
            shift
            ;;
        -r|--order)
            order=$2
            shift
            ;;
        -s|--strand)
            strand=$2
            shift
            ;;
        -a|--minqual)
            minqual=$2
            shift
            ;;
        -t|--feature)
            feature=$2
            shift
            ;;
        -m|--mode)
            mode=$2
            shift
            ;;
        -d|--add_attr)
            add_attr="--additional-attr=$2"
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
if [[ -z ${bam:-} ]]; then
    printf "\nERROR: --input argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${gtf:-} ]]; then
    printf "\nERROR: --gtf argument required\n"
    echo "$help_message"; exit 1
fi

# check files 
if [[ "${check:-}" = yes ]]; then
    if [[ ! -r $bam ]]; then
        printf "\nERROR: BAM cannot be read: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    elif [[ ! -r $gtf ]]; then
        printf "\nERROR: GTF cannot be read: %s/%s\n" $workdir $gtf
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
bam_base=$(basename "$bam")
gtf_base=$(basename "$gtf")

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    name=${bam%%.*}
    name=$(basename "$name")
fi

# set sample name to name if not provided
if [[ -z "${sm:-}" ]]; then
    sm="SAMPLE_NAME=$name"
fi

# set commands
htseq_cmd=("htseq-count"
           "-f $format"
           "-r $order"
           "--max-reads-in-buffer 1000000000"
           "-s $strand"
           "-a $minqual"
           "-t $feature"
           "-m $mode"
           "${add_attr:-}"
           "${bam_base}"
           "${gtf_base}"
           "> ${name}.htseq-counts.tsv")

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=20gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.htseq
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load modules
		module load anaconda3/personal
		source activate htseq

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log
		
		# copy resource files to scratch
		cp -L $workdir/$bam* . &>> $out_log
		cp -L $workdir/$gtf* . &>> $out_log

		# run htseq
		printf "\nRunning htseq-count:\n" >> $out_log 
		${htseq_cmd[@]} 2>> $out_log
		cp ${name}.htseq-counts.tsv $workdir/$outdir &>> $out_log

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
