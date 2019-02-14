#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
format='bigwig'
outdir=''
extend=200
quality=10
dedup='yes'
binsize=200
mode='log2'

# help message
help_message="
usage:
    bash $(basename "$0") [-options] -b <BAM>
purpose:
    # Wrapper to run deeptools bamCompare to generate tracks of difference in signal between two bam.
required arguments:
    -t|--test : BAM file for test sample
    -c|--cntrl : BAM file for control sample
    -F|--fasta : Reference FASTA file
optional arguments:
    -o|--outdir : output directory for track files (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output file (default = BAM filename)
    -f|--format : output format [bedgraph,bigwig] (default = bedgraph)
    -m|--mode : type of comparison [log2|ratio|subtract|mean] (default = 'log2')
    -g|--genome : genome version - used in trackline if provided (default = NULL)
    -e|--extend : bp extension for generating coverage track (default = 200)
    -q|--quality : min alignment quality to filter reads (default = 10)
                   (to disable simply set to 0)
    -d|--dedup : whether to filter duplicates [yes,no] (default = yes)
    -s|--binsize : binning for track output (default = 200)
additional info:
    # all paths should be relative to the current working directory
    # FASTA file is used to obtain a list of canonical chromosomes for filtering
      the BAM prior to generating tracks (regex for canonical = chr[0-9MXY]+)

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -t|--test)
            test_bam=$2
            shift
            ;;
        -c|--cntrl)
            control_bam=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
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
        -f|--format)
            format=$2
            shift
            ;;
        -m|--mode)
            mode=$2
            shift
            ;;
        -g|--genome)
            db=$2
            shift
            ;;
        -n|--name)
            name=$2
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
            echo "\nError: Illegal argument: %s %s" $1 $2
            echo "$help_message"; exit 1
        ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arguments
if [[ -z ${test_bam:-} ]]; then
    printf "\nERROR: no test BAM file provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: no FASTA file provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${control_bam:-} ]]; then
    printf "\nERROR: no control BAM file provided\n"
    echo "$help_message"; exit 1
fi

# if no name provided extract from BAM
if [[ -z "$name" ]]; then
    name=${test_bam%%.*}
fi

# set dedup argument
if [[ $dedup = yes ]]; then
    dedup_arg="-F 1024"
fi

# set format dependent arg
if [[ $format = 'bedgraph' ]]; then
	extension="bedGraph"
elif [[ $format = 'bigwig' ]]; then
    extension='bw'
else
    printf "\nERROR: format not recognised, must be either bedgraph or bigwig\n"
    echo "$help_message"; exit 1
fi
    
# get list of canonical chromosomes
chr_list=$(cat $fasta | grep -Eo "^>chr[0-9MXY]+\b" | \
           grep -Eo "chr[0-9XY]+"| tr "\n" " ") 

# set base
test_base=$(basename "$test_bam")
control_base=$(basename "$control_bam")

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set log names
scr_name=$(basename $0 .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# run job
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=10gb:ncpus=20
		#PBS -j oe
		#PBS -N $name.track
		#PBS -q med-bio
		#PBS -o $std_log
	
		# load modules
		module load samtools/1.2
		module load anaconda3/personal
		source activate deeptools-3.1.0    

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy resource files to scratch
		cp $workdir/$test_bam* . &>> $out_log
		cp $workdir/$control_bam* . &>> $out_log

		# filter bam prior to making track
		samtools view -b ${dedup_arg:-} -q $quality $test_base $chr_list > $name.test.bam
		samtools view -b ${dedup_arg:-} -q $quality $control_base $chr_list > $name.control.bam
		samtools index $name.test.bam
		samtools index $name.control.bam

		# generate coverage track
		bamCompare \
			--bamfile1 $name.test.bam \
			--bamfile2 $name.control.bam \
			--operation $mode \
			--extendReads $extend \
			--binSize $binsize \
			--scaleFactorsMethod readCount \
			--ignoreForNormalization chrX \
			--numberOfProcessors 20 \
			--outFileFormat $format \
			--outFileName $name.$extension &>> $out_log

		# copy output to outdir
		cp $name.$extension $workdir/$outdir/

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
