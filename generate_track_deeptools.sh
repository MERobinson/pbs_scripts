#!/bin/bash

# default arg
workdir=$PWD
trackdir=tracks
bamdir=bam
logdir=logs
genome=hg38
extend=200
quality=10
dedup=TRUE
binsize=200

# help message
help_message="
usage:
    bash $(basename "$0") [-wrtblgneh] -b <BAM>
purpose:
    # Wrapper to run deeptools bamCoverage to generate genome tracks.
required arguments:
    -i|--inbam : Input bam file
optional arguments:
    -w|--workdir : working directory - used as base dir for all input/output (default = pwd)
    -t|--trackdir : output directory for genome tracks (default = tracks)
    -b|--bamdir : directory to find bam files (default = bam)
    -l|--logdir : output directory for log files (default = logs)
    -g|--genome : genome version - used in trackline (default = hg38)
    -n|--name : name prefix for output file (default = BAM filename)
    -e|--extend : bp extension for generating coverage track (default = 200)
    -q|--quality : min alignment quality to filter reads (default = 10)
                   (to disable simply set to 0)
    -d|--dedup : flag to disable duplicate removal (default = on)
    -s|--binsize : binning for track output (default = 200)
example:
    bash $(basename "$0") -d -g hg38 -i REH_H3K27ac_ChIPseq.bam

"

# parse command line arg
while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
        -i|--inbam)
        bam=$2
        shift
        ;;
        -w|--workdir)
        workdir=$2
        shift
        ;;
        -t|--trackdir)
        trackdir=$2
        shift
        ;;
        -b|--bamdir)
        bamdir=$2
        shift
        ;;
        -l|--logdir)
        logdir=$2
        shift
        ;;
        -g|--genome)
        genome=$2
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
        dedup=FALSE
        shift
        ;;
        -s|--binsize)
        binsize=$2
        shift
        ;;
        *)
        echo "Error: Illegal argument"
        echo "$help_message"
        exit 1
        ;;
    esac
    shift
done

# check indir exists
if [[ ! -d "$workdir/$bamdir" ]]; then
    printf "\nERROR: input directory does not exist: %s/%s/%s\n" $workdir $bamdir
    echo "$help_message"; exit 1
fi

# check inbam argument is provided and file exists
if [[ -z "$bam" ]]; then
    printf "\nERROR: Input BAM file required\n"
    echo "$help_message"; exit 1
elif [[ ! -f "$workdir/$bamdir/$bam" ]]; then
    printf "\nERROR: input BAM file does not exist: %s/%s/%s\n" $workdir $bamdir $bam
    echo "$help_message"; exit 1
fi

# if no name provided extract from FASTQ filename
if [[ -z "$name" ]]; then
    name=${bam%%.*}
fi

# set dedup argument
if [[ $dedup = TRUE ]]; then
    dedup_arg="-F 1024"
fi

# set trackline text
trackline="track type=bedGraph name=$name visibility=2 windowingFunction=mean autoScale=off"
trackline=$trackline" smoothingWindow=10 viewLimits=0:20 maxHeightPixels=0:75:150"

# get list of canonical chromosomes
#chr_list=$(cat $workdir/$resdir/$fasta | grep -E "^>" | grep -Eo "chr[0-9MXY]+\b") 

# create required dirs
mkdir -p $workdir/$logdir
mkdir -p $workdir/$trackdir

# run job
JOBID=$(cat <<- EOS | qsub -N $name.bwa - 
	#!/bin/bash
	#PBS -l walltime=01:00:00
	#PBS -l select=1:mem=10gb:ncpus=20
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $workdir/$logdir/$name.generate_track_deeptools.log
	
	# load modules
	module load samtools/1.2
	module load anaconda/2.4.1
	module load deeptools/2.4.2
	source activate deeptools-2.4.2    

	# copy resource files to scratch
	cp $workdir/$bamdir/$bam* .

	# filter bam prior to making track
	samtools view -b $dedup_arg -q $quality $name.$genome.bam > $name.$genome.filt.bam

	# generate coverage track
	bamCoverage \
		--extendReads $extend \
		--binSize $binsize \
		--normalizeUsingRPKM \
		--ignoreForNormalization chrX \
		--numberOfProcessors 20 \
		--outFileFormat "bedgraph" \
		--bam $name.$genome.filt.bam \
		--outFileName $name.$genome.bedGraph
	echo $trackline | cat - $name.$genome.bedGraph > tmp.bedGraph
	mv tmp.bedGraph $workdir/$trackdir/$name.$genome.bedGraph

	ls -lhAR
	EOS
)