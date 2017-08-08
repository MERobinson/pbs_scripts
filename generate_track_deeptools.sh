#!/bin/bash

# default arg
workdir=$PWD
trackdir=tracks
bamdir=bam
resdir=resources
logdir=logs
genome=hg38
extend=200
quality=10
dedup=yes
binsize=200
fasta=genome.fa

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
    -r|--resdir : directory within workdir containing FASTA file (default = resources)
    -g|--genome : genome version - used in trackline (default = hg38)
    -n|--name : name prefix for output file (default = BAM filename)
    -e|--extend : bp extension for generating coverage track (default = 200)
    -q|--quality : min alignment quality to filter reads (default = 10)
                   (to disable simply set to 0)
    -d|--dedup : whether to filter duplicates [yes,no] (default = yes)
    -s|--binsize : binning for track output (default = 200)
    -f|--fasta : name of genome FASTA file within resdir (default = genome.fa)
example:
    bash $(basename "$0") -d -g hg38 -i REH_H3K27ac_ChIPseq.bam

"

# parse command line arg
while [[ $# -gt 1 ]]; do
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
        -r|--resdir)
        resdir=$2
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
        dedup=$2
        shift
        ;;
        -s|--binsize)
        binsize=$2
        shift
        ;;
        -f|--fasta)
        fasta=$2
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
if [[ $dedup = yes ]]; then
    dedup_arg="-F 1024"
fi

# set trackline text
trackline="track type=bedGraph name=$name visibility=2 windowingFunction=mean autoScale=off"
trackline=$trackline" smoothingWindow=10 viewLimits=0:50 maxHeightPixels=0:75:150"

# get list of canonical chromosomes
chr_list=$(cat $workdir/$resdir/$fasta | grep -Eo "^>chr[0-9MXY]+\b" | \
           grep -Eo "chr[0-9XY]+"| tr "\n" " ") 

# create required dirs
mkdir -p $workdir/$logdir
mkdir -p $workdir/$trackdir

# run job
JOBID=$(cat <<- EOS | qsub -N $name.bwa - 
	#!/bin/bash
	#PBS -l walltime=05:00:00
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
	samtools view -b $dedup_arg -q $quality $bam $chr_list > $name.$genome.filt.bam
    samtools index $name.$genome.filt.bam

    # temp - output idxstats
    samtools idxstats $bam
    samtools idxstats $name.$genome.filt.bam

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
	cat $name.$genome.bedGraph | grep -E "^chr[0-9XYM]+\b" > tmp.txt
	echo $trackline | cat - tmp.txt > $name.$genome.bedGraph
	
	cp $name.$genome.bedGraph $workdir/$trackdir
	ls -lhAR
	EOS
)
