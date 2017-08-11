#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'

# help message
help_message="
usage:
    bash $(basename "$0") [-options] -f1 <FASTQ1> -F <FASTA> -I <INDEX>
purpose:
    # Wrapper script to run STAR via PBS script
required arguments:
    -f1|--fastq1 : input FASTQ file
    -I|--index : STAR index directory
optional arguments:
    -f2|--fastq2 : mate FASTQ file if PE [FASTQ] (default = NULL)
    -n|--name : prefix for output files [string] (default = extracted from fastq)
    -o|--outdir : output directory [string] (default = ".")
    -l|--logdir : output dir for log files [string] (default = --outdir)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list of PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # example depenency list: 'afterok:123456,afterok:123457'
outputs:
    # outputs STAR results to --outdir & log file to --logdir 

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -f1|--fastq1)
            fq1=$2
            shift
            ;;
		-f2|--fastq2)
            fq2=$2
            shift
            ;;
        -I|--index)
            index=$2
            shift
            ;;
        -n|--name)
            name=$2
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
        --check)
            check=$2
            shift
            ;;
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arg
if [[ -z ${fq1:-} ]]; then
    printf "\nERROR: no --fastq1 argument provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${index:-} ]]; then
    printf "\nERROR: no --index argument provided\n"
    echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$fq1 ]]; then
        printf "\nERROR: FASTQ file is not readable: %s/%s\n" $workdir $fq1
        echo "$help_message"; exit 1
    fi
    if [[ ! -z ${fq2:-} ]] & [[ ! -r $workdir/$fq2 ]]; then
        printf "\nERROR: FASTQ file is not readable: %s/%s\n" $workdir $fq2
        echo "$help_message"; exit 1
    fi
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
    name=$(basename ${fq1})
    name=${name%%.*}
fi

# get basenames/prefix
fq1_base=$(basename "$fq1")
index_base=$(basename "$index")

# set PE args if mate provided
if [[ ! -z ${fq2:-} ]]; then 
    echo "Mate detected: running in PE mode"
    fq2_base=$(basename "$fq2")
    fq2_cp_arg="cp $workdir/$fq2 ."
fi

# set compression arg
fq_ext=${fq1##*.}
if [[ $fq_ext = "gz" ]]; then
    compress_arg="--readFilesCommand zcat"
fi

# setup output dir
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# run PBS script
jobid=$(cat <<- EOS | qsub -N $name.star - 
		#!/bin/bash
		#PBS -l walltime=20:00:00
		#PBS -l select=1:mem=40gb:ncpus=20
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.align_fastq_star.log
		${depend:-}
		
		# load modules
		module load samtools/1.2
		module load star/2.5.0
		module load intel-suite/2015
		module load gcc/5.4.0
		
		# copy files to scratch
		cp -rL $workdir/$index* .
		cp -L $workdir/$fq1 .
		$fq2_cp_arg

		# run STAR
		STAR --runMode alignReads \
			--runThreadN 20 \
			--genomeDir $index_base \
			--readFilesIn ${fq1_base} ${fq2_base:-} \
			--outFilterType BySJout \
			--outFilterMultimapNmax 20 \
			--alignSJoverhangMin 8 \
			--alignSJDBoverhangMin 1 \
			--outFilterMismatchNoverLmax 0.04 \
			--alignIntronMin 20 \
			--alignIntronMax 1000000 \
			--alignMatesGapMax 1000000 \
			--outSAMstrandField intonMotif \
			--quantMode GeneCounts \
			--readNameSeparator _ \
			--outFileNamePrefix ${name}_ \
			--outBAMsortingThreadN 20 \
			--outSAMtype BAM SortedByCoordinate \
			--outWigType wiggle \
			--outWigStrand Unstranded \
			${compress_arg:-}
	
		# index bam 
		samools index  ${name}_Aligned.sortedByCoord.out.bam 
		
		# copy output to outdir
		cp ${name}_* $workdir/$outdir

		ls -lhAR
		EOS
) 
echo "JOBID: $jobid"
