#!/bin/bash

workdir=$PWD
outdir=""

# parse arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam=$2
            shift
            ;;
        -c|--count)
            count=$2
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

# set n
n=$(echo "${count}*1000000" | bc)
n=${n%%.*}

# set name if not provided
if [[ -z ${name:-} ]]; then
    name=$(basename "${bam}")
    name=${name%%.*}
	name=$name.shuf.${count}M
fi

# get basenames/prefix
bam_base=$(basename "${bam}")

# setup output dir
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# run subsampling
jobid=$(cat <<- EOS | qsub -N $name.shuf -
		#/bin/bash
		#PBS -l walltime=02:00:00
		#PBS -l select=1:mem=10gb:ncpus=5
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.shuf_bam.log
		${depend:-}
		
		module load samtools/1.3.1
		
		# copy to scratch
		cp $workdir/$bam* .

		# shuf
		echo "shuffling"
		cat <(samtools view -H $bam_base) <(samtools view $bam_base | shuf -n $n) > $name.sam
		
		# convert
		head $name.sam
		samtools view -b $name.sam > tmp.bam
		head tmp.bam

		# reindex
		echo "sorting"
		samtools sort -@ 4 tmp.bam > $name.bam 
		echo "indexing"
		samtools index $name.bam

		# copy back 
		cp $name.bam* $workdir/$outdir

		ls -lhAR
	EOS
)

