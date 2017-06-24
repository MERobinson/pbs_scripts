#!/bin/bash

# default arg
WORKDIR=$PWD
RESDIR=resources
LOGDIR=logs
FASTA=genome.fa
GTF=transcriptome.gtf

# help message
USAGE="

bash $(basename "$0") [-wrfgh]

purpose:
	Simple wrapper to generate STAR index
required arguments:
        n/a
optional arguments:
	-w|--workdir : working directory - used as base dir for all input/output (default = pwd)
	-r|--resdir : dir containing [links to] ref FASTA and GTF resources (default = resources) 
        -f|--fasta : name of fasta file in resdir (default = genome.fa)
	-g|--gtf : name of gtf file in resdir (default = transcriptome.gtf)
	-l|--logdir : dir to output log file of run (default = logs)
        -h|--help : print current help message
example:
        bash $(basename "$0") -g lincRNA.gtf
additional info:
	> Recommended to setup a new dir per index containing a resources folder with links to 
	  all required resources (i.e. FASTA and GTF), then run script within this folder.

"

# parse command line arg
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-w|--workdir)
		WORKDIR=$2
		shift
		;;
		-r|--resdir)
		RESDIR=$2
		shift
		;;
		-f|--fasta)
		FASTA=$2
		shift
		;;
		-g|--gtf)
		GTF=$2
		shift
		;;
		-l|--logdir)
		LOGDIR=$2
		shift
		;;
		-h|--help)
		echo "$USAGE"
		exit 1
		;;
		*)
		echo "Error: Illegal argument"
		echo "$USAGE"
		exit 1
		;;
	esac
	shift
done

# check res files exist
if [[ ! -e "$WORKDIR/$RESDIR/$FASTA" ]]; then
	printf "\nERROR: Input FASTA file does not exist: %s/%s/%s\n" $WORKDIR $RESDIR $FASTA 
	echo "$USAGE"; exit 1
elif [[ ! -e "$WORKDIR/$RESDIR/$FASTA" ]]; then
	printf "\nERROR: Input GTF file does not exist: %s/%s/%s\n" $WORKDIR $RESDIR $GTF
	echo "$USAGE"; exit 1
fi

# create required dirs
mkdir -p $WORKDIR/$LOGDIR

JOBID=$(cat <<- EOS | qsub -N STAR_index - 
	#!/bin/bash
	#PBS -l walltime=40:00:00
	#PBS -l select=1:mem=40gb:ncpus=20
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $WORKDIR/$LOGDIR/generate_index_star.log
	
	# load modules
	module load star/2.5.0
	module load intel-suite/2015
	module load gcc/5.4.0
	
	# make temp dir and copy accross res
	mkdir outdir
	cp -L $WORKDIR/$RESDIR/$FASTA* .
	cp -L $WORKDIR/$RESDIR/$GTF* .

	# run STAR 
	STAR --runMode genomeGenerate \
		--runThreadN 20 \
		--genomeDir outdir \
		--genomeFastaFiles $FASTA \
		--sjdbGTFfile $GTF

	cp -r outdir/* $WORKDIR/ 
	ls -lhAR
	EOS
) 
