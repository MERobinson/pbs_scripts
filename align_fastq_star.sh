#!/bin/bash

# default arg
WORKDIR=$PWD
INDIR=fastq
GENOME=hg38
RESDIR=resources
TRACKDIR=tracks
BAMDIR=bam
QCDIR=qc
LOGDIR=logs
FASTA=genome.fa
INDEX=index
OUTDIR=star

# help message
USAGE="
usage:
	bash $(basename "$0") [-wrtbqlgnh] -f <FASTQ>
purpose:
	# Simple wrapper to run STAR alignment.
required arguments:
    -f|--fastq : comma separated list of FASTQ files
optional arguments:
	-m|--mate : comma separated list of read 2 FASTQ if paired end
	-w|--workdir : working directory - used as base dir for all input/output (default = pwd)
	-i|--indir : input directory containing FASTQ files (default = fastq)
	-o|--outdir : output directory to write STAR files to (default = star)
	-r|--resdir : directory containing index and ref FASTA resources (default = resources) 
	-t|--trackdir : output directory for genome tracks (default = tracks)
        -b|--bamdir : output directory for bam files (default = bam)
	-q|--qcdir : output directory for QC files (default = qc)
	-l|--logdir : output directory for log files (default = logs)
	-g|--genome : genome version (default = hg38)
        -n|--names : comma separated list of names for output files (default = FASTQ filename)
        -h|--help : print current help message
example:
        bash $(basename "$0") -g hg38 -f REH_H3K27ac_ChIPseq.fastq.gz
additional info:
	# if mates included, should be same length as fastq - leave empty slots for SE if necessary
	# e.g. -f sample1_R1.fq,sample2_R1.fq,sample3_R1.fq -m ,,sample3_R2.fq 

"

# parse command line arg
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-f|--fastq)
		FASTQ_LIST=$2
		shift
		;;
		-m|--mate)
		MATE_LIST=$2
		shift
		;;
		-w|--workdir)
		WORKDIR=$2
		shift
		;;
		-i|--indir)
		INDIR=$2
		shift
		;;
		-o|--outdir)
		OUTDIR=$2
		shift
		;;
		-r|--resdir)
		RESDIR=$2
		shift
		;;
		-t|--trackdir)
		TRACKDIR=$2
		shift
		;;
		-b|--bamdir)
		BAMDIR=$2
		shift
		;;
		-l|--logdir)
		LOGDIR=$2
		shift
		;;
		-q|--qcdir)
		QCDIR=$2
		shift
		;;
		-g|--genome)
		GENOME=$2
		shift
		;;
		-n|--names)
		NAMES_LIST=$2
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

# check FASTQ
if [[ -z "$FASTQ_LIST" ]]; then
	printf "\nERROR: FASTQ files required\n"
	echo "$USAGE"; exit 1
else
	IFS=',' read -r -a FASTQ_ARRAY <<< "$FASTQ_LIST"
fi

# check mates
if [[ ! -z "$MATE_LIST" ]]; then
	IFS=',' read -r -a MATE_ARRAY <<< "$MATE_LIST"
fi

# check names
if [[ -z "$NAMES_LIST" ]]; then
	NAMES_ARRAY=(${FASTQ_ARRAY[@]})
else
	IFS=',' read -r -a NAMES_ARRAY <<< "$NAMES_LIST"
fi

# strip fasta extension
FASTA_BASE=${FASTA%%.*}

# create required dirs
mkdir -p $WORKDIR/$LOGDIR
mkdir -p $WORKDIR/$QCDIR/fastqc
mkdir -p $WORKDIR/$QCDIR/metrics
mkdir -p $WORKDIR/$BAMDIR
mkdir -p $WORKDIR/$TRACKDIR
mkdir -p $WORKDIR/$OUTDIR

# for each SRA number, submit PBS job
for IDX in ${!FASTQ_ARRAY[@]}; do

	FASTQ=${FASTQ_ARRAY[$IDX]}
	MATE=${MATE_ARRAY[$IDX]}
	NAME=${NAMES_ARRAY[$IDX]}
	if [[ -z $NAMES_LIST ]]; then NAME=${NAME%%.*}; fi

	echo "FASTQ: $FASTQ, MATE: $MATE"

	# check if mate provided
	if [[ ! -z "$MATE" ]]; then
		PE_CP="cp $WORKDIR/$INDIR/$MATE ."
	fi

	JOBID=$(cat <<- EOS | qsub -N $NAME - 
		#!/bin/bash
		#PBS -l walltime=20:00:00
		#PBS -l select=1:mem=40gb:ncpus=20
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $WORKDIR/$LOGDIR/$NAME.align_fastq_star.log
		
		# load modules
		module load fastqc/0.11.2
		module load samtools/1.2
		module load java/jdk-8u66
		module load picard/2.6.0
		module load star/2.5.0
		module load intel-suite/2015
		module load gcc/5.4.0
		
		# make temp dir and copy accross index
		mkdir fastqc
		cp -rL $WORKDIR/$RESDIR/$INDEX .
		cp -L $WORKDIR/$RESDIR/$FASTA .
		cp $WORKDIR/$INDIR/$FASTQ .
		$PE_CP

		# run fastqc
		fastqc --noextract -o fastqc -t 20 $FASTQ $MATE
		cp fastqc/* $WORKDIR/$QCDIR/fastqc/
		
		# align with STAR
		STAR --runMode alignReads \
			--runThreadN 20 \
			--genomeDir $INDEX \
			--readFilesIn $FASTQ $MATE \
			--readFilesCommand zcat \
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
			--outFileNamePrefix ${NAME}_ \
			--outBAMsortingThreadN 20 \
			--outSAMtype BAM SortedByCoordinate \
			--outWigType wiggle \
			--outWigStrand Unstranded
	
		# index
		mv ${NAME}_Aligned.sortedByCoord.out.bam $NAME.$GENOME.bam
		samtools index $NAME.$GENOME.bam
		
		# copy results back
		cp $NAME.$GENOME.bam* $WORKDIR/$BAMDIR/
		cp ${NAME}_* $WORKDIR/$OUTDIR/

		# alignment metrics
		java -Xmx16g -jar /apps/picard/2.6.0/picard.jar CollectAlignmentSummaryMetrics \
			R=$FASTA \
			I=$NAME.$GENOME.bam \
			O=$NAME.alignment_summary_metrics
		cp $NAME.alignment_summary_metrics $WORKDIR/$QCDIR/metrics/
	
		## filter non-canon chr
		#cat tmp.bedGraph | grep -E "^track|^chr[0-9MXY]+\b" > $NAME.$GENOME.bedGraph
		#cp $NAME.$GENOME.bedGraph $WORKDIR/$TRACKDIR/

		ls -lhAR
		EOS
	) 
done
