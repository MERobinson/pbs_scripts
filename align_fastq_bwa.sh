#!/bin/bash

# default arg
WORKDIR=$PWD
INDIR=fastq
GENOME=hg38
RESDIR=resources
TRACKDIR=tracks
BAMDIR=bam
LOGDIR=logs
FASTA=genome.fa
INDEX=genome.fa
EXTEND=200
QCDIR=qc
PL="PLATFORM=ILLUMINA"

# help message
USAGE="
usage:
	bash $(basename "$0") [-wrtblgneh] -f <FASTQ>
purpose:
	# Align reads in FASTQ format with BWA and generate alignment metrics and track.
required arguments:
        -f|--fastq : FASTQ filename of read 1 
optional arguments:
	-m|--mate : FASTQ filename of read 2 - if paired end
	-w|--workdir : working directory - used as base dir for all input/output (default = pwd)
	-i|--indir : input directory containing FASTQ files (default = fastq)
	-r|--resdir : directory containing index and ref FASTA resources (default = resources) 
        -t|--trackdir : output directory for genome tracks (default = tracks)
        -b|--bamdir : output directory for bam files (default = bam)
	-q|--qcdir : output directory for qc metrics (default = qc)
	-l|--logdir : output directory for log files (default = logs)
	-g|--genome : genome version (default = hg38)
        -n|--name : name prefix for output file (default = FASTQ filename)
	-e|--extend : bp extension for generating coverage track (default = 200)
        -h|--help : print current help message
example:
        bash $(basename "$0") -g hg38 -f REH_H3K27ac_ChIPseq.fastq.gz
additional info:
	# Additional read group and sequence information can be provided to include in the
	  BAM header
	# Possible read group data to include:
		--id : read group ID (usually flow cell + lane)
		--lb : library name (name unique to each library)
		--pu : platform unit (usually flow cell + barcode + lane)
		--pl : sequencing platform sequenced on (e.g. ILLUMINA)
		--cn : sequencing centre sequencing was performed at 
		--dt : run date (Iso8601Date)
		--pi : predicted median insert size (e.g. 200)
		--pm : platform model - further discription of platform
	# by default will assume standard illumina format and extract ID and PU from readname
	# if reads do not have standard illumina name formats make sure to specify PU and ID 

"

# parse command line arg
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-f|--fastq)
		FASTQ=$2
		shift
		;;
		-m|--mate)
		MATE=$2
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
		-n|--name)
		NAME=$2
		shift
		;;
		-e|--extend)
		EXTEND=$2
		shift
		;;
		--id)
		ID="READ_GROUP_NAME=$2"
		shift
		;;
		--cn)
		CN="SEQUENCING_CENTER=$2"
		shift
		;;
		--dt)
		DT="RUN_DATE=$2"
		shift
		;;
		--lb)
		LB="LIBRARY_NAME=$2"
		shift
		;;
		--pi)
		PI="PREDICTED_INSERT_SIZE=$2"
		shift
		;;
		--pl)
		PL="PLATFORM=$2"
		shift
		;;
		--pm)
		PM="PLATFORM_MODEL=$2"
		shift
		;;
		--pu)
		PU="PLATFORM_UNIT=$2"
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

# check indir exists
if [[ ! -e "$WORKDIR/$INDIR" ]]; then
	printf "\nERROR: input directory does not exist: %s/%s/%s\n" $WORKDIR $INDIR $FASTQ
	echo "$USAGE"; exit 1
fi

# check FASTQ argument is provided and file exists
if [[ -z "$FASTQ" ]]; then
	printf "\nERROR: FASTQ file input argument required\n"
	echo "$USAGE"; exit 1
elif [[ ! -e "$WORKDIR/$INDIR/$FASTQ" ]]; then
	printf "\nERROR: input FASTQ file does not exist: %s/%s/%s\n" $WORKDIR $INDIR $FASTQ
	echo "$USAGE"; exit 1
fi

# check if mate exists and, if PE, set PE arguments
if [[ ! -z "$MATE" ]]; then
	if [[ ! -e "$WORKDIR/$INDIR/$MATE" ]]; then
		printf "\nERROR: input mate FASTQ does not exist: %s/%s/%s\n" $WORKDIR $INDIR $FASTQ
		echo "$USAGE"; exit 1
	fi
	PE_CP="cp $WORKDIR/$INDIR/$MATE ."
	PE_FQ2SAM_ARG="FASTQ2=$MATE"
	PE_BWAMEM_ARG="-p "
fi

# if no name provided extract from FASTQ filename
if [[ -z "$NAME" ]]; then
	NAME=${FASTQ%%.*}
fi

# check resource files exist
if [[ ! -e "$WORKDIR/$RESDIR/$FASTA" ]]; then
	printf "\nERROR: FASTA file does not exist: %s/%s/%s\n" $WORKDIR $RESDIR $FASTA
        echo "$USAGE"; exit 1
fi
if [[ ! -e "$WORKDIR/$RESDIR/$INDEX" ]]; then
	printf "\nERROR: Index file does not exist: %s/%s/%s\n" $WORKDIR $RESDIR $INDEX
	echo "$USAGE"; exit 1
fi

# strip fasta extension (to copy dict)
FASTA_BASE=${FASTA%%.*}

# get list of canonical chromosomes
CHR_LIST=$(cat $WORKDIR/$RESDIR/$FASTA | grep -E "^>" | grep -Eo "chr[0-9MXY]+\b") 

# extract read group info (if not provided)
if [[ $PL = "PLATFORM=ILLUMINA" ]]; then
	READ_NAME=$(gzip -dc "$WORKDIR/$INDIR/$FASTQ" | head -n 1)
	FCID=$(echo $READ_NAME | cut -d ":" -f 3) # flow cell ID
	FCLN=$(echo $READ_NAME | cut -d ":" -f 4) # flow cell lane
	RIDX=$(echo $READ_NAME | cut -d ":" -f 10) # barcode/index
	if [[ -z $ID ]]; then ID="READ_GROUP_NAME=$FCID.$FCLN"; fi
	if [[ -z $PU ]]; then PU="PLATFORM_UNIT=$FCID.$FCLN.$RIDX"; fi
fi

# create required dirs
mkdir -p $WORKDIR/$LOGDIR
mkdir -p $WORKDIR/$QCDIR/fastqc
mkdir -p $WORKDIR/$QCDIR/metrics
mkdir -p $WORKDIR/$BAMDIR
mkdir -p $WORKDIR/$TRACKDIR

# run job
JOBID=$(cat <<- EOS | qsub -N $NAME.bwa - 
	#!/bin/bash
	#PBS -l walltime=60:00:00
	#PBS -l select=1:mem=40gb:ncpus=20
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $WORKDIR/$LOGDIR/$NAME.align_fastq_bwa.log
	
	# load modules
	module load fastqc/0.11.2
	module load samtools/1.2
	module load java/jdk-8u66
	module load picard/2.6.0
	module load bio-bwa/0.7.10
	module load anaconda/2.4.1
	module load deeptools/2.4.2
	
	# copy resource files to scratch
	mkdir fastqc
	cp -L $WORKDIR/$RESDIR/$INDEX* .
	cp -L $WORKDIR/$RESDIR/$FASTA_BASE* .
	cp $WORKDIR/$INDIR/$FASTQ .
	$PE_CP

	# run fastqc
	fastqc --noextract -o fastqc -t 20 $FASTQ $MATE
	cp -r fastqc/* $WORKDIR/$QCDIR/fastqc/
	
	# covert fastq to ubam
	java -jar -Xmx32G -jar /apps/picard/2.6.0/picard.jar FastqToSam \
		FASTQ=$FASTQ \
		OUTPUT=$NAME.unaligned.bam \
		SAMPLE_NAME=$NAME \
		$PE_FQ2SAM_ARG $ID $LB $PL $PU $PM $CN $DT $PI

	# mark adapters
	java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkIlluminaAdapters \
		I=$NAME.unaligned.bam \
		O=$NAME.markadapters.bam \
		M=$NAME.mark_adapters_metrics
	cp $NAME.mark_adapters_metrics $WORKDIR/$QCDIR/metrics/

	# convert uBAM to interleaved fastq
	java -Xmx32G -jar /apps/picard/2.6.0/picard.jar SamToFastq \
		I=$NAME.markadapters.bam \
		FASTQ=$NAME.fq \
		CLIPPING_ATTRIBUTE=XT \
		CLIPPING_ACTION=2 \
		INTERLEAVE=true \
		NON_PF=true

	# align to genome with BWA
	bwa mem -M -t 20 $PE_BWAMEM_ARG\
		$INDEX \
		$NAME.fq \
		> $NAME.aligned.sam

	# merge uBAM and aligned
	java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MergeBamAlignment \
		R=$FASTA \
		UNMAPPED_BAM=$NAME.unaligned.bam \
		ALIGNED=$NAME.aligned.sam \
		O=$NAME.merge.bam \
		CREATE_INDEX=true \
		ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		ATTRIBUTES_TO_RETAIN=XS

	# mark dup
	java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkDuplicates \
		I=$NAME.merge.bam \
		O=$NAME.$GENOME.bam \
		M=$NAME.mark_duplicate_metrics \
		CREATE_INDEX=true
	cp $NAME.mark_duplicate_metrics $WORKDIR/$QCDIR/metrics/

	# alignment metrics
	java -Xmx16g -jar /apps/picard/2.6.0/picard.jar CollectAlignmentSummaryMetrics \
		R=$FASTA \
		I=$NAME.$GENOME.bam \
		O=$NAME.alignment_summary_metrics
	cp $NAME.alignment_summary_metrics $WORKDIR/$QCDIR/metrics/

	# copy final bam to outdir
	cp $NAME.$GENOME.bam* $WORKDIR/$BAMDIR/
	
	# filter bam prior to making track
	samtools view -b -F 1024 -q 10 $NAME.$GENOME.bam $CHR_LIST > $NAME.$GENOME.filt.bam

	# generate coverage track
	source activate deeptools-2.4.2
	bamCoverage \
		--extendReads $EXTEND \
		--binSize 200 \
		--normalizeUsingRPKM \
		--ignoreForNormalization chrX \
		--numberOfProcessors 20 \
		--outFileFormat "bedgraph" \
		--bam $NAME.$GENOME.filt.bam \
		--outFileName $NAME.$GENOME.bedGraph
	echo "track type=bedGraph name=$NAME visibility=2 windowingFunction=mean smoothingWindow=10 viewLimits=0:20 maxHeightPixels=0:75:150" | \
		cat - $NAME.$GENOME.bedGraph > tmp.bedGraph
	mv tmp.bedGraph $WORKDIR/$TRACKDIR/$NAME.$GENOME.bedGraph

	ls -lhAR
	EOS
) 
