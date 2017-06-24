#!/bin/bash

# default arg
WORKDIR=$PWD
GENOME=hg38
RESDIR=resources
TRACKDIR=tracks
BAMDIR=bam
LOGDIR=logs
FASTA=genome.fa
INDEX=genome.fa
SRADIR=sra
EXTEND=200

# help message
USAGE="bash $(basename "$0") [-wtblgneh] -s <SRA_accession>

purpose:
-- fetch SRA record, extract FASTQ and align data per sample accession number listed.
required arguments:
        -s|--sra : comma separated list of SRA run accession numbers
optional arguments:
	-w|--workdir : working directory - used as base dir for all input/output (default = pwd)
	-r|--resdir : directory containing [links to] index and ref FASTA resources (default = resources) 
        -t|--trackdir : output directory for genome tracks (default = tracks)
        -b|--bamdir : output directory for bam files (default = bam)
	-l|--logdir : output directory for log files (default = logs)
	-g|--genome : genome version (default = hg38)
        -n|--names : comma separated list of names for output files (default = SRA accession)
	-e|--extend : bp extension for generating coverage track (default = 300)
        -h|--help : print current help message
example:
        bash $(basename "$0") -g hs -s SRR3330011,SRR3330012 -n REH_H3K27ac_ChIPseq,REH_input_ChIPseq"

# parse command line arg
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-s|--sra)
		SRA_LIST=$2
		shift
		;;
		-w|--workdir)
		WORKDIR=$2
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
		-g|--genome)
		GENOME=$2
		shift
		;;
		-n|--names)
		NAMES_LIST=$2
		shift
		;;
		-e|--extend)
		EXTEND=$2
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

# check SRA accesion number
if [[ -z "$SRA_LIST" ]]; then
	echo "ERROR: SRA accession number required."; echo "$USAGE"; exit 1
else
	IFS=',' read -r -a SRA_ARRAY <<< "$SRA_LIST"
fi

# check sample names
if [[ -z "$NAMES_LIST" ]]; then
	NAMES_ARRAY=(${SRA_ARRAY[@]})
        echo "WARNING: No name argument given, using accession numbers."
else
	IFS=',' read -r -a NAMES_ARRAY <<< "$NAMES_LIST"
fi

# strip fasta extension
FASTA_BASE=${FASTA%%.*}

# create required dirs
mkdir -p $WORKDIR/$LOGDIR
mkdir -p $WORKDIR/qc/fastqc
mkdir -p $WORKDIR/qc/metrics
mkdir -p $WORKDIR/$BAMDIR
mkdir -p $WORKDIR/$TRACKDIR

# for each SRA number, submit PBS job
for IDX in ${!SRA_ARRAY[@]}; do

	SRA=${SRA_ARRAY[$IDX]}
	NAME=${NAMES_ARRAY[$IDX]}

	# fetch SRA info and check if SE or PE
	SRA_INFO=$(wget -qO- 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term='$SRA | \
		   grep $SRA)
	IFS="," read -r -a SRA_INFO <<< "$SRA_INFO"
	LIBTYPE=${SRA_INFO[15]}

	if [[ ! $LIBTYPE = SINGLE ]]; then
		PE_FQ2SAM_ARG="FASTQ2=${SRA}_2.fastq.gz"
		PE_BWAMEM_ARG="-p"
	fi

	JOBID=$(cat <<- EOS | qsub -N $NAME - 
		#!/bin/bash
		#PBS -l walltime=60:00:00
		#PBS -l select=1:mem=40gb:ncpus=20
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $WORKDIR/$LOGDIR/fetch_align_$SRA.log.txt
		
		# load modules
		module load sra-toolkit/2.8.1-3
		module load fastqc/0.11.2
		module load samtools/1.2
		module load java/jdk-8u66
		module load picard/2.6.0
		module load bio-bwa/0.7.10
		module load anaconda/2.4.1
		module load deeptools/2.4.2
		
		# make temp dir and copy accross index
		mkdir fastqc
		mkdir tmp
		cp -L $WORKDIR/$RESDIR/$INDEX* .
		cp -L $WORKDIR/$RESDIR/$FASTA_BASE* .

		# prefetch SRR and extract FASTQ
		if [[ ! -e $WORKDIR/$SRADIR/$SRA.sra ]]; then
			prefetch --max-size 100G $SRA
		fi
		cp $WORKDIR/$SRADIR/$SRA.sra .
		fastq-dump -v -I -W -B --skip-technical --outdir . \
		--gzip --split-files $SRA.sra
	
		# run fastqc
		fastqc --noextract -o fastqc -t 20 *.fastq.gz
		cp -r fastqc/* $WORKDIR/qc/fastqc/
		
		# covert fastq to ubam
		java -jar -Xmx32G -jar /apps/picard/2.6.0/picard.jar FastqToSam \
			FASTQ=${SRA}_1.fastq.gz \
			OUTPUT=$NAME.unaligned.bam \
			SAMPLE_NAME=$NAME $PE_FQ2SAM_ARG

		# mark adapters
		java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkIlluminaAdapters \
			I=$NAME.unaligned.bam \
			O=$NAME.markadapters.bam \
			M=$NAME.mark_adapters_metrics
		cp $NAME.mark_adapters_metrics $WORKDIR/qc/metrics/

		# convert uBAM to interleaved fastq
		java -Xmx32G -jar /apps/picard/2.6.0/picard.jar SamToFastq \
			I=$NAME.markadapters.bam \
			FASTQ=$NAME.fq \
			CLIPPING_ATTRIBUTE=XT \
			CLIPPING_ACTION=2 \
			INTERLEAVE=true \
			NON_PF=true

		# align to genome with BWA
		bwa mem -M -t 20 $PE_BWAMEM_ARG \
			$INDEX \
			$NAME.fq \
			>  $NAME.aligned.sam
	
		# merge uBAM and aligned
		java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MergeBamAlignment \
			R=$FASTA \
			UNMAPPED_BAM=$NAME.unaligned.bam \
			ALIGNED=$NAME.aligned.sam \
			O=$NAME.merge.bam \
			CREATE_INDEX=true \
			ADD_MATE_CIGAR=true \
			CLIP_ADAPTERS=true \
			CLIP_OVERLAPPING_READS=true \
			INCLUDE_SECONDARY_ALIGNMENTS=false \
			MAX_INSERTIONS_OR_DELETIONS=-1 \
			PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
			ATTRIBUTES_TO_RETAIN=XS

		# mark dup
		java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkDuplicates \
			I=$NAME.merge.bam \
			O=$NAME.$GENOME.bam \
			M=$NAME.mark_duplicate_metrics \
			CREATE_INDEX=true
		cp $NAME.mark_duplicate_metrics $WORKDIR/qc/metrics/

		# alignment metrics
		java -Xmx16g -jar /apps/picard/2.6.0/picard.jar CollectAlignmentSummaryMetrics \
			R=$FASTA \
			I=$NAME.$GENOME.bam \
			O=$NAME.alignment_summary_metrics
		cp $NAME.alignment_summary_metrics $WORKDIR/qc/metrics/
	
		# copy final bam to outdir
		cp $NAME.$GENOME.bam* $WORKDIR/$BAMDIR/
		
		# filter chromosomes to canonical


		# generate coverage track
		source activate deeptools-2.4.2
		bamCoverage \
			--extendReads $EXTEND \
			--binSize 200 \
			--normalizeUsingRPKM \
			--numberOfProcessors 20 \
			--outFileFormat "bedgraph" \
			--bam $NAME.$GENOME.bam \
			--outFileName $NAME.$GENOME.bedGraph
		echo "track type=bedGraph name=$NAME visibility=2 windowingFunction=mean smoothingWindow=10 viewLimits=0:20 maxHeightPixels=0:75:150" | \
			cat - $NAME.$GENOME.bedGraph > tmp.bedGraph
		
		# filter non-canon chr
		cat tmp.bedGraph | grep -E "^track|^chr[0-9MXY]+\b" > $NAME.$GENOME.bedGraph
		cp $NAME.$GENOME.bedGraph $WORKDIR/$TRACKDIR/

		ls -lhAR
		EOS
	) 
done
