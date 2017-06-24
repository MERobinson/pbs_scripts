#!/bin/bash

# default arg
WORKDIR=$PWD
PEAKDIR=peaks
BEDGDIR=bedgraph
BAMDIR=bam
GENOME=hs

# help message
USAGE="

usage:
	bash $(basename "$0") [-wiognh] -t <BAM>
purpose:
	Simple wrapper to call peaks with macs2
required arguments:
        -t|--test : test sample (aligned bam)
optional arguments:
	-c|--control : control sample (aligned bam)
	-w|--workdir : path of working directory to use (default = pwd)
	-b|--bamdir : directory in which to find input  bam files (default = bam)
	-bfg|--bedgdir : name of output directory for bedgraph (default = bedgraph)
	-p|--peakdir : name of output directory for peaks (default = peaks)
	-g|--genome : genome [hs,mm] (default = hs)
	-n|--name : name prefix for output files (deafult = name of test sample)
	--extsize : run macs2 with provided extension (default = unset)
	-h|--help : print current help message
example:
        bash $(basename "$0") -t K562_H3K27ac.bam -c K562_input.bam -g mm
additional info:
	> --bamdir should be an existing folder within --workdir
	> --outdir will be created, if not pre-existing
	> if --control is provided, will use as background for peak calling,
	  otherwise no control sample is used in macs
	> if not provided, name is extracted from test filename up to first period (.)

"

# parse command line arg
while [[ $# -gt 1 ]]; do
        key=$1
	case $key in
		-w|--workdir)
		WORKDIR=$2
		shift
		;;
                -i|--bamdir)
                BAMDIR=$2
                shift
                ;;
                -o|--outdir)
                OUTDIR=$2
                shift
                ;;
                -t|--test)
                TEST=$2
                shift
                ;;
		-c|--control)
		CONTROL=$2
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
		--extsize)
		EXTSIZE=$2
		shift
		;;
                -h|--help)
                echo "$USAGE"
                exit 1
                ;;
                *)
                printf "\nERROR: Illegal argument\n"
                echo "$USAGE"
                exit 1
                ;;
        esac
	shift
done

# check dir
if [[ ! -d "$WORKDIR/$BAMDIR" ]]; then
	printf "\nERROR: Input bam directory not found: $WORKDIR/$BAMDIR \n"
	echo "$USAGE"; exit 1
fi

# check required arguments
if [[ -z "$TEST" ]]; then
        printf "\nERROR: No test sample provided\n"
        echo "$USAGE"; exit 1
elif [[ ! -f $WORKDIR/$BAMDIR/$TEST ]]; then
	printf "\nERROR: Test file not found: $WORKDIR/$BAMDIR/$TEST \n"
	echo "$USAGE"; exit 1
fi

# get name if req
if [[ -z "$NAME" ]]; then
	NAME=${TEST%%.*}
fi

# set control arg
if [[ ! -z "$CONTROL" ]]; then
	CONTROL_ARG="-c $CONTROL"
	CONTROL_CP="cp $WORKDIR/$BAMDIR/$CONTROL* ."
	CONTROL_FILT="samtools view -b -F 1024 -q 10 $CONTROL $CHR_LIST > tmp.bam; mv tmp.bam $CONTROL"
fi

# create output dirs
mkdir -p $WORKDIR/$PEAKDIR
mkdir -p $WORKDIR/$BEDGDIR
mkdir -p $WORKDIR/logs
mkdir -p $WORKDIR/qc/fragsize

# run job
MACS_JOB=$(cat <<- EOS | qsub -N $NAME.$CHR.MT2 -
	#!/bin/bash
	#PBS -l walltime=10:00:00
	#PBS -l select=1:mem=10gb:ncpus=1
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $WORKDIR/logs/$NAME.macs2.callpeaks.log

	module load samtools/1.2
	module load macs/2.1.0

	cp $WORKDIR/$BAMDIR/$TEST* .
	$CONTROL_CP

	# filter bams to remove dup/lowqual/scaffold
	samtools view -b -F 1024 -q 10 $TEST $CHR_LIST > tmp.bam; mv tmp.bam $TEST
	$CONTROL_FILT

	# run predictd for qc
	macs2 predictd \
		-i $TEST \
		-g $GENOME \
		--outdir . \
		--rfile $NAME.predictd.R
	Rscript $NAME.predictd.R
	cp $NAME.predictd.R* $WORKDIR/qc/fragsize/

	# if extension size provided, run with no model
	if [[ ! -z "$EXTSIZE" ]]; then
		macs2 callpeak \
			-t $TEST \
			-g $GENOME \
			--outdir . \
			--nomodel \
			--extsize $EXTSIZE \
			-n $NAME \
			--bdg \
			$CONTROL_ARG
	# otherwise just run default macs2 callpeak
	else
		macs2 callpeak \
			-t $TEST \
			-g $GENOME \
			--outdir . \
			-n $NAME \
			--bdg \
			$CONTROL_ARG
	fi

	ls -lhAR
	
	cp ${NAME}_peaks.narrowPeak $WORKDIR/$PEAKDIR/
	cp ${NAME}_peaks.xls $WORKDIR/$PEAKDIR/
	cp ${NAME}_treat_pileup.bdg $WORKDIR/$BEDGDIR/
	cp ${NAME}_control_lambda.bdg $WORKDIR/$BEDGDIR/

	EOS
)
