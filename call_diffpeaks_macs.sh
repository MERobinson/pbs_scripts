#!/bin/bash

# default arg
WORKDIR=$PWD
PEAKDIR=peaks
BEDGDIR=bedgraph
OUTDIR=diff_peaks
LOGDIR=logs

# help message
USAGE="

usage:
	bash $(basename "$0") [-wbnh] -c1 <COND1> -c2 <COND2>
purpose:
	Simple wrapper to call differential peaks with macs2 bdgdiff
required arguments:
        -c1|--cond1 : prefix of macs2 call peak output files for condition 1
	-c2|--cond2 : prefix of macs2 call peak output files for condition 2
optional arguments:
	-w|--workdir : path of working directory to use (default = pwd)
	-b|--bedgdir : directory in which to find macs2 callpeak bedgraph files (default = bedgraph)
	-p|--peakdir : directory in which to find macs2 callpeak peak .xls files (default = peaks)
	-o|--outdir : output directory for results (default = diff_peaks)
	-l|--logdir : output directory for run logs (default = logs)
	-d1|--depth1 : minimum filtered tag count for cond1 (extracted from .xls)
	-d2|--depth2 : minimum filtered tag count for cond2 (extracted from .xls)
	-n|--name : name prefix for output files (deafult = <COND1>_vs_<COND2>)
	-h|--help : print current help message
example:
        bash $(basename "$0") -t K562_H3K27ac -c K562_input
additional info:
	> requires macs2 callpeak output files
	> --bedgdir should be an existing folder within --workdir
	> --bedgdir should contain macs2 callpeak output pileup and lambda bedgraph files
	> cond1 and cond2 prefixes should match names of macs2 callpeak output files
	> will attempt to extract read counts from macs2 callpeak .xls files by default

"

# parse command line arg
while [[ $# -gt 1 ]]; do
        key=$1
	case $key in
		-c1|--cond1)
		COND1=$2
		shift
		;;
                -c2|--cond2)
                COND2=$2
                shift
                ;;
                -w|--workdir)
                WORKDIR=$2
                shift
                ;;
                -b|--bedgdir)
                BEDGDIR=$2
                shift
                ;;
		-o|--outdir)
		OUTDIR=$2
		shift
		;;
		-l|--logdir)
		LOGDIR=$2
		shift
		;;
		-d1|--depth1)
		DEPTH1=$2
		shift
		;;
                -d2|--depth2)
                DEPTH2=$2
                shift
                ;;
                -n|--name)
                NAME=$2
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
if [[ ! -d "$WORKDIR/$BEDGDIR" ]]; then
	printf "\nERROR: Input bedgraph directory not found: $WORKDIR/$BEDGDIR \n"
	echo "$USAGE"; exit 1
fi

# check required arguments
if [[ -z "$COND1" ]] || [[ -z "$COND2" ]]; then
        printf "\nERROR: Requires both test and control prefixes to be specified\n"
        echo "$USAGE"; exit 1
elif [[ ! -f $WORKDIR/$BEDGDIR/${COND1}_treat_pileup.bdg ]]; then
	printf "\nERROR: Test file not found: $WORKDIR/$BEDGDIR/${COND1}_treat_pileup.bdg \n"
	echo "$USAGE"; exit 1
elif [[ ! -f $WORKDIR/$BEDGDIR/${COND1}_control_lambda.bdg ]]; then
	printf "\nERROR: Test file not found: $WORKDIR/$BEDGDIR/${COND1}_control_lambda.bdg \n"
	echo "$USAGE"; exit 1
elif [[ ! -f $WORKDIR/$BEDGDIR/${COND2}_treat_pileup.bdg ]]; then
        printf "\nERROR: Test file not found: $WORKDIR/$BEDGDIR/${COND2}_treat_pileup.bdg \n"
        echo "$USAGE"; exit 1
elif [[ ! -f $WORKDIR/$BEDGDIR/${COND2}_control_lambda.bdg ]]; then
        printf "\nERROR: Test file not found: $WORKDIR/$BEDGDIR/${COND2}_control_lambda.bdg \n"
        echo "$USAGE"; exit 1
fi

# get name if req
if [[ -z "$NAME" ]]; then NAME=${COND1}_vs_${COND2}; fi

# if no depth provided, extract from peaks file
GREP="^# tags after filtering in (treatment|control)"
if [[ -z "$DEPTH1" ]]; then
	DEPTHS=($(grep -E "$GREP" "$WORKDIR/$PEAKDIR/${COND1}_peaks.xls" | grep -Eo "[0-9]+"))
	if [[ ${DEPTHS[0]} > ${DEPTHS[1]} ]]; then DEPTH1=${DEPTHS[1]}; else DEPTH1=${DEPTHS[0]}; fi
	echo "WARNING: No read depth provided for cond1, using extracted value: $DEPTH1"
fi
if [[ -z "$DEPTH2" ]]; then
	DEPTHS=($(grep -E "$GREP" "$WORKDIR/$PEAKDIR/${COND2}_peaks.xls" | grep -Eo "[0-9]+"))
	if [[ ${DEPTHS[0]} > ${DEPTHS[1]} ]]; then DEPTH2=${DEPTHS[1]}; else DEPTH2=${DEPTHS[0]}; fi
	echo "WARNING: No read depth provided for cond2, using extracted value: $DEPTH2"
fi

# create output dirs
mkdir -p $WORKDIR/$LOGDIR
mkdir -p $WORKDIR/$OUTDIR

# run job
MACS_JOB=$(cat <<- EOS | qsub -N $NAME.MACS -
	#!/bin/bash
	#PBS -l walltime=10:00:00
	#PBS -l select=1:mem=10gb:ncpus=1
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $WORKDIR/$LOGDIR/$NAME.macs2.diffpeaks.log

	module load samtools/1.2
	module load macs/2.1.0

	cp $WORKDIR/$BEDGDIR/$COND1* .
	cp $WORKDIR/$BEDGDIR/$COND2* .

	macs2 bdgdiff \
		--t1 ${COND1}_treat_pileup.bdg \
		--c1 ${COND1}_control_lambda.bdg \
		--t2 ${COND2}_treat_pileup.bdg \
		--c2 ${COND2}_control_lambda.bdg \
		--d1 $DEPTH1 \
		--d2 $DEPTH2 \
		--o-prefix $NAME

	cp $NAME* $WORKDIR/$OUTDIR/	

	ls -lhAR	
	EOS
)
