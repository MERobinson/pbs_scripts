#!/bin/bash

# default arg
workdir=$PWD
outdir=macs
bamdir=bam
logdir=logs
qcdir=qc
genome=hs

# help message
help_message="

usage:
    bash $(basename "$0") [-cwbogne] -t <BAM>
purpose:
    Wrapper to call peaks with macs2
required arguments:
    -t|--test : test sample (aligned bam)
optional arguments:
    -c|--control : control sample (aligned bam)
    -w|--workdir : path of working directory for all input/output (default = pwd)
    -b|--bamdir : folder within workdir containing input bam (default = bam)
    -o|--outdir : name of output directory for macs files (default = macs)
    -l|--logdir : name of output directory for log files (default = logs)
    -q|--qcdir : name of output directory for QC files (default = qc)
    -g|--genome : genome [hs,mm,<numeric>] (default = hs)
    -n|--name : name prefix for output files (deafult = test sample name)
    -e|--extsize : run macs2 with provided extension (default = unset)
    -v|--validate : whether to check files prior to running [yes,no] (default = yes)
    -d|--depend : list of dependencies for PBS script (e.g. afterok:012345,afterok:012346)
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
        -t|--test)
        test=$2
        shift
        ;;
        -c|--control)
        control=$2
        shift
        ;;
        -w|--workdir)
        workdir=$2
        shift
        ;;
        -b|--bamdir)
        bamdir=$2
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
        -q|--qcdir)
        qcdir=$2
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
        -e|--extsize)
        extsize=$2
		shift
		;;
        -v|--validate)
        validate=$2
        shift
        ;;
        -d|--depend)
        depend=$2
        shift
        ;;
        *)
        printf "\nERROR: Illegal argument\n"
        echo "$help_message"; exit 1
        ;;
    esac
	shift
done

# check required arguments
if [[ -z "$test" ]]; then
    printf "\nERROR: No test sample provided\n"
    echo "$help_message"; exit 1
fi

if [[ $validate = yes ]]; then
    # check dir
    if [[ ! -d "$workdir/$bamdir" ]]; then
        printf "\nERROR: Input bam directory not found: $workdir/$bamdir \n"
        echo "$help_message"; exit 1
    fi

    # check bam
    if [[ ! -f $workdir/$bamdir/$test ]]; then
        printf "\nERROR: Test file not found: $workdir/$bamdir/$test \n"
        echo "$help_message"; exit 1
    fi
fi

# get name if nec
if [[ -z "$name" ]]; then
    name=${test%%.*}
fi

# set control arg
if [[ ! -z "$control" ]]; then
    control_macs_arg="-c $control"
    control_cp_arg="cp $workdir/$bamdir/$control* ."
    control_filter_arg="samtools view -b -F 1024 -q 10 $control $chr_list > tmp.bam"
    control_filter_arg="${control_filter_arg}; mv tmp.bam $control"
fi

# set depend arg
if [[ ! -z ${depend:-} ]]; then
    depend="#PBS -W depend=$depend"
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir/fragsize

# run job
jobid=$(cat <<- EOS | qsub -N $name.macs -
	#!/bin/bash
	#PBS -l walltime=10:00:00
	#PBS -l select=1:mem=10gb:ncpus=4
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $workdir/$logdir/$name.call_peaks_macs.log
	$depend

	module load samtools/1.2
	module load macs/2.1.0

	cp $workdir/$bamdir/$test* .
	$control_cp_arg

	# filter bams to remove dup/lowqual/scaffold
	samtools view -b -F 1024 -q 10 $test $chr_list > tmp.bam; mv tmp.bam $test
	$control_filter_arg

	# run predictd for qc
	macs2 predictd \
		-i $test \
		-g $genome \
		--outdir . \
		--rfile $name.predictd.R
	Rscript $name.predictd.R
	cp $name.predictd.R* $workdir/$qcdir/fragsize/

	# if extension size provided, run with no model
	if [[ ! -z "$extsize" ]]; then
		macs2 callpeak \
			-t $test \
			-g $genome \
			--outdir . \
			--nomodel \
			--extsize $extsize \
			-n $name \
			--bdg \
			$control_macs_arg
	# otherwise just run default macs2 callpeak
	else
		macs2 callpeak \
			-t $test \
			-g $genome \
			--outdir . \
			-n $name \
			--bdg \
			$control_macs_arg
	fi

	cp ${name}_peaks.narrowPeak $workdir/$outdir/
	cp ${name}_peaks.xls $workdir/$outdir/
	cp ${name}_treat_pileup.bdg $workdir/$outdir/
	cp ${name}_control_lambda.bdg $workdir/$outdir/

	ls -lhAR
	EOS
)
echo "JOBID: $jobid"
