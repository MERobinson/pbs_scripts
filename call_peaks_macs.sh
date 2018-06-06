#!/bin/bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
genome='hs'
check='yes'

# help message
help_message="
Wrapper to call peaks with MACS2 

usage:
    bash $(basename "$0") [-options] -t <BAM>
required arguments:
    -t|--test : test sample [BAM]
    -F|--fasta : whole genome fasta file
optional arguments:
    -c|--control : control sample (aligned bam)
    -o|--outdir : output directory for macs files (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -q|--qcdir : output directory for qc files (default = ---outdir)
    -g|--genome : genome size [hs,mm,<numeric>] (default = hs)
    -n|--name : name prefix for output files (deafult = test sample name)
    -e|--extsize : run macs2 with provided extension (default = unset)
    --check : whether to check files prior to running [yes,no] (default = yes)
    --depend : list of dependencies for PBS script (e.g. afterok:012345,afterok:012346)
additional info:
    # all paths should be given relative to working directory
    # --outdir will be created, if not pre-existing
    # if --control is provided, will be used as background for peak calling,
      otherwise no control sample is used in macs
    # if not provided, name is extracted from test filename up to first period (.)

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -t|--test)
            test=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -c|--control)
            control=$2
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
        --check)
            check=$2
            shift
            ;;
        -d|--depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        *)
            printf "\nERROR: Unrecognised argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
	shift
done

# set logdir and qcdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
if [[ -z ${qcdir:-} ]]; then
    qcdir=$outdir
fi

# check required arguments
if [[ -z "${test:-}" ]]; then
    printf "\nERROR: No test BAM provided\n"
    echo "$help_message"; exit 1
elif [[ -z "${fasta:-}" ]]; then
    printf "\nERROR: No FASTA file provided\n"
    echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r "$workdir/$test" ]]; then
        printf "\nERROR: Input bam cannot be read: %s/%s\n" $workdir $test
        echo "$help_message"; exit 1
    elif [[ ! -z ${control:-} ]] & [[ ! -r $workdir/${control:-} ]]; then
        printf "\nERROR: Control bam cannot be read: %s/%s\n" $workdir ${control:-}
        echo "$help_message"; exit 1
    elif [[ ! -r "$workdir/$fasta" ]]; then
        printf "\nERROR: FASTA cannot be read: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    fi
fi

# set name if not provided
if [[ -z ${name:-} ]]; then
    name=$(basename $test)
    name=${name%%.*}
fi

# set basename
test_base=$(basename "$test")

# get canonical chr from fasta
chr_list=$(cat $fasta | grep -Eo "^>chr[0-9MXY]+\b" | \
           grep -Eo "chr[0-9XY]+"| tr "\n" " ")

# set control arg
if [[ ! -z "${control:-}" ]]; then
    control_base=$(basename $control)
    control_macs="-c $control_base"
    control_cp="cp $workdir/$control* ."
    control_filter="samtools view -b -F 1024 -q 10 $control_base $chr_list > tmp.bam"
    control_filter="${control_filter}; mv tmp.bam $control_base; samtools index $control_base"
fi

# set extsize arg
if [[ ! -z "${extsize:-}" ]]; then
	extsize_macs="--nomodel --extsize $extsize"
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# run job
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=10gb:ncpus=4
		#PBS -j oe
		#PBS -N $name.macs
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.call_peaks_macs.log
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# load modules
		module load samtools/1.2 &>> $out_log
		module load macs/2.1.0 &>> $out_log

		# copy to scratch
		cp $workdir/$test* . &>> $out_log
		${control_cp:-} 

		# filter bams to remove dup/lowqual/scaffold
		printf "\nFiltering input bams:\n" >> $out_log
		samtools view -b -F 1024 -q 10 $test_base $chr_list > tmp.bam
		mv tmp.bam $test_base; samtools index $test_base
		${control_filter:-}

		# run predictd for qc
		printf "\nPredicting insert size:\n" >> $out_log
		macs2 predictd \
			-i $test_base \
			-g $genome \
			--outdir . \
			--rfile $name.predictd.R &>> $out_log
		Rscript $name.predictd.R &>> $out_log
		cp $name.predictd.R* $workdir/$qcdir/ &>> $out_log

		# call peaks
		printf "\nCalling peaks:\n" >> $out_log
		macs2 callpeak \
			-t $test_base \
			-g $genome \
			--outdir . \
			-n $name \
			--bdg \
			${control_macs:-} \
			${extsize_macs:-} &>> $out_log

		cp ${name}_peaks.narrowPeak $workdir/$outdir/ &>> $out_log
		cp ${name}_peaks.xls $workdir/$outdir/ &>> $out_log
		cp ${name}_treat_pileup.bdg $workdir/$outdir/ &>> $out_log
		cp ${name}_control_lambda.bdg $workdir/$outdir/ &>> $out_log

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR &>> $out_log
		cp $out_log $workdir/$logdir/
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
