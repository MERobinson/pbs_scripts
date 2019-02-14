#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default args
workdir=$PWD
outdir=''
check='yes'

# help message
help_message="
Wrapper to run Picard's CollectAlignmentSummaryMetrics.

usage:
    $(basename "$0") [-options] -b <bam,files>
required arguments:
    -b|--bam : BAM file to analyse
    -F|--fasta : reference FASTA
    -BI|--bait_intervals : bait intervals file
    -TI|--target_intervals : targets interval file
optional arguments:
    -o|--outdir : output directory for metric files (default = pwd)
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output files (default = bam filename)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling

"

# parse arguments
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam=$2
            shift
            ;;
		-F|--fasta)
            fasta=$2
            shift
            ;;
        -BI|--bait_intervals)
            bait_intervals=$2
            shift
            ;;
        -TI|--target_intervals)
            target_intervals=$2
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
        -n|--name)
            name=$2
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
            printf "ERROR: Undefined argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
	shift
done

# set logdir & create outdirs
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
mkdir -p $workdir/$logdir
mkdir -p $workdir/$outdir

# check req arg
if [[ -z ${bam:-} ]]; then
	printf "\nERROR: --bam argument required\n"
	echo "$help_message"; exit 1
elif [[ -z ${bait_intervals:-} ]]; then
    printf "\nERROR: --bait_intervals argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${target_intervals:-} ]]; then
    printf "\nERROR: --target_intervals argument required\n"
    echo "$help_message"; exit 1
fi
if [[ ${check:-} = yes ]]; then
    if [[ ! -r ${workdir}/${bam} ]]; then
        printf "\nERROR: BAM file is not readable: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    elif [[ ! -r ${workdir}/${bait_intervals} ]]; then
        printf "\nERROR: Bait intervals file is not readable: %s/%s\n" $workdir $bait_intervals
        echo "$help_message"; exit 1
    elif [[ ! -r ${workdir}/${target_intervals} ]]; then
        printf "\nERROR: Target intervals file is not readable: %s/%s\n" $workdir $target_intervals
        echo "$help_message"; exit 1
    elif [[ -n ${fasta:-} ]] && [[ ! -r ${workdir}/${fasta:-} ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    fi
fi

# set name if not given
if [[ -z ${name:-} ]]; then
    name=$(basename ${bam})
    name=${name%%.*}
fi

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# set commands
hs_cmd=("java -Xmx16G -jar /apps/picard/2.6.0/picard.jar CollectHsMetrics"
        "I=$(basename ${bam})"
        "O=${name}.hs_metrics"
        "TARGET_INTERVALS=$(basename ${target_intervals})"
        "BAIT_INTERVALS=$(basename ${bait_intervals})")
if [[ -n ${fasta:-} ]]; then
    hs_cmd+=("R=$(basename ${fasta})")
    fasta_cp="cp $workdir/${fasta%.*}* . &>> $out_log"
fi

# run job
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:ncpus=1:mem=16gb
		#PBS -j oe
		#PBS -q med-bio
		#PBS -N HS_met_$name
		#PBS -o $std_log
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# load modules
		module load java/jdk-8u66 &>> $out_log
		module load picard/2.6.0 &>> $out_log
		module load samtools/1.2 &>> $out_log
	
		# copy to scratch
		cp -r $workdir/$bam* . &>> $out_log
		cp -r $workdir/$bait_intervals . &>> $out_log
		cp -r $workdir/$target_intervals . &>> $out_log
		${fasta_cp:-}

		# run
		printf "\nCollecting HS metrics:\n" >> $out_log
		${hs_cmd[@]} &>> $out_log

		# copy metrics to outdir
		cp $name.hs_metrics $workdir/$outdir/ &>> $out_log
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR >> $out_log
		ls -lhAR
		cp $out_log $workdir/$logdir/
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
