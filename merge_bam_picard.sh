#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'
date=`date '+%Y-%m-%d %H:%M:%S'`
scr_name=$(basename "$0")

# help message
help_message="
Merge BAM files from separate runs and generate metrics for merged alignments.

usage:
    bash $scr_name [-options] -b <bam,list>
required arguments:
    -b|--bam_list : comma separated list of input BAM files to merge
    -F|--fasta : whole genome FASTA
optional arguments:
    -n|--name : name prefix for output files (default = extracted from first bam)
    -o|--outdir : outut directory for merged BAM (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --depend : comma-sep list of job dependencies - format afterok:<jobid>,afterok:<jobid>
    --check : whether to check files prior to running job [yes,no] (default = yes)
additional_info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # bam/log/qc output directories inherit from --outdir unless specified

"

# parse command line arg
while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
        -b|--bam_list)
            bam_list=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -n|--name)
            name=$2
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
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required arg provided
if [[ -z "${bam_list:-}" ]]; then
    printf "\nERROR: --bam_list argument required\n"
    echo "$help_message"; exit 1
else
    IFS="," read -r -a bam_array <<< "$bam_list"
    for bam in ${bam_array[@]}; do
        if [[ $check = "yes" ]] && [[ -z $workdir/$bam ]]; then
            printf "\nERROR: BAM file not readable: %s/%s\n" $workdir $bam
            echo "$help_message"; exit 1
        else
            bam_base=$(basename "$bam")
            bam_mdup_arg="${bam_mdup_arg:-}I=$bam_base "
            bam_cp_arg="${bam_cp_arg:-}cp $workdir/$bam* .; "
        fi
    done
    bam_cp_arg=${bam_cp_arg%;*}
fi
if [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
else
    if [[ $check = 'yes' ]] && [[ ! -r $workdir/$fasta ]]; then 
        printf "\nERROR: BAM file not readable: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    else
        fasta_prefix=${fasta%%.*}
        fasta_base=$(basename "$fasta")
    fi
fi

# if no name provided extract from first bam
if [[ -z "$name" ]]; then
    name=${bam_array[0]%%.*}".merged"
fi

# set output directories
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
if [[ -z "${qcdir:-}" ]]; then
    qcdir=$outdir
fi

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir

# set commands
merge_cmd=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar MergeSamFiles"
            "${megre_sam_input}"
            "OUTPUT=${name}.merge.bam")
mdup_cmd=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar MarkDuplicates"
            "INPUT=${name}.bam"
            "OUTPUT=${name}.mdup.bam" 
            "METRICS_FILE=${name}.mark_dup_metrics.txt"
            "CREATE_INDEX=true")
mdup2_command=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar MarkDuplicates"
                "$name.bam O=tmp M=$name.mark_duplicate_metrics CREATE_INDEX=true")
casm_command=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar CollectAlignmentSummaryMetrics"
                "R=$fasta_base I=$name.bam O=$name.alignment_summary_metrics")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$workdir/$logdir/$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=18gb:ncpus=12
		#PBS -j oe
		#PBS -N $name.merge
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# load modules
		module load samtools/1.2
		module load java/jdk-8u66
		module load picard/2.6.0

		# copy accross bam and fasta
		cp $workdir/$fasta_prefix* . &>> $out_log
		$bam_cp_arg &>> $out_log

		# merge & mark dup
		printf "\nMerging BAM and re-marking duplicates\n" &>> $out_log 
		${mdup_command[@]} &>> $out_log
		cp $name.mark_duplicate_metrics $workdir/$qcdir/ &>> $out_log

		# alignment metrics
		printf "\nCollecting alignment metrics\n" >> $out_log
		${casm_command[@]} &>> $out_log
		cp $name.alignment_summary_metrics $workdir/$qcdir/ &>> $out_log

		# copy final bam to outdir
		samtools index $name.bam &>> $out_log
		cp $name.bam* $workdir/$outdir/ >> $out_log

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` &>> $out_log
		ls -lhAR &>> $out_log
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
