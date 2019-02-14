#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='y'
date=`date '+%Y-%m-%d %H:%M:%S'`

# help message
help_message="
Merge BAM files from separate runs and generate metrics for merged alignments.

usage:
    bash $(basename "$0") [-options] -i <bam,list>
required arguments:
    -i|--input : comma separated list of BAM files to merge
    -F|--fasta : whole genome FASTA
optional arguments:
    -n|--name : name prefix for output files (default = extracted from first bam)
    -o|--outdir : outut directory (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --depend : list of job dependencies to pass to PBS job script
    --check : whether to check files prior to running job [y|n] (default = y)
additional_info:
    # all paths should be relative to working directory
    # check and depend options useful for job scheduling
      e.g. --depend afterok:123456,afterok:123457
    # log/qc output directories inherit from outdir unless specified

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--input)
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
if [[ -z ${bam_list:-} ]]; then
    printf "\nERROR: --input argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
fi
    
# split bam list to array
IFS="," read -r -a bam_array <<< "$bam_list"
    
# check files
if [[ $check = 'y' ]]; then
    if [[ ! -r $workdir/$fasta ]]; then 
        printf "\nERROR: FASTA file not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    fi
    for bam in ${bam_array[@]}; do
        if [[ ! -r $workdir/$bam ]]; then
            printf "\nERROR: BAM file not readable: %s/%s\n" $workdir $bam
            echo "$help_message"; exit 1
        fi
    done
fi

# get basenames/set merge inputs
fasta_prefix=${fasta%%.*}
fasta_base=$(basename "$fasta")
for bam in ${bam_array[@]}; do
    bam_base=$(basename "$bam")
    bam_merge_arg="${bam_merge_arg:-}I=$bam_base "
    bam_cp_arg="${bam_cp_arg:-}cp $workdir/$bam* .; "
done
bam_cp_arg=${bam_cp_arg%;*}

# if no name provided extract from first bam
if [[ -z ${name:-} ]]; then
    name=${bam_array[0]%%.*}
    name=$(basename "$name")
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
merge_command=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar MergeSamFiles"
                "${bam_merge_arg}"
                "OUTPUT=${name}.merge.bam")
mdup_command=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar MarkDuplicates"
                "INPUT=${name}.bam"
                "OUTPUT=${name}.mdup.bam"
                "METRICS_FILE=${name}.mark_duplicate_metrics"
                "VALIDATION_STRINGENCY=SILENT"
                "OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500")
sort_command=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar SortSam"
                "INPUT=${name}.mdup.bam"
                "OUTPUT=${name}.sorted.bam"
                "SORT_ORDER=coordinate"
                "CREATE_INDEX=false"
                "CREATE_MD5_FILE=false")
fixtags_command=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar SetNmAndUqTags"
                "INPUT=${name}.sorted.bam"
                "OUTPUT=${name}.bam"
                "REFERENCE_SEQUENCE=${fasta_base}"
                "CREATE_INDEX=true"
                "CREATE_MD5_FILE=true")
casm_command=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar CollectAlignmentSummaryMetrics"
                "R=${fasta_base}"
                "I=${name}.bam"
                "O=${name}.alignment_summary_metrics")

# set log file names
scr_name=$(basename "$0" .sh)
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
		${merge_command[@]} &>> $out_log
		mv ${name}.merge.bam ${name}.bam
		${mdup_command[@]} &>> $out_log
		cp ${name}.mark_duplicate_metrics $workdir/$qcdir &>> $out_log

		# sort and fix tags
		printf "\nSorting and fixing tags\n" &>> $out_log
		${sort_command[@]} &>> $out_log
		${fixtags_command[@]} &>> $out_log

		# alignment metrics
		printf "\nCollecting alignment metrics\n" >> $out_log
		${casm_command[@]} &>> $out_log
		cp ${name}.alignment_summary_metrics $workdir/$qcdir &>> $out_log

		# copy final bam to outdir
		cp ${name}.bam* $workdir/$outdir &>> $out_log
		cp ${name}.bai $workdir/$outdir &>> $out_log

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
