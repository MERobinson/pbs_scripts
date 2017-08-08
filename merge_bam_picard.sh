#!/bin/bash

# default arg
workdir=$PWD
bamdir=bam
resdir=resources
logdir=logs
qcdir=qc
fasta=genome.fa
validate=yes

# help message
help_message="
usage:
    bash $(basename "$0") [-wrtblgneh] -f <FASTQ>
purpose:
    # Merge BAM files from separate runs and generate metrics for merged alignments.
required arguments:
    -i|--bamlist : comma separated list of input BAM files to merge
optional arguments:
    -w|--workdir : working directory - used as base dir for all input/output (default = pwd)
    -b|--bamdir : directory within workdir containing BAM files (default = bam)
    -r|--resdir : directory within workdir containing resource files (default = resources)
    -q|--qcdir : output directory for qc metrics (default = qc)
    -l|--logdir : output directory for log files (default = logs)
    -n|--name : name prefix for output file (default = FASTQ filename)
    -d|--depend : comma-sep list of job dependencies - format afterok:<jobid>,afterok:<jobid>
    -f|--fasta : name of FASTA reference genome file within resdir (default = genome.fa)
    -v|--validate : whether to check files prior to running job [yes,no] (default = yes)
example:
    bash $(basename "$0") -g hg38 -i REH_H3K27ac_run1.bam,REH_H3K27ac_run2.bam
additional_info:
    # validate option useful for scheduling future jobs when files don't exist yet
"

# parse command line arg
while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
        -i|--bamlist)
        bam_list=$2
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
        -r|--resdir)
        resdir=$2
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
        -n|--name)
        name=$2
        shift
        ;;
        -f|--fasta)
        fasta=$2
        shift
        ;;
        -d|--depend)
        depend_list=$2
        shift
        ;;
        -v|--validate)
        validate=$2
        shift
        ;;
        *)
        echo "Error: Illegal argument"
        echo "$help_message"
        exit 1
        ;;
    esac
    shift
done

# check BAM arg provided
if [[ -z "$bam_list" ]]; then
    printf "\nERROR: -i|--bamlist argument required\n"
    echo "$help_message"; exit 1
else
    IFS="," read -r -a bam_array <<< "$bam_list"
    for bam in ${bam_array[@]}; do
        bam_mdup_arg="${bam_mdup_arg}I=$bam "
        bam_cp_arg="${bam_cp_arg}cp -t . $workdir/$bamdir/$bam*; "
    done
fi

# validate files
if [[ "$validate" = yes ]]; then
    # check indir exists
    if [[ ! -d "$workdir/$bamdir" ]]; then
        printf "\nERROR: input directory does not exist: %s/%s/%s\n" $workdir $bamdir
        echo "$help_message"; exit 1
    fi

    # check resource files
    if [[ ! -f "$workdir/$resdir/$fasta" ]]; then
        printf "\nERROR: FASTA file not found: %s/%s/%s\n" $workdir $resdir $fasta
    else
        fasta_base=${fasta%%.*}
    fi

    # check BAM files exist
    for bam in ${bam_array[@]}; do
        if [[ ! -f "$workdir/$bamdir/$bam" ]]; then
            printf "\nERROR: input BAM file does not exist: %s/%s/%s\n" $workdir $bamdir $bam
            echo "$help_message"; exit 1
        fi
    done
fi

# if no name provided extract from first bam
if [[ -z "$name" ]]; then
    name=${bam_array[0]%%.*}
fi

# set dependency argument
if [[ ! -z $depend_list ]]; then
    depend="#PBS -W depend=$depend_list"
fi

# create required dirs
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir/metrics

# run job
jobid=$(cat <<- EOS | qsub -N $name.merge - 
	#!/bin/bash
	#PBS -l walltime=10:00:00
	#PBS -l select=1:mem=30gb:ncpus=8
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $workdir/$logdir/$name.merge_bam_picard.log
	$depend

	# load modules
	module load samtools/1.2
	module load java/jdk-8u66
	module load picard/2.6.0

	# copy accross bam and fasta
	cp $workdir/$resdir/$fasta_base* .
	$bam_cp_arg

	# merge/mark dup
	java -Xmx28G -jar /apps/picard/2.6.0/picard.jar MarkDuplicates \
		$bam_mdup_arg \
		O=$name.bam \
		M=$name.mark_duplicate_metrics \
		CREATE_INDEX=true
	cp $name.mark_duplicate_metrics $workdir/$qcdir/metrics/

	# alignment metrics
	java -Xmx16g -jar /apps/picard/2.6.0/picard.jar CollectAlignmentSummaryMetrics \
		R=$fasta \
		I=$name.bam \
		O=$name.alignment_summary_metrics
	cp $name.alignment_summary_metrics $workdir/$qcdir/metrics/

	# copy final bam to outdir
	samtools index $name.bam
	cp $name.bam* $workdir/$bamdir/

	ls -lhAR
	EOS
)
echo "JOBID: $jobid"
