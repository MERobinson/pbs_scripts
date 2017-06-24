#!/usr/bin/env bash

set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
fqdir=fastq
resdir=resources
logdir=logs
scriptdir=scripts
delim=","
genome=hg38
fq1_idx=1
fq2_idx=2
sm_idx=3
fasta=genome.fa
index=genome.fa
log_file=$(basename $0 .sh)_$(date "+%Y-%m-%d").log
verbose=n

# help message
help_message="
usage:
    bash $(basename "$0") [wfrdg] -s <SAMPLE_INFO>
purpose:
    # parse a sample info file and run ChIPseq analysis pipeline for each sample
    # by default, will map reads with BWA MEM, call peaks and generate metrics & tracks
required arguments:
    -i|--sampleinfo : a delimited text file containing sample information (see below)
optional arguments:
    -w|--workdir : working directory (default = pwd)
    -f|--fqdir : directory within workdir containing FASTQ files (default = fastq)
    -r|--resdir : directory within workdir containing resource files (default = resources)
    -l|--logdir : directory within workdir to output log files (default = logs)
    -s|--scriptdir : directory within workdir containing pipeline scripts (default = scripts)
    -d|--delim : delimiter used in sample info file (default = ',')
    -g|--genome : genome version used (default = hg38)
additional info:
    # sample info file must contain FASTQ names and sample names
    # peak calling will only be conducted if relevant control sample name in -cs_idx column
    # possible columns to include in sample_info:
        - Filename of FASTQ R1 [required]
        - Filename of FASTQ R2 [required - leave blank column for SE]
        - Sample name [required]
        - Control sample name [optional]
        - Read group ID [optional]
        - Platform unit [optional]
        - Library name [optional]
        - Sequencing centre [optional]
        - Platform [default = ILLUMINA]
        - Platform model [optional]
        - Run date [optional]
        - Predicted median insert size [optional] 
    # to include any of the above in BAM header, provide column indexes as detailed below
    # following resource files must be present in resdir:
        - Whole genome sequence [FASTA and corresponding disctionary file]
        - BWA index files [all index components required]
    # to alter expected resource filenames, see resource filename arguments
column index arguments:
    -f1_idx : index of column containing FASTQ R1 filenames (default = 1)
    -f2_idx : index of column containing FASTQ R2 filenames (default = 2)
    -sm_idx : index of column containing sample names (default = 3)
    -id_idx : index of column containing read group ID (default = null)
    -pu_idx : index of column containing platform unit (default = null)
    -lb_idx : index of column containing library name (default = null)
    -cn_idx : index of column containing sequencing centre (default = null)
    -pl_idx : index of column containing platform (default = null)
    -pi_idx : index of column containing median insert size (deafult = null)
    -pm_idx : index of column containing platform model (default = null)
    -dt_idx : index of column containing run date (default = null)
resource filename arguments:
    -ix|--index_base : basename of index files in resdir (default = genome.fa)
    -fa|--fasta : name of FASTA file in resdir (default = genome.fa) 

"

while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
        -i|--sampleinfo)
            sample_info=$2
            shift
            ;;
        -w|--workdir)
            workdir=$2
            shift
            ;;
        -f|--fqdir)
            fqdir=$2
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
        -d|--delim)
            delim=$2
            shift
            ;;
        -g|--genome)
            genome=$2
            shift
            ;;
        -f1_idx)
            fq1_idx=$2
            shift
            ;;
        -f2_idx)
            fq2_idx=$2
            shift
            ;;
        -sm_idx)
            sm_idx=$2
            shift
            ;;
        -lb_idx) 
            lb_idx=$2
            shift
            ;;
        -pu_idx)
            pu_idx=$2
            shift
            ;;
        -pl_idx)
            pl_idx=$2
            shift
            ;;
        -cn_idx) 
            cn_idx=$2
            shift
            ;;
        -pi_idx)
            pi_idx=$2
            shift
            ;;
        -pm_idx)
            pm_idx=$2
            shift
            ;;
        -dt_idx)
            dt_idx=$2
            shift
            ;;
        -ix|--index_base)
            index=$2
            shift
            ;;
        -fa|--fasta)
            fasta=$2
            shift
            ;;
        *) 
            echo "Error: Illegal argument"
            echo "$help_message"; exit 1
	    ;;
     esac
shift
done

# check required arg
if [[ -z "${sample_info:-}" ]]; then
	printf "\nERROR: Sample info file not provided.\n"
	echo "$help_message"; exit 1
fi
if [[ ! -r "$workdir/$sample_info" ]]; then
	printf "\nERROR: Input sample file doesnt exist/isnt readable: %s/%s\n" $workdir $sample_info
	echo "$help_message"; exit 1
fi
if [[ ! -d "$workdir/$resdir" ]]; then
	printf "\nERROR: Resources directory doesnt exist: %s/%s\n" $workdir $resdir
	echo "$help_message"; exit 1
fi
if [[ !	-d "$workdir/$fqdir" ]]; then
    printf "\nERROR: FASTQ directory doesnt exist: %s/%s\n" $workdir $fqdir
    echo "$help_message"; exit 1
fi
if [[ !	-r "$workdir/$resdir/$index" ]]; then
    printf "\nERROR: Index file doesnt exist/isnt readable: %s/%s/%s\n" $workdir $resdir $index
    echo "$help_message"; exit 1
fi
if [[ ! -r "$workdir/$resdir/$fasta" ]]; then
    printf "\nERROR: FASTA file doesnt exist/isnt readable: %s/%s/%s\n" $workdir $resdir $fasta
    echo "$help_message"; exit 1
fi

# function to extract fields from sample info file
function extract_info {
	awk -F $delim -v v1=$sm_idx -v v2=$sample_name -v v3=$1 \
	'$v1 == v2 {print $v3}' $workdir/$sample_info
}

# get unique sample names
samples=$(tail -n +2 "$sample_info" | cut -d "$delim" -f "$sm_idx" | sort | uniq)
echo "Sample list: ${samples[@]}" > $log_file

# loop through sample names
for sample_name in ${samples[@]}; do

    # get all FASTQ files for current sample -> array
    fq1_array=($(extract_info $fq1_idx))
    fq2_array=($(extract_info $fq2_idx))

    # do same for optional header info
    if [[ ! -z "${id_idx:-}" ]]; then id_array=($(extract_info $id_idx)); fi
    if [[ ! -z "${pu_idx:-}" ]]; then pu_array=($(extract_info $pu_idx)); fi
    if [[ ! -z "${lb_idx:-}" ]]; then lb_array=($(extract_info $lb_idx)); fi
    if [[ ! -z "${cn_idx:-}" ]]; then cn_array=($(extract_info $cn_idx)); fi
    if [[ ! -z "${pl_idx:-}" ]]; then pl_array=($(extract_info $pl_idx)); fi
    if [[ ! -z "${pi_idx:-}" ]]; then pi_array=($(extract_info $pi_idx)); fi
    if [[ ! -z "${pm_idx:-}" ]]; then pm_array=($(extract_info $pm_idx)); fi
    if [[ ! -z "${dt_idx:-}" ]]; then dt_array=($(extract_info $dt_idx)); fi

    # loop through each FASTQ1
    for idx in ${!fq1_array[@]}; do
	
        # assemble command for alignment script with any provided optional args
        command="bash $workdir/$scriptdir/align_fastq_bwa.sh -f ${fq1_array[$idx]} -n $sample_name"
        
        if [[ ! -z "${fq2_array[$idx]:-}" ]]; then command="${command} --mate ${fq2_array[$idx]} "; fi 
        if [[ ! -z "${id_array[$idx]:-}" ]]; then command="${command} --id ${id_array[$idx]} "; fi
        if [[ ! -z "${pu_array[$idx]:-}" ]]; then command="${command} --pu ${pu_array[$idx]} "; fi
        if [[ ! -z "${lb_array[$idx]:-}" ]]; then command="${command} --lb ${lb_array[$idx]} "; fi
       	if [[ ! -z "${cn_array[$idx]:-}" ]]; then command="${command} --cn ${cn_array[$idx]} "; fi
       	if [[ ! -z "${pl_array[$idx]:-}" ]]; then command="${command} --pl ${pl_array[$idx]} "; fi
       	if [[ ! -z "${pi_array[$idx]:-}" ]]; then command="${command} --pi ${pi_array[$idx]} "; fi
       	if [[ ! -z "${pm_array[$idx]:-}" ]]; then command="${command} --pm ${pm_array[$idx]} "; fi
        if [[ ! -z "${dt_array[$idx]:-}" ]]; then command="${command} --dt ${dt_array[$idx]} "; fi
	
	    # run alignment script
        echo "$command" > $log_file
	    $command

    done
done


