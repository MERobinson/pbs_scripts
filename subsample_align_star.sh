#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
fqdir=fastq
bamdir=bam
logdir=logs
scriptdir=scripts
log_file=$(basename $0 .sh)_$(date "+%Y-%m-%d").log

# help message
help_message="
usage:
    bash $(basename "$0") [wfrdg] -f <FASTQ1> -r <FASTQ2>
purpose:
    # Takes a pair of FASTQ files, subsamples reads and re-runs STAR alignment
required arguments:
    -f1|--fastq1 : FASTQ R1
    -f2|--fastq2 : FASTQ R2
optional arguments:
    -w|--workdir : working directory (default = pwd)
    -o|--outdir : output directory for star files files (default = star)
    -f|--fqdir : directory within workdir containing FASTQ files (default = fastq)
    -r|--resdir : directory within workdir containing resource files (default = resources)
    -l|--logdir : directory within workdir to output log files (default = logs)
    -s|--scriptdir : directory within workdir containing pipeline scripts (default = scripts)

"

while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
        -f1|--fastq1)
            fq1=$2
            shift
            ;;
        -f2|--fastq2)
            fq2=$2
            shift
            ;;
        -w|--workdir)
            workdir=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
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
        -s|--scriptdir)
            scriptdir=$2
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
if [[ -z "${fq1:-}" ]] | [[ -z "${fq2:-}" ]]; then
	printf "\nERROR: FASTQ files not provided.\n"
	echo "$help_message"; exit 1
else
    if [[ ! -r $workdir/$fqdir/$fq1 ]] | [[ ! -r $workdir/$fqdir/$fq2 ]]; then
        printf "\nERROR: FASTQ files not readable\n"
        echo "$help_message"; exit 1
    fi
fi

# get name
if [[ -z ${name:-} ]]; then
    name=${fq1%%.*}
fi

# get total count
module load samtools
total_count=$(awk '{s++}END{print s/4}' "$workdir/$fqdir/$fq1")

# subsample in increments of 2.5M and run alignment
i=2500000
while [ $i -lt $total_count ]; do 
    
    j=$[$i / 100000]

    paste $workdir/$fqdir/$fq1 $workdir/$fqdir/$fq2 |\
    awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' |\
    shuf |\ 
    head $i | sed 's/\t\t/\n/g' |\
    awk '{print $1 > "$workdir/$fqdir/${name}_sub$j.R1.fq"; print $2 > "$workdir/$fqdir/${name}_sub$j.R2.fq"}'

    star="bash $workdir/$scriptdir/align_fastq_star.sh -f ${name}_sub$j.R1.fq -m ${name}_sub$j.R2.fq"
    star="${star} -n ${name}_sub$j"
    $star

    i=$[$i + 2500000]
    exit
done
    
