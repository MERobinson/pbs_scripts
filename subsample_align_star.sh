#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=sub
fqdir=fastq
logdir=logs
resdir=resources
scriptdir=scripts
index=index
log_file=$(basename $0 .sh)_$(date "+%Y-%m-%d").log
start_val=2500000
increment=2500000

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
    -n|--name : name prefix for output files (default = extracted from fq1)
    -o|--outdir : output directory for star files files (default = sub)
    -f|--fqdir : directory within workdir containing FASTQ files (default = fastq)
    -r|--resdir : directory within workdir containing resource files (default = resources)
    -l|--logdir : directory within workdir to output log files (default = logs)
    -s|--scriptdir : directory within workdir containing pipeline scripts (default = scripts)
    -I|--index : name of STAR index (deafault = 'index')
    --start : starting value for subsampling (default = 2,500,000)
    --increment : step size for subsampling (default = 2,500,000)
    --max : maximum value for subsampling (default = read number of input FASTQ)

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
        -n|--name)
            name=$2
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
        -I|--index)
            index=$2
            shift
            ;;
        --start)
            start_val=$2
            shift
            ;;
        --increment)
            increment=$2
            shift
            ;;
        --max)
            max_val=$2
            shift
            ;;
        *) 
            printf "\nERROR: Illegal argument: %s %s\n" $1 $2
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
    if [[ ! -r $workdir/$fq1 ]] | [[ ! -r $workdir/$fq2 ]]; then
        printf "\nERROR: FASTQ files not readable\n"
        echo "$help_message"; exit 1
    fi
fi

# get name
if [[ -z ${name:-} ]]; then
    name=$(basename $fq1)
    name=${name%%.*}
fi

# get total count
if [[ -z ${max_val:-} ]]; then
    module load samtools
    max_val=$(awk '{s++}END{print s/4}' "$workdir/$fq1")
fi

# subsample in increments of 2.5M and run alignment
i=$start_val
while [ $i -le $max_val ]; do 

    # set name to M
    j=$(bc <<< "scale=1; $i / 1000000")

    # subsample
    sub="bash $scriptdir/subsample_fastq.sh -f1 $fq1 -f2 $fq2 -s $i -n ${name}_${j}M"
    sub="${sub} --outdir $fqdir --logdir $logdir"
    $sub > tmp.log

    # check return/get jobid
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: subsampling script failed\n"; rm tmp.log; exit 1
    else
        jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
        rm tmp.log
    fi

    # run star
    star="bash $scriptdir/align_fastq_star.sh -f1 $fqdir/${name}_${j}M.R1.fq"
    star="${star} -f2 $fqdir/${name}_${j}M.R2.fq -n ${name}_${j}M --depend afterok:$jobid"
    star="${star} -I $resdir/$index --check no --outdir $outdir"
    $star

    i=$[$i + $increment]
   
done 
