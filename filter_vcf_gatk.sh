#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'

# help message
help_message="

usage:
    bash $(basename "$0") [-options] -v <VCF>
purpose:
    Wrapper to filter VCF files with GATK v4.0 SelectVariants
required arguments:
    -v|--vcf : VCF file to filter [VCF]
optional arguments:
    -o|--outdir : output directory [string] (default = '.')
    -l|--logdir : output directory for log files [string] (default = --outdir) 
    -n|--name : prefix for output files [string] (deafult = extracted from bam)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list pf PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -v|--vcf)
            vcf=$2
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
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arg
if [[ -z "${vcf:-}" ]]; then
    printf "\nERROR: No --vcf argument provided\n"
    echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$vcf ]]; then
        printf "\nERROR: VCF file is not readable: %s/%s\n" $workdir $vcf
        echo "$help_message"; exit 1
    fi
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
    name=$(basename "$vcf")
    name=${name%%.*}
fi

# get basenames/prefixes
vcf_base=$(basename "$vcf")
fasta_base=$(basename "$fasta")
fasta_prefix=${fasta%%.*}

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# run job 
jobid=$(cat <<- EOS | qsub -N $name.filter_vcf_gatk -
		#!/bin/bash
		#PBS -l walltime=20:00:00
		#PBS -l select=1:mem=20gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.filter_vcf_gatk.log
		${depend:-}

		# load required modules
		module load java/jdk-8u66
		module load gatk/3.6

		# copy required inputs to scratch
		cp $workdir/$vcf* .
		cp $workdir/$fasta_prefix* .

		# run gatk filtering 
		java -Xmx32G -jar /apps/gatk/3.6/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $fasta_base \
			-o $name.vcf \
			--variant $vcf_base \
			--filterExpression "QSS / AD > 20" \
			--filterName "average QSS < 20" \
			--filterExpression "! vc.getGenotype(\"${name}\").isHomRef()" \
			--filterName "is hom ref" 

		cp $name.vcf* $workdir/$outdir/
		
		ls -lhRA
		EOS
)
echo "JOBID: $jobid"
