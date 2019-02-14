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
Wrapper to call structural variants with Delly2

usage:
    bash $scr_name [-options] -b <bam> -v <vartype>
required arguments:
    -b|--bam : comma separated list of bam files [BAM]
    -v|--vartype : variant type to analyse [INS,DEL,INV,BND,DUP]
    -F|--fasta : whole genome fasta [FASTA]
optional arguments:
    -o|--outdir : output directory [string] (default = '.')
    -ld|--logdir : output directory for log files [string] (default = --outdir)
    -E|--exclude : regions to exclude (default = NULL) 
    -n|--name : prefix for output files [string] (deafult = extracted from bam)
    --vcf : VCF file for re-genotyping [string] (default = NULL)
    --delly_arg : any additional arg to pass to delly [string] (default = NULL)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list pf PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # exclusion regions should be provided as ...
    # --delly_arg should be provided as quoted list of arguments,
      and are passed as-is to delly. e.g. --delly_arg '-q 20 --noindels'
    # see Delly2 man for additional info: https://github.com/dellytools/delly

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam_list=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -ld|--logdir)
            logdir=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -E|--exclude)
            exclude=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -v|--vartype)
            vartype=$2
            shift
            ;;
        --vcf)
            vcf=$2
            shift
            ;;
        --delly_arg)
            delly_arg=$2
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
if [[ -z "${bam_list:-}" ]]; then
    printf "\nERROR: No --bam argument provided\n"
    echo "$help_message"; exit 1
elif [[ -z "${vartype:-}" ]]; then
    printf "\nERROR: No --vartype argument provided\n"
    echo "$help_message"; exit 1
elif [[ -z "${fasta:-}" ]]; then
    printf "\nERROR: No --fasta argument provided\n"
    echo "$help_message"; exit 1  
fi

# split bam list
IFS="," read -r -a bam_array <<< "$bam_list"

# check files unless flagged
if [[ $check = yes ]]; then
    for bam in ${bam_array[@]}; do
        if [[ ! -r $workdir/$bam ]]; then
            printf "\nERROR: BAM file is not readable: %s/%s\n" $workdir $bam
            echo "$help_message"; exit 1
        fi
    done
    if [[ ! -r $workdir/$fasta ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    fi
fi

# check variant type 
vartype_list="INS DEL INV BND DUP"
if [[ ! "$vartype_list" =~ (^|[[:space:]])$vartype($|[[:space:]]) ]]; then
    printf "\nERROR: input variant type not recognised: %s\n" $vartype
    echo "$help_message"; exit 1
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
    name=$(basename ${bam_array[0]})
    name=${name%%.*}
fi

# set arguments for bam files
for bam in ${bam_array[@]}; do
    bam_cp="${bam_cp:-}; cp -rL $workdir/$bam* ."
    bam_base=$(basename "$bam")
    bam_arg="${bam_arg:-} $bam_base"
done
bam_cp=${bam_cp#*; }
bam_arg=${bam_arg#* }

# get basename/prefix of fasta
fasta_base=$(basename "$fasta")
fasta_prefix=${fasta%%.*}

# set exclusion call arguments if provided
if [[ ! -z ${exclude:-} ]]; then
	excl_cp="cp -rL $workdir/$exclude ."
    excl_base=$(basename "$exclude")
    excl_arg="-x $excl_base"
fi

# set vcf args if provided
if [[ ! -z ${vcf:-} ]]; then
    vcf_cp="cp $workdir/$vcf* ."
    vcf_base=$(basename "$vcf")
    vcf_arg="-v $vcf_base"
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set commands
call_command=("delly call -t $vartype -g $fasta_base -o $name.bcf ${excl_arg:-}"
              "${vcf_arg:-} ${delly_arg:-} $bam_arg")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=46gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.delly
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load required modules
		module load bcftools/1.2
		module load anaconda3/personal
		source activate delly

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy required inputs to scratch
		cp -rL $workdir/$fasta_prefix* .
		${excl_cp:-} &>> $out_log
		${vcf_cp:-} &>> $out_log
		${bam_cp:-} &>> $out_log

		# run delly calling
		${call_command[@]} &>> $out_log

		cp $name.bcf* $workdir/$outdir/ &>> $out_log

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhRA >> $out_log
		cp $out_log $workdir/$logdir
		EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
