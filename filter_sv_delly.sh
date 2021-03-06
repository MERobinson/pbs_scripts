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
Wrapper to filter structural variants with Delly2

usage:
    bash $scr_name [-options] -b <BCF> -v <VARTYPE>
required arguments:
    -b|--bcf : structural variants file [BCF]
    -v|--vartype : variant type [INS,DEL,INV,BND,DUP]
optional arguments:
    -s|--sample_info : sample info file for somatic variants [tsv] (default = NULL)
    -o|--outdir : output directory [string] (default = '.')
    -l|--logdir : output directory for log files [string] (default = --outdir) 
    -n|--name : prefix for output files [string] (deafult = extracted from bam)
    --delly_arg : any additional arg to pass to delly [string] (default = NULL)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list pf PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # filtering is run in germline mode unless --sample_info argument is provided
    # sample info file should be a tab separated text file with sample ID in
      first column and either 'tumor' or 'control' in second column.
    # --delly_arg should be provided as quoted list of arguments,
      and are passed as-is to delly. e.g. --delly_arg '-q 20 --noindels'
    # see Delly2 man for additional info: https://github.com/dellytools/delly

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bcf)
            bcf=$2
            shift
            ;;
        -v|--vartype)
            vartype=$2
            shift
            ;; 
        -s|--sample_info)
            sample_info=$2
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
if [[ -z "${bcf:-}" ]]; then
    printf "\nERROR: No --bcf argument provided\n"
    echo "$help_message"; exit 1
elif [[ -z "${vartype:-}" ]]; then
    printf "\nERROR: No --vartype argument provided\n"
    echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$bcf ]]; then
        printf "\nERROR: BCF file is not readable: %s/%s\n" $workdir $bcf
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
    name=$(basename "$bcf")
    name=${name%%.*}
fi

# get basenames/prefixes
bcf_base=$(basename "$bcf")

# if control provided setup somatic filtering args
if [[ ! -z ${sample_info:-} ]]; then
    echo "Filtering in somatic mode"
    filter_type="somatic"
    si_cp="cp $workdir/$sample_info ."
    sample_info=$(basename "$sample_info")
    si_arg="-s $sample_info"
else
    echo "Filtering in germline mode"
    filter_type="germline"
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set commands
filter_command=("delly filter -t $vartype -f $filter_type -o $name.bcf"
                "${si_arg:-} $bcf_base")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=18gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.filt
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load required modules
		module load bcftools/1.2
		module load anaconda3/personal
		source activate delly

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy required inputs to scratch
		cp $workdir/$bcf* . &>> $out_log
		${si_cp:-} &>> $out_log

		# run delly
		${filter_command[@]} &>> $out_log

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
