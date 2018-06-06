#!/usr/bin/env bah
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'
date=`date '+%Y-%m-%d %H:%M:%S'`
scr_name=$(basename "$0")
merge='id'
format='b'
force_sm='no'
print_head='no'
missing='yes'

# help message
help_message="
Wrapper to merge a list of BCF/VCF with bcftools.

usage:
    bash $scr_name [-options] -b <FASTQ>
required arguments:
    -b|--bcf_list : comma separated list of BCF files 
optional arguments:
    -m|--merge : how to merge records [none|snps|indels|both|all|id] (default = id)
    -O|--out-type : output format [b|u|z|v] (default = b)
    -n|--name : name prefix for output files (default = extracted from first bcf)
    -o|--outdir : output directory for merged bcf (default = PWD)
    -ld|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [yes|no] (default = yes)
    --depend : list of PBS dependencies (default = NULL)
    --force-samples : whether to create unique sample names if dulicated [yes|no] (default = no)
    --print-header : print merged header only [yes|no] (default = no)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log output directories inherit from --outdir unless specified

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bcf_list)
            bcf_list=$2
            shift
            ;;
        -m|--merge)
            merge=$2
            shift
            ;;
        -O|--output-type)
            format=$2
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
        -ld|--logdir)
            logdir=$2
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
        --force-samples)
            force_sm=$2
            shift
            ;;
        --print-header)
            print_head=$2
            shift
            ;;
        *)
            echo "Error: Illegal argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z "${bcf_list:-}" ]]; then
    printf "\nERROR: --bcf_list argument required\n"
    echo "$help_message"; exit 1
else
    IFS="," read -r -a bcf_array <<< $bcf_list
fi

# check files
if [[ "${check:-}" = yes ]]; then
    for bcf in ${bcf_array[@]}; do
        if [[ ! -r $workdir/$bcf ]]; then
            printf "\nERROR: BCF file cannot be read: %s/%s\n" $workdir $bcf
            echo "$help_message"; exit 1
        fi
    done
fi

# set output directories
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set basenames
bcf_cp=$(printf "cp $workdir/%s* .; " ${bcf_array[@]})
bcf_cp=${bcf_cp%;*}
bcf_arg=$(printf "%s " ${bcf_array[@]##*/})

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    name="$(basename ${bcf_array[0]%%.*}).merge"
fi

# set format arg
if [[ $format = 'b' ]]; then
    format_arg="-O $format"; ext="bcf"
elif [[ $format = 'u' ]]; then
    format_arg="-O $format"; ext="bcf"
elif [[ $format = 'z' ]]; then
    format_arg="-O $format"; ext="vcf.gz"
elif [[ $format = 'v' ]]; then
    format_arg="-O $format"; ext="vcf"
else
    printf "\nERROR: format not recognised\n"
    echo "$help_message"; exit 1
fi

# parse optional program arg
if [[ $force_sm = 'yes' ]]; then
    force_arg="--force-samples"
fi
if [[ $print_head = 'yes' ]]; then
    head_arg="--print-header"
fi

# set commands
merge_command=("bcftools merge -o $name.$ext --merge $merge $format_arg ${force_arg:-}"
               "${head_arg:-} --threads 8 $bcf_arg")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=20gb:ncpus=8
		#PBS -j oe
		#PBS -N $name.merge
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load modules
		module load bcftools

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log
		
		# copy resource files to scratch
		$bcf_cp &>> $out_log

		# alignment metrics
		${merge_command[@]} &>> $out_log

		bcftools index $name.$ext
		cp $name.$ext* $workdir/$outdir/ &>> $out_log
 
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR &>> $out_log 
		cp $out_log $workdir/$logdir/
	EOS
) 
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
