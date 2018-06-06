#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'
name="germline_sv"
date=`date '+%Y-%m-%d %H:%M:%S'`
scr_name=$(basename "$0")

# help message
help_message="
Wrapper to merge structural variants (by variant) with Delly2

usage:
    bash $(basename "$0") [-options] -s <bcf,list> -v <vartype>
required arguments:
    -s|--sv_list : comma separated list of sv files to merge [BCF]
    -v|--vartype : variant type to analyse [INS,DEL,INV,BND,DUP]
optional arguments:
    -o|--outdir : output directory [string] (default = '.')
    -ld|--logdir : output directory for log files [string] (default = --outdir)
    -n|--name : prefix for output files [string] (deafult = 'germline_sv')
    --max : maximum sv size [integer] (default = 1000000)
    --min : minimum sv size [integer] (default = 0)
    --offset : maximum breakpoint offset [integer] (default = 1000)
    --overlap : minimum recipricol overlap [float] (default = 0.8)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list pf PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # see Delly2 man for additional info: https://github.com/dellytools/delly

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -s|--sv_list)
            sv_list=$2
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
        -n|--name)
            name=$2
            shift
            ;;
        -v|--vartype)
            vartype=$2
            shift
            ;;
        --max)
            max_arg="-n $2"
            shift
            ;;
        --min)
            min_arg="-m $2"
            shift
            ;;
        --overlap)
            ol_arg="-r $2"
            shift
            ;;
        --offset)
            bp_arg="-b $2"
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
if [[ -z "${sv_list:-}" ]]; then
    printf "\nERROR: No --sv_list argument provided\n"
    echo "$help_message"; exit 1
else
    IFS="," read -ra sv_array <<< "$sv_list"
fi
if [[ -z "${vartype:-}" ]]; then
    printf "\nERROR: No --vartype argument provided\n"
    echo "$help_message"; exit 1 
fi

# check files unless flagged
if [[ $check = yes ]]; then
    for sv_file in "${sv_array[@]}"; do
        if [[ ! -r $workdir/$sv_file ]]; then
            printf "\nERROR: SV file is not readable: %s/%s\n" $workdir $sv_file
            echo "$help_message"; exit 1
        fi
    done
fi

# check variant type 
vartype_list="INS DEL INV BND DUP"
if [[ ! "$vartype_list" =~ (^|[[:space:]])$vartype($|[[:space:]]) ]]; then
    printf "\nERROR: input variant type not recognised: %s\n" $vartype
    echo "$help_message"; exit 1
fi

# set input arguments for bcf files
for sv_file in "${sv_array[@]}"; do
    bcf_cp="${bcf_cp:-}; cp $workdir/$sv_file* ."
    sv_file=$(basename "$sv_file")
    bcf_arg="${bcf_arg:-} $sv_file" 
done
bcf_cp=${bcf_cp#*;}
bcf_arg=${bcf_arg#* }

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set commands
merge_command=("delly merge -t $vartype -o $name.bcf ${bp_arg:-} ${ol_arg:-}"
               "${max_arg:-} ${min_arg:-} $bcf_arg")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write script 
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=24gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.merge
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load required modules
		module load bcftools/1.2
		module load anaconda3/personal
		source activate delly

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy inputs to scratch
		$bcf_cp &>> $out_log

		# run delly for each variant type
		${merge_command[@]} &>> $out_log

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
