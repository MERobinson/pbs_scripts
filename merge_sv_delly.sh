#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'
name="germline_sv"

# help message
help_message="

usage:
    bash $(basename "$0") [-options] -s <SV,LIST> -v <VAR>
purpose:
    Wrapper to merge structural variants with Delly2
required arguments:
    -s|--sv_list : list of sv files to merge [BCF]
    -v|--vartype : variant type to analyse [INS,DEL,INV,BND,DUP]
optional arguments:
    -o|--outdir : output directory [string] (default = '.')
    -l|--logdir : output directory for log files [string] (default = --outdir)
    -n|--name : prefix for output files [string] (deafult = 'germline_sv')
    --delly_arg : any additional arg to pass to delly [string] (default = NULL)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list pf PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # --delly_arg should be provided as quoted list of arguments,
      and are passed as-is to delly. e.g. --delly_arg '-q 20 --noindels'
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
        -l|--logdir)
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
if [[ -z "${sv_list:-}" ]]; then
    printf "\nERROR: No --sv_list argument provided\n"
    echo "$help_message"; exit 1
else
    IFS="," read -ra sv_array <<< "$sv_list"
fi
if [[ -z "${vartype:-}" ]]; then
    printf "\nERROR: No --var argument provided\n"
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

# run job 
jobid=$(cat <<- EOS | qsub -N $name.merge_sv_delly -
		#!/bin/bash
		#PBS -l walltime=20:00:00
		#PBS -l select=1:mem=40gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.merge_sv_delly.log
		${depend:-}

		# load required modules
		module load bcftools/1.2
		module load anaconda3/personal
		source activate delly

		# copy required inputs to scratch
		$bcf_cp

		# run delly for each variant type
		delly merge -t $vartype \
			-o $name.bcf \
			-m 500 -b 500 -r 0.5 \
			${delly_arg:-} $bcf_arg

		cp $name.bcf* $workdir/$outdir/
		
		ls -lhRA
		EOS
)
echo "JOBID: $jobid"
