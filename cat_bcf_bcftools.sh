#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='yes'
format='bcf'
date=`date '+%Y-%m-%d %H:%M:%S'`
scr_name=$(basename "$0")

# help message
help_message="
Wrapper to GATK v3 CatVariants - to concatenate VCFs.

usage:
    bash $scr_name [-options] -v <list,of,VCF> -F <FASTA>
required arguments:
    -v|--vcf_list : comma separated list of VCF/BCF files to cat
optional arguments:
    -o|--outdir : output directory for concatenated VCF (default = PWD)
    -f|--format : output format [vcf|bcf] (default = BCF)
    -ld|--logdir : output directory for log files (default = --outdir)
    -n|--name : name of output file [string] (default = extracted from first VCF)
    --check : whether to check input files [yes,no] (default = 'yes')
    --depend : dependency list to pass to PBS script (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # depend arguments should have PBS format, e.g. 'afterok:123456,afterok:123457'

"

# parse arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -v|--vcf_list)
            vcf_list=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -f|--format)
            format=$2
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
        --check)
            check=$2
            shift
            ;;
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided\n"
            echo "$help_message"; exit 2
            ;;
    esac
    shift
done

# set outdirs
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arg
if [[ -z ${vcf_list:-} ]]; then
	printf "\nERROR: no vcf argument provided\n"
	echo "$help_message"; exit 2
else
    IFS="," read -r -a vcf_array <<< "$vcf_list"
fi

# check files unless flagged
if [[ $check = yes ]]; then
    for vcf in ${vcf_array[@]}; do
        if [[ -z $workdir/${vcf} ]]; then
            printf "\nERROR: VCF file is not readable: $workdir/$vcf\n"
            echo "$help_message"; exit 2
        fi
    done
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
    name="${vcf_array[0]%%.*}.cat"
fi

# set arg for VCF list copying and command input
cp_arg=$(printf "cp $workdir/%s* . ;" ${vcf_array[@]})

# set format
if [[ $format = "vcf" ]]; then
    format_arg="-O v"; ext="vcf"
elif [[ $format = "bcf" ]]; then
    format_arg="-O b"; ext="bcf"
else
    printf "\nERROR: format argument not recognised\n"
    echo "$help_message"; exit 1
fi 

# setup output directories
mkdir -p $workdir/$logdir
mkdir -p $workdir/$outdir

# set commands
cat_command=("bcftools concat -o $name.$ext ${format_arg} --threads 8 ${vcf_array[@]##*/}")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# cat variants
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=18gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.cat
		##PBS -q med-bio
		#PBS -o $std_log
		${depend:-}		

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# load modules
		module load bcftools

		# copy files to scratch
		$cp_arg &>> $out_log

		# call
		${cat_command[@]} &>> $out_log

		bcftools index $name.$ext
		cp $name.$ext* $workdir/$outdir/  &>> $out_log
		
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		printf "JOBID: %s\n" \${PBS_JOBID:-} >> $out_log
		ls -lhAR >> $out_log
		cp $out_log $workdir/$logdir/
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
