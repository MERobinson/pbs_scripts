#!/usr/bin/env bash
set -o pipefail
set -o nounset 

# default arg
workdir=$PWD
outdir=''
check='on'

# help message
help_message="
Wrapper to call differential footprints from DNase data with Wellington

usage:
    bash $(basename "$0") [-options] -r <BED> -b1 <BAM> -b2 <BAM>
required arguments:
    -r|--regions : accessible regions to profile [BED]
    -t|--treatment : aligned reads for treatment/test [BAM]
    -c|--control : aligned reads for control [BAM] 
optional arguments:
    -o|--outdir : output directory name (default = PWD)
    -n|--name : name of comparison & outdir subfolder (default = <tname>_vs_<cname>)
    -tn|--tname : treatment name (default = extracted from treatment)
    -cn|--cname : control name (defailt = extracted from control)
    --fdr : FDR cutoff (default = 0.01)
    --fdrlimit : minimum pval to consider sig for FDR (default = -20)
    --logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on,off] (default = on)
    --depend : list of PBS dependencies (default = NULL)
additional info:

"

# parse command line arguments
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-r|--regions)
		    regions=$2
		    shift
		    ;;
		-t|--treatment)
		    treatment=$2
		    shift
		    ;;
        -c|--control)
            control=$2
            shift
            ;;
		-o|--outdir)
		    outdir=$2
		    shift
		    ;;
        -n|--name)
            name=$2
            shift
            ;;
		-tn|--tname)
		    tname=$2
		    shift
		    ;;
        -cn|--cname)
            cname=$2
            shift
            ;;
        --fdr)
            fdr="-fdr $2"
            shift
            ;;
        --fdrlimit)
            fdrlimit="-fdrlimit $2"
            shift
            ;;
        --logdir)
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
		*)
		    printf "ERROR: Unrecognised argument: %s %s" $1 $2
		    echo "$help_message"; exit 1
		    ;;
	esac
	shift
done

# check required args
if [[ -z ${regions:-} ]]; then
    printf "\nERROR: --region argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${treatment:-} ]]; then
    printf "\nERROR: --treatment argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${control:-} ]]; then
    printf "\nERROR: --control argument required\n"
    echo "$help_message"; exit 1
fi

# set optional args
if [[ -z ${tname:-} ]]; then
    tname=$(basename "${treatment}")
    tname=${tname%%.*}
fi
if [[ -z ${cname:-} ]]; then
    cname=$(basename ${control})
    cname=${cname%%.*}
fi
if [[ -z ${name:-} ]]; then
    name="${tname}_vs_${cname}"
fi
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# check files
if [[ ${check} = 'on' ]]; then
    if [[ ! -r ${workdir}/${regions} ]]; then
        printf "\nERROR: input file cannot be read: %s/%s\n" $workdir $regions
        echo "$help_message"; exit 1
    elif [[ ! -r ${workdir}/${treatment} ]]; then
        printf "\nERROR: bam file cannot be read: %s/%s\n" $workdir $treatment
        echo "$help_message"; exit 1
    elif [[ ! -r ${workdir}/${control} ]]; then
        printf "\nERROR: bam file cannot be read: %s/%s\n" $workdir $control
        echo "$help_message"; exit 1
    fi
fi

# set commands
pydnase_cmd=("wellington_bootstrap.py -A ${fdr_limit:-} -p 12" 
             "${tname}.srt.bam ${cname}.srt.bam preprocessed.bed"
             "${name}/${tname}_specific_fp.txt"
             "${name}/${cname}_specific_fp.txt")

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$workdir/$logdir/$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:ncpus=12:mem=24gb
		#PBS -j oe
		#PBS -N fp.$name
		#PBS -q med-bio
		#PBS -o ${std_log}
		${depend:-}

		# load modules
		module load anaconda3/personal
		source activate pydnase

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy input & index to scratch
		cp -L ${workdir}/${treatment}* .
		cp -L ${workdir}/${control}* .
		cp -L ${workdir}/${regions}* .

		# pre process bed
		cut -f 1-3 $(basename $regions) > preprocessed.bed 

		# pre process bam
		samtools sort -@ 12 -o ${tname}.srt.bam $(basename ${treatment})
		samtools index ${tname}.srt.bam
		samtools sort -@ 12 -o ${cname}.srt.bam $(basename ${control})
		samtools index ${cname}.srt.bam 

		# run wellington
		mkdir -p ${name}
		${pydnase_cmd[@]} &>> $out_log
		cp -r ${name} $workdir/$outdir/

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR &>> $out_log
		ls -lhAR 
	EOS
)
echo "$script" > $pbs_log

# submit job, echo id and exit
jobid=$(qsub "$pbs_log")
echo "JOBID: $jobid"
exit 0
