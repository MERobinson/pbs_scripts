#!/user/bin/env bash
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
    bash $(basename ${0}) [-options] -b <BAM>
purpose:
    Wrapper to run dnase_wig_tracks.py from pyDNase.
required arguments:
    -b|--bam : Input reads [BAM]
    -r|--regions : Input regions [BED]
optional arguments:
    -o|--outdir : output directory for track files (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output file (default = BAM filename)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL) 
additional info:
    # all paths should be relative to the current working directory

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam=$2
            shift
            ;;
        -r|--regions)
            regions=$2
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
            printf "\nError: Illegal argument: %s %s" $1 $2
            echo "$help_message"; exit 1
        ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arguments
if [[ -z ${bam:-} ]]; then
    printf "\nERROR: --bam argument is required.\n"
    echo "$help_message"; exit 1
elif [[ -z ${regions:-} ]]; then
    printf "\nERROR: --regions argument is required.\n"
    echo "$help_message"; exit 1
fi

# check input files
if [[ ${check} = 'on' ]]; then
    if [[ ! -r ${workdir}/${bam} ]]; then
        printf "\nERROR: BAM file cannot be read: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    elif [[ ! -r ${workdir}/${regions} ]]; then
        printf "\nERROR: BED file cannot be read: %s/%s\n" $workdir $regions
        echo "$help_message"; exit 1
    fi
fi

# if no name provided extract from BAM
if [[ -z "${name:-}" ]]; then
    name=$(basename ${bam})
    name=${name%%.*}
fi

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set log file names
scr_name=$(basename $0 .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$workdir/$logdir/$name.$scr_name.out.log

# compile commands
pydnase_cmd=("dnase_wig_tracks.py preprocessed.bed"
             "${name}.srt.bam ${name}_f.wig ${name}_r.wig")

# run job
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=4gb:ncpus=1
		#PBS -j oe
		#PBS -N ${name}.dnasewig
		#PBS -q med-bio
		#PBS -o ${std_log}
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > ${out_log}

		# load modules
		module load anaconda3/personal
		source activate pydnase

		# copy resource files to scratch
        cp ${workdir}/${bam}* . &>> ${out_log}
        cp ${workdir}/${regions}* . &>> ${out_log}

		# preprocess bed
		cut -f 1-3 $(basename $regions) > preprocessed.bed
		samtools sort -o ${name}.srt.bam $(basename $bam)
		samtools index ${name}.srt.bam

		# copy resource files to scratch
		cp ${workdir}/${bam}* . &>> ${out_log}
		cp ${workdir}/${regions}* . &>> ${out_log}

		# generate coverage track
		${pydnase_cmd[@]} &>> ${out_log}

		cp ${name}_f.wig ${workdir}/${outdir}/ &>> ${out_log}
		cp ${name}_r.wig ${workdir}/${outdir}/ &>> ${out_log}

		# temp
		mv ${name}.srt.bam $(basename $bam)
		mv ${name}.srt.bam.bai $(basename $bam).bai
		cp $(basename $bam)* ${workdir}/${outdir}/

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> ${out_log}
		ls -lhAR >> ${out_log}
		ls -lhAR
	EOS
)
echo "${script}" > ${pbs_log}

# submit job
jobid=$(qsub "${pbs_log}")

# echo jobid and exit
echo "JOBID: ${jobid}"
exit 0
