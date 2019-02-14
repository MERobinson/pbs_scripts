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

# default ChIA-PET2 args
linkerA='GTTGGATAAG'
linkerB='GTTGGAATGT'
start_phase=1
mode=0
err=0
keepempty=0
thread=16
short=0
mapq=30
length=15

# help message
help_message="
Wrapper to run ChIA-PET2 pipeline from aligning reads to calling interactions

usage:
    bash $scr_name [-options] -fq1 <FASTQ1> -fq2 <FASTQ2>
required arguments:
    -fq1|--fastq1 : FASTQ filename of read 1
    -fq2|--fastq2 : FASTQ filename of read 2
    -C|--chrom : chromosome sizes file
    -I|--index : BWA index prefix 
optional arguments:
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -ld|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL)
optional ChIA-PET2 arguments:
    --linkerA : sequence of linker A (default = GTTGGATAAG)
    --linkerB : sequence of linker B (default = GTTGGAATGT)
    --start : start phase [1-8] (default=1):
         1=trim, 2=map, 3=PETs, 4=peaks, 5=interactions, 6=qc, 7=stats, 8=phase
    --mode : 0=A/B, 1=bridge, 2=enzyme [0,1,2] (default = 0)
    --err : maximum mismatches in linker sequence (default = 0)
    --keepempty : no. of linker-empty reads to keep [0,1,2] (default = 0)
    --thread : number of threads (default = 16)
    --short : flag for if reads are < 70bp [0,1] (default = 0) 
    --mapq : mapq cutoff (default = 30)
    --length : min. length of reads after trimming (default = 15)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified
    # if multiple runs, provide commma seperated list of fastqs to merge pre-alignment

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -fq1|--fastq1)
            fq1=$2
            shift
            ;;
        -fq2|--fastq2)
            fq2=$2
            shift
            ;;
        -C|--chrom)
            chrom_size=$2
            shift
            ;;
        -I|--index)
            index=$2
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
        --check)
            check=$2
            shift
            ;;
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        -A|--linkerA)
            linkerA="$2"
            shift
            ;;
        -B|--linkerB)
            linkerB="$2"
            shift
            ;;
        --start)
            start_phase="$2"
            shift
            ;;
        --mode)
            mode="$2"
            shift
            ;;
        --err)
            err="$2"
            shift
            ;;
        --keepempty)
            keepempty="$2"
            shift
            ;;
        --thread)
            thread="$2"
            shift
            ;;
        --short)
            short="$2"
            shift
            ;;
        --mapq)
            mapq="$2"
            shift
            ;;
        --length)
            length="$2"
            shift
            ;;
        *)
            printf "Error: Illegal argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z ${fq1:-} ]]; then
    printf "\nERROR: --fastq1 argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${fq2:-} ]]; then
    printf "\nERROR: --fastq2 argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${chrom_size:-} ]]; then
    printf "nERROR: --chrom argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${index:-} ]]; then
    printf "\nERROR: --index argument required\n"
    echo "$help_message"; exit 1
fi

# check files 
if [[ ${check:-} = yes ]]; then
    if [[ ! -r ${chrom_size:-} ]]; then
        printf "\nERROR: Chromosome sizes file cannot be read: %s/%s\n" $workdir $chrom_size
        echo "$help_message"; exit 1
    elif [[ ! -e ${workdir:-}/${index:-} ]]; then
        printf "\nERROR: Index files do not exist: %s/%s\n" $workdir $index
        echo "$help_message"; exit 1
    fi
fi

# if multiple fastq set concat args
IFS=',' read -r -a fq1_array <<< "$fq1"
IFS=',' read -r -a fq2_array <<< "$fq2"
if [[ ${#fq1_array[@]} > 1 ]] && [[ ${#fq1_array[@]} = ${#fq2_array[@]} ]]; then
    fq1_merge="cat"; fq2_merge="cat"
    fq_ext=${fq1_array[0]#*.}
    for idx in ${!fq1_array[@]}; do
        if [[ $check = 'yes' ]] && [[ ! -r ${fq1_array[idx]} ]]; then
            printf "\nERROR: FASTQ r1 cannot be read: %s/%s\n" $workdir ${fq1_array[idx]}
            echo "$help_message"; exit 1
        elif [[ $check = 'yes' ]] && [[ ! -r ${fq2_array[idx]} ]]; then
            printf "\nERROR: FASTQ r2 cannot be read: %s/%s\n" $workdir ${fq2_array[idx]}
            echo "$help_message"; exit 1
        fi
        fq1_merge="${fq1_merge} $(basename ${fq1_array[idx]})"
        fq2_merge="${fq2_merge} $(basename ${fq2_array[idx]})"
        fq1_cp="${fq1_cp:-}; cp $workdir/${fq1_array[idx]} ."
        fq2_cp="${fq2_cp:-}; cp $workdir/${fq2_array[idx]} ."
    done
    fq1_merge="${fq1_merge} > fastq1.$fq_ext"; fq2_merge="${fq2_merge} > fastq2.$fq_ext"
    fq1_cp=${fq1_cp#*;}; fq2_cp=${fq2_cp#*;}
    fq1="fastq1.$fq_ext"; fq2="fastq2.$fq_ext"

else
    if [[ $check = 'yes' ]] &&[[ ! -r ${fq1} ]]; then
        printf "\nERROR: FASTQ r1 cannot be read: %s/%s\n" $workdir $fq1
        echo "$help_message"; exit 1
    elif [[ $check = 'yes' ]] && [[ ! -r ${fq2} ]]; then
        printf "\nERROR: FASTQ r2 cannot be read: %s/%s\n" $workdir $fq2
        echo "$help_message"; exit 1
    fi
    fq1_cp="cp $workdir/$fq1 ."; fq2_cp="cp $workdir/$fq2 ."
fi

# set output directories
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set basenames
fq1_base=$(basename "$fq1")
fq2_base=$(basename "$fq2")
chrom_base=$(basename "$chrom_size")
index_base=$(basename "$index")

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    name=${fq1_base%%.*}
fi

# set commands
chiapet_command=("ChIA-PET2 -f $fq1_base -r $fq2_base -s $start_phase"
                 "-A $linkerA -B $linkerB"
                 "-b $chrom_base -g $index_base"
                 "-o tmp -n $name"
                 "-t $thread -m $mode"
                 "-e $err -k $keepempty"
                 "-d $short -Q $mapq -l $length")

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=48:00:00
		#PBS -l select=1:mem=64gb:ncpus=$thread
		#PBS -j oe
		#PBS -N $name.chiapet
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load modules
		module load bio-bwa/0.7.5a
		module load samtools/1.2
		module load macs/2.1.0
		module load bedtools/2.25
		module load R/3.3.2
		export PATH=/work/mrobinso/tools/ChIA-PET2/bin/:$PATH

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log
	
		# temp
		mkdir tmp
		cp -rL $workdir/$outdir/$name* tmp/
	
		# copy resource files to scratch
		mkdir tmp
		cp -L $workdir/$index* . &>> $out_log
		cp -L $workdir/$chrom_size . &>> $out_log
		${fq1_cp} &>> $out_log
		${fq2_cp} &>> $out_log

		# concat fastq if nec
		${fq1_merge:-}
		${fq2_merge:-}

		# run chiapet2 
		printf "\nRunning ChIA-PET2:\n" >> $out_log 
		${chiapet_command[@]} &>> $out_log

		# copy final bam to outdir
		cp -r tmp/* $workdir/$outdir/ &>> $out_log
 
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR &>> $out_log
		ls -lhAR 
		cp $out_log $workdir/$logdir/
        
	EOS
) 
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
