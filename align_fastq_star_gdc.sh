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
Wrapper to align reads with STAR following the GDC/ICGC settings.

usage:
    bash $(basename "$0") [-options] -f1 <FASTQ1> -F <FASTA> -I <INDEX>
required arguments:
    -fq1|--fastq1 : input FASTQ file
    -I|--index : STAR index directory
    -F|--fasta : GTF file of transcripts
optional arguments:
    -fq2|--fastq2 : mate FASTQ file if PE [FASTQ] (default = NULL)
    -n|--name : prefix for output files [string] (default = extracted from fastq)
    -o|--outdir : output directory [string] (default = '.')
    -l|--logdir : output dir for log files [string] (default = --outdir)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list of PBS dependencies [string] (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend arguments are used for job scheduling/pipelines
    # example depenency list: 'afterok:123456,afterok:123467'

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
        -I|--index)
            index=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
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
        -l|--logdir)
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
if [[ -z ${fq1:-} ]]; then
    printf "\nERROR: no --fastq1 argument provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${index:-} ]]; then
    printf "\nERROR: no --index argument provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: no --fasta argument provided\n"
    echo "$help_message"; exit 1
fi    

# get sample name if not provided
if [[ -z ${name:-} ]]; then
    name=$(basename "$fq1")
    name=${name%%.*}
fi

# set compression arg
fq_ext=${fq1##*.}
if [[ $fq_ext = "gz" ]]; then
    compress_arg="--readFilesCommand zcat"
    fq_ext="fastq.gz"
else
    compress_arg="--readFilesCommand zcat"
    fq_ext="fastq"
fi

# check if multiple fastq per read
IFS=',' read -r -a fq1_array <<< $fq1
if [[ ${#fq1_array[@]} > 1 ]]; then
    fq1_merge="cat ${fq1_array[@]} > $name.R1_merge.$fq_ext"
    fq1_cp=$(printf "cp $workdir/%s .;" ${fq1_array[@]})
    fq1=$name.R1_merge.$fq_ext
else
    fq1_cp="cp $workdir/$fq1 ."
fi
if [[ -n ${fq2:-} ]]; then
    echo "Mate input - running in PE mode"
    IFS=',' read -r -a fq2_array <<< $fq2
    if [[ ${#fq2_array} > 1 ]]; then
        fq2_merge="cat ${fq2_array[@]} > $name.R2_merge.$fq_ext"
        fq2_cp=$(printf "cp $workdir/%s .;" ${fq2_array[@]})
        fq2=$name.R2_merge.$fq_ext
    else
        fq2_cp="cp $workdir/$fq2 ."
    fi
else
    echo "No mate input - running in SE mode"
fi

# check files unless flagged
if [[ $check = yes ]]; then
    for idx in  "${!fq1_array[@]}"; do
        if [[ ! -r $workdir/${fq1_array[$idx]} ]]; then
            printf "\nERROR: FASTQ file is not readable: %s/%s\n" $workdir ${fq1_array[$idx]}
            echo "$help_message"; exit 1
        elif [[ ! -z ${fq2:-} ]] & [[ ! -r $workdir/${fq2_array[$idx]} ]]; then
            printf "\nERROR: FASTQ file is not readable: %s/%s\n" $workdir ${fq2_array[$idx]}
            echo "$help_message"; exit 1
        fi
    done
    if [[ ! -r $workdir/$index ]]; then
        printf "\nERROR: Index folder is not readable: %s/%s\n" $workdir $index
        echo "$help_message"; exit 1
    elif [[ ! -r $workdir/$gtf ]]; then
        printf "\nERROR: GTF file is not readable: %s/%s\n" $workdir $gtf
        echo "$help_message"; exit 1
    fi
fi

# get basenames/prefix
fq1_base=$(basename "$fq1")
if [[ -n ${fq2:-} ]]; then fq2_base=$(basename "$fq2"); fi
index_base=$(basename "$index")
fasta_base=$(basename "$fasta")

# setup output dir
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set commands
star_first_pass=("STAR --runMode alignReads"
                "--runThreadN 20"
                "--genomeDir ${index_base}"
                "--readFilesIn ${fq1_base} ${fq2_base:-}"
                "--outFilterMultimapScoreRange 1"
                "--outFilterMultimapNmax 20"
                "--outFilterMismatchNmax 10"
                "--alignIntronMax 500000"
                "--alignMatesGapMax 1000000"
                "--sjdbScore 2"
                "--alignSJDBoverhangMin 1"
                "--outFilterMatchNminOverLread 0.33"
                "--outFilterScoreMinOverLread 0.33"
                "--sjdbOverhang 100"
                "--outSAMstrandField intronMotif"
                "--outSAMtype None"
                "--outSAMmode None"
                "${compress_arg:-}")
index_generation=("STAR --runMode genomeGenerate"
                  "--genomeDir first_pass"
                  "--genomeFastaFiles ${fasta_base}"
                  "--sjdbOverhang 100"
                  "--runThreadN 20"
                  "--sjdbFileChrStartEnd SJ.out.tab")
star_second_pass=("STAR --runMode alignReads"
                  "--genomeDir first_pass"
                  "--readFilesIn ${fq1_base} ${fq2_base:-}"
                  "--runThreadN 20"
                  "--outFilterMultimapScoreRange 1"
                  "--outFilterMultimapNmax 20"
                  "--outFilterMismatchNmax 10"
                  "--alignIntronMax 500000"
                  "--alignMatesGapMax 1000000"
                  "--sjdbScore 2"
                  "--alignSJDBoverhangMin 1"
                  "--outFilterMatchNminOverLread 0.33"
                  "--outFilterScoreMinOverLread 0.33"
                  "--sjdbOverhang 100"
                  "--outSAMstrandField intronMotif"
                  "--outSAMattributes NH HI NM MD AS XS"
                  "--outSAMunmapped Within"
                  "--outSAMtype BAM SortedByCoordinate"
                  "--outFileNamePrefix ${name}.star."
                  "${compress_arg:-}")

# set log file names
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# run PBS script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=40gb:ncpus=20
		#PBS -j oe
		#PBS -N $name.star
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# load modules
		module load samtools/1.2  &>> $out_log
		module load star/2.5.0  &>> $out_log
		module load intel-suite/2015  &>> $out_log
		module load gcc/5.4.0  &>> $out_log
		
		# copy files to scratch
		cp -rL $workdir/$index . &>> $out_log
		cp -rL $workdir/$fasta . &>> $out_log
		${fq1_cp:-} &>> $out_log
		${fq2_cp:-} &>> $out_log

		# merge fastq if multiple
		${fq1_merge:-}
		${fq2_merge:-}

		# run STAR
		mkdir first_pass
		${star_first_pass[@]} &>> $out_log
		${index_generation[@]} &>> $out_log
		${star_second_pass[@]} &>> $out_log		

		# index bam
		mv ${name}.star.Aligned.sortedByCoord.out.bam $name.star.bam 
		samtools index ${name}.star.bam
		
		# copy output to outdir
		cp ${name}.star.* $workdir/$outdir &>> $out_log

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
