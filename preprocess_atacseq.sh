#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='on'
scr_name=$(basename "$0")
chrlist='hs'

# help message
help_message="

Takes an aligned BAM and filters & shifts tags.

usage:
    bash $scr_name [-options] -i <BAM>
required arguments:
    -i|--input : input BAM file [BAM]
    -bl|--blacklist : input blacklist regions [BED]
optional arguments:
    -chr|--chrlist : chr to include - either species or comma sep list of chr
                     [mm|hs|list] (default = 'hs')
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on|off] (default = on)
    --depend : list of PBS dependencies (default = NULL)
additional info:
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified
    # if --chrlist = mm/hs, chr included = chr[0-9XY]+

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--input)
            input=$2
            shift
            ;;
        -bl|--blacklist)
            blacklist=$2
            shift
            ;;
        -chr|--chrlist)
            chrlist=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -q|--qcdir)
            qcdir=$2
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
            printf "ERROR: Unrecognised argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z ${input:-} ]]; then
    printf "\nERROR: --input argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${blacklist:-} ]]; then
    printf "\nERROR: --blacklist argument required\n"
    echo "$help_message"; exit 1
fi

# check files 
if [[ "${check:-}" = on ]]; then
    if [[ ! -r ${workdir}/${input} ]]; then
        printf "\nERROR: Input BAM cannot be read: %s/%s\n" $workdir $input
        echo "$help_message"; exit 1
    elif [[ ! -r ${workdir}/${blacklist} ]]; then
        printf "\nERROR: Input blacklist BED file cannot be read\n" $workdir $blacklist
        echo "$help_message"; exit 1
    fi
fi

# set output directories
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
if [[ -z "${qcdir:-}" ]]; then
    qcdir=$outdir
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    name=$(basename ${input})
    name=${name%%.*}
fi

# set chr list
if [[ ${chrlist} = 'hs' ]]; then
    chr_array=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" 
               "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" 
               "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")
elif [[ ${chrlist} = 'mm' ]]; then
    chr_array=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9"
               "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19"
               "chrX" "chrY")
else
    IFS=',' read -r -a chr_array <<< ${chrlist}
fi

# set commands
chrfilt_cmd=("samtools view -q 10 -o ${name}.chrfilt.bam ${name}.sort.bam ${chr_array[@]}")
blfilt_cmd=("bedtools intersect -v -a ${name}.chrfilt.bam"
            "-b $(basename ${blacklist}) > ${name}.blfilt.bam")
filt_cmd=("samtools view -F 1548 -u ${name}.blfilt.bam |"
           "samtools sort -@ 8 -n -o ${name}.dedup.bam -")
fixmate_cmd=("samtools fixmate -r ${name}.dedup.bam ${name}.fix.bam")
bedpe_cmd=("bedtools bamtobed -bedpe -mate1 -i ${name}.fix.bam | gzip -nc > ${name}.bedpe.gz")
shift_cmd=("zcat -f ${name}.bedpe.gz |"
           "awk 'BEGIN {OFS = \"\t\"}"
           "{if (\$9 == \"+\") {\$2 = \$2 + 4}"
           "else if (\$9 == \"-\") {\$3 = \$3 - 5}"
           "if (\$10 == \"+\") {\$5 = \$5 + 4}"
           "else if (\$10 == \"-\") {\$6 = \$6 - 5} print \$0}' |"
           "gzip -nc > ${name}.tn5.bedpe.gz")
bedtota_cmd=("zcat ${name}.tn5.bedpe.gz |"
             "awk 'BEGIN{OFS=\"\t\"}"
             "{printf \"%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n\","
             "\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' |"
             "gzip -nc > ${name}.tn5.tagAlign.gz")
spr_cmd=("zcat ${name}.tn5.bedpe.gz |"
         "shuf | split -d -l \${nlines} - ${name}.sub")
spr1tota_cmd=("awk 'BEGIN{OFS=\"\t\"}"
              "{printf \"%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n\","
              "\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' ${name}.sub00 |"
              "gzip -nc > ${name}.sub1.tn5.tagAlign.gz")
spr2tota_cmd=("awk 'BEGIN{OFS=\"\t\"}"
              "{printf \"%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n\","
              "\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' ${name}.sub01 |"
              "gzip -nc > ${name}.sub2.tn5.tagAlign.gz")

# set log file names
scr_name=$(basename ${0} .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=48:00:00
		#PBS -l select=1:mem=18gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.atac
		#PBS -o $std_log
		${depend:-}

		# load modules
		module load anaconda3/personal
		source activate sambedtools

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
		
		# copy files to scratch
		cp -L ${workdir}/${input%.*}* .
		cp -L ${workdir}/${blacklist} .

		# make sure sorted and indexed
		samtools sort -o ${name}.sort.bam $(basename ${input})
		samtools index ${name}.sort.bam

		# filter
		echo "Filtering BAM"
		${chrfilt_cmd[@]}
		${blfilt_cmd[@]}
		${filt_cmd[@]}

		# fix mate
		echo "Fixing mates"
		${fixmate_cmd[@]}

		# bedpe conversion
		echo "Converting to BEDPE"
		${bedpe_cmd[@]} 

		# tag shifting
		echo "Tn5 shifting BEDPE reads"
		${shift_cmd[@]}
		cp ${name}.tn5.bedpe.gz $workdir/$outdir/

		# tag conversion
		echo "Converting filtered, shifted BEDPE to tagAlign"
		${bedtota_cmd[@]} 
		cp ${name}.tn5.tagAlign.gz $workdir/$outdir/

		# subsampling bedpe
		echo "Subsampling BEDPE to 2 SPR" 
		nlines=\$(zcat ${name}.tn5.bedpe.gz | wc -l )
		nlines=\$(( (nlines + 1) / 2 ))
		${spr_cmd[@]}
		${spr1tota_cmd[@]}
		cp ${name}.sub1.tn5.tagAlign.gz $workdir/$outdir/  
		${spr2tota_cmd[@]} 
		cp ${name}.sub2.tn5.tagAlign.gz $workdir/$outdir/ 

		# sort and index bam
		samtools sort -o ${name}.filt.bam ${name}.fix.bam
		samtools index ${name}.filt.bam
		samtools flagstat ${name}.filt.bam > ${name}.filt.flagstats.txt
		cp ${name}.filt.bam* $workdir/$outdir/
		cp ${name}.filt.flagstats.txt $workdir/$qcdir/

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
		ls -lhAR 
	EOS
) 
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0 
