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

# help message
help_message="

Wrapper to align reads with BWA

usage:
    bash $scr_name [-options] -fq1 <FASTQ>
purpose:
    # Align reads in FASTQ format with BWA and generate alignment metrics.
required arguments:
    -fq1|--fastq1 : FASTQ filename of read 1 
    --fasta : reference FASTA used in index generation
    --index : BWA index prefix to align against 
optional arguments:
    -fq2|--fastq2 : FASTQ filename of read 2 (if paired end)
    -o|--outdir : output directory for bam files (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output files (default = FASTQ filename)
    -e|--extend : bp extension for generating coverage track (default = 200)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL)
bam header arguments:
    --sm : sample name (default = name)
    --id : read group ID (usually flow cell + lane)
    --lb : library name (name unique to each library)
    --pu : platform unit (usually flow cell + barcode + lane)
    --pl : sequencing platform sequenced on, valid options:
           [ILLUMINA,PACBIO,IONTORRENT,ONT,CAPILLARY,LS454,SOLID,HELICOS]
    --cn : sequencing centre sequencing was performed at 
    --dt : run date (Iso8601Date)
    --pi : predicted median insert size (e.g. 200)
    --pm : platform model - further discription of platform
additional info:
    # qc/log output directories inherit from --outdir unless specified
    # any of the bam header arguments provided will be included in the header
    # if platform is ILLUMINA - will automatically extract PU and ID from read name 
    # check and depend args useful for jobs scheduling

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -fq1|--fastq1)
            fq1=$2
            shift
            ;;
        --fasta)
            fasta=$2
            shift
            ;;
        --index)
            index=$2
            shift
            ;;
        -fq2|--fastq2)
            fq2=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -q|--qcdir)
            logdir=$2
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
        -e|--extend)
            extend=$2
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        --depend)
            depend=$2
            shift
            ;;
        --sm)
            sm="SAMPLE_NAME=$2"
            shift
            ;;
        --id)
            id="READ_GROUP_NAME=$2"
            shift
            ;;
        --cn)
            cn="SEQUENCING_CENTER=$2"
            shift
            ;;
        --dt)
            dt="RUN_DATE=$2"
            shift
            ;;
        --lb)
            lb="LIBRARY_NAME=$2"
            shift
            ;;
        --pi)
            pi="PREDICTED_INSERT_SIZE=$2"
            shift
            ;;
        --pl)
            pl="PLATFORM=$2"
            shift
            ;;
        --pm)
            pm="PLATFORM_MODEL=$2"
            shift
            ;;
        --pu)
            pu="PLATFORM_UNIT=$2"
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
if [[ -z "${fq1:-}" ]]; then
    printf "\nERROR: --fastq1 argument required\n"
    echo "$help_message"; exit 1
elif [[ -z "${fasta:-}" ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
elif [[ -z "${index:-}" ]]; then
    printf "\nERROR: --index argument required\n"
    echo "$help_message"; exit 1
fi

# check files 
if [[ "${check:-}" = yes ]]; then
    if [[ ! -r "${workdir:-}/${fq1:-}" ]]; then
        printf "\nERROR: FASTQ r1 cannot be read: %s/%s\n" $workdir $fq1
        echo "$help_message"; exit 1
    elif [[ ! -z ${fq2:-} ]] & [[ ! -r "${workdir:-}/${fq2:-}" ]]; then
        printf "\nERROR: FASTQ r2 cannot be read: %s/%s\n" $workdir $fq2
        echo "$help_message"; exit 1
    elif [[ ! -r "${workdir:-}/${fasta:-}" ]]; then
        printf "\nERROR: FASTA file cannot be read: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    elif [[ ! -e "${workdir:-}/${index:-}" ]]; then
        printf "\nERROR: Index file does not exist: %s/%s\n" $workdir $index
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
if [[ -z "${trackdir:-}" ]]; then
    trackdir=$outdir
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir/fastqc
mkdir -p $workdir/$qcdir/metrics
mkdir -p $workdir/$trackdir

# set basenames
fq1_base=$(basename "$fq1")
fasta_base=$(basename "$fasta")
fasta_prefix=${fasta%.*}
index_base=$(basename "$index")

# set PE/SE arguments
if [[ ! -z "${fq2:-}" ]]; then
    fq2_cp="cp $workdir/$fq2 ."
    fq2_base=$(basename "$fq2")
    fq2_fq2sam="FASTQ2=$fq2_base"
    fq2_bowtie="-1 $fq1_base -2 $fq2_base"
else
    fq2_bowtie="$fq1_base"
fi

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    name=${fq1_base%%.*}
fi

# set sample name to name if not provided
if [[ -z "${sm:-}" ]]; then
    sm="SAMPLE_NAME=$name"
fi

# extract read group info (if not provided)
if [[ ${pl:-} = "PLATFORM=ILLUMINA" ]]; then
    read_name=$(gzip -dc "$workdir/$fq1" | head -n 1)
    fcid=$(echo $read_name | cut -d ":" -f 3) # flow cell ID
    fcln=$(echo $read_name | cut -d ":" -f 4) # flow cell lane
    ridx=$(echo $read_name | cut -d ":" -f 10) # barcode/index
    if [[ -z ${id:-} ]]; then id="READ_GROUP_NAME=$fcid.$fcln"; fi
    if [[ -z ${pu:-} ]]; then pu="PLATFORM_UNIT=$fcid.$fcln.$ridx"; fi
fi

# set commands
fastqc_command=("fastqc --noextract -o fastqc -t 20 $fq1_base ${fq2_base:-}")
fq2sam_command=("java -jar -Xmx32G -jar /apps/picard/2.6.0/picard.jar FastqToSam" 
                "FASTQ=$fq1_base" 
                "OUTPUT=$name.unaligned.bam" 
                "${fq2_fq2sam:-} ${id:-} ${lb:-} ${pl:-}" 
                "${pu:-} ${pm:-} ${cn:-} ${dt:-} ${pi:-} ${sm:-}")
madapt_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkIlluminaAdapters"
                "I=$name.unaligned.bam"
                "O=$name.markadapters.bam"
                "M=$name.mark_adapters_metrics")
sam2fq_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar SamToFastq"
                "I=$name.markadapters.bam" 
                "FASTQ=$name.fq" 
                "CLIPPING_ATTRIBUTE=XT" 
                "CLIPPING_ACTION=2" 
                "INTERLEAVE=true" 
                "NON_PF=true")
align_command=("bowtie $index_base $fq2_bowtie --sam --best"
                "--threads 20 > $name.sam")
merge_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MergeBamAlignment"
                "R=$fasta_base" 
                "UNMAPPED_BAM=$name.unaligned.bam" 
                "ALIGNED=$name.aligned.sam"
                "O=$name.merge.bam"
                "CREATE_INDEX=true"
                "ADD_MATE_CIGAR=true"
                "CLIP_ADAPTERS=false"
                "CLIP_OVERLAPPING_READS=true"
                "INCLUDE_SECONDARY_ALIGNMENTS=true" 
                "MAX_INSERTIONS_OR_DELETIONS=-1"
                "PRIMARY_ALIGNMENT_STRATEGY=MostDistant" 
                "ATTRIBUTES_TO_RETAIN=XS")
mdup_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkDuplicates"
            	"I=$name.merge.bam"
            	"O=$name.bam"
            	"M=$name.mark_duplicate_metrics"
            	"CREATE_INDEX=true")
alstat_command=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar CollectAlignmentSummaryMetrics"
            	"R=$fasta_base"
            	"I=$name.bam"
            	"O=$name.alignment_summary_metrics")


# set logfile name
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
out_log=$workdir/$logdir/$name.$scr_name.out.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log

# run job
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=30:00:00
		#PBS -l select=1:mem=40gb:ncpus=20
		#PBS -j oe
		#PBS -N $name.bwa
		#PBS -q med-bio
		#PBS -o $std_log

		# load modules
		module load fastqc/0.11.2
		module load samtools/1.2
		module load java/jdk-8u66
		module load picard/2.6.0
		module load bowtie/1.1.1
        
		printf "\nSTART: %s\n" date

		# copy resource files to scratch
		mkdir -p fastqc
		cp -L $workdir/$index* .
		cp -L $workdir/$fasta_base* .
		cp $workdir/$fq1 .
		${fq2_cp:-}

		# run fastqc
		printf "\nRunning FASTQC\n" >> $out_log 
		${fastqc_command[@]} >> $out_log
		cp -r fastqc/* $workdir/$qcdir/fastqc/

		# covert fastq to ubam
		printf "\nConverting to uBAM\n" >> $out_log
		${fq2sam_command[@]} >> $out_log

		# mark adapters
		printf "\nMarking adapters\n" >> $out_log
		${madapt_command[@]} >> $out_log
		cp $name.mark_adapters_metrics $workdir/$qcdir/metrics/

		# convert uBAM to interleaved fastq
		printf "\nConverting to FASTQ\n" >> $out_log
		${sam2fq_command[@]} >> $out_log

		# align to genome with Bowtie1
		printf "\nAligning to genome\n" >> $out_log
		${align_command[@]}

		# merge uBAM and aligned
		printf "\nMerging aligned with uBAM\n" >> $out_log
		${merge_command[@]} >> $out_log

		# mark dup
		printf "\nMarking duplicates\n" >> $out_log
		${mdup_command[@]} >> $out_log
		cp $name.mark_duplicate_metrics $workdir/$qcdir/metrics/

		# alignment metrics
		printf "\nCollecting alignment metrics\n" >> $out_log
		${alstat_command[@]} >> $out_log
		cp $name.alignment_summary_metrics $workdir/$qcdir/metrics/

		# copy final bam to outdir
		samtools index $name.bam
		cp $name.bam* $workdir/$outdir/
 
		printf "\nEND: %s\n" date
		ls -lhAR
	EOS
) 
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0 
