#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='on'
short='off'
scr_name=$(basename "$0")

# help message
help_message="
Wrapper to align reads with BWA

usage:
    bash $scr_name [-options] -fq1 <FASTQ>
required arguments:
    -fq1|--fastq1 : FASTQ filename of read 1 
    -F|--fasta : reference FASTA used in index generation
    -I|--index : BWA index prefix to align against 
optional arguments:
    -fq2|--fastq2 : FASTQ filename of read 2 (if paired end)
    -s|--short : align in short read SE mode [on|off] (default = off)
    -n|--name : name prefix for output files (default = FASTQ filename)
    -o|--outdir : output directory for bam files (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files [on,off] (default = on)
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
    # all paths should be relative to working directory
    # check and depend options used for job scheduling
    # log/qc output directories inherit from --outdir unless specified
    # any of the bam header arguments provided will be included in the header
    # if platform is ILLUMINA - will automatically extract PU and ID from read name 

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -fq1|--fastq1)
            fq1=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -I|--index)
            index=$2
            shift
            ;;
        -fq2|--fastq2)
            fq2=$2
            shift
            ;;
        -s|--short)
            short=$2
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
            printf "ERROR: Unrecognised argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z ${fq1:-} ]]; then
    printf "\nERROR: --fastq1 argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${index:-} ]]; then
    printf "\nERROR: --index argument required\n"
    echo "$help_message"; exit 1
fi

# check files 
if [[ "${check:-}" = on ]]; then
    if [[ ! -r $fq1 ]]; then
        printf "\nERROR: FASTQ r1 cannot be read: %s/%s\n" $workdir $fq1
        echo "$help_message"; exit 1
    elif [[ -n "${fq2:-}" ]] && [[ ! -r "${fq2:-}" ]]; then
        printf "\nERROR: FASTQ r2 cannot be read: %s/%s\n" $workdir ${fq2:-}
        echo "$help_message"; exit 1
    elif [[ ! -r $fasta ]]; then
        printf "\nERROR: FASTA file cannot be read: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    elif [[ ! -e $index.bwt ]]; then
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

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir

# set basenames
fq1_base=$(basename "$fq1")
fasta_base=$(basename "$fasta")
fasta_prefix=${fasta%.*}
index_base=$(basename "$index")

# set PE/SE arguments
if [[ -n "${fq2:-}" ]]; then
    fq2_cp="cp $workdir/$fq2 ."
    fq2_base=$(basename "$fq2")
    fq2_fq2sam="FASTQ2=$fq2_base"
    fq2_bwamem="-p"
fi

# extract filename prefix if not provided
if [[ -z "${name:-}" ]]; then
    echo "WARNING: name not set, extracting from: $fq1_base"
    name=${fq1_base%%.*}
fi

# set sample name to name if not provided
if [[ -z "${sm:-}" ]]; then
    sm="SAMPLE_NAME=$name"
fi

# extract read group info (if not provided)
if [[ ${pl:-} = "PLATFORM=ILLUMINA" ]]; then
    IFS=":" read -r -a read_name <<< $(gzip -dc $fq1 | head -n 1)
    fcid=${read_name[2]} # flow cell ID
    fcln=${read_name[3]} # flow cell lane
    ridx=${read_name[9]} # barcode/index
    if [[ -z ${id:-} ]]; then id="READ_GROUP_NAME=$fcid.$fcln"; fi
    if [[ -z ${pu:-} ]]; then pu="PLATFORM_UNIT=$fcid.$fcln.$ridx"; fi
fi

# set commands
fastqc_command=("fastqc --noextract --dir tmp -o fastqc -t 20 ${fq1_base} ${fq2_base:-}")
fq2sam_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar FastqToSam" 
                "FASTQ=${fq1_base}" 
                "OUTPUT=${name}.bam"
                "TMP_DIR=tmp"
                "${fq2_fq2sam:-} ${id:-} ${lb:-} ${pl:-}" 
                "${pu:-} ${pm:-} ${cn:-} ${dt:-} ${pi:-} ${sm:-}")
madapt_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkIlluminaAdapters"
                "I=${name}.bam"
                "O=${name}.madapt.bam"
                "M=${name}.mark_adapters_metrics"
                "TMP_DIR=tmp")
sam2fq_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar SamToFastq"
                "I=${name}.madapt.bam" 
                "FASTQ=${name}.unaligned.fq" 
                "CLIPPING_ATTRIBUTE=XT" 
                "CLIPPING_ACTION=2" 
                "INTERLEAVE=true" 
                "NON_PF=true"
                "TMP_DIR=tmp")
if [[ $short = 'on' ]]; then
    sai_command=("bwa aln ${index_base} ${name}.unaligned.fq > ${name}.unaligned.sai")
    align_command=("bwa samse ${index_base} ${name}.unaligned.sai"
                   "${name}.unaligned.fq > ${name}.aligned.sam")
else
    align_command=("bwa mem -M -t 20 ${fq2_bwamem:-} ${index_base}"
                   "${name}.unaligned.fq > ${name}.aligned.sam")
fi
merge_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MergeBamAlignment"
                "R=${fasta_base}" 
                "UNMAPPED_BAM=${name}.madapt.bam" 
                "ALIGNED=${name}.aligned.sam"
                "O=${name}.bam"
                "SORT_ORDER=coordinate"
                "CREATE_INDEX=true"
                "ADD_MATE_CIGAR=true"
                "CLIP_ADAPTERS=false"
                "CLIP_OVERLAPPING_READS=true"
                "INCLUDE_SECONDARY_ALIGNMENTS=true" 
                "MAX_INSERTIONS_OR_DELETIONS=-1"
                "PRIMARY_ALIGNMENT_STRATEGY=MostDistant" 
                "ATTRIBUTES_TO_RETAIN=XS"
                "TMP_DIR=tmp")
mdup_command=("/apps/sambamba/0.6.5/bin/sambamba_v0.6.5 markdup" 
            	"--nthreads=20"
            	"--tmpdir=tmp"
                "${name}.bam ${name}.mdup.bam")
alstat_command=("java -Xmx16g -jar /apps/picard/2.6.0/picard.jar CollectAlignmentSummaryMetrics"
            	"R=${fasta_base}"
            	"I=${name}.bam"
            	"O=${name}.alignment_summary_metrics"
                "TMP_DIR=tmp")


# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$workdir/$logdir/$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=32gb:ncpus=20
		#PBS -j oe
		#PBS -N $name.bwa
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load modules
		module load fastqc
		module load samtools
		module load java
		module load picard
		module load bio-bwa
		module load sambamba

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log
		
		# copy resource files to scratch
		mkdir -p fastqc &>> $out_log
		mkdir -p tmp
		cp -L $workdir/$index* . &>> $out_log
		cp -L $workdir/$fasta_prefix* . &>> $out_log
		cp $workdir/$fq1 . &>> $out_log
		${fq2_cp:-}

		# run fastqc
		printf "\nRunning FASTQC:\n" >> $out_log 
		${fastqc_command[@]} &>> $out_log
		cp -r fastqc/*fastqc.zip $workdir/$qcdir

		# covert fastq to ubam
		printf "\nConverting to uBAM:\n" >> $out_log
		${fq2sam_command[@]} &>> $out_log

		# mark adapters
		printf "\nMarking adapters:\n" >> $out_log
		${madapt_command[@]} &>> $out_log
		cp -r $name.mark_adapters_metrics $workdir/$qcdir &>> $out_log

		# convert uBAM to interleaved fastq
		printf "\nConverting to FASTQ:\n" >> $out_log
		${sam2fq_command[@]} &>> $out_log

		# align to genome with bwa
		printf "\nAligning to genome:\n" >> $out_log
		${sai_command[@]:-}
		${align_command[@]}

		# merge uBAM and aligned
		printf "\nMerging aligned with uBAM:\n" >> $out_log
		${merge_command[@]} &>> $out_log

		# mark dup
		printf "\nMarking duplicates:\n" >> $out_log
		${mdup_command[@]} &>> $out_log

		# sort & index final bam & copy to outdir
		samtools sort -o ${name}.bam ${name}.mdup.bam
		samtools index $name.bam &>> $out_log
		cp -r $name.bam* $workdir/$outdir/ &>> $out_log

		# alignment metrics
		printf "\nCollecting alignment metrics:\n" >> $out_log
		${alstat_command[@]} &>> $out_log
		cp -r $name.alignment_summary_metrics $workdir/$qcdir &>> $out_log
 
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR &>> $out_log
		ls -lhAR 
		cp -r $out_log $workdir/$logdir/
        
	EOS
) 
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0 
