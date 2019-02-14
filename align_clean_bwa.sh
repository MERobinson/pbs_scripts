#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='y'
date=`date '+%Y-%m-%d %H:%M:%S'`
bwa_version='0.7.15'

# default header args
pl="PLATFORM=ILLUMINA"

# help message
help_message="
Wrapper to convert FASTQ to cleaned, aligned BAM for WGS/WES data & generate metrics

usage:
    bash $(basename $0) [-options] -fq1 <FASTQ_R1> -fq2 <FASTQ_R2>
required arguments:
    -fq1|--fastq1 : FASTQ R1 (single sample/run)
    -F|--fasta : reference FASTA used in index generation
    -I|--index : BWA index prefix
optional arguments:
    -fq2|--fastq2 : FASTQ R2 (single sample/run)
    -o|--outdir : output directory (default = pwd)
    -n|--name : name prefix for output files (deafult = FASTQ filename)
    -q|--qcdir : output directory for qc files (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --check : whether to check input files prior to run [y|n] (default = y)
    --depend : list of PBS dependencies to pass to PBS script (default = NULL)
bam header arguments:
    --sm : sample anme (default = name)
    --id : read group ID (usually flow cell + lane)
    --lb : library name (name unique to each library)
    --pu : platform unit (usually flow cell + barcode + lane)
    --pl : sequencing platform (default = ILLUMINA), valid options;
           [ILLUMINA,PACBIO,IONTORRENT,ONT,CAPILLARY,LS454,SOLID,HELICOS]
    --cn : sequencing centre
    --dt : run date (Iso8601Date)
    --pi : predicted median insert size
    --pm : platform model (further descrition of platform)
additional info:
    # all paths should be relative to working directory
    # check and depend options useful for job scheduling
      e.g. --depend afterok:123456,afterok:123457
    # log/qc output directories inherited from outdir unless specified
    # any bam header arguments specified will be added to output bam header
    # if platform = ILLUMINA (default), will extract PU and ID from read name

"
# parse command line arguments
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
        -F|--fasta)
            fasta=$2
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
            printf "\nERROR: Unrecognised argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
	esac
	shift
done

# check required arguments
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
if [[ ${check} = 'y' ]]; then
    if [[ ! -r ${fq1} ]]; then
        printf "ERROR: file is not readable: %s/%s" $workdir ${fq1}
        echo "$help_message"; exit 1
    elif [[ -n ${fq2:-} && ! -r ${fq2:-} ]]; then
        printf "ERROR: file is not readable: %s/%s" $workdir ${fq2}
        echo "$help_message"; exit 1
    elif [[ ! ${fq1} =~ .(fq|fastq)(.gz)? ]]; then
        printf "ERROR: file extension not recognised: %s/%s" $workdir ${fq1}
        echo "$help_message"; exit 1
    elif [[ -n ${fq2:-} && ! ${fq2:-} =~ .(fq|fastq)(.gz)? ]]; then
        printf "ERROR: file extension not recognised: %s/%s" $workdir ${fq2}
        echo "$help_message"; exit 1
    elif [[ ! -r ${fasta} ]]; then
        printf "ERROR: FASTA is not readable: %s/%s" $workdir ${fasta}
        echo "$help_message"; exit 1
    elif [[ ! -r ${index}.alt ]]; then
        printf "ERROR: BWA ALT index is not readable: %s/%s" $workdir ${index}.alt
        echo "$help_message"; exit 1
    fi
fi

# check/set name
if [[ -z ${name:-} ]]; then
	name=$(basename $fq1)
    name=${name%%.*}
fi
if [[ -z ${sm:-} ]]; then
    sm="SAMPLE_NAME=$name"
fi

# check/set outdirs
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
if [[ -z "${qcdir:-}" ]]; then
    qcdir=$outdir
fi
mkdir -p $logdir
mkdir -p $qcdir/fastqc
mkdir -p $qcdir/metrics

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
    fq2_bwamem="-p"
    fq2_qc_mv="mv fastqc/${fq2_base%%.*}_fastqc.zip $workdir/$qcdir/fastqc"
fi

# extract read group info
if [[ ${pl:-} = "PLATFORM=ILLUMINA" ]]; then
    IFS=":" read -r -a read_name <<< $(gzip -dc $fq1 | head -n 1)
    fcid=${read_name[2]} # flow cell ID
    fcln=${read_name[3]} # flow cell lane
    ridx=$(echo ${read_name[9]} | grep -Eo "[CGAT]+") # barcode/index
    if [[ -z ${id:-} ]]; then id="READ_GROUP_NAME=$fcid.$fcln"; fi
    if [[ -z ${pu:-} ]]; then pu="PLATFORM_UNIT=$fcid.$fcln.$ridx"; fi
fi

# set commands
fastqc_command=("fastqc --noextract -o fastqc -t 20 $fq1_base ${fq2_base:-}")
fq2sam_command=("java -jar -Xmx32G -jar /apps/picard/2.6.0/picard.jar FastqToSam"
                "FASTQ=$fq1_base"
                "OUTPUT=$name.unaligned.bam"
                "${fq2_fq2sam:-} ${sm:-}"
                "${id:-} ${lb:-} ${pl:-} ${pu:-}"
                "${pm:-} ${cn:-} ${dt:-} ${pi:-}")
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
align_command=("bwa mem -M -t 20 ${fq2_bwamem:-} $index_base $name.fq > $name.aligned.sam")
merge_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MergeBamAlignment"
                "R=${fasta_base}"
                "UNMAPPED_BAM=${name}.unaligned.bam"
                "ALIGNED=${name}.aligned.sam"
                "O=${name}.bam"
                "VALIDATION_STRINGENCY=SILENT"
				"SORT_ORDER=unsorted"
                "ALIGNED_READS_ONLY=false"
                "ADD_MATE_CIGAR=true"
                "CLIP_ADAPTERS=false"
                "CLIP_OVERLAPPING_READS=true"
                "INCLUDE_SECONDARY_ALIGNMENTS=true"
                "MAX_INSERTIONS_OR_DELETIONS=-1"
                "PRIMARY_ALIGNMENT_STRATEGY=MostDistant"
                "ATTRIBUTES_TO_RETAIN=X0"
                "PROGRAM_RECORD_ID=bwamem"
                "PROGRAM_GROUP_VERSION=0.7.15"
                "PROGRAM_GROUP_COMMAND_LINE=\"${align_command[@]}\""
                "PROGRAM_GROUP_NAME=bwamem"
                "UNMAP_CONTAMINANT_READS=true")
mdup_command=("java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkDuplicates"
                "I=${name}.merge.bam"
                "O=tmp.bam"
                "M=${name}.mark_duplicate_metrics"
				"OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"
				"ASSUME_SORT_ORDER=queryname")

# set log file names
scr_name=$(basename $0 .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=72:00:00
		#PBS -l select=1:mem=32gb:ncpus=20
		#PBS -j oe
		#PBS -N $name.align
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		# load modules
		module load picard/2.6.0
		module load java/jdk-8u66
		module load bio-bwa/$bwa_version
		module load fastqc/0.11.5
		module load samtools/1.2

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy to scratch
		mkdir fastqc &>> $out_log
		cp -L $workdir/$fq1 . &>> $out_log
		${fq2_cp:-}
		cp -L $workdir/${fasta%%.*}* . &>> $out_log
		cp -L $workdir/$index* . &>> $out_log

		# run fastqc
		printf "\nRunning FASTQC:\n" >> $out_log 
		${fastqc_command[@]} &>> $out_log
		mv fastqc/${fq1_base%%.*}_fastqc.zip $workdir/$qcdir/fastqc &>> $out_log
		${fq2_mv:-}

		# covert fastq to ubam
		printf "\nConverting to uBAM:\n" >> $out_log
		${fq2sam_command[@]} &>> $out_log

		# mark adapters
		printf "\nMarking adapters:\n" >> $out_log
		${madapt_command[@]} &>> $out_log
		cp $name.mark_adapters_metrics $workdir/$qcdir/metrics/ &>> $out_log

		# convert uBAM to interleaved fastq
		printf "\nConverting to FASTQ:\n" >> $out_log
		${sam2fq_command[@]} &>> $out_log

		# align to genome with bwa
		printf "\nAligning to genome:\n" >> $out_log
		${align_command[@]}

		# merge uBAM and aligned
		printf "\nMerging aligned with uBAM:\n" >> $out_log
		${merge_command[@]} &>> $out_log

		# mark dup
		printf "\nMarking duplicates:\n" >> $out_log
		${mdup_command[@]} &>> $out_log
		cp $name.mark_duplicate_metrics $workdir/$qcdir/metrics/ &>> $out_log

		# copy final bam to outdir
		samtools index $name.bam &>> $out_log
		cp $name.bam* $workdir/$outdir/ &>> $out_log

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
