#!/bin/bash

# default arg
workdir=$PWD
fqdir=fastq
resdir=resources
trackdir=tracks
bamdir=bam
logdir=logs
qcdir=qc
fasta=genome.fa
index=genome.fa
extend=200
pl="PLATFORM=ILLUMINA"
validate=yes

# help message
help_message="
usage:
    bash $(basename "$0") [-wrtblgneh] -f <FASTQ>
purpose:
    # Align reads in FASTQ format with BWA and generate alignment metrics.
required arguments:
    -f|--fastq : FASTQ filename of read 1 
optional arguments:
    -m|--mate : FASTQ filename of read 2 - if paired end
    -w|--workdir : working directory - used as base dir for all input/output (default = pwd)
    -i|--fqdir : input directory containing FASTQ files (default = fastq)
    -r|--resdir : directory containing index and ref FASTA resources (default = resources) 
    -t|--trackdir : output directory for genome tracks (default = tracks)
    -b|--bamdir : output directory for bam files (default = bam)
    -q|--qcdir : output directory for qc metrics (default = qc)
    -l|--logdir : output directory for log files (default = logs)
    -n|--name : name prefix for output files (default = FASTQ filename)
    -e|--extend : bp extension for generating coverage track (default = 200)
    -v|--validate : whether to validate existance of input files [yes,no] (default = yes) 
example:
    bash $(basename "$0") -g hg38 -f REH_H3K27ac_ChIPseq.fastq.gz
additional info:
    # Additional read group information can be provided to include in BAM header
    # Possible read group data to include:
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
    # if platform is ILLUMINA - will automatically extract PU and ID from read name 
    # validate argument is useful for scheduling future jobs when files dont currently exist

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -f|--fastq)
        fastq=$2
        shift
        ;;
        -m|--mate)
        mate=$2
        shift
        ;;
        -w|--workdir)
        workdir=$2
        shift
        ;;
        -i|--fqdir)
        fqdir=$2
        shift
        ;;
        -r|--resdir)
        resdir=$2
        shift
        ;;
        -b|--bamdir)
        bamdir=$2
        shift
        ;;
        -l|--logdir)
        logdir=$2
        shift
        ;;
        -q|--qcdir)
        qcdir=$2
        shift
        ;;
        -n|--name)
        name=$2
        shift
        ;;
        -v|--validate)
        validate=$2
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
        echo "Error: Illegal argument"
        echo "$help_message"
        exit 1
        ;;
    esac
    shift
done

# check required argument
if [[ -z "$fastq" ]]; then
    printf "\nERROR: FASTQ argument required\n"
    echo "$help_message"; exit 1
fi

if [[ "$validate" = yes ]]; then
    # check indir exists
    if [[ ! -d "$workdir/$fqdir" ]]; then
        printf "\nERROR: input directory does not exist: %s/%s/%s\n" $workdir $fqdir
        echo "$help_message"; exit 1
    fi

    # check FASTQ exists
    if [[ ! -f "$workdir/$fqdir/$fastq" ]]; then
        printf "\nERROR: input FASTQ file does not exist: %s/%s/%s\n" $workdir $fqdir $fastq
        echo "$help_message"; exit 1
    fi

    # check mate (if set)
    if [[ ! -z "$mate" ]] && [[ ! -f "$workdir/$fqdir/$mate" ]]; then
        printf "\nERROR: input mate FASTQ does not exist: %s/%s/%s\n" $workdir $fqdir $fastq
        echo "$help_message"; exit 1
    fi

    # check resource files exist
    if [[ ! -e "$workdir/$resdir/$fasta" ]]; then
        printf "\nERROR: FASTA file does not exist: %s/%s/%s\n" $workdir $resdir $fasta
        echo "$help_message"; exit 1
    fi
    if [[ ! -e "$workdir/$resdir/$index" ]]; then
        printf "\nERROR: Index file does not exist: %s/%s/%s\n" $workdir $resdir $index
        echo "$help_message"; exit 1
    fi
fi

# if mate - set PE args
if [[ ! -z "$mate" ]]; then
    pe_cp_arg="cp $workdir/$fqdir/$mate ."
    pe_fq2sam_arg="FASTQ2=$mate"
    pe_bwamem_arg="-p "
fi

# extract filename prefix if not provided
if [[ -z "$name" ]]; then
    name=${fastq%%.*}
fi

# set sample name to filename if not provided
if [[ -z "$sm" ]]; then
    sm="SAMPLE_NAME=$name"
fi

# strip fasta extension (to copy dict)
fastq_base=${fasta%%.*}

# extract read group info (if not provided)
if [[ $pl = "PLATFORM=ILLUMINA" ]]; then
    read_name=$(gzip -dc "$workdir/$fqdir/$fastq" | head -n 1)
    fcid=$(echo $read_name | cut -d ":" -f 3) # flow cell ID
    fcln=$(echo $read_name | cut -d ":" -f 4) # flow cell lane
    ridx=$(echo $read_name | cut -d ":" -f 10) # barcode/index
    if [[ -z $id ]]; then id="READ_GROUP_NAME=$fcid.$fcln"; fi
    if [[ -z $pu ]]; then pu="PLATFORM_UNIT=$fcid.$fcln.$ridx"; fi
fi

# create required dirs
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir/fastqc
mkdir -p $workdir/$qcdir/metrics
mkdir -p $workdir/$bamdir
mkdir -p $workdir/$trackdir

# run job
jobid=$(cat <<- EOS | qsub -N $name.bwa - 
	#!/bin/bash
	#PBS -l walltime=60:00:00
	#PBS -l select=1:mem=40gb:ncpus=20
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $workdir/$logdir/$name.align_fastq_bwa.log
	
	# load modules
	module load fastqc/0.11.2
	module load samtools/1.2
	module load java/jdk-8u66
	module load picard/2.6.0
	module load bio-bwa/0.7.10
        
	# copy resource files to scratch
	mkdir -p fastqc
	cp -L $workdir/$resdir/$index* .
	cp -L $workdir/$resdir/$fasta_base* .
	cp $workdir/$fqdir/$fastq .
	$pe_cp_arg

	# run fastqc
	fastqc --noextract -o fastqc -t 20 $fastq $mate
	cp -r fastqc/* $workdir/$qcdir/fastqc/
	
	# covert fastq to ubam
	java -jar -Xmx32G -jar /apps/picard/2.6.0/picard.jar FastqToSam \
		FASTQ=$fastq \
		OUTPUT=$name.unaligned.bam \
		$pe_fq2sam_arg $id $lb $pl $pu $pm $cn $dt $pi $sm

	# mark adapters
	java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkIlluminaAdapters \
		I=$name.unaligned.bam \
		O=$name.markadapters.bam \
		M=$name.mark_adapters_metrics
	cp $name.mark_adapters_metrics $workdir/$qcdir/metrics/

	# convert uBAM to interleaved fastq
	java -Xmx32G -jar /apps/picard/2.6.0/picard.jar SamToFastq \
		I=$name.markadapters.bam \
		FASTQ=$name.fq \
		CLIPPING_ATTRIBUTE=XT \
		CLIPPING_ACTION=2 \
		INTERLEAVE=true \
		NON_PF=true

	# align to genome with BWA
	bwa mem -M -t 20 $pe_bwamem_arg\
		$index \
		$name.fq \
		> $name.aligned.sam

	# merge uBAM and aligned
	java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MergeBamAlignment \
		R=$fasta \
		UNMAPPED_BAM=$name.unaligned.bam \
		ALIGNED=$name.aligned.sam \
		O=$name.merge.bam \
		CREATE_INDEX=true \
		ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		ATTRIBUTES_TO_RETAIN=XS

	# mark dup
	java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkDuplicates \
		I=$name.merge.bam \
		O=$name.bam \
		M=$name.mark_duplicate_metrics \
		CREATE_INDEX=true
	cp $name.mark_duplicate_metrics $workdir/$qcdir/metrics/

	# alignment metrics
	java -Xmx16g -jar /apps/picard/2.6.0/picard.jar CollectAlignmentSummaryMetrics \
		R=$fasta \
		I=$name.bam \
		O=$name.alignment_summary_metrics
	cp $name.alignment_summary_metrics $workdir/$qcdir/metrics/

	# copy final bam to outdir
	samtools index $name.bam
	cp $name.bam* $workdir/$bamdir/

	ls -lhAR
	EOS
)
