#!/bin/bash

# default arg
BASEDIR=$WORK
WORKDIR=$PWD
GENOME=mm10
DELIM=","
PLATFORM="ILLUMINA"

# default sample info fields
FQ1_IDX=1
FQ2_IDX=2
SMID_IDX=3
FCID_IDX=4
FCLN_IDX=5
BCID_IDX=6
CSID_IDX=7

# help message
USAGE="$(basename "$0") [-gbwdhv] -s <sample_info>

purpose:
--pipeline to align, clean, generate metrics, and call mut in ChIP-seq samples
required arguments:
	-si|--sample_info : delimited file containing sample info,
			    # see sample info format spec below
optional arguments:
	-g|--genome : genome version [hg38,mm10] (default = mm10)
	-b|--basedir : home directory (default = \$WORK)
	-w|--workdir : working directory to output all files (default = \$PWD)
	-d|--delim : delimiter used in sample info file	(default = ",")
	-sl|--sample_list : comma separated list of sample names - specifies a subset  
			    to analyse from sample info (default = all samples)
	-sc|--seq_centre : optional sequencing centre flag for readgroup (default = none)
	-pl|--platform : platform [ILLUMINA,SOLID,LS454,HELICOS,PACBIO] (default = ILLUMINA)
	-po|--pon : panel of normals (default = none)
	-v|--verbose : print additional information to STDOUT [0,1] (default = 0)
	-h|--help : print this help message and	exit
example:
	$(basename "$0") -g hg38 -si sample_info.csv

sample info format:
--the sample info file can be any delimited text file with the following columns:
	-FASTQ_read_1 - file name of forward reads in FASTQ format (incl extension) [required]
	-FASTQ_read_2 - file name of reverse reads in FASTQ format (incl extension) [optional]
	-SMID - unique name for each biological samples [optional]
	        - used as output file prefix
		- if not provided, will use FASTQ filename + assay type
	-FCID - flow cell ID for each run [optional]
	-FCLN - lane of flow cell for each library sample [optional]
	-INID - ID of the sequencing instrument [optional]
	-BCID - barcode index for each library sample [optional]
	-CSID - control sample ID [required for variant calling]
		- SMID of control to compare against in mut calling
		- if none provided, no variant calling will be run
		- control sample CSID fields should be left blank
--If not provided, FCID, FCLN, PLID and BCID will be extracted from the
  readname assuming standard illumina format, i.e.:
  @<INID>:<run_number>:<FCID>:<FCLN>:<tile>:<x-pos>:<y-pos> <read>:<filtered>:<control>:<BCID>"

# parse arg
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-si|--sample_info)
		SAMPLE_INFO=$2
		shift
		;;
		-g|--genome)
		GENOME=$2
		shift
		;;
		-b|--basedir)
		BASEDIR=$2
		shift
		;;
		-w|--workdir)
		WORKDIR=$2
		shift
		;;
		-d|--delim)
		DELIM=$2
		shift
		;;
		-sl|--sample_list)
		SAMPLE_LIST=$2
		shift
		;;
		-sc|--seq_centre)
		SEQ_CENTRE=$2
		shift
		;;
		-pl|--platform)
		PLATFORM=$2
		shift
		;;
		-po|--pon)
		PON=$2
		shift
		;;
		-v|--verbose)
		VERBOSITY=$2
		shift
		;;
		-h|--help)
		echo "$USAGE"
		exit 1
		;;
		*)
		echo "Error: Undefined argument provided"
		echo "$USAGE"
		exit 1
		;;
	esac
	shift
done

# check required arg
if [[ -z $SAMPLE_INFO ]]; then
	echo "Error: no sample info provided"
	echo "$USAGE"; exit 1
else
	SAMPLE_INFO=$(realpath $WORKDIR/$SAMPLE_INFO)
fi
if [[ ! -e $SAMPLE_INFO ]]; then
	echo "Sample info file not found, input path: $SAMPLE_INFO"
	echo "$USAGE"; exit 1
fi

# get unique sample names
if [[ -z $SAMPLE_LIST ]]; then
	SAMPLES=$(tail -n +2 $SAMPLE_INFO | cut -d $DELIM -f $SMID_IDX | sort | uniq)
else
	IFS=$DELIM read -r -a SAMPLES <<< "$SAMPLE_LIST"
fi

# set genome resources
if [[ $GENOME = mm10 ]]; then
	DBSNP=$BASEDIR/Resources/Mus_musculus/gatk_mm10bundle/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
	REFFA=$BASEDIR/Resources/Mus_musculus/gatk_mm10bundle/Mus_musculus_mm10
	INDEX=$BASEDIR/Resources/Mus_musculus/gatk_mm10bundle/bwa_index/genome.fa
elif [[ $GENOME = hg38 ]]; then
	DBSNP=$BASEDIR/Resources/Homo_sapiens/gatk_hg38bundle/Homo_sapiens_assembly38.dbsnp.vcf.gz
	COSMIC=$BASEDIR/Resources/Homo_sapiens/gatk_hg38bundle/cosmic_coding_noncoding.vcf.gz
	REFFA=$BASEDIR/Resources/Homo_sapiens/gatk_hg38bundle/Homo_sapiens_assembly38
	INDEX=$BASEDIR/Resources/Homo_sapiens/gatk_hg38bundle/bwa_index/hg38bundle
else
	echo "Genome version not recognised"; echo "$USAGE"; exit 1
fi
REFFA_BASE=$(basename "$REFFA")
DBSNP_BASE=$(basename "$DBSNP")
INDEX_BASE=$(basename "$INDEX")

# check if PON and COSMIC provided
if [[ ! -z $PON ]]; then
	PONARG="-PON $PON"
fi
if [[ ! -z $COSMIC ]]; then
	COSMIC=$(realpath "$COSMIC")
	COSMICBASE=$(basename "$COSMIC")
	COSMICARG="--cosmic $COSMICBASE"
fi

# setup output directories
mkdir -p $WORKDIR/vcf
mkdir -p $WORKDIR/bam
mkdir -p $WORKDIR/logs
mkdir -p $WORKDIR/qc/fastqc
mkdir -p $WORKDIR/qc/metrics

# parse sample info and run
unset VAR_DEPEND
for SMID in ${SAMPLES[@]}; do
	FASTQ1=$(awk -F $DELIM -v smid_idx="$SMID_IDX" -v smid="$SMID" -v field="$FQ1_IDX" \
		'$smid_idx==smid {print $field}' "$SAMPLE_INFO")
	FASTQ2=$(awk -F $DELIM -v smid_idx="$SMID_IDX" -v smid="$SMID" -v field="$FQ2_IDX" \
		'$smid_idx==smid {print $field}' "$SAMPLE_INFO")
	FCIDS=$(awk -F $DELIM -v smid_idx="$SMID_IDX" -v smid="$SMID" -v field="$FCID_IDX" \
		'$smid_idx==smid {print $field}' "$SAMPLE_INFO")
	FCLNS=$(awk -F $DELIM -v smid_idx="$SMID_IDX" -v smid="$SMID" -v field="$FCLN_IDX" \
                '$smid_idx==smid {print $field}' "$SAMPLE_INFO")
	BCIDS=$(awk -F $DELIM -v smid_idx="$SMID_IDX" -v smid="$SMID" -v field="$BCID_IDX" \
                '$smid_idx==smid {print $field}' "$SAMPLE_INFO")
	CSIDS=$(awk -F $DELIM -v smid_idx="$SMID_IDX" -v smid="$SMID" -v field="$CSID_IDX" \
                '$smid_idx==smid {print $field}' "$SAMPLE_INFO")
	
	for IDX in ${!FASTQ1[@]}; do
		FQ1=${FASTQ1[$IDX]}
		FQ2=${FASTQ2[$IDX]}
		FCID=${FCIDS[$IDX]}
		FCLN=${FCLNS[$IDX]}
		BCID=${BCIDS[$IDX]}
		CSID=${CSIDS[$IDX]}

		# check if paired end
		if [[ -z "$FQ2" ]]; then
			PE=0
		else
			PE=1
			PE_FQ2SAM_ARG="FASTQ2=$FQ2"
			PE_BWAMEM_ARG="-p"
		fi
		
		# if readgroup info not in sample info file, extract from read name
		READNAME=$(gzip -dc $WORKDIR/raw_data/$FQ1 | head -n 1)
		if [[ -z $FCID ]]; then
			if [[ $VERBOSITY > 0 ]]; then echo "Extracting FCID from read name"; fi
			FCID=$(echo $READNAME | cut -d ":" -f 3)
		fi
		if [[ -z $FCLN ]]; then
			if [[ $VERBOSITY > 0 ]]; then echo "Extracting FCLN from read name"; fi
			FCLN=$(echo $READNAME | cut -d ":" -f 4)
		fi
		if [[ -z $BCID ]]; then
			if [[ $VERBOSITY > 0 ]]; then echo "Extracting BCID from read name"; fi
			BCID=$(echo $READNAME | cut -d ":" -f 10)
		fi
		RGID=$FCID.$FCLN # read group ID
		PU=$FCID.$FCLN.$BCID # platform unit
		
		# output if verbose
		if [[ $VERBOSITY > 0 ]]; then
			echo "Sample Info parsed for sample $SMID:
				-FASTQ r1 = $FQ1
				-FASTQ r2 = $FQ2
				-Flow cell ID = $FCID
				-Flow cell lane = $FCLN
				-Barcode = $BCID
				-Control condition = $CSID"
		fi

		# reset variables
		unset BAM_LIST
		unset DEPEND

		if [[ ! -e $WORKDIR/bam/$SMID.$RGID.bam ]]; then
			# align and clean each readgroup
			ALIGNJOB=$(cat <<- EOS | qsub -N $PU.aln -
				#!/bin/bash
				#PBS -l select=1:mem=40gb:ncpus=20
				#PBS -l walltime=50:00:00
				#PBS -j oe
				#PBS -q med-bio
				#PBS -o $WORKDIR/logs/$SMID.$PU.alignment.log.txt
				module load picard/2.6.0
				module load java/jdk-8u66
				module load bio-bwa/0.7.10
				module load fastqc/0.11.2
		
				# copy data to scratch
				cp $WORKDIR/raw_data/$FQ1 .
				if [[ $PE = 1 ]]; then cp $WORKDIR/raw_data/$FQ2 .; fi
				cp -rL $INDEX* .
				cp -rL $REFFA* .
		
				# FASTQC
				fastqc --noextract $FQ1 $FQ2
				cp *fastqc.zip $WORKDIR/qc/fastqc/
		
				# convert FASTQ to SAM
				java -Xmx32G -jar /apps/picard/2.6.0/picard.jar FastqToSam \
					FASTQ=$FQ1 \
					OUTPUT=$SMID.$RGID.unaligned.bam \
					READ_GROUP_NAME=$RGID \
					SAMPLE_NAME=$SMID \
					LIBRARY_NAME=$FCID.$BC \
					PLATFORM_UNIT=$PU \
					PLATFORM=$PLATFORM $PE_FQ2SAM_ARG
		
				# mark adapters
				java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkIlluminaAdapters \
					I=$SMID.$RGID.unaligned.bam \
					O=$SMID.$RGID.markadapters.bam \
					M=$SMID.$RGID.markadapters.metrics.txt \
					TMP_DIR=./picard_tmp/
				cp $SMID.$RGID.markadapters.metrics.txt $WORKDIR/qc/metrics/
		
				# convert uBAM to interleaved FASTQ
				java -Xmx32G -jar /apps/picard/2.6.0/picard.jar SamToFastq \
					I=$SMID.$RGID.markadapters.bam \
					FASTQ=$SMID.$RGID.interleaved.fq \
					CLIPPING_ATTRIBUTE=XT \
					CLIPPING_ACTION=2 \
					INTERLEAVE=true \
					NON_PF=true \
					TMP_DIR=picard_tmp
		
				# align
				bwa mem -M -t 20 $PE_BWAMEM_ARG \
					$INDEX_BASE \
					$SMID.$RGID.interleaved.fq > \
					$SMID.$RGID.aligned.sam
		
				# merge uBAM and aligned
				java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MergeBamAlignment \
					R=$REFFA_BASE.fa \
					UNMAPPED_BAM=$SMID.$RGID.markadapters.bam \
					ALIGNED_BAM=$SMID.$RGID.aligned.sam \
					O=$SMID.$RGID.merged.bam \
					CREATE_INDEX=true \
					ADD_MATE_CIGAR=true \
					CLIP_ADAPTERS=false \
					CLIP_OVERLAPPING_READS=true \
					INCLUDE_SECONDARY_ALIGNMENTS=true \
					MAX_INSERTIONS_OR_DELETIONS=-1 \
					PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
					ATTRIBUTES_TO_RETAIN=XS \
					TMP_DIR=picard_tmp
		
 				# mark duplicates
				java -Xmx32G -jar /apps/picard/2.6.0/picard.jar MarkDuplicates \
					I=$SMID.$RGID.merged.bam \
					O=$SMID.$RGID.bam \
					M=$SMID.$RGID.markduplicates.metrics.txt \
					TMP_DIR=./picard_tmp \
					CREATE_INDEX=true
				cp $SMID.$RGID.markduplicates.metrics.txt $WORKDIR/qc/metrics/
		
				# alignment metrics
				java -Xmx32G -jar /apps/picard/2.6.0/picard.jar \
					CollectAlignmentSummaryMetrics \
					R=$REFFA_BASE.fa \
					I=$SMID.$RGID.bam \
					O=$SMID.$RGID.alignmentsummary.metrics.txt
				cp $SMID.$RGID.alignmentsummary.metrics.txt $WORKDIR/qc/metrics/
		
				# copy clean BAM back
				cp $SMID.$RGID.bam* $WORKDIR/bam/
				
				ls -lhAR
				EOS
			)
		fi
		
		# record job ID & bam name for merge job
		DEPEND="$DEPEND,afterok:$ALIGNJOB"
		BAM_LIST="$BAM_LIST $SMID.$RGID.bam"
	done
	
	# remove leading comma/space
	DEPEND=${DEPEND#*,}
	BAM_LIST=${BAM_LIST#* }
	
	# merge all read groups per sample
	MERGEJOB=$(cat <<- EOS | qsub -N $SMID.merge -
		#!/bin/bash
		#PBS -l walltime=20:00:00
		#PBS -l select=1:mem=20gb:ncpus=1
		#PBS -j oe
		#PBS -W depend=$DEPEND
		#PBS -q med-bio
		#PBS -o $WORKDIR/logs/$SMID.mergebam.runinfo.txt
		module load samtools/1.2
		
		cp $WORKDIR/bam/$SMID* .
		
		if [[ ${#FASTQ1[@]} > 1 ]]; then
			samtools merge $SMID.bam $BAM_LIST
		else
			mv $BAM_LIST $SMID.bam
		fi

		samtools index $SMID.bam
		cp $SMID.bam* $WORKDIR/bam/ 
		
		ls -lhAR
		EOS
	)

	# add all merge jobs to dependency list
	VAR_DEPEND="$VAR_DEPEND,afterok:$MERGEJOB"
done

# remove leading comma
VAR_DEPEND=${VAR_DEPEND#*,}

# call variants
for SMID in ${SAMPLES[@]}; do
        
	CSID=$(awk -F $DELIM -v smid_idx="$SMID_IDX" -v smid="$SMID" -v field="$CSID_IDX" \
		'$smid_idx==smid {print $field}' "$SAMPLE_INFO" | uniq)

	if [[ -z $CSID ]]; then
		continue
	fi
	
	# Peak calling
	if [[ -e $WORKDIR/peaks/$SMID.macs2.narrowPeak ]]; then
		echo "$SMID.macs2.narrowPeak already exsists, peak calling not run"
	else
		MACSJOB=$(cat <<- EOS | qsub -N $SMID.MACS - 
			#!/bin/bash
			#PBS -l walltime=10:00:00
			#PBS -l select=1:mem=20gb:ncpus=1
			#PBS -j oe
			#PBS -W depend=afterok:$MERGEJOB
			#PBS -q med-bio
			#PBS -o $WORKDIR/logs/$SMID.peakCalling.runinfo.txt
			module load samtools/1.2
			module load macs/2.1.0
			
			REGEX="^altd <- \\([0-9]\\)"
			FRAGL=$(cat ${SMID}_predictd.R | grep "^altd <- ()" 

	# Variant calling
	if [[ -e $WORKDIR/vcf/$SMID.mutect2.vcf ]]; then
		echo "$SMID.mutect2.vcf already exists, variant calling not run"
	else
		MT2JOB=$(cat <<- EOS | qsub -N $SMID.MT2 -
			#!/bin/bash
			#PBS -l walltime=40:00:00
			#PBS -l select=1:mem=40gb:ncpus=1
			#PBS -j oe
			#PBS -W depend=afterok:$MERGEJOB
			#PBS -q med-bio
			#PBS -o $WORKDIR/logs/$SMID.variantCalling.runinfo.txt
			module load picard/2.6.0
			module load gatk/3.6
			module load java/jdk-8u66
			module load samtools/1.2

			# copy control and test samples and ref accross
			cp $WORKDIR/bam/$SMID.bam* .
			cp $WORKDIR/bam/$CSID.bam* .
			cp $WORKDIR/vcf/$PON* .
			cp -rL $REFFA* .
			cp -rL $DBSNP* .
			cp -rL $COSMIC* .

			# call
			java -jar /apps/gatk/3.6/GenomeAnalysisTK.jar \
				-T MuTect2 \
				-R $REFFA_BASE.fa \
				-I:normal $CSID.bam \
				-I:tumor $SMID.bam \
				--dbsnp $DBSNP_BASE \
				-o $SMID.mutect2.vcf $PONARG $COSMICARG

			ls -lhAR
			cp $SMID.mutect2.vcf* $WORKDIR/vcf/
			EOS
		)
	fi
done
