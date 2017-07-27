#!/bin/bash

# default arguments
SPN=5
FQ1=1
FQ2=2
INDIR=$PWD/raw_data
PHASE=1
DELIM=","

# help message
USAGE="$(basename "$0") [-frn] -s <sample_info.txt> 
-- script to parse sample info file and run WGS analysis pipeline for each sample
-- for each unique sample name, will call wgs_pipeline.sh to map and clean all fastq prior to aggregating

required arguments:
	-s|--samplefile : path to delimited txt file with sample info
optional arguments:
	-p|--phase : which analysis phase to run (default = 1)
			[ 1 = mapping and cleaning bams,
			  2 = merging maps per sample and recalibrating base scores ]
	-i|--indir : directory containing FASTQ (defualt = pwd/raw_data)
	-n|--name : column containing sample names (default = 5)
	-f|--forward : column containing forwards fastq filenames (default = 2)
	-r|--reverse : column containing reverse fastq filenames (default = 3)
	-l|--list : list of sample names (in sample_info.txt) to run pipeline on,
		    if not provided, all sample names are used (default = none)
	-d|--delim : delimiter for txt file (default = ',')"

# parse command line arguments
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-s|--sampleinfo)
		SAMPLE_INFO=$2
		shift
		;;
		-p|--phase)
		PHASE=$2
		shift
		;;
		-i|--indir)
		INDIR=$2
		shift
		;;
		-f|--forward)
		FQ1=$2
		shift
		;;
		-r|--reverse)
		FQ2=$2
		shift
		;;
		-n|--name)
		SPN=$2
		shift
		;;
		-l|--list)
		SAMPLE_LIST=$2
		shift
		;;
		-d|--delim)
		DELIM=$2
		shift
		;;
		*)
		echo "Error: Illegal argument"
		echo "$USAGE"
		exit 1
		;;
	esac
	shift
done

# check input file
if [[ -z "$SAMPLE_INFO" ]] || [[ ! -e "$SAMPLE_INFO" ]]; then
	echo "Could not find sample info file: $SAMPLE_INFO"
	echo "$USAGE"
	exit 1
fi

# get unique sample names
if [[ -z $SAMPLE_LIST ]]; then
	SAMPLES=$(tail -n +2 $SAMPLE_INFO | cut -d $DELIM -f $SPN | sort | uniq)
else
	IFS=$DELIM read -r -a SAMPLES <<< "$SAMPLE_LIST"
fi

# extract fastq names from sample info and run alignment/cleaning jobs
if [[ $PHASE = 1 ]]; then
	for SAMPLE_ID in ${SAMPLES[@]}; do
		echo $SAMPLE_ID
		FASTQ1=$(tail -n +2 $SAMPLE_INFO | cut -d $DELIM -f $FQ1 | grep -E "^${SAMPLE_ID}_")
		FASTQ2=$(tail -n +2 $SAMPLE_INFO | cut -d $DELIM -f $FQ2 | grep -E "^${SAMPLE_ID}_")
		FASTQ1=${FASTQ1//$'\n'/.fastq.gz,"$INDIR"\/}
		FASTQ2=${FASTQ2//$'\n'/.fastq.gz,"$INDIR"\/}
		FASTQ1=$FASTQ1".fastq.gz"
		FASTQ2=$FASTQ2".fastq.gz"
		FASTQ1=${FASTQ1/#/$INDIR\/}
		FASTQ2=${FASTQ2/#/$INDIR\/}
		echo "processing sample $SAMPLE_ID with command:
			bash /work/mrobinso/Scripts/HPC/NF/align_clean.sh -n $SAMPLE_ID -f $FASTQ1 -r $FASTQ2" \
			> logs/$SAMPLE_ID.run_wgs.log.txt
		bash /work/mrobinso/Scripts/HPC/NF/wgs_align_clean.sh \
			-n $SAMPLE_ID -f $FASTQ1 -r $FASTQ2 \
			>> logs/$SAMPLE_ID.run_wgs.log.txt
	done
fi

# merge and recal seperately mapped readgroups from each sample
if [[ $PHASE = 2 ]]; then
	for SAMPLE_ID in ${SAMPLES[@]}; do
		echo "merging and recalibrating $SAMPLE_ID with command:
			bash /work/mrobinso/Scripts/HPC/NF/wgs_merge_recal.sh -s $SAMPLE_ID" \
			>> logs/$SAMPLE_ID.run_wgs.log.txt
		bash /work/mrobinso/Scripts/HPC/NF/wgs_merge_recal.sh -s $SAMPLE_ID \
			>> logs/$SAMPLE_ID.run_wgs.log.txt
	done
fi

# call variants
if [[ $PHASE = 3 ]]; then
	for SAMPLE_ID in ${SAMPLES[@]}; do
		if [[ $SAMPLE_ID =~ "C" ]]; then
			C_LIST=$C_LIST,$SAMPLE_ID.recal.bam
		else
			T_LIST=$T_LIST,$SAMPLE_ID.recal.bam
		fi
	done
	C_LIST=${C_LIST#*,}
	T_LIST=${T_LIST#*,}
	bash /work/mrobinso/Scripts/HPC/NF/wgs_call_variants.sh \
		-t $T_LIST -c $C_LIST \
		-p AK_WGS_mutect2_PON.vcf 
fi
