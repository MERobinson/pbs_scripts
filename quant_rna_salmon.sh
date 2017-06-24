#!/bin/bash 

# default arguments
WORKDIR=$PWD
OUTDIR=salmon_out
INDIR=fastq
INDEX=Homo_sapiens.GRCh38.cds.all_index
LIBTYPE=ISR

# help message
USAGE="

usage: bash $(basename "$0") [-bioxlh] -s <SAMPLE_NAME>
purpose:
	Simple wrapper to run salmon on RNAseq data
optional arguments:
	-b|--basedir : the base directory to find input and output dir (default = pwd)
	-i|--indir : name of directory (in basedir) containing FASTQ for analysis (default = fastq)
	-o|--outdir : output directory name (will be created in basedir) (default = salmon_out)
	-x|--index : index filename (should be located in basedir) (default = Homo_sapiens.GRCh38.cds.all_index)
	-l|--libtype : library types string to pass to salmon (default = ISR)
	-h|--help : display this help message
example:
	bash $(basename "$0") --basedir $WORK/experiment1/
expectations:
	-expects FASTQ filenames in standard illumina format, i.e:
		<SAMPLE_NAME>_<LANE>_<READ1/2>_001.fastq.gz

"

# parse command line arguments
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-b|--basedir)
		WORKDIR=$2
		shift
		;;
		-i|--indir)
		INDIR=$2
		shift
		;;
		-o|--outdir)
		OUTDIR=$2
		shift
		;;
		-x|--index)
		INDEX=$2
		shift
		;;
		-l|--libtype)
		LIBTYPE=$2
		shift
		;;
		-h|--help)
		echo "$USAGE"; exit 1
		;;
		*)
		echo "Error: Illegal command line argument"
		echo "$USAGE"
		exit 1
		;;
	esac
	shift
done

# create outdir if needed
mkdir -p $WORKDIR/$OUTDIR/logs

# check args
if [[ ! -e $WORKDIR/$INDIR ]]; then
	echo "ERROR: Input directory not found"; echo "$USAGE"; exit 1
elif [[ ! -e $WORKDIR/$INDEX ]]; then
	echo "ERROR: Index not found"; echo "$USAGE"; exit
fi

# list R1 FASTQ files
FASTQ_LIST=($(ls $WORKDIR/$INDIR/*_R1_*.fastq*))

# check files
if [[ -z "$FASTQ_LIST" ]]; then
	echo "ERROR: No FASTQ files identified in $WORKDIR/$INDIR"
	echo "$USAGE"; exit 1
fi

# loop through FASTQ and run salmon script
for FASTQ in ${FASTQ_LIST[@]}; do

	# extract sample names
	FASTQ=$(basename "$FASTQ")
	FILE_NAME=${FASTQ%%_R1_*}
	SAMPLE_NAME=${FASTQ%%_L00*}

	# run salmon
	JOBID=$(cat <<- EOS | qsub -N $SAMPLE_NAME.salmon -
		#!/bin/bash
		#PBS -j oe
		#PBS -l walltime=10:00:00
		#PBS -l select=1:ncpus=8:mem=20gb
		#PBS -o $WORKDIR/$OUTDIR/logs/${SAMPLE_NAME}.salmonQuant.log.txt

		# load modules
		module load salmon/0.8.2

		# copy FASTQ and INDEX files to scratch directory
		cp $WORKDIR/$INDIR/*.fastq* .
		cp -rL $WORKDIR/$INDEX* .

		# run Salmon
		salmon quant \
			-i $INDEX \
			-l $LIBTYPE \
			-p 8 \
			--seqBias \
			--gcBias \
			-1 ${FILE_NAME}_R1_001.fastq.gz \
			-2 ${FILE_NAME}_R2_001.fastq.gz \
			-o $SAMPLE_NAME
	
		# copy results to output directory
		cp -r  $SAMPLE_NAME/ $WORKDIR/$OUTDIR/
		EOS
	)
done
