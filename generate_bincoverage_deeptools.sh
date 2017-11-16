#!/bin/bash

# default arg
workdir=$PWD
outdir=''
quality=10
dedup='yes'
extend=200
binsize=1000

# help message
help_message="
usage:
    bash $(basename "$0") [-wrtblgneh] -b <BAM>
purpose:
    # Wrapper to run deeptools multiBamSumarry to generate genome coverage data.
required arguments:
    -b|--bam : comma seperated list of bam files
optional arguments:
    -o|--outdir : out directory for .npz file (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output file (default = BAM filename)
    -r|--region : region to analyse (e.g. 'chr1' 'chr1:1000-100000') (default = NULL)
    -e|--extend : length to extend reads (bp) (default = 200)
    -q|--quality : min alignment quality to filter reads (default = 10)
                   (to disable simply set to 0)
    -d|--dedup : whether to filter duplicates [yes,no] (default = yes)
    -s|--binsize : binning for track output (default = 1000)

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam_list=$2
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
        -n|--name)
            name=$2
            shift
            ;;
        -r|--region)
            region=$2
            shift
            ;;
        -e|--extend)
            extend=$2
            shift
            ;;
        -q|--quality)
            quality=$2
            shift
            ;;
        -d|--dedup)
            dedup=$2
            shift
            ;;
        -s|--binsize)
            binsize=$2
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

# check bam argument is provided and split
if [[ -z "${bam_list:-}" ]]; then
    printf "\nERROR: Input BAM file required\n"
    echo "$help_message"; exit 1
else
    IFS=',' read -r -a bam_array <<< $bam_list
fi

# check bam files exist and add to args
for bam in ${bam_array[@]}; do
    if [[ ! -f "$bam" ]]; then
        printf "\nERROR: input BAM file does not exist: %s/%s/%s\n" $bam
        echo "$help_message"; exit 1
    else
        bam_base=$(basename "$bam")
        bam_arg="${bam_arg:-} $bam_base"
        bam_cp="${bam_cp:-}; cp $workdir/$bam* ."
    fi
done
bam_arg=${bam_arg#* }
bam_cp=${bam_cp#*; }

# if no name provided extract from first bam
if [[ -z "$name" ]]; then
    name=$(basename "${bam_array[1]}")
    name=${name%%.*}
fi

# set optional arg
if [[ $dedup = yes ]]; then
    dedup_arg="--samFlagExclude 1024"
fi
if [[ ! -z ${binsize:-} ]]; then
    bs_arg="--binSize $binsize"
fi
if [[ ! -z ${region:-} ]]; then
    region_arg="--region $region"
fi
if [[ ! -z ${extend:-} ]]; then
    extend_arg="--extendReads $extend"
fi
if [[ ! -z ${quality:-} ]]; then
    qual_arg="--minMappingQuality $quality"
fi

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# run job
jobid=$(cat <<- EOS | qsub -N $name.multibam - 
	#!/bin/bash
	#PBS -l walltime=05:00:00
	#PBS -l select=1:mem=10gb:ncpus=20
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $workdir/$logdir/$name.generate_bincoverage_deeptools.log
	
	# load modules
	module load samtools/1.2
	module load anaconda/2.4.1
	module load deeptools/2.4.2
	source activate deeptools-2.4.2    

	# copy to scratch
	$bam_cp

	# generate coverage
	multiBamSummary bins \
		--bamfiles $bam_arg \
		--outFileName $name.npz \
		--numberOfProcessors 20 \
		${region_arg:-} \
		${extend_arg:-} \
        ${qual_arg:-} \
        ${dedup_arg:-} \
        ${bs_arg:-} \
	
	cp $name.npz $workdir/$outdir
	ls -lhAR
	EOS
)
