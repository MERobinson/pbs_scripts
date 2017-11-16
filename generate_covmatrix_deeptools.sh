#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''

# help message
help_message="
usage:
    bash $(basename "$0") [-options] -b <bigWig,files> -r <regions>
purpose:
    # Wrapper to run deeptools computeMatrix to generate a coverage matrix.
required arguments:
    -b|--bigwig : comma separated list of bigwig files [bigWig]
    -r|--regions : regions to copute coverage for [BED,GTF]
    -m|--mode : runmode - [scale-regions,reference-point]
optional arguments:
    -o|--outdir : output directory for track files (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output file (default = BAM filename)
    --upstream : bp upstream to include [integer] (default = 0)
    --downstream : bp downstream to include [integer] (default = 0)
    --skip_zeros : whether to skip zero score regions [yes,no] (default = no)
    --ref_point : feature site for reference-point [TSS,TES,center] (default = TSS) 
additional info:
    # all paths should be relative to the current working directory
    # coverage can be calculated for a fixed region (reference-point) or scaled accross the
      entire region (scale-regions) mode

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bigwig)
            bw_list=$2
            shift
            ;;
        -r|--regions)
            region_list=$2
            shift
            ;;
        -m|--mode)
            mode=$2
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
        --upstream)
            upstream="-b $2"
            shift
            ;;
        --downstream)
            downstream="-a $2"
            shift
            ;;
        --skip_zeros)
            skip_zeros="--skipZeros"
            shift
            ;;
        --ref_point)
            ref_point="--referencePoint $2"
            shift
            ;;
        *)
            echo "\nError: Illegal argument: %s %s" $1 $2
            echo "$help_message"; exit 1
        ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arguments
if [[ -z ${bw_list:-} ]]; then
    printf "\nERROR: no BAM file provided\n"
    echo "$help_message"; exit 1
elif [[ -z ${region_list:-} ]]; then
    printf "\nERROR: no FASTA file provided\n"
    echo "$help_message"; exit 1
fi

# split lists
IFS="," read -r -a bw_array <<< "$bw_list"
IFS="," read -r -a rg_array <<< "$region_list"

# check files and add to arg
bw_arg="-S"; rg_arg="-R"
for bw in ${bw_array[@]}; do
    if [[ ! -r $workdir/$bw ]]; then
        printf "\nERROR: bigWig file not readable: %s\n" $bw
        echo "$help_message"; exit 1
    else
        bw_cp="${bw_cp:-}; cp $workdir/$bw* ."
        bw_arg="${bw_arg:-} $(basename $bw)"
    fi
done
for rg in ${rg_array[@]}; do
    if [[ ! -r $workdir/$rg ]]; then
        printf "\nERROR: BED file not readable: %s\n" $rg
        echo "$help_message"; exit 1
    else
        rg_cp="${rg_cp:-}; cp $workdir/$rg* ."
        rg_arg="${rg_arg:-} $(basename $rg)"
    fi
done

# strip leading sperators
bw_cp=${bw_cp#*; }; rg_cp=${rg_cp#*; }

# if no name provided extract from BAM
if [[ -z "$name" ]]; then
    name=${bam%%.*}
fi

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# run job
jobid=$(cat <<- EOS | qsub -N $name.covmat - 
		#!/bin/bash
		#PBS -l walltime=02:00:00
		#PBS -l select=1:mem=20gb:ncpus=20
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.generate_covmatrix_deeptools.log
	
		# load modules
		module load samtools/1.2
		module load anaconda/2.4.1
		module load deeptools/2.4.2
		source activate deeptools-2.4.2    

		# copy resource files to scratch
		${bw_cp:-}
		${rg_cp:-}

		# generate coverage track
		computeMatrix ${mode} \
			${bw_arg} \
			${rg_arg} \
			${upstream:-} \
			${downstream:-} \
			${ref_point:-} \
			${skip_zeros:-} \
			--numberOfProcessors 20 \
			-o $name.cov_matrix.gz \
			--outFileSortedRegions $name.sorted_regions.bed

		cp $name.cov_matrix.gz $workdir/$outdir/
		cp $name.sorted_regions.bed $workdir/$outdir/
		ls -lhAR
	EOS
)
