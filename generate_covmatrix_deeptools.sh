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
    -r|--regions : regions to compute coverage for [BED,GTF]
    -m|--mode : runmode - [scale-regions,reference-point]
optional arguments:
    -o|--outdir : output directory for track files (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output file (default = BAM filename)
    --upstream : bp upstream to include [integer] (default = 0)
    --downstream : bp downstream to include [integer] (default = 0)
    --skip_zeros : whether to skip zero score regions [yes,no] (default = no)
    --ref_point : feature site for reference-point [TSS,TES,center] (default = TSS) 
    --sample_labels : comma sep list of labels matching samples (default = NULL)
    --depend : PBS dependencies to pass to job (default = NULL)
additional info:
    # all paths should be relative to the current working directory
    # coverage can be calculated for a fixed region (reference-point) or scaled accross the
      entire region (scale-regions) mode
    # multiple region files can be provided as comma separated list for grouped coverage
    # dependencies should be in form 'afterok:123456,afterok:123457'

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
        --sample_labels)
            sample_labels="$2"
            shift
            ;;
        --depend)
            depend="#PBS -W depend=$2"
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
    printf "\nERROR: no regions provided\n"
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

# if no name provided extract from BAM
if [[ -z ${name:-} ]]; then
    name=$(basename ${bw_array[0]})
    name=${name%%.*}.cov
fi

# sort labels if provided
if [[ -n ${sample_labels:-} ]]; then
    IFS="," read -r -a lab_array <<< "$sample_labels"
    lab_arg="--samplesLabel ${lab_array[@]}"
else
    lab_arg="--smartLabels"
fi

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set log filenames
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# set commands
bw_cp="${bw_cp#*; } &>> $out_log"
rg_cp="${rg_cp#*; } &>> $out_log"
cov_cmd=("computeMatrix ${mode} --numberOfProcessors 20 ${bw_arg} ${rg_arg}"
         "-o $name.cov_matrix.gz --outFileSortedRegions $name.sorted_regions.bed"
         "${upstream:-} ${downstream:-} ${ref_point:-} ${skip_zeros:-} ${lab_arg:-}")

# write job script
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=20gb:ncpus=20
		#PBS -j oe
		#PBS -N $name.covmat
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}
	
		# load modules
		module load samtools/1.2
		module load anaconda3/personal
		source activate deeptools-3.1.0

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy resource files to scratch
		${bw_cp:-}
		${rg_cp:-}

		# generate coverage track
		${cov_cmd[@]} &>> $out_log

		cp $name.cov_matrix.gz $workdir/$outdir/ &>> $out_log
		cp $name.sorted_regions.bed $workdir/$outdir/ &>> $out_log

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR >> $out_log
		
		cp $out_log $workdir/$logdir/
		ls -lhAR
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
