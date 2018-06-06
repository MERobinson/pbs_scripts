#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
format='png'
outliers='yes'
skipzeros='yes'
plotnumbers='yes'
colormap='RdBu'
check='yes'

# help message
help_message="
Wrapper to run DeepTools plotCorrelation and plotPCA.

usage:
    bash $(basename "$0") [-options] -i <bincov_matrix.npz>
required arguments:
    -i|--input : input bin coverage matrix generated with multiBamSummary 
optional arguments:
    -o|--outdir : out directory for plots (default = '.')
    -l|--logdir : output directory for log files (default = --outdir)
    -n|--name : name prefix for output files (default = --input filename)
    --labels : altenative labels for samples, comma sep (default = NULL)
    --title : plot title (default = NULL)
    --format : plot format [png|eps|pdf|svg] (default = png)
    --matrix : filename to save corr matrix (default = NULL)
    --min : minimum corr value for scale (default = NULL)
    --max : maximum corr value for scale (default = NULL)
    --plotnumbers: whether to plot corr values [yes|no] (default = yes)
    --outliers : whether to remove outlier bins [yes|no] (default = yes)
    --skipzeros : whether to skip bins with zero counts [yes|no] (default = yes)
    --colormap : any name color map from matplotlib (default = RdBu)
    --check : whether to check input files prior to running [yes|no] (default = yes) 
    --depend : pbs dependencies (default = NULL)

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--input)
            input=$2
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
        --labels)
            labels="$2"
            shift
            ;;
        --title)
            title="--plotTitle $2"
            shift
            ;;
        --format)
            format="$2"
            shift
            ;;
        --outliers)
            outliers=$2
            shift
            ;;
        --matrix)
            matrix="--outFileCorMatrix $2"
            shift
            ;;
        --min)
            min="-min $2"
            shift
            ;;
        --max)
            max="-max $2"
            shift
            ;;
        --plotnumbers)
            plotnumers="$2"
            shift
            ;;
        --colormap)
            colormap="--colorMap $2"
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
        *)
            printf "Error: Unrecognised argument: %s %s" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check bam argument is provided and split
if [[ -z "${input:-}" ]]; then
    printf "\nERROR: --input argument required\n"
    echo "$help_message"; exit 1
fi

# check bam files exist and add to args
if [[ $check = 'yes' ]] && [[ ! -r $workdir/$input ]]; then
    printf "\nERROR: input file cannot be read: %s/%s" $workdir $input
    echo "$help_message"; exit 1
fi

# if no name provided extract from first bam
if [[ -z ${name:-} ]]; then
    name=$(basename "$input")
    name=${name%%.*}
fi

# set optional arg
if [[ $outliers = yes ]]; then
    outliers_arg="--removeOutliers"
fi
if [[ -n ${labels:-} ]]; then
    labels=$(echo $labels | sed 's/,/ /g')
    labels="--labels $labels"
fi
if [[ $plotnumbers = 'yes' ]]; then
    plot_num_arg="--plotNumbers"
fi
if [[ $skipzeros = 'yes' ]]; then
    skipzeros_arg="--skipZeros"
fi

# set input basename
input_base=$(basename $input)

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set log files
scr_name=$(basename "$0" .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$name.$scr_name.out.log

# set commands
plot_cor_cmd=("plotCorrelation -in $input_base --whatToPlot heatmap --corMethod spearman"
              "--plotFile $name.corr_hm.$format --plotFileFormat $format"
              "--colorMap $colormap ${outlier_arg:-} ${plot_num_arg:-} ${labels:-}"
              "${title:-} ${min:-} ${max:-}  ${matrix:-} ${skipzeros_arg:-}")
plot_pca_cmd=("plotPCA -in $input_base --plotFile $name.pca.$format"
              "--plotFileFormat $format ${title:-}") 

# run job
script=$(cat <<- EOS
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=12gb:ncpus=1
		#PBS -j oe
		#PBS -N $name.plotqc
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}
	
		# load modules
		module load samtools/1.2
		module load anaconda
		source activate my_deeptools-2.4.2    

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# copy to scratch
		cp $workdir/$input . &>> $out_log

		# plot corr heatmap
		${plot_cor_cmd[@]} &>> $out_log

		# plot pca
		${plot_pca_cmd[@]} &>> $out_log
	
		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log

		cp $name* $workdir/$outdir/ &>> $out_log
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
