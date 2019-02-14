#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
sep='_'
outdir=''
quality=10
binsize=200
check='yes'
scriptdir='scripts'
scr_name=$(basename "$0")

# help message
help_message="
Run EpicSeg to segment chromatin by histone marks.

usage:
    bash $scr_name [-options] -m <list,of,BAM>
required arguments:
    -B|--bed : bed file of chromsome regions to analyse
    -ns|--n_states : number of states to segment into
    # one of either --counts or --bam_list
    -b|--bam_list : comma seperated list of aligned BAM (histone marks)
    -c|--counts_file : pre-computed counts file from epicseg getcounts
    # one of either --marks or --field and --sep args
    -m|--marks : comma sep list of which mark each bam specifies [string]
    -f|--field : field from filename to extract mark for each bam [interger]
# optional arguments:
    -s|--sep : field seperator in filename for extracting mark [string] (default = '_')
    -o|--outdir : output directory (default = PWD)
    -ld|--logdir : output directory for log files (default = --outdir)
    -sd|--scriptdir: directory containing epicseg.R script (default = 'scripts')
    -n|--name : name prefix for output file (default = BAM filename)
    -g|--genome : genome version - used in trackline if provided (default = NULL)
    -q|--quality : min alignment quality to filter reads (default = 10)
                   (to disable simply set to 0)
    -bs|--binsize : binning for segmentation (default = 200)
    --check : whether to check input files [yes,no] (default = yes)
    --depend : list of PBS dependencies (default = NULL) 
additional info:
    # all paths should be relative to the current working directory
    # each input bam requires matching mark name - e.g. H3K4me3,
      these can be extracted from filename with --sep and --field args,
      or manually provided with --marks arg. 

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -ns|--n_states)
            n_states=$2
            shift
            ;;
        -b|--bam_list)
            bam_list=$2
            shift
            ;;
        -c|--counts_file)
            counts_file=$2
            shift
            ;;
        -B|--bed)
            bed=$2
            shift
            ;;
        -m|--marks)
            marks=$2
            shift
            ;;
        -f|--field)
            field=$2
            shift
            ;;
        -s|--sep)
            sep=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -ld|--logdir)
            logdir=$2
            shift
            ;;
        -sd|--scriptdir)
            scriptdir=$2
            shift
            ;;
        -g|--genome)
            db=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -q|--quality)
            quality=$2
            shift
            ;;
        -bs|--binsize)
            binsize=$2
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
            printf "\nError: Illegal argument: %s %s" $1 $2
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
if [[ -z ${n_states:-} ]]; then
    printf "\nERROR: --n_states must be specified\n"
    echo "$help_message"; exit 1
fi
if [[ -z ${bam_list:-} ]] & [[ -z ${counts_file:-} ]]; then
    printf "\nERROR: Either --counts_file or --bam_list must be specified\n"
    echo "$help_message"; exit 1
elif [[ -n ${bam_list:-} ]]; then
    IFS="," read -r -a bam_array <<< $bam_list
fi
if [[ -z ${bed:-} ]]; then
    printf "\nERROR: no BED file provided\n"
    echo "$help_message"; exit 1
else
    if [[ $check = "yes" ]] & [[ ! -r $bed ]]; then
        printf "\nERROR: BED file isnt readable: %s/%s\n" $workdir $bed
        echo "$help_message"; exit 1
    else
        bed_base=$(basename $bed)
    fi
fi

# if no marks provided, extract from filenames
if [[ -n ${counts_file:-} ]]; then
    :
elif [[ -n "${marks:-}" ]]; then
    IFS="," read -r -a marks_array <<< $marks
elif [[ -z ${field:-} ]]; then
    printf "\nERROR: either --marks or --field arguments need specifying\n"
    echo "$help_message"; exit 1
else
    marks_array=($(for f in ${bam_array[@]}; do echo $(basename "${f%.*}") | \
					cut -d "$sep" -f "$field"; done))
fi

# set marks args
for idx in "${!bam_array[@]}"; do
    bam_base=$(basename ${bam_array[idx]})
    if [[ $idx = 0 ]]; then 
        bam_cp="cp ${bam_array[idx]}* ."
        marks_arg="-m ${marks_array[idx]}:$bam_base"
    else
        bam_cp="${bam_cp}; cp ${bam_array[idx]}* ."
        marks_arg="${marks_arg} -m ${marks_array[idx]}:$bam_base"
    fi
done

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set log file names
scr_name=${scr_name%.*}
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log
out_log=$workdir/$logdir/$name.$scr_name.out.log

# compile commands
if [[ -z ${counts_file:-} ]]; then
	counts_command=("Rscript $workdir/$scriptdir/epicseg.R getcounts $marks_arg -r $bed_base"
                	"--target $name.counts.rda --nthreads 8 --mapq 10")
	counts_base="$name.counts.rda"
elif [[ $check = "yes" ]] & [[ ! -r $counts_file ]]; then
	printf "\nERROR: counts file is not readable: %s/%s\n" $workdir $counts_file
	echo "$help_message"; exit 1
else
	counts_command=("cp $workdir/$counts_file .")
	counts_base=$(basename "$counts_file")
fi
bed_command=("if [[ -e ${name}_refined_regions.bed ]]; then"
			"bed_file=${name}_refined_regions.bed; else"
			"bed_file=$bed_base; fi")
seg_command=("Rscript $workdir/$scriptdir/epicseg.R segment --counts $counts_base --nthreads 8"
                "--regions \$bed_file --nstates $n_states --prefix $name.seg$n_states.")

# run job
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=10gb:ncpus=8
		#PBS -j oe
		#PBS -N $name.seg
		#PBS -q med-bio
		#PBS -o $std_log
		${depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# load modules
		module load samtools/1.2
		module load anaconda3/personal    

		# copy resource files to scratch
		${bam_cp:-}
		cp ${workdir}/${bed} .

		# get counts
		${counts_command[@]} &>> $out_log

		# check if reformatted 
		${bed_command[@]} &>> $out_log

		# segment
		${seg_command[@]} &>> $out_log

		cp ${name}* $workdir/$outdir &>> $out_log

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
		ls -lhAR >> $out_log
	EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo jobid and exit
echo "JOBID: $jobid"
exit 0
