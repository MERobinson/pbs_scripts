#!/bin/bash

# default arg
workdir=$PWD
macsdir=macs
outdir=macs
logdir=logs

# help message
help_message="

usage:
    bash $(basename "$0") [-wbnh] -c1 <cond1> -c2 <cond2>
purpose:
    Simple wrapper to call differential peaks with macs2 bdgdiff
required arguments:
    -c1|--cond1 : prefix of macs2 call peak output files for condition 1
    -c2|--cond2 : prefix of macs2 call peak output files for condition 2
optional arguments:
    -w|--workdir : path of working directory to use (default = pwd)
    -m|--macsdir : directory in which to find macs2 callpeak peak .xls files (default = peaks)
    -o|--outdir : output directory for results (default = diff_peaks)
    -l|--logdir : output directory for run logs (default = logs)
    -d1|--depth1 : minimum filtered tag count for cond1 (extracted from .xls)
    -d2|--depth2 : minimum filtered tag count for cond2 (extracted from .xls)
    -n|--name : name prefix for output files (deafult = <cond1>_vs_<cond2>)
    -d|--depend : list of dependencies to pass to PBS job (default = none)
    -v|--validate : whether to check files [yes,no] - useful if scheduling jobs (default = yes)
example:
    bash $(basename "$0") -t K562_H3K27ac -c K562_input
additional info:
    > macs2 callpeak needs to be run prior to this script
    > --macsdir should contain macs2 callpeak output - pileup & lambda bdg and xls
    > cond1 and cond2 prefixes should match names of macs2 callpeak output files
    > will attempt to extract read counts from macs2 callpeak .xls files by default

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -c1|--cond1)
        cond1=$2
        shift
        ;;
        -c2|--cond2)
        cond2=$2
        shift
        ;;
        -w|--workdir)
        workdir=$2
        shift
        ;;
        -m|--macsdir)
        macsdir=$2
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
        -d1|--depth1)
        depth1=$2
        shift
        ;;
        -d2|--depth2)
        depth2=$2
        shift
        ;;
        -n|--name)
        name=$2
        shift
        ;;
        *)
        printf "\nERROR: Illegal argument\n"
        echo "$help_message"
        exit 1
        ;;
    esac
    shift
done

# check required arguments
if [[ -z "$cond1" ]] || [[ -z "$cond2" ]]; then
    printf "\nERROR: Requires both test and control prefixes to be specified\n"
    echo "$help_message"; exit 1
fi

# function to check files
function check_file {
    local filename=$1
    if [[ ! -f $filename ]]; then
        printf "\nERROR: input file not found: %s\n" $filename
        echo $help_message; exit 2
    fi
}

# check input files exist
if [[ $validate = yes ]]; then
    check_file "$workdir/$macsdir/${cond1}_treat_pileup.bdg"
    check_file "$workdir/$macsdir/${cond1}_control_lamda.bdg"
    check_file "$workdir/$macsdir/${cond2}_treat_pileup.bdg"
    check_file "$workdir/$macsdir/${cond1}_control_lamda.bdg"
fi

# get name if req
if [[ -z "$name" ]]; then name=${cond1}_vs_${cond2}; fi

# if no depth provided, extract from peaks file
pattern="\"^# tags after filtering in (treatment|control)\""
if [[ -z "$depth1" ]]; then
	depth_arg_1="depths=(\$(grep -E "$pattern" "${cond1}_peaks.xls" | grep -Eo "[0-9]+"))"
	depth_arg_1="$depth_arg_1; if [[ \${depths[0]} -gt \${depths[1]} ]]"
    depth_arg_1="$depth_arg_1; then depth1=\${depths[1]}; else depth1=\${depths[0]}; fi"
else
    depth_arg="depth1=$depth1"
fi
if [[ -z "$depth2" ]]; then
    depth_arg_2="depths=(\$(grep -E "$pattern" "${cond2}_peaks.xls" | grep -Eo "[0-9]+"))"
    depth_arg_2="$depth_arg_2; if [[ \${depths[0]} -gt \${depths[1]} ]]" 
    depth_arg_2="$depth_arg_2; then depth2=\${depths[1]}; else depth2=\${depths[0]}; fi"
else
    depth_arg="depth2=$depth2"
fi

# handle dependencies
if [[ ! -z $depend ]]; then depend="#PBS -W depend=$depend"; fi

# create output dirs
mkdir -p $workdir/$logdir
mkdir -p $workdir/$outdir

# run job
jobid=$(cat <<- EOS | qsub -N $name.macsd -
	#!/bin/bash
	#PBS -l walltime=10:00:00
	#PBS -l select=1:mem=10gb:ncpus=1
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $workdir/$logdir/$name.call_diffpeaks_macs.log
	$depend

	module load samtools/1.2
	module load macs/2.1.0

	cp $workdir/$macsdir/$cond1* .
	cp $workdir/$macsdir/$cond2* .

	$depth_arg_1
	$depth_arg_2

    echo \$depth1; echo \$depth2

	macs2 bdgdiff \
		--t1 ${cond1}_treat_pileup.bdg \
		--c1 ${cond1}_control_lambda.bdg \
		--t2 ${cond2}_treat_pileup.bdg \
		--c2 ${cond2}_control_lambda.bdg \
		--d1 \$depth1 \
		--d2 \$depth2 \
		--o-prefix $name

	cp $name* $workdir/$outdir/

	ls -lhAR
	EOS
)
