#!/bin/bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
genome='hg38'
check='y'
force='n'

# help message
help_message="
Wrapper to call peaks with epic2 

Usage:
    bash $(basename "$0") [-options] -t <BAM>
required arguments:
    -t|--treatment : treatment sample(s) [BAM|SAM|BED]
optional arguments:
    -c|--control : control sample(s) [BAM|SAM|BED]
    -o|--outdir : output dir for macs files (default = '.')
    -l|--logdir : output dir for log files (default = --outdir)
    -q|--qcdir : output dir for qc files (default = ---outdir)
    -g|--genome : genome version [hg38,hg19,etc] (default = hg38)
    -n|--name : prefix for output files (default = treatment sample name)
    -f|--fragsize : reads are extended by half fragsize (default = 150)
    --fdr : FDR cutoff [float] (default = 0.05)
    --gaps : gaps allowed [integer] (default = 3)
    --force : whether to run even if output already exists [y|n] (default = n)
    --check : whether to check files prior to running [y|n] (default = y)
    --depend : list of dependencies for PBS script (e.g. afterok:012345,afterok:012346)
additional info:
    # all paths should be given relative to working directory
    # --outdir will be created, if not pre-existing
    # if --control is provided, will be used as background for peak calling,
      otherwise no control sample is used in macs
    # if not provided, name is extracted from treatment filename up to first period
    # treatment/control can be given as single inputs or matched comma separated lists

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -t|--treatment)
            treatment=$2
            shift
            ;;
        -c|--control)
            control=$2
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
        -q|--qcdir)
            qcdir=$2
            shift
            ;;
        -g|--genome)
            genome=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        -f|--fragsize)
            fragsize=$2
		    shift
		    ;;
        --fdr)
            fdr=$2
            shift
            ;;
        --gaps)
            gaps=$2
            shift
            ;;
        --force)
            force=$2
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
            printf "\nERROR: Unrecognised argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
	shift
done

# set logdir and qcdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
if [[ -z ${qcdir:-} ]]; then
    qcdir=$outdir
fi

# check required arguments
if [[ -z "${treatment:-}" ]]; then
    printf "\nERROR: No treatment BAM provided\n"
    echo "$help_message"; exit 1
fi

# split input lists to arrays
IFS="," read -r -a t_array <<< $treatment
if [[ -n ${control:-} ]]; then
    IFS="," read -r -a c_array <<< $control
else
    echo "WARNING: no control samples input"
fi

# set name if not provided
if [[ -z ${name:-} ]]; then
    name_array=(${t_array[@]/%.*/})
    name_array=(${name_array[@]/#*\//})
else
	IFS="," read -r -a name_array <<< $name
fi

# parse optional args
if [[ -n ${fragsize:-} ]]; then
    fragsize="--fragment-size $fragsize"
fi
if [[ -n ${fdr:-} ]]; then
    fdr="-fdr $fdr"
fi
if [[ -n ${gaps:-} ]]; then
    gaps="--gaps-allowed $gaps"
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir

for idx in ${!t_array[@]}; do

	treatment=${t_array[$idx]}
	control=${c_array[$idx]:-}
	name=${name_array[$idx]}

    # check input
    if [[ $check = 'y' ]] && [[ ! -r $workdir/$treatment ]]; then
        printf "WARNING: treatment BAM cannot be read - skipping: %s/%s\n" $workdir $treatment
    elif [[ -r $outdir/${name}_peaks.xls ]] && [[ $force != 'y' ]]; then
        printf "WARNING: peaks already exist - skipping: %s/%s_peaks.xls\n" $workdir $name
    else
		# set control arg
		if [[ -n ${control:-} ]] && [[ ${control:-} != 'none' ]]; then
			if [[ $check = 'y' ]] && [[ ! -r $workdir/$control ]]; then
				printf "\nERROR: control BAM cannot be read: %s/%s\n" $workdir $control
				echo "$help_message"; exit 1
			fi
			cntrl_cmd="--control $(basename ${control})"
			cntrl_cp="cp ${workdir}/${control}* ."
        else
            echo "WARNING: no control sample input."
        fi

        # set commands
        epic_cmd=("epic2 --treatment $(basename ${treatment})"
                  "--genome $genome"
                  "${cntrl_cmd:-}"
                  "${fragsize:-}"
                  "${gaps:-}"
                  "${fdr:-}")

		# set log file names
		std_log=$workdir/$logdir/$name.$(basename "$0" .sh).std.log
		pbs_log=$workdir/$logdir/$name.$(basename "$0" .sh).pbs.log
		out_log=$workdir/$logdir/$name.$(basename "$0" .sh).out.log

		# run job
		script=$(cat <<- EOS
				#!/bin/bash
				#PBS -l walltime=24:00:00
				#PBS -l select=1:mem=10gb:ncpus=4
				#PBS -j oe
				#PBS -N epic.$name
				#PBS -q med-bio
				#PBS -o $std_log
				${depend:-}

				printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

				# load modules
				module load samtools/1.2 &>> ${out_log}
				module load anaconda3/personal &>> ${out_log}
				source activate epic2 &>> ${out_log}

				# copy to scratch
				cp ${workdir}/${treatment}* . &>> ${out_log}
				${cntrl_cp:-} 

				# call peaks
				printf "\nCalling peaks:\n" >> $out_log
				${epic_cmd[@]} > ${name}.epic2_peaks.bed 2> $out_log
				cp ${name}.epic2_peaks.bed $workdir/$outdir/

				printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
				ls -lhAR &>> $out_log
			EOS
		)
		echo "$script" > $pbs_log

		# submit job
		jobid=$(qsub "$pbs_log")

		# echo job id and exit
		echo "JOBID: $jobid"
	fi
done
exit 0
