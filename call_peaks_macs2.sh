#!/bin/bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
genome='hs'
check='on'
force='off'
bdg='off'
broad='off'
nomodel='off'
summits='off'
format='BAM'

# help message
help_message="
Wrapper to call peaks with MACS2 

Usage:
    bash $(basename "$0") [-options] -t <BAM>
required arguments:
    -t|--test : test sample(s) [BAM]
    -F|--fasta : genome seq [FASTA]
optional arguments:
    -c|--control : control sample(s) [BAM]
    -o|--outdir : output dir for macs files (default = '.')
    -l|--logdir : output dir for log files (default = --outdir)
    -q|--qcdir : output dir for qc files (default = ---outdir)
    -g|--genome : genome size [hs,mm,<numeric>] (default = hs)
    -n|--name : prefix for output files (default = test sample name)
    -e|--extsize : run macs2 with provided extension (default = unset)
    -f|--format : format of test file [BAM|BED|...] (default = 'BAM')
    --shift : shift size [INT] (default = NULL)
    --broad : whether to run in broad peak mode [on|off] (default = off)
    --bdg : whether to generate bedgraph tracks [on|off] (default = off)
    --nomodel : whether to run in nomodel mode [on|off] (default = off)
    --pval : pvalue cutoff [INT] (default NULL)
    --summits : whether to call summits [on|off] (default = off)
    --force : whether to re-run if output already exist [on|off] (default = off)
    --check : whether to check files prior to running [on|off] (default = on)
    --depend : list of dependencies for PBS script (e.g. afterok:012345,afterok:012346)
additional info:
    # all paths should be given relative to working directory
    # --outdir will be created, if not pre-existing
    # if --control is provided, will be used as background for peak calling,
      otherwise no control sample is used in macs
    # if not provided, name is extracted from test filename up to first period
    # test/control can be given as single inputs or matched comma separated lists

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -t|--test)
            test=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
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
        -e|--extsize)
            extsize=$2
		    shift
		    ;;
        -f|--format)
            format="$2"
            shift
            ;;
        --nomodel)
            nomodel=$2
            shift
            ;;
        --shift)
            shift_arg="--shift $2"
            shift
            ;;
        --pval)
            pval_arg="--pvalue $2"
            shift
            ;;
        --summits)
            summits=$2
            shift
            ;;
        --broad)
            broad=$2
            shift
            ;;
        --bdg)
            bdg=$2
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
if [[ -z "${test:-}" ]]; then
    printf "\nERROR: No test BAM provided\n"
    echo "$help_message"; exit 1
elif [[ -z "${fasta:-}" ]]; then
    printf "\nERROR: No FASTA file provided\n"
    echo "$help_message"; exit 1
elif [[ ${check} = 'on' ]] && [[ ! -r $workdir/$fasta ]]; then
    printf "\nERROR: FASTA cannot be read: %s/%s\n" $workdir $fasta
    echo "$help_message"; exit 1
fi

# check optional arg
if [[ ${bdg} = 'on' ]]; then
    bdg_arg="--bdg --SPMR"
fi
if [[ ${summits} = 'on' ]]; then
    summits_arg="--call-summits"
fi
if [[ ${broad} = 'on' ]]; then
    broad_arg="--broad"
fi
if [[ ${nomodel} = 'on' ]]; then
    model_arg="--nomodel"
fi
if [[ -n ${extsize:-} ]]; then
    ext_arg="--extsize ${extsize}"
    model_arg="--nomodel"
fi

# split input lists to arrays
IFS="," read -r -a test_array <<< $test
if [[ -n ${control:-} ]]; then
    IFS="," read -r -a cntrl_array <<< $control
else
    echo "WARNING: no control samples input"
fi

# set name if not provided
if [[ -z ${name:-} ]]; then
    name_array=(${test_array[@]/%.*/})
    name_array=(${name_array[@]/#*\//})
else
	IFS="," read -r -a name_array <<< $name
fi

# get canonical chr from fasta
chr_list=$(cat $fasta | grep -Eo "^>chr[0-9MXY]+\b" | \
           grep -Eo "chr[0-9XY]+"| tr "\n" " ")

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir

for idx in ${!test_array[@]}; do

	test=${test_array[$idx]}
	control=${cntrl_array[$idx]:-}
	name=${name_array[$idx]}

    # check input
    if [[ $check = 'on' ]] && [[ ! -r $workdir/$test ]]; then
        printf "WARNING: test BAM cannot be read - skipping: %s/%s\n" $workdir $test
    elif [[ -r $outdir/${name}_peaks.xls ]] && [[ $force != 'on' ]]; then
        printf "WARNING: peaks already exist - skipping: %s/%s/%s_peaks.xls\n" $workdir $outdir $name
    else
		# set control arg
		if [[ -n ${control:-} ]] && [[ ${control:-} != 'none' ]]; then
			if [[ $check = 'on' ]] && [[ ! -r $workdir/$control ]]; then
				printf "\nERROR: control BAM cannot be read: %s/%s\n" $workdir $control
				echo "$help_message"; exit 1
			fi
			cntrl_arg="-c $(basename ${control})"
			cntrl_cp="cp ${workdir}/${control}* ."
			cntrl_filt=("samtools view -b -F 1024 -q 10"
						"$(basename ${control}) ${chr_list} > tmp.bam;"
						"mv tmp.bam $(basename ${control});"
						"samtools index $(basename ${control})")
		else
            cntrl_arg="--nolambda"
            echo "WARNING: no control - running with --nolambda."
        fi

        # set commands
        predictd_cmd=("macs2 predictd"
                      "-i $(basename ${test})"
                      "-g ${genome}"
                      "--outdir ."
                      "--format $format"
                      "--rfile ${name}.predictd.R")
        callpeaks_cmd=("macs2 callpeak"
                       "-t $(basename ${test})"
                       "-g $genome"
                       "--outdir ."
                       "-n $name"
                       "--keep-dup all"
                       "--format ${format}"
                       "${bdg_arg:-}"
                       "${cntrl_arg:-}"
                       "${ext_arg:-}"
                       "${broad_arg:-}"
                       "${pval_arg:-}"
                       "${model_arg:-}"
                       "${shift_arg:-}"
                       "${summits_arg:-}")

		# set log file names
		std_log=$workdir/$logdir/$name.$(basename "$0" .sh).std.log
		pbs_log=$workdir/$logdir/$name.$(basename "$0" .sh).pbs.log
		out_log=$name.$(basename "$0" .sh).out.log

		# run job
		script=$(cat <<- EOS
				#!/bin/bash
				#PBS -l walltime=24:00:00
				#PBS -l select=1:mem=10gb:ncpus=1
				#PBS -j oe
				#PBS -N macs.$name
				#PBS -o $std_log
				${depend:-}

				printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

				# load modules
				module load samtools/1.2 &>> ${out_log}
				module load anaconda3/personal &>> ${out_log}
				source activate macs2 &>> ${out_log}

				# copy to scratch
				cp ${workdir}/${test}* . &>> ${out_log}
				${cntrl_cp:-} 

				# filter bams to remove dup/lowqual/scaffold
				if [[ $format = "BAM" ]]; then 
					printf "\nFiltering input bams:\n" >> ${out_log}
					samtools view -b -F 1024 -q 10 $(basename ${test}) ${chr_list} > tmp.bam
					mv tmp.bam $(basename ${test}); samtools index $(basename ${test})
					${cntrl_filt[@]:-}
				fi

				# run predictd for qc
				printf "\nPredicting insert size:\n" >> ${out_log}
				${predictd_cmd[@]} &>> ${out_log}
				Rscript ${name}.predictd.R &>> $out_log
				cp ${name}.predictd.R* $workdir/$qcdir/ &>> $out_log

				# call peaks
				printf "\nCalling peaks:\n" >> $out_log
				${callpeaks_cmd[@]} &>> $out_log

				cp ${name}_peaks* $workdir/$outdir/ &>> $out_log
				cp ${name}_*.bdg $workdir/$outdir/ &>> $out_log

				printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` >> $out_log
				ls -lhAR &>> $out_log
				cp $out_log $workdir/$logdir/
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
