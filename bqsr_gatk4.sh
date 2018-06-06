#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
check='y'
date=`date '+%Y-%m-%d %H:%M:%S'`

# help message
help_message="
Perform base quality score recalibration following GATK4 protocol. 

usage:
    bash $(basename $0) [-options] -i <bam>
required arguments:
    -i|--input : aligned & cleaned BAM file to recalibrate
    -F|--fasta : whole genome FASTA
    -I|--intervals : either a reference dict file or interval_list file
    -SNP|--dbsnp : dbsnp vcf
    -IND|--indels : vcf of known indels
optional arguments:
    -n|--name : name prefix for output files (default = extracted from first bam)
    -o|--outdir : outut directory (default = PWD)
    -q|--qcdir : output directory for qc metrics (default = --outdir)
    -l|--logdir : output directory for log files (default = --outdir)
    --depend : list of job dependencies to pass to PBS job script
    --check : whether to check files prior to running job [y|n] (default = y)
additional_info:
    # all paths should be relative to working directory
    # check and depend options useful for job scheduling
      e.g. --depend afterok:123456,afterok:123457
    # log/qc output directories inherit from outdir unless specified
    # 

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--input)
            bam=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -I|--intervals)
            intervals=$2
            shift
            ;;
        -SNP|--dbsnp)
            dbsnp=$2
            shift
            ;;
        -IND|--indels)
            indels=$2
            shift
            ;;
        -n|--name)
            name=$2
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
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required arg provided
if [[ -z ${bam:-} ]]; then
    printf "\nERROR: --input argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${intervals:-} ]]; then
    printf "\nERROR: --intervals argument required\n"
    echo "$help_message"; exit 1
elif [[ -z ${dbsnp:-} ]]; then
	printf "\nERROR: --dbsnp argument required\n"
	echo "$help_message"; exit 1
elif [[ -z ${indels:-} ]]; then
	printf "\nERROR: --indels argument required\n"
	echo "$help_message"; exit 1
fi 
    
# check files
if [[ $check = 'y' ]]; then
    if [[ ! -r $workdir/$fasta ]]; then 
        printf "\nERROR: FASTA file not readable: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    elif [[ ! -r $workdir/$bam ]]; then
        printf "\nERROR: BAM file not readable: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    elif [[ ! -r $workdir/$intervals ]]; then
        printf "\nERROR: intervals file is not readable: %s/%s\n" $workdir $intervals
        echo "$help_message"; exit 1
    elif [[ ! -r $workdir/$dbsnp ]]; then
		printf "\nERROR: dbsnp file is not readable: %s/%s\n" $workdir $dbsnp
		echo "$help_message"; exit 1
	elif [[ ! -r $workdir/$indels ]]; then
		printf "\nERROR: indels file is not readable: %s/%s\n" $workdir $indels
		echo "$help_message"; exit 1
	fi
fi

# get basenames/set merge inputs
fasta_prefix=${fasta%%.*}
fasta_base=$(basename "$fasta")
bam_base=$(basename "$bam")
indels_base=$(basename "$indels")
dbsnp_base=$(basename "$dbsnp")

# if no name provided extract from first bam
if [[ -z ${name:-} ]]; then
    name=${bam_base%%.*}
fi

# set/create output directories
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
if [[ -z "${qcdir:-}" ]]; then
    qcdir=$outdir
fi
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir

# determine interval type
if [[ ${intervals##*.} = 'dict' ]]; then
    printf "\nScattering intervals from dict file: %s\n" $intervals
    sn_array=($(grep ^@SQ $intervals | cut -f 2 | grep -Po "SN:\K.+$"))
    ln_array=($(grep ^@SQ $intervals | cut -f 3 | grep -Po "LN:\K[0-9]+$"))
    int_array=(); curr_sum=0
    for idx in "${!sn_array[@]}"; do
        curr_sum=$(($curr_sum + ${ln_array[$idx]}))
        if [[ $idx = 0 ]]; then
            sn_list="-L ${sn_array[$idx]}"
        elif [[ ${curr_sum} -gt 500000000 ]]; then
            int_array+=("$sn_list")
            sn_list="-L ${sn_array[$idx]}"; curr_sum=0
        else
            sn_list="${sn_list} -L ${sn_array[$idx]}"
        fi
    done
    int_array+=("$sn_list")
elif [[ ${intervals##*.} = 'interval_list' ]]; then
    echo "ERROR: script does not currently handle interval_list files"; exit 1
    #printf "\nScattering intervals from interval_list file: %s\n" $intervals
fi

## tmp - just scatter over standard chromosomes
#chr_array=($(grep ^@SQ $intervals | cut -f 2 | grep -E "^SN:" | grep -Eo "chr[0-9MYX]{1,2}$"))

# run base recal report generation scattered over intervals
for idx in ${!int_array[@]}; do
    
    gather_cp="${gather_cp:-}; cp $workdir/$outdir/${name}.${idx}.recal_report* ."
    gather_input="${gather_input:-}-I=${name}.${idx}.recal_report "

	# skip if exists
	if [[ -r "${outdir}/${name}.${idx}.recal_report" ]]; then
        continue
	fi

    # set intervals
    intervals=${int_array[$idx]}

	# set log file names
	std_log="${workdir}/${logdir}/${name}.${idx}.recal_rep.std.log"
	pbs_log="${workdir}/${logdir}/${name}.${idx}.recal_rep.pbs.log"
	out_log="${workdir}/${logdir}/${name}.${idx}.recal_rep.out.log"

	# write job script
	script=$(cat <<- EOS 
			#!/bin/bash
			#PBS -l walltime=24:00:00
			#PBS -l select=1:mem=12gb:ncpus=1
			#PBS -j oe
			#PBS -N bqsr.${name}.${idx}
			#PBS -q med-bio
			#PBS -o ${std_log}
			${depend:-}

			printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

			# load modules
			module load samtools/1.2
			module load java/jdk-8u66
			module load picard/2.6.0
			module load gatk/4.0

			# copy accross bam and fasta
			cp -L $workdir/$fasta_prefix* . &>> $out_log
			cp -L $workdir/$bam* . &>> $out_log
			cp -L $workdir/$indels* . &>> $out_log
			cp -L $workdir/$dbsnp* . &>> $out_log

			# gen bqsr report
			JAVA_OPTS="-Xmx12G" gatk BaseRecalibrator \
				-R ${fasta_base} \
				-I ${bam_base} \
				${intervals} \
				--use-original-qualities \
				-O ${name}.${idx}.recal_report \
				--known-sites ${dbsnp_base} \
				--known-sites ${indels_base} &>> $out_log

			# copy output back
			cp ${name}.${idx}.recal_report* ${workdir}/${outdir} &>> $out_log

			printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` &>> $out_log
			ls -lhAR &>> $out_log
			ls -lhAR
		EOS
	)
	echo "$script" > $pbs_log

	# submit job
	jobid=$(qsub "$pbs_log")

	# add to gather dependencies & input
	gather_depend="${gather_depend:-},afterok:$jobid"
done

# remove leading delim
if [[ -n ${gather_depend:-} ]]; then
    gather_depend="#PBS -W depend=${gather_depend#*,}"
fi
gather_cp=${gather_cp#*;}

# gather recal reports
if [[ -r "${workdir}/${outdir}/${name}.recal_report" ]]; then
	:
else
	std_log="${workdir}/${logdir}/${name}.gather_recal_rep.std.log"
	pbs_log="${workdir}/${logdir}/${name}.gather_recal_rep.pbs.log"
	out_log="${workdir}/${logdir}/${name}.gather_recal_rep.out.log"
	script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=12gb:ncpus=1
		#PBS -j oe
		#PBS -N gather_rep_${name}
		#PBS -q med-bio
		#PBS -o $std_log
		${gather_depend:-}

		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

		# load modules
		module load samtools/1.2
		module load java/jdk-8u66
		module load picard/2.6.0
		module load gatk/4.0

		# copy reports to scratch
		${gather_cp} &>> $out_log

		# gen bqsr report
		JAVA_OPTS="-Xmx12G" gatk GatherBQSRReports \
			${gather_input} \
			-O ${name}.recal_report &>> $out_log

		# copy output back
		cp ${name}.recal_report* $workdir/$outdir &>> $out_log

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` &>> $out_log
		ls -lhAR &>> $out_log
		ls -lhAR
		EOS
	)
	echo "$script" > $pbs_log
	jobid=$(qsub "$pbs_log")
	depend="#PBS -W depend=afterok:${jobid}"
fi

# reset dependencies
unset gather_depend
unset gather_cp

# apply recal scattered over intervals
for idx in ${!int_array[@]}; do

    gather_cp="${gather_cp:-}; cp ${workdir}/${outdir}/${name}.${idx}.recal.bam* ."
    gather_input="${gather_input:-}INPUT=${name}.${idx}.recal.bam "

    # skip if exists
    if [[ -r "${outdir}/${name}.${idx}.recal.bam" ]]; then
        continue
    fi

    # set log file names
    std_log="${workdir}/${logdir}/${name}.${idx}.recal.std.log"
    pbs_log="${workdir}/${logdir}/${name}.${idx}.recal.pbs.log"
    out_log="${workdir}/${logdir}/${name}.${idx}.recal.out.log"

    # set cmd
    apply_bqsr=("JAVA_OPTS=\"-Xmx12G\" gatk ApplyBQSR"
                "-R ${fasta_base}"
                "-I ${bam_base}" 
                "${int_array[$idx]}" 
                "-O ${name}.${idx}.recal.bam" 
                "-bqsr ${name}.recal_report" 
                "--use-original-qualities" 
                "--emit-original-quals"
                "--create-output-bam-index" 
                "--create-output-bam-md5" 
                "--add-output-sam-program-record")

    # write job script
    script=$(cat <<- EOS 
			#!/bin/bash
			#PBS -l walltime=24:00:00
			#PBS -l select=1:mem=12gb:ncpus=1
			#PBS -j oe
			#PBS -N apply.${name}.${idx}
			#PBS -q med-bio
			#PBS -o $std_log
			${depend:-}

			printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

			# load modules
			module load samtools/1.2
			module load java/jdk-8u66
			module load picard/2.6.0
			module load gatk/4.0

			# copy accross bam and fasta
			cp -L $workdir/$fasta_prefix* . &>> $out_log
			cp -L $workdir/$bam* . &>> $out_log
			cp -L $workdir/$outdir/${name}.recal_report . &>> $out_log

			# gen bqsr report
			${apply_bqsr[@]} &>> $out_log

			# copy output back
			cp ${name}.${idx}.recal.bam* $workdir/$outdir &>> $out_log

			printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` &>> $out_log
			ls -lhAR &>> $out_log
			ls -lhAR
		EOS
	)
    echo "$script" > $pbs_log

    # submit job
    jobid=$(qsub "$pbs_log")

    # add to gather dependencies & input
    gather_depend="${gather_depend:-},afterok:$jobid"
done

# remove leading delim
if [[ -n ${gather_depend:-} ]]; then
    gather_depend="#PBS -W depend=${gather_depend#*,}"
fi
gather_cp=${gather_cp#*;}

# gather bam
if [[ -r "${workdir}/${outdir}/${name}.recal_report" ]]; then
    echo "Final recalibrated bam already exists: %s/%s/%s" $workdir $outdir $name.recal.bam
else
    std_log="${workdir}/${logdir}/${name}.gather_bam.std.log"
    pbs_log="${workdir}/${logdir}/${name}.gather_bam.pbs.log"
    out_log="${workdir}/${logdir}/${name}.gather_bam.out.log"
	script=$(cat <<- EOS 
			#!/bin/bash
			#PBS -l walltime=24:00:00
			#PBS -l select=1:mem=12gb:ncpus=1
			#PBS -j oe
			#PBS -N gather_bam_${name}
			#PBS -q med-bio
			#PBS -o $std_log
			${gather_depend:-}

			printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` > $out_log

			# load modules
			module load samtools/1.2
			module load java/jdk-8u66
			module load picard/2.6.0
			module load gatk/4.0

			# copy reports to scratch
			${gather_cp} &>> $out_log
			cp $workdir/$bam* . &>> $out_log

			# gen bqsr report
			java -Xmx12G -jar /apps/picard/2.6.0/picard.jar GatherBamFiles \
				${gather_input} \
				INPUT=$bam_base \
				OUTPUT=${name}.recal.bam \
				CREATE_INDEX=true \
				CREATE_MD5=true &>> $out_log

			# copy output back
			cp ${name}.recal.bam* $workdir/$outdir &>> $out_log
			cp ${name}.recal.bai $workdir/$outdir &>> $out_log

			printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\` &>> $out_log
			ls -lhAR &>> $out_log
			ls -lhAR
		EOS
	)
	echo "$script" > $pbs_log
	jobid=$(qsub "$pbs_log")
	echo "JOBID: $jobid"
fi
exit 0
