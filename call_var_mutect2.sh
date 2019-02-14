#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=""
check="yes"

# help message
help_message="
Wrapper to run Mutect2 (GATK v3).

usage:
    bash $(basename "$0") [-options] -t <TUMOR_BAM> -f <FASTA>
required arguments:
    -t|--tumor : tumor sample BAM file
    -f|--fasta : whole genome FASTA file
optional arguments:
    -c|--control : control sample alignments [BAM] (default = NULL)
    --intervals : intervals file [BED,interval] (default = NULL)
    --chr : chromosome to restrict analysis to [string] (default = NULL)
    --dbsnp : dbSNP file [VCF] (default = NULL)
    --cosmic : cosmic file [VCF] (default = NULL)
    --pon : PON file [VCF] (default = NULL)
    -o|--outdir : output directory [string] (default = '.')
    -l|--logdir : output directory for log files [string] (default = --outdir)
    -n|--name : prefix to give output files [string] (default = extracted from bam)
    --check : whether to check input files [yes,no] (default = 'yes')
additional info:
    # check argument is useful for scheduling jobs where input files 
      may not presently exist
    # control, intervals, dbsnp, cosmic and pon arguments are all optional but
      if provided, they will be passed to mutect2
output:
    # creates logs and variants directories in outdir containing 
      run log and output VCF files respectively
"

# parse arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -t|--tumor)
            tumor=$2
            shift
            ;;
        -f|--fasta)
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
        -n|--name)
            name=$2
            shift
            ;;
        --check)
            check=$2
            shift
            ;;
        --intervals)
            intervals=$2
            shift
            ;;
        --chr)
            chr=$2
            shift
            ;;
        --dbsnp)
            dbsnp=$2
            shift
            ;;
        --cosmic)
            cosmic=$2
            shift
            ;;
        --pon)
            pon=$2
            shift
            ;;
        *)
            printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 2
            ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arg
if [[ -z ${tumor:-} ]]; then
	printf "\nERROR: no tumor argument provided\n"
	echo "$help_message"; exit 2
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: no FASTA argument provided\n"
    echo "$help_message"; exit 2
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$tumor ]]; then
        printf "\nERROR: tumor file is not readable: %s/%s\n" $workdir $tumor
        echo "$help_message"; exit 2
    elif [[ ! -r $workdir/$fasta ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 2
    fi
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
	name=$(basename ${tumor})
    name=${name%%.*}
fi

# get basenames/prefix
tumor_base=$(basename "$tumor")
fasta_base=$(basename "$fasta")
fasta_prefix=${fasta%%.*}

# setup output directories
mkdir -p $workdir/$outdir/logs
mkdir -p $workdir/$outdir/variants

# set optional arguments
if [[ ! -z ${control:-} ]]; then
   if [[ $check = yes ]] && [[ ! -r $workdir/$control ]]; then
        printf "\nERROR: control file is not readable: %s/%s\n" $workdir $control
        echo "$help_message"; exit 2
    else
        control_base=$(basename $control)
        control_cp="cp $workdir/$control* ."
        control_arg="-I:normal $control_base"
    fi
fi
if [[ ! -z ${intervals:-} ]]; then
   if [[ $check = yes ]] && [[ ! -r $workdir/$intervals ]]; then
        printf "\nERROR: intervals file is not readable: %s/%s\n" $workdir $intervals
        echo "$help_message"; exit 2
    else
        intervals_base=$(basename $intervals)
        intervals_cp="cp $workdir/$intervals* ."
        intervals_arg="-L $intervals_base"
    fi
fi
if [[ ! -z ${dbsnp:-} ]]; then
    if [[ $check = yes ]] && [[ ! -r $workdir/$dbsnp ]]; then
        printf "\nERROR: dbSNP file is not readable: %s/%s\n" $workdir $dbsnp
        echo "$help_message"; exit 2
    else
        dbsnp_base=$(basename $dbsnp)
        dbsnp_cp="cp $workdir/$dbsnp* ."
        dbsnp_arg="--dbsnp $dbsnp_base"
    fi
fi
if [[ ! -z ${cosmic:-} ]]; then
    if [[ $check = yes ]] && [[ ! -r $workdir/$cosmic ]]; then
        printf "\nERROR: cosmic file is not readable: %s/%s\n" $workdir $cosmic
        echo "$help_message"; exit 2
    else
        cosmic_base=$(basename $cosmic)
        cosmic_cp="cp $workdir/$cosmic* ."
        cosmic_arg="--cosmic $cosmic_base"
    fi
fi
if [[ ! -z ${pon:-} ]]; then
    if [[ $check = yes ]] && [[ ! -r $workdir/$pon ]]; then
        printf "\nERROR: PON file is not readable: %s/%s\n" $workdir $pon
        echo "$help_message"; exit 2
    else
        pon_base=$(basename $pon)
        pon_cp="cp $workdir/$pon* ."
        pon_arg="--normal_panel $pon_base"
    fi
fi
if [[ ! -z ${chr:-} ]]; then
    chr_arg="-L $chr"
fi


# call variants
jobid=$(cat <<- EOS | qsub -N $name.var -
		#!/usr/bin/env bash
		#PBS -l walltime=40:00:00
		#PBS -l select=1:mem=40gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.call_var_mutect2.log
		
		# load modules
		module load picard/2.6.0
		module load gatk/3.6
		module load java/jdk-8u66
		module load samtools/1.2

		# copy files to scratch
		cp $workdir/$tumor* .
		cp -rL $workdir/$fasta_prefix* .
		${control_cp:-}
		${intervals_cp:-}
		${dbsnp_cp:-}
		${cosmic_cp:-}
		${pon_cp:-}

		# call
		java -jar /apps/gatk/3.6/GenomeAnalysisTK.jar \
			-T MuTect2 \
			-R $fasta_base \
			-I:tumor $tumor_base \
			-o $name.mutect2.vcf \
			${control_arg:-} \
			${intervals_arg:-} \
			${chr_arg:-} \
			${dbsnp_arg:-} \
			${cosmic_arg:-} \
			${pon_arg:-}

		cp $name.mutect2.vcf* $workdir/$outdir/
		ls -lhAR
	EOS
)
