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
usage:
    bash $(basename "$0") [-options] -t <sample1.T.bam> -f <genome.fa>
purpose:
    # simple wrapper to run mutect2 (GATK v3)
required arguments:
    -t|--tumor : tumor sample BAM file
    -f|--fasta : whole genome FASTA file
optional arguments:
    -c|--control : control sample alignments [BAM] (default = NULL)
    --intervals : intervals file [BED,interval] (default = NULL)
    --dbsnp : dbSNP file [VCF] (default = NULL)
    --cosmic : cosmic file [VCF] (default = NULL)
    --pon : PON file [VCF] (default = NULL)
    -o|--outdir : output directory to create within workdir [string] (default = '')
    -n|--name : prefix to give output files [string] (default = extracted from bam)
    --check : whether to check input files [yes,no] (default = 'yes')
additional info:
    # check argument is useful for scheduling futre jobs where input files 
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
            printf "\nERROR: Undefined argument provided\n"
            echo "$help_message"; exit 2
            ;;
    esac
    shift
done

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
        control_arg="-I:control $control_base"
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
        cosmic_arg="--dbsnp $cosmic_base"
    fi
fi
if [[ ! -z ${pon:-} ]]; then
    if [[ $check = yes ]] && [[ ! -r $workdir/$pon ]]; then
        printf "\nERROR: PON file is not readable: %s/%s\n" $workdir $pon
        echo "$help_message"; exit 2
    else
        pon_base=$(basename $pon)
        pon_cp="cp $workdir/$pon* ."
        pon_arg="--dbsnp $pon_base"
    fi
fi


# call variants
jobid=$(cat <<- EOS | qsub -N $name.var -
		#!/usr/bin/env bash
		#PBS -l walltime=40:00:00
		#PBS -l select=1:mem=40gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$outdir/logs/$name.call_var_mutect2.log
		
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
			${control_arg:-} ${intervals_arg:-} \
			${dbsnp_arg:-} ${cosmic_arg:-} ${pon_arg:-}

		cp $name.mutect2.vcf* $workdir/$outdir/variants/
		ls -lhAR
	EOS
)
