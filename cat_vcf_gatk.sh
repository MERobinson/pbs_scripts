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
    bash $(basename "$0") [-options] -v <list,of,vcf> -f <genome.fa>
purpose:
    # simple wrapper to concatenate VCF files (GATK v3)
required arguments:
    -v|--vcf : comma separated list of VCF files to cat
    -f|--fasta : whole genome FASTA file
optional arguments:
    -o|--outdir : output directory to create within workdir [string] (default = '')
    -n|--name : name of output file [string] (default = 'combined.vcf')
    --check : whether to check input files [yes,no] (default = 'yes')
    --depend : dependency list to pass to PBS script - e.g. "afterok:1920830,afterok:1920831"
additional info:
    # check argument is useful for scheduling future jobs where input files 
      may not presently exist
output:
    # creates logs and variants directories in outdir containing 
      run log and output VCF files respectively

"

# parse arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -v|--vcf)
            vcf_list=$2
            shift
            ;;
        -f|--fasta)
            fasta=$2
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
        --depend)
            depend="#PBS -W depend=$2"
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
if [[ -z ${vcf_list:-} ]]; then
	printf "\nERROR: no vcf argument provided\n"
	echo "$help_message"; exit 2
else
    IFS="," read -r -a vcf_array <<< "$vcf_list"
fi
if [[ -z ${fasta:-} ]]; then
    printf "\nERROR: no fasta argument provided\n"
    echo "$help_message"; exit 2
fi

# check files unless flagged
if [[ $check = yes ]]; then
    for vcf in ${vcf_array[@]}; do
        if [[ -z $workdir/${vcf} ]]; then
            printf "\nERROR: VCF file is not readable: $workdir/$vcf\n"
            echo "$help_message"; exit 2
        fi
    done
    if [[ ! -r $workdir/$fasta ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 2
    fi
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
    name="combined"
fi

# get basename/prefix for fasta
fasta_base=$(basename "$fasta")
fasta_prefix=${fasta%%.*}

# set arg for VCF list copying and command input
cp_arg=$(printf "cp $workdir/%s* . ;" ${vcf_array[@]})
v_arg=$(printf -- "-V %s " ${vcf_array[@]##*/})

# setup output directories
mkdir -p $workdir/$outdir/logs
mkdir -p $workdir/$outdir/variants

# cat variants
jobid=$(cat <<- EOS | qsub -N $name.cat -
		#!/bin/bash
		#PBS -l walltime=08:00:00
		#PBS -l select=1:mem=20gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$outdir/logs/$name.cat_vcf_gatk.log
		${depend:-}		

		# load modules
		module load picard/2.6.0
		module load gatk/3.6
		module load java/jdk-8u66
		module load samtools/1.2

		# copy files to scratch
		$cp_arg
		cp -rL $workdir/$fasta_prefix* .

		# call
		java -cp /apps/gatk/3.6/GenomeAnalysisTK.jar \
			org.broadinstitute.gatk.tools.CatVariants \
			-R $fasta_base \
			$v_arg \
			-out $name.vcf

		cp $name.vcf* $workdir/$outdir/
		
		printf "DATE: %s\n" $(date "+%Y-%m-%d")
		printf "JOBID: %s\n" '${PBS_JOBID:-}'
		ls -lhAR
	EOS
)
echo "JOBID: $jobid"
