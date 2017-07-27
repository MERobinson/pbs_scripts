#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=""
check="yes"
gvcf="no"

# help message
help_message="
usage:
    bash $(basename "$0") [-options] -b <sample1.bam> -f <genome.fa>
purpose:
    # simple wrapper to run HaplotypeCaller (GATK v3)
required arguments:
    -b|--bam : aligned and preprocessed BAM file
    -f|--fasta : whole genome FASTA file
optional arguments:
    -o|--outdir : output directory to create within workdir [string] (default = '')
    -n|--name : prefix to give output files [string] (default = extracted from bam)
    --gvcf : whether to run in GVCF mode [yes,no] (default = 'no')
    --intervals : intervals file [BED,interval] (default = NULL)
    --dbsnp : dbSNP file [VCF] (default = NULL)
    --check : whether to check input files [yes,no] (default = 'yes')
    --depend : list of dependencies for PBS script (e.g. afterok:012345,afterok:012346)
additional info:
    # check argument is useful for scheduling future jobs where input files 
      may not presently exist
    # interval, dbsnp and gvcf arguments are all optional but
      if provided, they will be passed to HaplotypeCaller
output:
    # creates logs and variants directories in outdir containing 
      run log and output VCF files respectively
"

# parse arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bam)
            bam=$2
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
        --gvcf)
            gvcf=$2
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
if [[ -z ${bam:-} ]]; then
	printf "\nERROR: no bam argument provided\n"
	echo "$help_message"; exit 2
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: no FASTA argument provided\n"
    echo "$help_message"; exit 2
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$bam ]]; then
        printf "\nERROR: bam file is not readable: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 2
    elif [[ ! -r $workdir/$fasta ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 2
    fi
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
	name=$(basename ${bam})
    name=${name%%.*}
fi

# get basenames/prefix
bam_base=$(basename "$bam")
fasta_base=$(basename "$fasta")
fasta_prefix=${fasta%%.*}

# setup output directories
mkdir -p $workdir/$outdir/logs
mkdir -p $workdir/$outdir/variants

# gvcf options
if [[ $gvcf = "yes" ]]; then
    vcf_ext="g.vcf"
    gvcf_arg="-ERC GVCF"
else
    vcf_ext="vcf"
fi

# deal with potential inputs
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

# call variants
jobid=$(cat <<- EOS | qsub -N $name.var -
		#!/usr/bin/env bash
		#PBS -l walltime=60:00:00
		#PBS -l select=1:mem=40gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$outdir/logs/$name.call_var_hc.$vcf_ext.log
		${depend:-}

		# load modules
		module load picard/2.6.0
		module load gatk/3.6
		module load java/jdk-8u66
		module load samtools/1.2

		# copy files to scratch
		cp $workdir/$bam* .
		cp -rL $workdir/$fasta_prefix* .
		${intervals_cp:-}
		${dbsnp_cp:-}

		# call
		java -jar /apps/gatk/3.6/GenomeAnalysisTK.jar \
			-T HaplotypeCaller \
			-R $fasta_base \
			-I $bam_base \
			-o $name.hc.$vcf_ext \
			${gvcf_arg:-} \
			${intervals_arg:-} \
			${dbsnp_arg:-}

		cp $name.hc.$vcf_ext* $workdir/$outdir/variants/
		
        printf "DATE: %s\n" $(date "+%Y-%m-%d")
        printf "JOBID: %s\n" '${PBS_JOBID:-}'
        ls -lhAR
	EOS
)
echo "JOBID: $jobid"
