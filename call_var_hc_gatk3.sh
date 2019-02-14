#!/usr/bin/env bash
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
    -o|--outdir : output directory for VCF files [string] (default = '.')
    -n|--name : prefix to give output files [string] (default = extracted from bam)
    -l|--logdir : output dir for log files [string] (default = --outdir value)
    --gvcf : whether to run in GVCF mode [yes,no] (default = 'no')
    --intervals : intervals file [BED,interval] (default = NULL)
    --chr : chromosome to restrict analysis to [string] (default = NULL)
    --dbsnp : dbSNP file [VCF] (default = NULL)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list of dependencies for PBS script [string] (e.g. afterok:012345,afterok:012346)
additional info:
    # outdir/logdir paths should be given relative to current working dir, 
      and are created if not already present 
    # check and depend arguments are used for job scheduling/pieplines
    # intervals, chr, dbsnp and gvcf arguments are all optional but
      if provided, they will be passed to HaplotypeCaller
output:
    # outputs VCF and index file to outdir & run log to logdir
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
        -l|--logdir)
            logdir=$2
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
        --chr)
            chr=$2
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
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# set logdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi

# check required arg
if [[ -z ${bam:-} ]]; then
	printf "\nERROR: no bam argument provided\n"
	echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: no FASTA argument provided\n"
    echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$bam ]]; then
        printf "\nERROR: bam file is not readable: %s/%s\n" $workdir $bam
        echo "$help_message"; exit 1
    elif [[ ! -r $workdir/$fasta ]]; then
        printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
        echo "$help_message"; exit 1
    fi
fi

# get sample name if not provided
if [[ -z ${name:-} ]]; then
	name=$(basename ${bam})
    name=${name%%.*}
fi

# get basenames/prefix
bam_base=$(basename "$bam")
fasta_prefix=${fasta%%.*}
fasta=$(basename "$fasta")

# setup output dir
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

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
if [[ ! -z ${chr:-} ]]; then
    chr_arg="-L $chr"
fi

# call variants
jobid=$(cat <<- EOS | qsub -N $name.var -
		#!/usr/bin/env bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=40gb:ncpus=8
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.call_var_hc.$vcf_ext.log
		${depend:-}

		# load modules
		module load picard/2.6.0
		module load gatk/3.7
		module load java
		module load samtools/1.2

		# copy files to scratch
		cp $workdir/$bam* .
		cp -rL $workdir/$fasta_prefix* .
		${intervals_cp:-}
		${dbsnp_cp:-}

		# call
		java -jar /apps/gatk/3.7/GenomeAnalysisTK.jar \
			-T HaplotypeCaller \
			-R $fasta \
			-I $bam_base \
			-o $name.hc.$vcf_ext \
			-nct 1 \
			${chr_arg:-} \
			${gvcf_arg:-} \
			${intervals_arg:-} \
			${dbsnp_arg:-}

		cp $name.hc.$vcf_ext* $workdir/$outdir
		
        ls -lhAR
	EOS
)
echo "JOBID: $jobid"
