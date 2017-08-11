#!/user/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=""
check="yes"
name="germline_var"

# help message
help_message="
usage:
    bash $(basename $0) [-options] -g <gvcf,list> -f <genome.fa>
purpose:
    # simple wrapper to run GenotypeGVCFs (GATK v3)
required arguments:
    -g|--gvcf : comma separated list of GVCF files
    -f|--fasta : whole genome FASTA used for generating GVCF
optional aruments:
    -o|--outdir : output directory for VCF files (default = '.')
    -n|--name : prefix to give output files [string] (default = 'germline_var')
    -l|--logir : output dir for log files [string] (default = --outdir value)
    --chr : chromosome to restrict analysis to [string] (default = NULL)
    --intervals : intervals file (passed to GATK as intervals arg)
    --check : whether to check input files [yes|no] (default = 'yes')
    --depend : list of job dependencies to pass to PBS script (default = NULL)
output:
    # outputs VCF and index in outdir & run log in logdir
additional info:
    # outdir/logdir paths should be given relative to current working dir, 
      and are created if not already present 
    # check and depend arguments are used for job scheduling/pieplines
    # intervals and chr arguments are optional but if provided,
      they will be passed to GenotypeGVCFs

"

# parse arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -g|--gvcf)
            gvcf_list=$2
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
        --check)
            check=$2
            shift
            ;;
        --depend)
            depend=$2
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
        --depend)
            depend="#PBS -W depend=$2"
            shift
            ;;
        *)
            printf "\nERROR: undefined argument provided: %s\n" $1
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
if [[ -z ${gvcf_list:-} ]]; then
    printf "\nERROR: no GVCF files provided\n"
    echo "$help_message"; exit 1
else
    IFS=',' read -ra gvcf_array <<< "$gvcf_list"
fi  
if [[ -z ${fasta:-} ]]; then
    printf "\nERROR: no FASTA file provided\n"
    echo "$help_message"; exit 1
fi

# check files unless flagged
if [[ $check = yes ]]; then
    if [[ ! -r $workdir/$fasta ]]; then
        printf "\nERROR: FASTA file not found/not readable: %s/%s" $workdir $fasta
        echo "$help_message"; exit 1
    fi
    for gvcf in ${gvcf_array[@]}; do
        if [[ ! -r $workdir/$gvcf ]]; then
            printf "\nERROR: GVCF file not found/not readable: %s/%s" $workdir $gvcf
            echo "$help_message"; exit 1
        fi
    done
fi

# fasta prefix (for copying dict/idx)
fasta_prefix=${fasta%%.*}
fasta=$(basename "$fasta")

# seyup output dir
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# set args for each gcvf
for gvcf in ${gvcf_array[@]}; do
    gvcf_base=$(basename "$gvcf")
    gvcf_gatk_arg="${gvcf_gatk_arg:-} --variant $gvcf_base"
    gvcf_cp_arg="${gvcf_cp_arg:-}; cp $workdir/$gvcf* ."
done
gvcf_gatk_arg=${gvcf_gatk_arg#* }
gvcf_cp_arg=${gvcf_cp_arg#*; }

# optional intervals arg
if [[ ! -z ${intervals:-} ]]; then
    if [[ $check = yes ]] && [[ ! -r $workdir/$intervals ]]; then
        printf "\nERROR: intervals file not found/not readable: %s/%s\n" $workdir $intervals
        echo "$help_message"; exit 2
    else
        intervals_base=$(basename $intervals)
        intervals_cp_arg="cp $workdir/$intervals* ."
        intervals_gatk_arg="-L $intervals_base"
    fi
fi
if [[ ! -z ${chr:-} ]]; then
    chr_gatk_arg="-L $chr"
fi

# run pbs script
jobid=$(cat <<- EOS | qsub -N $name.ggvcf -
		#!/usr/bin/env bash
		#PBS -l walltime=10:00:00
		#PBS -l select=1:mem=40gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/$logdir/$name.genotype_gvcf_gatk.log
		${depend:-}

		# load modules
		module load java/jdk-8u66
		module load gatk/3.6
		module load samtools/1.2
		module load picard/2.6.0

		# copy files to scratch
		$gvcf_cp_arg
		${intervals_cp_arg:-}
		cp -rL $workdir/$fasta_prefix* .

		java -jar /apps/gatk/3.6/GenomeAnalysisTK.jar \
			-T GenotypeGVCFs \
			-R $fasta \
			-o $name.vcf \
			$gvcf_gatk_arg \
			${intervals_gatk_arg:-} \
			${chr_gatk_arg:-}

		cp $name.vcf* $workdir/$outdir
	EOS
)
echo "JOBID: $jobid"
