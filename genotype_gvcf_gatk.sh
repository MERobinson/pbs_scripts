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
    -o|--outdir : output directory relative to current directory for output (default = NULL)
    -n|--name : prefix to give output files [string] (default = 'germline_var')
    -l|--logir : output directory for log file (if different to outdir) (default = NULL)
    --check : whether to check arguments prior to running PBS script [yes|no] (default = 'yes')
    --depend : list of job dependencies to pass to PBS script (default = NULL)
output:
    # outputs VCF and index in outdir & run log in logdir

"

# parse arguments
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
        *)
            printf "\nERROR: undefined argument provided\n"
            echo "$help_message"
            ;;
    esac
    shift
done

# set scriptdir
if [[ -z ${scriptdir:-} ]]; then
    scriptdir=$outdir
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

# fasta file basename
fasta_base=${fasta%%.*}
fasta=$(basename "$fasta")

# set args for each gcvf
for gvcf in ${gvcf_array[@]}; do
    gvcf_base=$(basename "$gvcf")
    gatk_arg="${gatk_arg:-} --variant $gvcf_base"
    cp_arg="${cp_arg:-}; cp $workdir/$gvcf* ."
done
gatk_arg=${gatk_arg#* }
cp_arg=${cp_arg#*; }

# run pbs script
jobid=$(cat <<- EOS | qsub -N $name.ggvcf -
		#!/bin/bash
		#PBS -j oe
		#PBS -l walltime=05:00:00
		#PBS -l select=1:mem=40gb:ncpus=1
		#PBS -o $workdir/$outdir/$name.genotype_gvcf_gatk.log

		# load modules
		module load java/jdk-8u66
		module load gatk/3.6
		module load samtools/1.2
		module load picard/2.6.0

		# copy files to scratch
		$cp_arg
		cp -rL $workdir/$fasta_base* .

		java -jar /apps/gatk/3.6/GenomeAnalysisTK.jar \
			-T GenotypeGVCFs \
			-R $fasta \
			-o $name.g.vcf \
			$gatk_arg

		cp $name.g.vcf* $workdir/$outdir
	EOS
)
echo "JOBID: $jobid"
