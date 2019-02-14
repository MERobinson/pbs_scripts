#!/bin/bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''

# help message
help_message="

Simple wrapper to build Bowtie2 index with bowtie2-build

usage:
    bash $(basename "$0") [-options] -F <FASTA>
required arguments:
    -f|--fasta : path to input fasta to index
optional arguments:
    -n|--name : name prefix for index files (default = fa name)
    -o|--outdir : output directory for index files (default = pwd)
    -l|--logdir : output directory for log files (default = outdir)
    -s|--species : species name - added to sequence dictionary if provided
    -g|--genome : genome version - added to sequence dictionary if provided
additional info:
    # will generate fasta idx and dictionary if not already present

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -f|--fasta)
            fasta=$2
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
        -s|--species)
            species="SPECIES=$2"
            shift
            ;;
        -g|--genome)
            genome="GENOME_ASSEMBLY=$2"
            shift
            ;;
        *)
            echo "Error: Illegal argument"
            echo "$help_message"; exit 1
            ;;
    esac
    shift
done

# check required argument
if [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
fi

# check res files exist
if [[ ! -r ${workdir}/${fasta} ]]; then
	printf "\nERROR: Input fasta file does not exist: %s/%s\n" $workdir $fasta 
	echo "$help_message"; exit 1
fi

# create required outdir
if [[ -z ${logdir:-} ]]; then
    logdir=$outdir
fi
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# get basename/prefix
fasta_base=$(basename ${fasta})
fasta_prefix=${fasta%.fa}

# set outfile name
if [[ -z ${name:-} ]]; then
    name=${fasta_base%.fa}
fi

# check if index & dict exist
if [[ ! -e ${workdir}/${fasta_prefix}.fa.fai ]]; then
    echo "WARNING: no fai file detected, will generate."
    index_cmd=("samtools faidx ${fasta_base}")
fi
if [[ ! -e ${workdir}/${fasta_prefix}.dict ]]; then
    echo "WARNING: no dict file detected, will generate."
    dict_cmd=("java -jar /apps/picard/2.6.0/picard.jar CreateSequenceDictionary"
              "REFERENCE=${fasta_base}"
              "OUTPUT=${fasta_base%.fa}.dict"
              "${species:-} ${genome:-}")
fi

# set logfiles
scr_name=$(basename $0 .sh)
std_log=$workdir/$logdir/$name.$scr_name.std.log
pbs_log=$workdir/$logdir/$name.$scr_name.pbs.log

# run job
script=$(cat <<- EOS 
		#!/bin/bash
		#PBS -l walltime=24:00:00
		#PBS -l select=1:mem=24gb:ncpus=8
		#PBS -j oe
		#PBS -N bwa_index
		#PBS -q med-bio
		#PBS -o ${std_log}
	
		# load modules
		module load samtools/1.2
		module load java/jdk-8u66
		module load picard/2.6.0
		module load bowtie2/2.2.9
	
		printf "\nSTART: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`

		# copy fasta to scratch
		cp -L ${workdir}/${fasta_prefix}* .

		# create seq dict if req
		${dict_cmd[@]:-}

		# create fasta index if req
		${index_cmd[@]:-}

		# create BWA index
		bowtie2-build --threads 8 ${fasta_base} $name 
		
		# copy index files to outdir
		cp ${name}* $workdir/$outdir
		cp ${fasta_base}* $workdir/$outdir

		printf "\nEND: %s %s\n" \`date '+%Y-%m-%d %H:%M:%S'\`
		ls -lhAR

		EOS
)
echo "$script" > $pbs_log

# submit job
jobid=$(qsub "$pbs_log")

# echo job id and exit
echo "JOBID: $jobid"
exit 0
