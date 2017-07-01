#!/bin/bash

# default arg
workdir=$PWD
resdir=resources
logdir=logs
fasta=genome.fa

# help message
help_message="

bash $(basename "$0") [-wrfgh]

purpose:
	Simple wrapper to generate BWA index
required arguments:
    n/a
optional arguments:
    -w|--workdir : working directory - used as base dir for all input/output (default = pwd) 
    -f|--fasta : name of fasta file in workdir (default = genome.fa)
    -l|--logdir : output directory for log files (default = logs)
    -s|--species : species name - added to sequence dictionary if provided
    -g|--genome : genome version - added to sequence dictionary if provided
example:
    bash $(basename "$0") -g lincRNA.gtf
additional info:
    > Recommended to setup a new dir per index containing a link to the genome
      FASTA file for indexing - picard CreateSequenceDictionary, samtools faidx and
      bwa index will be run to create all necessary genome resources.

"

# parse command line arg
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-w|--workdir)
		workdir=$2
		shift
		;;
		-f|--fasta)
		fasta=$2
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
		echo "$help_message"
		exit 1
		;;
	esac
	shift
done

# check res files exist
if [[ ! -e "$workdir/$fasta" ]]; then
	printf "\nERROR: Input fasta file does not exist: %s/%s\n" $workdir $fasta 
	echo "$help_message"; exit 1
fi

# create required dirs
mkdir -p $workdir/$logdir

# basename
fasta_base=${fasta%%.*}

# run job
jobid=$(cat <<- EOS | qsub -N BWA_index - 
	#!/bin/bash
	#PBS -l walltime=40:00:00
	#PBS -l select=1:mem=20gb:ncpus=20
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $workdir/$logdir/generate_index_bwa.log
	
	# load modules
	module load samtools/1.2
	module load java/jdk-8u66
	module load picard/2.6.0
	module load bio-bwa/0.7.15
	
	# copy accross fasta
	cp -L $workdir/$fasta* .
	mkdir -p tmp

	# create seq dict
	java -jar /apps/picard/2.6.0/picard.jar CreateSequenceDictionary \
		REFERENCE=$fasta \
		OUTPUT=$fasta_base.dict \
		$species $genome
	cp $fasta_base.dict $workdir

	# create fasta index
	samtools faidx $fasta
	cp $fasta.fai $workdir

	# create BWA index
	bwa index -a bwtsw $fasta 
	cp $fasta* $workdir

	EOS
)
