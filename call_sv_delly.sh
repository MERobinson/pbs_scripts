#!/bin/bash

# default arg
WORKDIR=$PWD
OUTDIR=variants
INDIR=bam
GENOME=hg38
RESDIR=resources
FASTA=genome.fa
EXCLUDE=exclude.tsv
TOOLDIR=tools

# help message
USAGE="

usage:
	bash $(basename "$0") [-wirofenh] -t <BAM> -c <BAM>
purpose:
	Simple wrapper to call structural variants with Delly
required arguments:
        -t|--tumor : test sample (aligned bam)
	-c|--control : control sample (aligned bam)
optional arguments:
	-w|--workdir : path of working directory to use (default = pwd)
	-i|--indir : directory in which to find input bam files (default = bam)
	-r|--resdir : directory in which to find resource files (default = resources)
	-o|--outdir : name for output directory (default = variants)
	-f|--fasta : name of whole genome fasta file within resdir (default = genome.fa)
	-e|--excl : name of exclusion regions to mask (default = exclude.tsv) 
	-n|--name : name prefix for output files (deafult = extracted from tumor sample)
	-v|--vartypes : comma sep. list of variant types [INS,DEL,INV,BND,DUP] (default = all)
	--plot : whether to summarise and plot SVs with svprops (see additional info)
	--tooldir : dir containing additional tools (default = tools)
	-h|--help : print current help message
example:
        bash $(basename "$0") -t patient1_T_WGS.bam -c patient1_C_WGS.bam
inputs:
	> Tumor and Control BAM files (PE, aligned BAM) [required]
	> Whole genome fasta [required]
	> Regions to exclude [optional]
additional info:
	> --indir should be an existing folder within --workdir
	> --resdir should be an existing folder within --workdir
	> --resdir should contain links to required genome resources for analysis
	> --outdir will be created, if not pre-existing
	> name is extracted from test filename up to first period (.)
	> if wanting to run with no region exclusion set --excl NULL
	> --plot requires installation of svprops:
	  1) clone and install from https://github.com/dellytools/svprops
	  2) link svprops sampleprops and svprops.R to --tooldir
	> see Delly man for additional info: https://github.com/dellytools/delly

"

# parse command line arg
while [[ $# -gt 1 ]]; do
        key=$1
	case $key in
		-t|--tumor)
		TUMOR=$2
		shift
		;;
		-c|--control)
		CONTROL=$2
		shift
		;;
		-w|--workdir)
		WORKDIR=$2
		shift
		;;
		-i|--indir)
		INDIR=$2
		shift
		;;
		-r|--resdir)
		RESDIR=$2
		shift
		;;
		-o|--outdir)
		OUTDIR=$2
		shift
		;;
		-f|--fasta)
		FASTA=$2
		shift
		;;
		-e|--excl)
		EXCLUDE=$2
		shift
		;;
		-n|--name)
		NAME=$2
		shift
		;;
		-v|--vartypes)
		VARTYPES=$2
		shift
		;;
		--plot)
		PLOT=1
		shift
		;;
		--tooldir)
		TOOLDIR=$2
		shift
		;;
		-h|--help)
		echo "$USAGE"
		exit 1
		;;
		*)
		printf "\nERROR: Illegal argument\n"
		echo "$USAGE"
		exit 1
		;;
	esac
	shift
done

# check dir
if [[ ! -d $WORKDIR/$INDIR ]]; then
	printf "\nERROR: Input directory not found: $WORKDIR/$INDIR \n"
	echo "$USAGE"; exit 1
elif [[ ! -d $WORKDIR/$RESDIR ]]; then
	printf "\nERROR: Resource directory not found: $WORKDIR/$RESDIR \n"
fi

# check required arguments
if [[ -z "$TUMOR" ]]; then
        printf "\nERROR: No tumor sample provided\n"
        echo "$USAGE"; exit 1
elif [[ ! -e $WORKDIR/$INDIR/$TEST ]]; then
	printf "\nERROR: Tumor file not found: $WORKDIR/$INDIR/$TUMOR \n"
	echo "$USAGE"; exit 1
elif [[ -z "$CONTROL" ]]; then
	printf "\nERROR: No control sample provided\n"
	echo "$USAGE"; exit 1
elif [[ ! -e $WORKDIR/$INDIR/$TEST ]]; then
	printf "\nERROR: Control file not found: $WORKDIR/$INDIR/$CONTOL \n"
	echo "$USAGE"; exit 1
elif [[ ! -e $WORKDIR/$RESDIR/$FASTA ]]; then
	printf "\nERROR: FASTA file not found: $WORKDIR/$RESDIR/$FASTA \n"  
fi

# get name if req
if [[ -z "$NAME" ]]; then
	NAME=${TUMOR%%.*}
fi
TUMOR_NAME=${TUMOR%%.*}
CNTRL_NAME=${CONTROL%%.*}

# split fasta filename
FASTA_BASE=${FASTA%%.*}

# set exclusion argument if non-null
if [[ ! $EXCLUDE = "NULL" ]]; then
	EXCLUDE_ARG="-x $EXCLUDE"
	EXCLUDE_CP="cp -rL $WORKDIR/$RESDIR/$EXCLUDE ."
fi

# check tools if plotting flagged
if [[ $PLOT = 1 ]]; then
	if [[ ! -e $WORKDIR/$INDIR/$TOOLDIR/svprops ]] || \ 
	   [[ ! -e $WORKDIR/$INDIR/$TOOLDIR/sampleprops ]] || \
	   [[ ! -e $WORKDIR/$INDIR/$TOOLDIR/svprops.R ]]; then
		printf "\nPlotting flagged but required tools not detected in $WORKDIR/$INDIR/$TOOLDIR\n"
	fi
	mkdir -p $WORKDIR/qc/svprops
fi

# set sv type array
if [[ -z $VARTYPES ]]; then
	VAR_ARRAY=('DEL' 'INS' 'INV' 'DUP' 'BND')
else
	IFS="," read -r -a VAR_ARRAY <<< "$VARTYPES"
fi

# create output dirs
mkdir -p $WORKDIR/$OUTDIR
mkdir -p $WORKDIR/logs

# run job per sv type
for VAR_TYPE in ${VAR_ARRAY[@]}; do
	JOBID=$(cat <<- EOS | qsub -N $NAME.$VAR_TYPE.DELLY -
		#!/bin/bash
		#PBS -l walltime=60:00:00
		#PBS -l select=1:mem=20gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $WORKDIR/logs/$NAME.$VAR_TYPE.delly.log

		# load required modules
		module load bcftools/1.2
		module load anaconda3/personal
		source activate delly

		# copy required inputs to scratch
		cp $WORKDIR/$INDIR/$TUMOR* .
		cp $WORKDIR/$INDIR/$CONTROL* .
		cp -rL $WORKDIR/$RESDIR/$FASTA_BASE* .
		$EXCLUDE_CP

		# create sample file for filtering
		printf "%s\ttumor\n" $TUMOR_NAME > sample_file.tsv
		printf "%s\tcontrol\n" $CNTRL_NAME >> sample_file.tsv
		cat sample_file.tsv
	
		# run delly for each variant type
		delly call -t $VAR_TYPE \
			-g $FASTA \
			-o $NAME.$VAR_TYPE.delly.unfilt.bcf \
			$EXCLUDE_ARG $TUMOR $CONTROL
		delly filter -t $VAR_TYPE \
			-p \
			-f somatic \
			-o $NAME.$VAR_TYPE.delly.bcf \
			-a 0.1 \
			-s sample_file.tsv \
			$NAME.$VAR_TYPE.delly.unfilt.bcf

		cp $NAME.$VAR_TYPE.delly.bcf $WORKDIR/$OUTDIR/
		
		ls -lhRA
		EOS
	)
	DEPEND="$DEPEND,afterok:$JOBID"
done

# remove leading comma from dependency
DEPEND=${DEPEND#*,}

# list bcf
BCF_LIST=$(printf " $NAME.%s.delly.bcf" "${VAR_ARRAY[@]}")
BCF_LIST=${BCF_LIST#* }

#  merge variants
JOBID=$(cat <<- EOS | qsub -N $NAME.merge -
	#!/bin/bash
        #PBS -l walltime=10:00:00
        #PBS -l select=1:mem=5gb:ncpus=1
        #PBS -j oe
        #PBS -q med-bio
        #PBS -o $WORKDIR/logs/$NAME.merge.delly.log
	#PBS -W depend=$DEPEND

        # load required modules
        module load R/3.3.2
        module load bcftools/1.2
        module load anaconda3/personal
        source activate delly

	cp $WORKDIR/$OUTDIR/$NAME* .

	bcftools concat -a -O b -o $NAME.delly.bcf $BCF_LIST

	# plot if flagged
	if [[ "$PLOT" = 1 ]]; then
		$WORKDIR/$TOOLDIR/svprops $NAME.delly.bcf > $NAME.delly.tab
		Rscript $WORKDIR/$TOOLDIR/svprops.R $NAME.delly.tab 
		cp $NAME.delly.tab* $WORKDIR/qc/svprops/
	fi

	ls -lhAR	
	EOS
)
