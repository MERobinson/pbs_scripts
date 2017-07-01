#!/bin/bash

# default arg
workdir=$PWD
outdir=variants
bamdir=bam
resdir=resources
scriptdir=scripts
logdir=logs
fasta=genome.fa
exclude=exclude.tsv

# help message
help_message="

usage:
    bash $(basename "$0") [-wborsfenv] -t <BAM> -c <BAM>
purpose:
    wrapper script to call structural variants with Delly
required arguments:
    -t|--tumor : tumor sample (aligned bam)
    -c|--control : control sample (aligned bam)
optional arguments:
    -w|--workdir : path of working directory to use (default = pwd)
    -b|--bamdir : directory within workdir containing input bam files (default = bam)
    -o|--outdir : directory to output variant BCF files (default = variants)
    -r|--resdir : directory within workdir containing resource files (default = resources)
    -s|--scriptdir: directory within workdir containing scripts (default = scripts)
    -l|--logdir : directory to output log files (default = logs)
    -f|--fasta : name of whole genome fasta file within resdir (default = genome.fa)
    -e|--excl : name of exclusion file within resdir (default = exclude.tsv) 
    -n|--name : name prefix for output files (deafult = tumor BAM filename)
    -v|--vartypes : comma sep. list of variant types to analyse
                    [INS,DEL,INV,BND,DUP] (default = all)
    --plot : whether to summarise and plot SVs with svprops (see additional info)
example:
    bash $(basename "$0") -t patient1_tumor.bam -c patient1_control.bam
inputs:
    > Tumor and control BAM files (PE, aligned BAM) [required]
    > Whole genome fasta [required]
    > Regions to exclude [optional]
additional info:
    > outdir will be created in workdir if not pre-existing
    > name is extracted from test filename up to first period (.)
    > if wanting to run with no region exclusion set --excl NULL
    > --plot requires installation of svprops:
      1) clone and install from https://github.com/dellytools/svprops
      2) link svprops sampleprops and svprops.R to --scriptdir
    > see Delly man for additional info: https://github.com/dellytools/delly

"

# parse command line arg
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -t|--tumor)
        tumor=$2
        shift
        ;;
        -c|--control)
        control=$2
        shift
        ;;
        -w|--workdir)
        workdir=$2
        shift
        ;;
        -b|--bamdir)
        bamdir=$2
        shift
        ;;
        -o|--outdir)
        outdir=$2
        shift
        ;;
        -r|--resdir)
        resdir=$2
        shift
        ;;
        -s|--scriptdir)
        scriptdir=$2
        shift
        ;;
        -l|--logdir)
        logdir=$2
        shift
        ;;
        -f|--fasta)
        fasta=$2
        shift
        ;;
        -e|--excl)
        exclude=$2
        shift
        ;;
        -n|--name)
        name=$2
        shift
        ;;
        -v|--vartypes)
        vartypes=$2
        shift
        ;;
        --plot)
        plot=1
        shift
        ;;
        *)
        printf "\nERROR: Illegal argument\n"
        echo "$help_message"; exit 1
        ;;
    esac
    shift
done

# check dir
if [[ ! -d $workdir/$bamdir ]]; then
    printf "\nERROR: Input directory not found: $workdir/$bamdir \n"
    echo "$help_message"; exit 1
elif [[ ! -d $workdir/$resdir ]]; then
    printf "\nERROR: Resource directory not found: $workdir/$resdir \n"
    echo "$help_message"; exit 1
fi

# check required arguments
if [[ -z "$tumor" ]]; then
    printf "\nERROR: No tumor sample provided\n"
    echo "$help_message"; exit 1
elif [[ ! -e $workdir/$bamdir/$tumor ]]; then
    printf "\nERROR: Tumor file not found: $workdir/$bamdir/$tumor \n"
    echo "$help_message"; exit 1
elif [[ -z "$control" ]]; then
    printf "\nERROR: No control sample provided\n"
    echo "$help_message"; exit 1
elif [[ ! -e $workdir/$bamdir/$control ]]; then
    printf "\nERROR: Control file not found: $workdir/$bamdir/$control \n"
    echo "$help_message"; exit 1
elif [[ ! -e $workdir/$resdir/$fasta ]]; then
    printf "\nERROR: fasta file not found: $workdir/$resdir/$fasta \n"
    echo "$help_message"; exit 1  
fi

# get name if req
if [[ -z "$name" ]]; then
    name=${tumor%%.*}
fi
tumor_name=${tumor%%.*}
control_name=${control%%.*}

# split fasta filename
fasta_base=${fasta%%.*}

# set exclusion argument if non-null
if [[ ! $exclude = "NULL" ]]; then
	exclude_arg="-x $exclude"
	exclude_cp="cp -rL $workdir/$resdir/$exclude ."
fi

# check tools if plotting flagged
if [[ $plot = 1 ]]; then
    if [[ ! -e $workdir/$scriptdir/svprops ]] || \ 
       [[ ! -e $workdir/$scriptdir/sampleprops ]] || \
       [[ ! -e $workdir/$scriptdir/svprops.R ]]; then
        printf "\nPlotting flagged but required tools not detected in %s/%s\n" $workdir $scriptdir
        echo "$help_message"; exit 1
    fi
    mkdir -p $workdir/qc/svprops
fi

# set sv type array
if [[ -z $vartypes ]]; then
	var_array=('DEL' 'INS' 'INV' 'DUP' 'BND')
else
	IFS="," read -r -a var_array <<< "$vartypes"
fi

# create output dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# run job per sv type
for var_type in ${var_array[@]}; do
	jobid=$(cat <<- EOS | qsub -N $name.$var_type.delly -
		#!/bin/bash
		#PBS -l walltime=60:00:00
		#PBS -l select=1:mem=50gb:ncpus=1
		#PBS -j oe
		#PBS -q med-bio
		#PBS -o $workdir/logs/$name.$var_type.delly.log

		# load required modules
		module load bcftools/1.2
		module load anaconda3/personal
		source activate delly

		# copy required inputs to scratch
		cp $workdir/$bamdir/$tumor* .
		cp $workdir/$bamdir/$control* .
		cp -rL $workdir/$resdir/$fasta_base* .
		$exclude_cp

		# create sample file for filtering
		printf "%s\ttumor\n" $tumor_name > sample_file.tsv
		printf "%s\tcontrol\n" $control_name >> sample_file.tsv
		cat sample_file.tsv
	
		# run delly for each variant type
		delly call -t $var_type \
			-g $fasta \
			-o $name.$var_type.delly.unfilt.bcf \
			$exclude_arg $tumor $control
		delly filter -t $var_type \
			-p \
			-f somatic \
			-o $name.$var_type.delly.bcf \
			-a 0.1 \
			-s sample_file.tsv \
			$name.$var_type.delly.unfilt.bcf

		cp $name.$var_type.delly.bcf* $workdir/$outdir/
		
		ls -lhRA
		EOS
	)
	depend="$depend,afterok:$jobid"
done

# remove leading comma from dependency
depend=${depend#*,}

# list bcf
bcf_list=$(printf " $name.%s.delly.bcf" "${var_array[@]}")
bcf_list=${bcf_list#* }

#  merge variants
jobid=$(cat <<- EOS | qsub -N $name.merge -
	#!/bin/bash
	#PBS -l walltime=10:00:00
	#PBS -l select=1:mem=10gb:ncpus=1
	#PBS -j oe
	#PBS -q med-bio
	#PBS -o $workdir/logs/$name.merge.delly.log
	#PBS -W depend=$depend

	# load required modules
	module load R/3.3.2
	module load bcftools/1.2
	module load anaconda3/personal
	source activate delly

	cp $workdir/$outdir/$name* .

	bcftools concat -a -O b -o $name.delly.bcf $bcf_list

	# plot if flagged
	if [[ "$plot" = 1 ]]; then
		$workdir/$scriptdir/svprops $name.delly.bcf > $name.delly.tab
		Rscript $workdir/$scriptdir/svprops.R $name.delly.tab 
		cp $name.delly.tab* $workdir/qc/svprops/
	fi

	ls -lhAR	
	EOS
)
