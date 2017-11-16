#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=''
bamdir='bam'
logdir='logs'
qcdir='qc'
fqdir='fastq'
scriptdir='scripts'
trackdir='tracks'

# help message
help_message="

Pipeline to fetch SRA, extract FASTQ, align with BWA MEM and call peaks with MACS2.

usage:
    bash $(basename "$0") [-options] -s <SRA_accession>
required arguments:
    -s|--sra : comma separated list of SRA run accession numbers
    -F|--fasta : whole genome fasta file 
    -I|--index : BWA index file base for alignment
    -g|--genome : genome version [hg19,hg38,mm9,mm10]
optional arguments:
    -bd|--bamdir : output directory for bam files (default = bam)
    -ld|--logdir : output directory for log files (default = logs)
    -fd|--fqdir : output directory for fastq files (default = fastq)
    -qd|--qcdir : output directory for qc files (default = qc)
    -sd|--scriptdir : directory containing scripts required (default = scripts)
    -td|--trackdir : output directory for genome tracks (default = tracks)
    -n|--names : comma separated list of names for output files (default = SRA accession)
    -e|--extend : bp extension for generating coverage track (default = 200)
additional info:
    # all paths should be relative to working directory
    # genome argument is used for macs2 genome size argument and filenames/headers
    # required scripts = fetch_fastq_sra.sh, align_fastq_bwa.sh & generate_tracks_deeptools.sh 
    # either SRA run IDs or GSE sample IDs of the form SRRxxxxxx or GSMxxxxxx can be supplied
    # run info will be retrieved for all runs of a given sample and merged post-alignment

"

# parse command line arg
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-s|--sra)
		    sra_list=$2
		    shift
		    ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -I|--index)
            index=$2
            shift
            ;;
		-o|--outdir)
            outdir=$2
            shift
            ;;
		-bd|--bamdir)
		    bamdir=$2
		    shift
		    ;;
		-ld|--logdir)
		    logdir=$2
		    shift
		    ;;
        -qd|--qcdir)
            qcdir=$2
            shift
            ;;
        -fd|--fqdir)
            fqdir=$2
            shift
            ;;
        -sd|--scriptdir)
            scriptdir=$2
            shift
            ;;
        -td|--trackdir)
            trackdir=$2
            shift
            ;;
		-g|--genome)
		    genome=$2
		    shift
		    ;;
		-n|--names)
		    names_list=$2
		    shift
		    ;;
		-e|--extend)
		    extend=$2
		    shift
		    ;;
		*)
		    printf "\nERROR: Undefined argument provided: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
	esac
	shift
done

# check/parse required arguments
if [[ -z ${sra_list:-} ]]; then
    printf "\nERROR: --sra argument required.\n"
    echo "$help_message"; exit 1
else
    IFS=',' read -r -a sra_array <<< "$sra_list"
fi
if [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required.\n"
    echo "$help_message"; exit 1
elif [[ ! -r $workdir/$fasta ]]; then
    printf "\nERROR: Input FASTA file cannot be read: %s/%s\n" $workdir $fasta
    echo "$help_message"; exit 1
else
    fasta_base=$(basename $fasta)
    fasta_prefix=${fasta%.*}
fi
if [[ -z ${index:-} ]]; then
    printf "\nERROR: --index argument required.\n"
    echo "$help_message"; exit 1
elif [[ ! -r "$workdir/$index.ann" ]]; then
    printf "\nERROR: Input index files cannot be read: %s/%s\n" $workdir $index
    echo "$help_message"; exit 1
else
    index_base=$(basename $index)
fi
if [[ -z ${genome:-} ]]; then
    printf "\nERROR: --genome argument required.\n"
    echo "$help_message"; exit 1
fi

# check sample names
if [[ -z "$names_list" ]]; then
	names_array=(${sra_array[@]})
    echo "WARNING: No name argument given, using accession numbers."
else
	IFS=',' read -r -a names_array <<< "$names_list"
fi

# check genome
if [[ -n ${genome:-} ]]; then
    genome_arg="--genome $genome"
fi

# create required dirs
mkdir -p $workdir/$logdir
mkdir -p $workdir/$qcdir/fastqc
mkdir -p $workdir/$qcdir/metrics
mkdir -p $workdir/$bamdir
mkdir -p $workdir/$trackdir

# for each SRA number, submit PBS job
source ~/.bashrc
for idx in ${!sra_array[@]}; do

	sra=${sra_array[$idx]}
	name=${names_array[$idx]}

    printf "\nSRA: %s; NAME: %s\n" $sra $name

	# fetch SRA run info
    $(esearch -db sra -query $sra | efetch -format runinfo \
         > $workdir/$logdir/$name.run_info.csv)

    # fetch and align per run
    while IFS=$'\n' read line; do

        if [[ ! $line = SRR* ]]; then continue; fi

	    IFS="," read -r -a srr_info <<< "$line"
        srr=${srr_info[0]}
        libtype=${srr_info[15]}
        platform=${srr_info[18]}
	    if [[ ! $libtype = SINGLE ]]; then
		    fq2_arg="-fq2 $fqdir/${srr}_2.fastq.gz"
	    fi

        # get fastq
        sra_call="bash $scriptdir/fetch_fastq_sra.sh --sra $srr"
        sra_call="${sra_call} --outdir $fqdir --logdir $logdir"
        $sra_call &> tmp.log
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: FASTQ fetching failed for sample: %s\n" $srr
            echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            depend="afterok:$jobid"
        fi

        # align
        bwa_call="bash $scriptdir/align_fastq_bwa.sh -fq1 $fqdir/${srr}_1.fastq.gz ${fq2_arg:-}"
        bwa_call="${bwa_call} --fasta $fasta --index $index --outdir $bamdir --qcdir $qcdir"
        bwa_call="${bwa_call} --logdir $logdir --depend $depend --check no --name $srr"
        $bwa_call &> tmp.log
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: Alignment failed for sample: %s\n" $sra
            echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            merge_depend="${merge_depend:-},afterok:$jobid"
            merge_arg="${merge_arg:-},$bamdir/$srr.bam"
        fi
    done < $workdir/$logdir/$name.run_info.csv

    # merge runs
    merge_depend=${merge_depend#*,}; merge_arg=${merge_arg#*,}
    merge_call="bash $scriptdir/merge_bam_picard.sh --bam_list $merge_arg --fasta $fasta"
    merge_call="${merge_call} --name $name --bamdir $bamdir --qcdir $qcdir --logdir $logdir"
    merge_call="${merge_call} --depend $merge_depend --check no"
    $merge_call &> tmp.log
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: Merging failed for sample: $name\n"
        echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
    else
        jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
        depend="afterok:$jobid"
    fi

    # generate track
    track_call="bash $scriptdir/generate_track_deeptools.sh --bam $bamdir/$name.bam"
    track_call="${track_call} --outdir $trackdir --fasta $fasta --logdir $logdir"
    track_call="${track_call} --name $name --depend $depend --check no ${genome_arg:-}"
    $track_call &> tmp.log
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: Track generation failed for sample: %s\n" $sra
        echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
    fi
    rm tmp.log
    
done

