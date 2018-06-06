#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir='star'
logdir='logs'
fqdir='fastq'
scriptdir='scripts'
trackdir='tracks'

# help message
help_message="

Pipeline to fetch SRA, extract FASTQ and align and quantify reads with STAR2.

usage:
    bash $(basename "$0") [-options] -s <SRA_accession>
required arguments:
    -s|--sra : comma separated list of SRA run accession numbers
    -I|--index : STAR index directory for alignment
    -F|--fasta : FASTA file corresponding to index
    -g|--genome : genome version [hg19,hg38,mm9,mm10]
optional arguments:
    -o|--outdir : output directory for star files (default = star)
    -ld|--logdir : output directory for log files (default = logs)
    -fd|--fqdir : output directory for fastq files (default = fastq)
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
        -I|--index)
            index=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
		-o|--outdir)
            outdir=$2
            shift
            ;;
		-ld|--logdir)
		    logdir=$2
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
if [[ -z ${index:-} ]]; then
    printf "\nERROR: --index argument required.\n"
    echo "$help_message"; exit 1
elif [[ ! -r "${workdir}/${index}" ]]; then
    printf "\nERROR: Input index files cannot be read: %s/%s\n" $workdir $index
    echo "$help_message"; exit 1
else
    index_base=$(basename $index)
fi
if [[ -z ${genome:-} ]]; then
    printf "\nERROR: --genome argument required.\n"
    echo "$help_message"; exit 1
elif [[ -z ${fasta:-} ]]; then
    printf "\nERROR: --fasta argument required.\n"
    echo "$help_message"; exit 1
elif [[ ! -r "${workdir}/${fasta}" ]]; then
    printf "\nERROR: FASTA index files cannot be read: %s/%s\n" $workdir $fasta
    echo "$help_message"; exit 1
fi

# check sample names
if [[ -z ${names_list:-} ]]; then
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
mkdir -p $workdir/$outdir
mkdir -p $workdir/$trackdir

# for each SRA number, submit PBS job
source ~/.bashrc
for idx in ${!sra_array[@]}; do

    unset depend
    unset star_fq1
    unset star_fq2

	sra=${sra_array[$idx]}
	name=${names_array[$idx]}

    if [[ -r $workdir/$outdir/$name.bam ]]; then
        printf "\nBAM already exists: %s/%s/%s\n" $workdir $outdir $name.bam
        continue
    else
        printf "\nProcessing sample: %s, name: %s\n" $sra $name
    fi

	# fetch SRA run info
    if [[ ! -r $logdir/$name.run_info.csv ]]; then
        $(esearch -db sra -query $sra | efetch -format runinfo \
            > $logdir/$name.run_info.csv)
    fi

    n_runs=$(wc -l $workdir/$logdir/$name.run_info.csv | cut -f1 -d' ') 
    if [[ $n_runs < 2 ]]; then
        echo "WARNING: no runs found for sample $sra"
        continue
    fi

    # fetch fastq per run
    while IFS=$'\n' read line; do

	    IFS="," read -r -a srr_info <<< "$line"
        srr=${srr_info[0]}
        libtype=${srr_info[15]}
        platform=${srr_info[18]}

        printf "\tprocessing run: %s\n" $srr

        # get fastq
        if [[ -r $workdir/$fqdir/${srr}_1.fastq.gz ]]; then
            echo "FASTQ already exists for $srr, skipping fetch"
        else
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
        fi

        star_fq1="${star_fq1:-},$fqdir/${srr}_1.fastq.gz"
        star_fq2="${star_fq2:-},$fqdir/${srr}_2.fastq.gz"

    done <<< $(tail -n +2 $workdir/$logdir/$name.run_info.csv) 

    star_fq1=${star_fq1#*,}
    star_fq2=${star_fq2#*,}

    # run STAR
    if [[ -r $workdir/$outdir/$name.bam ]]; then
        echo "BAM file already exists for $name, skipping alignment"
    else
        star_call="bash $scriptdir/align_fastq_star_gdc.sh -fq1 $star_fq1"
        star_call="${star_call} --index $index --fasta $fasta --outdir $outdir"
        star_call="${star_call} --logdir $logdir --check no --name $name"
        if [[ -n ${depend:-} ]]; then star_call="${star_call} --depend $depend"; fi
        if [[ ! $libtype = SINGLE ]]; then star_call="${star_call} -fq2 $star_fq2"; fi
        $star_call &> tmp.log
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: STAR failed for sample: %s\n" $sra
            echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            depend="afterok:$jobid"
        fi
    fi

    # generate track
    track_call="bash $scriptdir/generate_track_deeptools.sh --bam $outdir/$name.bam"
    track_call="${track_call} --fasta $fasta --outdir $trackdir --logdir $logdir"
    track_call="${track_call} --name $name --check no ${genome_arg:-}"
    if [[ -n ${depend:-} ]]; then track_call="${track_call} --depend $depend"; fi
    $track_call &> tmp.log
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: Track generation failed for sample: %s\n" $sra
        echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
    fi
    rm tmp.log

done


