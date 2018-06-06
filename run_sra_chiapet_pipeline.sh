#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir='chiapet'
fqdir='fastq'
logdir='logs'
scriptdir='scripts'

# help message
help_message="

Pipeline to fetch SRA, align and call interactions with ChIA-PET2.

usage:
    bash $(basename "$0") [-options] -s <SRA_accession>
required arguments:
    -s|--sra : comma separated list of SRA run accession numbers
    -C|--chrom : chromosome sizes file
    -I|--index : BWA index file base 
optional arguments:
    -o|--outdir : output directory for ChIA-PET2 output (default = chiapet)
    -fd|--fqdir : output directory for fastq files (default = fastq)
    -ld|--logdir : output directory for log files (default = logs)
    -sd|--scriptdir : directory containing scripts required (default = scripts)
    -n|--names : comma separated list of names for output files (default = SRA accession)
additional info:
    # all paths should be relative to working directory
    # required scripts = fetch_fastq_sra.sh & call_interactions_chiapet2.sh 
    # either SRA run IDs or GSE sample IDs of the form SRRxxxxxx or GSMxxxxxx can be supplied

"

# parse command line arg
while [[ $# -gt 1 ]]; do
	key=$1
	case $key in
		-s|--sra)
		    sra_list=$2
		    shift
		    ;;
        -C|--chrom)
            chrom_sizes=$2
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
        -fd|--fqdir)
            fqdir=$2
            shift
            ;;
		-ld|--logdir)
		    logdir=$2
		    shift
		    ;;
        -sd|--scriptdir)
            scriptdir=$2
            shift
            ;;
		-n|--names)
		    names_list=$2
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
if [[ -z ${chrom_sizes:-} ]]; then
    printf "\nERROR: --chrom argument required.\n"
    echo "$help_message"; exit 1
elif [[ ! -r $chrom_sizes ]]; then
    printf "\nERROR: Input chromosome sizes file cannot be read: %s/%s\n" $workdir $chrom_sizes
    echo "$help_message"; exit 1
elif [[ -z ${index:-} ]]; then
    printf "\nERROR: --index argument required.\n"
    echo "$help_message"; exit 1
elif [[ ! -e $index ]]; then
    printf "\nERROR: Input index files cannot be read: %s/%s\n" $workdir $index
    echo "$help_message"; exit 1
fi

# check sample names
if [[ -z ${names_list:-} ]]; then
	names_array=(${sra_array[@]})
    echo "WARNING: No name argument given, using accession numbers."
else
	IFS=',' read -r -a names_array <<< "$names_list"
fi

# create required dirs
mkdir -p $workdir/$outdir
mkdir -p $workdir/$logdir

# for each SRA number, submit fetch job
source ~/.bashrc
for idx in ${!sra_array[@]}; do

    unset depend
    unset fq1_arg
    unset fq2_arg

	sra=${sra_array[$idx]}
	name=${names_array[$idx]}

	# fetch SRA run info
    if [[ ! -e $workdir/$logdir/$name.run_info.csv ]]; then
        $(esearch -db sra -query $sra | efetch -format runinfo \
          > $workdir/$logdir/$name.run_info.csv)
    fi

    n_runs=$(wc -l $workdir/$logdir/$name.run_info.csv | cut -f1 -d' ') 
    if [[ $n_runs < 2 ]]; then
        echo "WARNING: no runs found for sample $sra"
        continue
    fi

    # fetch and align per run
    while IFS='' read line; do

        # skip header line
        if [[ ! $line == SRR* ]]; then
            continue
        fi
        
        # split line and extract info
	    IFS="," read -r -a srr_info <<< "$line"
        srr=${srr_info[0]}
        libtype=${srr_info[15]}
        platform=${srr_info[18]}
	    if [[ $libtype = SINGLE ]]; then
		    printf "\nERROR: libtype for %s = $libtype, only paired end supported" $srr $libtype
            continue
	    else
            printf "\tprocessing run: %s\n" $srr
        fi

        # get fastq
        if [[ ! -e $workdir/$fqdir/${srr}_1.fastq.gz ]]; then
            sra_call=("bash $scriptdir/fetch_fastq_sra.sh --sra $srr"
                      "--outdir $fqdir --logdir $logdir")
            ${sra_call[@]} &> tmp.log

            # check return/get jobid
            if [[ $? -ne 0 ]]; then
                printf "\nERROR: FASTQ fetching failed for sample: %s\n" $srr
                echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
            else
                jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
                depend="${depend:-},afterok:$jobid"
            fi
        fi

        fq1_arg="${fq1_arg:-},$fqdir/${srr}_1.fastq.gz"
        fq2_arg="${fq2_arg:-},$fqdir/${srr}_2.fastq.gz"
        
    done < $workdir/$logdir/$name.run_info.csv

    # trim leading commas
    if [[ -n ${depend:-} ]]; then depend=${depend#*,}; fi
    fq1_arg=${fq1_arg#*,}; fq2_arg=${fq2_arg#*,}

    # run chiapet2
    chia_call=("bash $scriptdir/call_interactions_chiapet2.sh"
               "-fq1 $fq1_arg -fq2 $fq2_arg"
               "--linkerA GTTGGATAAGATATCGC"
               "--linkerB GTTGGAATGTATATCGC"
               "--chrom $chrom_sizes --index $index"
               "--outdir $outdir --logdir $logdir"
               "--name $name --short 1")
    if [[ -n ${depend:-} ]]; then 
        echo $depend
        chia_call="${chia_call[@]} --depend $depend --check no"
    fi
    ${chia_call[@]} &> tmp.log

    # check return
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: ChIA-PET2 call failed for sample: %s\n" $sra
        echo "stderr: "; cat tmp.log; rm tmp.log; exit 1
    else
        rm tmp.log
    fi

done
