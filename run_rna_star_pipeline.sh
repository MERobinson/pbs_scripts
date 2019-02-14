#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
outdir=$workdir
fqdir='fastq'
qcdir='qc'
logdir='logs'
scriptdir='scripts'
trackdir='tracks'
chimeric='on'
delim=','
fq1_idx=1
fq2_idx=2
sm_idx=3

# help message
help_message="
Run RNA-seq analysis pipeline. Steps: 
    1) align & quantify with STAR
    2) generate alignment metrics
    3) generate genome tracks

usage:
    bash $(basename "$0") [-options] -s <SAMPLE_INFO>
purpose:
    # parse a sample info file and run ChIPseq analysis pipeline for each sample
    # by default, will map reads with BWA MEM, call peaks and generate metrics & tracks
required arguments:
    -i|--sample_info : a delimited text file containing sample information (see below)
    -F|--fasta : path to whole genome fasta file
    -I|--index : path to STAR index to align against
optional arguments:
    -o|--outdir : output directory for STAR files (default = pwd)
    -f|--fqdir : directory containing FASTQ files (default = fastq)
    -l|--logdir : output directory for log files (default = logs)
    -q|--qcdir : output directory for metric files (default = qc)
    -s|--scriptdir : directory containing pipeline scripts (default = scripts)
    -t|--trackdir : output directory for track files (default = tracks)
    -d|--delim : delimiter used in sample info file (default = ',')
    --chimeric : whether to detect chimeric reads [on|off] (default = on)
    --sjdb : a comma sep list of SJ.out.tab files to pass to STAR with
             --sjdbFileChrStartEnd argument (default = NULL)
additional info:
    # sample info file must minimally contain FASTQ names and sample names
column index arguments:
    -f1_idx : index of column containing FASTQ R1 filenames (default = 1)
    -f2_idx : index of column containing FASTQ R2 filenames (default = 2)
    -sm_idx : index of column containing sample names (default = 3)

"

while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
        -i|--sample_info)
            sample_info=$2
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
        -f|--fqdir)
            fqdir=$2
            shift
            ;;
        -l|--logdir)
            logdir=$2
            shift
            ;;
        -q|--qcdir)
            qcdir=$2
            shift
            ;;
        -t|--trackdir)
            trackdir=$2
            shift
            ;;
        --chimeric)
            chimeric=$2
            shift
            ;;
        --sjdb)
            sjdb=$2
            shift
            ;;
        -d|--delim)
            delim=$2
            shift
            ;;
        -f1_idx)
            fq1_idx=$2
            shift
            ;;
        -f2_idx)
            fq2_idx=$2
            shift
            ;;
        -sm_idx)
            sm_idx=$2
            shift
            ;;
        *) 
            printf "\nError: Unrecognised argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
            ;;
    esac
shift
done

# create any output dirs
mkdir -p $outdir
mkdir -p $logdir
mkdir -p $qcdir
mkdir -p $trackdir

# set log file
log_file=$logdir/$(basename $0 .sh).$(date "+%Y-%m-%d").log

# check required arg
if [[ -z "${sample_info:-}" ]]; then
	printf "\nERROR: --sample_info argument is required\n"
	echo "$help_message"; exit 1
elif [[ ! -r "$workdir/$sample_info" ]]; then
	printf "\nERROR: Input sample file is not readable: %s/%s\n" $workdir $sample_info
	echo "$help_message"; exit 1
elif [[ -z "${fasta:-}" ]]; then
    printf "\nERROR: --fasta argument required\n"
    echo "$help_message"; exit 1
elif [[ ! -r "$workdir/$fasta" ]]; then
    printf "\nERROR: FASTA file is not readable: %s/%s\n" $workdir $fasta
    echo "$help_message"; exit 1
elif [[ -z "${index:-}" ]]; then
    printf "\nERROR: --index argument is required\n"
    echo "$help_message"; exit 1
elif [[ ! -r "$workdir/$index" ]]; then
    printf "\nERROR: STAR index file is not readable: %s/%s\n" $workdir $index
    echo "$help_message"; exit 1
fi

# function to extract fields from sample info file
function extract_info {
	awk -F $delim -v v1=$sm_idx -v v2=$sample_name -v v3=$1 \
	'$v1 == v2 {print $v3}' $workdir/$sample_info
}

# function to print elements of array with delim
function join_by { local IFS="$1"; shift; echo "$*"; }

# check if chimeric on
if [[ $chimeric = 'on' ]]; then
    chimeric_arg="--chimeric on"
fi

# check sjdb
if [[ -n ${sjdb:-} ]]; then
    sjdb_arg="--sjdb $sjdb"
fi

# get unique sample names
samples=$(tail -n +2 "$sample_info" | cut -d "$delim" -f "$sm_idx" | sort | uniq)
printf "\nSample list:\n" > $log_file
printf "\t%s\n" ${samples[@]} >> $log_file

# loop through sample names
for sample_name in ${samples[@]}; do
    
    echo "Processing $sample_name"
    printf "\nProcessing %s:\n" $sample_name >> $log_file

    if [[ -r $outdir/$sample_name.star.bam ]]; then
        printf "STAR alignments already exists for %s:\n\t%s/%s.star.bam\n" \
            $sample_name $outdir $sample_name
    else
        # get all FASTQ files for current sample -> array
        fq1_array=($(extract_info $fq1_idx))
        fq2_array=($(extract_info $fq2_idx))

        # check if pairs
        if [[ ${#fq2_array} = ${#fq1_array} ]]; then
            pe=1; printf "\tRunning in PE mode\n" >> $log_file
        else
            printf "\tRunning in SE mode\n" >> $log_file
        fi
    
        # print list of fastq
        printf "\tInput FASTQ files:\n" >> $log_file
        printf "\t\t%s\n" ${fq1_array[@]} >> $log_file
        if [[ ${pe:-} ]]; then printf "\t\t%s\n" ${fq1_array[@]} >> $log_file; fi

        # set fastq arguments by joining arrays into comma sep lists
        fq1_arg=$(printf ",$fqdir/%s" "${fq1_array[@]}")
        fq1_arg=${fq1_arg#*,}
        if [[ ${pe:-} ]]; then 
            fq2_arg=$(printf ",$fqdir/%s" "${fq2_array[@]}")
            fq2_arg="--fastq2 ${fq2_arg#*,}"
        fi

        # assemble command for alignment script
        align_cmd=("bash $scriptdir/align_fastq_star.sh"
                   "--fastq1 $fq1_arg"
                   "--name $sample_name"
                   "--index $index"
                   "--outdir $outdir"
                   "--logdir $logdir"
                   "${fq2_arg:-}"
                   "${chimeric_arg:-}"
                   "${sjdb_arg:-}")
	
	    # run alignment script
        printf "\tAligning with command:\n" $sample_name >> $log_file
        echo ${align_cmd[@]} >> $log_file
	    ${align_cmd[@]} > tmp.log
       
        # check return value and add jobid/filename to merge input lists
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: alignment failed for sample %s. STDERR output:\n" $sample_name
            cat tmp.log; rm tmp.log; exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            depend="--depend afterok:${jobid}"
        fi
    fi
    if [[ -r $qcdir/metrics/${sample_name}.alignment_summary_metrics ]]; then
        printf "Metrics already exist for %s:\n\t%s/metrics/%s.alignment_summary_metrics\n" \
            $sample_name $qcdir $sample_name
    else
        # generate metrics command
        metrics="bash $scriptdir/collect_alignment_metrics.sh -b $outdir/$sample_name.star.bam"
        metrics="${metrics} -F $fasta -n $sample_name --check no ${depend:-}"
        metrics="${metrics} --outdir $qcdir/metrics --logdir $logdir"
    
        # run metrics
        printf "\tCollecting metrics command:\n" >> $log_file
        printf "\t\t%s %s \\ \n" $metrics >> $log_file
        $metrics > tmp.log

        # check return/jobids
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: metrics job failed for sample %s. STDERR output:\n" $sample_name
            cat tmp.log; rm tmp.log; exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            depend="--depend afterok:${jobid}"
        fi
    fi
    if [[ -r $trackdir/$sample_name.bedGraph ]]; then
        printf "Track already exists for %s:\n\t%s/%s.bedGraph\n" $sample_name $trackdir $sample_name
    else
        # generate tracks
        tracks="bash $scriptdir/generate_track_deeptools.sh --bam $outdir/$sample_name.star.bam"
        tracks="${tracks} --outdir $trackdir --logdir $logdir --fasta $fasta ${depend:-} --check no"
        $tracks
    fi 
done
rm tmp.log
