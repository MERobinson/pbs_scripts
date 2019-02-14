#!/usr/bin/env bash
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
logdir=logs
vardir=variants
scriptdir=scripts
tmpdir=tmp
delim=","
log_file=$(basename $0 .sh)_$(date "+%Y-%m-%d").log

# sample info index defaults
sm_idx=1
bam_idx=2

# help message
help_message="
usage:
    bash $(basename "$0") -s <SAMPLE_INFO>
purpose:
    Run germline SV calling with Delly2
required arguments:
    -i|--info : a delimited text file containing sample information (see below)
    -F|--fasta : whole genome FASTA file
    -E|--exclude : regions to exclude from SV calling (see Delly2 man for details)
optional arguments:
    -l|--logdir : output directory for log files (default = logs)
    -v|--vardir : output directory for variant files (default = variants)
    -s|--scriptdir : input directory containing pipeline scripts (default = scripts)
    -t|--tmpdir : temporary directory for intermediate files (default = tmp)
    --delim : delimiter used in sample info file (default = ',')
    --sm_idx : column index (1-based) containing sample names [integer] (default = '1')
    --bam_idx: column index (1-based) containing BAM path [integer] (default = '2')
additional info:
    # sample info file must contain sample names and the relative path to mapped BAM files
    # by default expects a comma separated file with sample names in first column
      and BAM file in second column but these can be altered with --*idx and --delim args
    # scriptdir must contain scripts called by the pipeline: 
        > call_sv_delly.sh
        > merge_sv_delly.sh
        > filter_sv_delly.sh
    # all paths should be given relative to current working directory
outputs:
    # outputs called structural variants to --vardir named <sample_name>.<sv_type>.delly.bcf
    # outputs run logs to --logdir
    # intermediate files are stored in tmpdir

"
# parse arguments
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--info)
            sample_info=$2
            shift
            ;;
        -F|--fasta)
            fasta=$2
            shift
            ;;
        -E|--exclude)
            exclude=$2
            shift
            ;;
        -l|--logdir)
            logdir=$2
            shift
            ;;
        -v|--vardir)
            vardir=$2
            shift
            ;;
        -s|--scriptdir)
            scriptdir=$2
            shift
            ;;
        -t|--tmpdir)
            tmpdir=$2
            shift
            ;;
        -d|--delim)
            delim=$2
            shift
            ;;
        --sm_idx)
            fq1_idx=$2
            shift
            ;;
        --bam_idx)
            fq2_idx=$2
            shift
            ;;
        *)
            printf "\nERROR: Illegal argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
        ;;
    esac
shift
done

# inputs
if [[ -z "${sample_info:-}" ]]; then
    printf "\nERROR: Sample info file not provided.\n"
    echo "$help_message"; exit 1
elif [[ -z "${fasta:-}" ]]; then
    printf "\nERROR: FASTA file not provided.\n"
    echo "$help_message"; exit 1
elif [[ ! -r "$workdir/$sample_info" ]]; then
    printf "\nERROR: Input sample file not readable: %s/%s\n" $workdir $sample_info
    echo "$help_message"; exit 1
elif [[ ! -d "$workdir/$scriptdir" ]]; then
    printf "\nERROR: Script directory doesnt exist: %s/%s\n" $workdir $scriptdir
    echo "$help_message"; exit 1
elif [[ ! -r "$workdir/$fasta" ]]; then
    printf "\nERROR: FASTA file doesnt not readable: %s/%s\n" $workdir $fasta
    echo "$help_message"; exit 1
elif [[ ! -r "$workdir/$scriptdir/call_sv_delly.sh" ]]; then
    printf "\nERROR: script not found: %s/%s/call_sv_delly.sh\n" $workdir $scriptdir
    echo "$help_message"; exit 1
elif [[ ! -r "$workdir/$scriptdir/merge_sv_delly.sh" ]]; then
    printf "\nERROR: script not found: %s/%s/merge_sv_delly.sh\n" $workdir $scriptdir
    echo "$help_message"; exit 1
elif [[ ! -r "$workdir/$scriptdir/filter_sv_delly.sh" ]]; then
    printf "\nERROR: script not found: %s/%s/filter_sv_delly.sh\n" $workdir $scriptdir
    echo "$help_message"; exit 1
fi

# function to extract fields from sample info file
function extract_info {
    awk -F $delim -v v1=$sm_idx -v v2=$sample_name -v v3=$1 \
    '$v1 == v2 {print $v3}' $workdir/$sample_info
}

# get unique sample names
samples=($(cat "${sample_info}" | cut -d "${delim}" -f "${sm_idx}" | sort | uniq))
printf "Input sample list:\n" > $logdir/$log_file
printf "\t%s\n" ${samples[@]} >> $logdir/$log_file

# create temporary directory
mkdir -p $workdir/$tmpdir

# run separately for each SV type
echo "Calling germline structural variants"
cat_arg=''
cat_depend=''
for vartype in INS DEL INV BND DUP; do

    # discover var
    echo "Calling ${vartype}..."
    merge_arg=''
    merge_depend=''
    for sample_name in ${samples[@]}; do

        # get bam file path
        bam=$(extract_info "${bam_idx}")

        # assemble command
        delly_call="bash $scriptdir/call_sv_delly.sh --test $bam --vartype $vartype"
        delly_call="${delly_call} --outdir $tmpdir --logdir $logdir"
        delly_call="${delly_call} --name $sample_name.$vartype.c1"
        delly_call="${delly_call} --fasta $fasta --exclude $exclude --check no"

        # run delly to discover sv
        printf "\tCalling %s variants for sample %s with command:\n" \
                $vartype $sample_name >> $logdir/$log_file
        printf "\t\t%s %s \\ \n" $delly_call >> $logdir/$log_file
        $delly_call > $tmpdir/tmp.log

        # check return value
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: call_sv_delly failed for sample: %s\n" $sample_name
            echo "stderr: "; cat $tmpdir/tmp.log
            exit 1
        else
            jobid=$(cat $tmpdir/tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            merge_arg=${merge_arg}",$tmpdir/$sample_name.$vartype.c1.bcf"
            merge_depend=${merge_depend}",afterok:${jobid}"
        fi
    done
    merge_arg=${merge_arg#*,}
    merge_depend=${merge_depend#*,}

    # merge var accross samples
    merge_call="bash $scriptdir/merge_sv_delly.sh --sv_list $merge_arg --vartype $vartype"
    merge_call="${merge_call} --outdir $tmpdir --logdir $logdir --name germline_sv.$vartype.m1"
    merge_call="${merge_call} --depend $merge_depend --check no"     
    printf "\tMerging %s variants with command:\n" \
           $vartype >> $logdir/$log_file
    printf "\t\t%s %s \\ \n" $merge_call >> $logdir/$log_file
    $merge_call > $tmpdir/tmp.log

    # check return value
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: merge_sv_delly failed for variant: %s\n" $vartype
        echo "stderr: "; cat $tmpdir/tmp.log
        exit 1
    else
        jobid=$(cat $tmpdir/tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
        recall_depend="afterok:${jobid}"
    fi

    # joint genotyping of sv
    merge_arg=''
    merge_depend=''
    for sample_name in ${samples[@]}; do
        
        # get bam file path
        bam=$(extract_info "${bam_idx}")

        # re-genotype each sample with merged sv file
        delly_call="bash $scriptdir/call_sv_delly.sh --test $bam --vartype $vartype"
        delly_call="${delly_call} --outdir $tmpdir --logdir $logdir" 
        delly_call="${delly_call} --name $sample_name.$vartype.c2"
        delly_call="${delly_call} --fasta $fasta --exclude $exclude"
        delly_call="${delly_call} --depend $recall_depend --check no"
        delly_call="${delly_call} --vcf $tmpdir/germline_sv.$vartype.m1.bcf"
        printf "\tRe-calling %s variants for sample %s with command:\n" \
                $vartype $sample_name >> $logdir/$log_file
        printf "\t\t%s %s \\ \n" $delly_call >> $logdir/$log_file
        $delly_call > $tmpdir/tmp.log

        # check return value
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: call_sv_delly failed for sample: %s\n" $sample_name
            echo "stderr: "; cat $tmpdir/tmp.log
            exit 1
        else
            jobid=$(cat $tmpdir/tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            merge_arg="${merge_arg},$tmpdir/$sample_name.$vartype.c2.bcf"
            merge_depend="${merge_depend},afterok:${jobid}"
        fi
    done
    merge_arg=${merge_arg#*,}
    merge_depend=${merge_depend#*,}

    # merge re-called sv accross samples
    merge_call="bash $scriptdir/merge_sv_delly.sh --sv_list $merge_arg --vartype $vartype"
    merge_call="${merge_call} --outdir $tmpdir --logdir $logdir --name germline_sv.$vartype.m2"
    merge_call="${merge_call} --depend $merge_depend --check no" 
    printf "\tMerging %s variants with command:\n" $vartype >> $logdir/$log_file
    printf "\t\t%s %s \\ \n" $merge_call >> $logdir/$log_file
    $merge_call > $tmpdir/tmp.log

    # check return value
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: merge_sv_delly failed for SV: %s\n" $vartype
        echo "stderr: "; cat $tmpdir/tmp.log
        exit 1
    else
        jobid=$(cat $tmpdir/tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
        filter_depend="afterok:${jobid}"
    fi

    # filter final call set
    filter_call="bash $scriptdir/filter_sv_delly.sh --vartype $vartype" 
    filter_call="${filter_call} --bcf $tmpdir/germline_sv.$vartype.m2.bcf"
    filter_call="${filter_call} --name germline_sv.$vartype.delly"
    filter_call="${filter_call} --depend $filter_depend --check no"
    filter_call="${filter_call} --outdir $vardir --logdir $logdir" 
    printf "\tFiltering %s callset with command:\n" $vartype >> $logdir/$log_file
    printf "\t\t%s %s \\ \n" $filter_call >> $logdir/$log_file
    $filter_call > $tmpdir/tmp.log

done
