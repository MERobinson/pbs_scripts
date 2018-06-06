#!/usr/bin/env bash
set -o pipefail

# default arg
workdir=$PWD
fqdir=fastq
resdir=resources
logdir=logs
scriptdir=scripts
delim=","
genome_size=hs
fq1_idx=1
fq2_idx=2
sm_idx=3
fasta=genome.fa
index=genome.fa
log_file=$(basename $0 .sh)_$(date "+%Y-%m-%d").log
extract='no'

# help message
help_message="
usage:
    bash $(basename "$0") [-options] -i <sample_info.csv>
purpose:
    # parse a sample info file and run ChIPseq analysis pipeline for each sample
    # by default, will map reads with BWA MEM, call peaks and generate metrics & tracks
required arguments:
    -i|--sampleinfo : a delimited text file containing sample information (see below)
optional arguments:
    -w|--workdir : working directory (default = pwd)
    -f|--fqdir : directory within workdir containing FASTQ files (default = fastq)
    -r|--resdir : directory within workdir containing resource files (default = resources)
    -l|--logdir : directory within workdir to output log files (default = logs)
    -s|--scriptdir : directory within workdir containing pipeline scripts (default = scripts)
    -d|--delim : delimiter used in sample info file (default = ',')
    -g|--genome_size : genome size argument passed to macs (default = hs)
    --extract : whether to extract ID and PU from read name [yes,no] (default = no)
additional info:
    # sample info file must minimally contain FASTQ names and sample names
    # additional columns may be included as listed below:
        - Filename of FASTQ R1 [required]
        - Filename of FASTQ R2 [required - leave blank column for SE]
        - Sample name [required]
        - Control sample name [optional]
        - Read group ID [optional]
        - Platform unit [optional]
        - Library name [optional]
        - Sequencing centre [optional]
        - Platform [default = ILLUMINA]
        - Platform model [optional]
        - Run date [optional]
        - Predicted median insert size [optional] 
    # to include any of the above - provide column indexes (see column index arguments)
    # read group ID and PU can be extracted from read names with '--extract yes' argument 
      if read names are standard illumina format
    # NB: peak calling is ONLY conducted if control sample index (cs_idx) is provided,
          and corresponding field is not empty.
    # following resource files must be present in resdir:
        - Whole genome sequence [FASTA and corresponding dictionary file]
        - BWA index files [all index components required]
    # to alter expected resource filenames, see resource filename arguments
column index arguments:
    -f1_idx : index of column containing FASTQ R1 filenames (default = 1)
    -f2_idx : index of column containing FASTQ R2 filenames (default = 2)
    -sm_idx : index of column containing sample names (default = 3)
    -cs_idx : index of column containing sample name of control (default = null) 
    -id_idx : index of column containing read group ID (default = null)
    -pu_idx : index of column containing platform unit (default = null)
    -lb_idx : index of column containing library name (default = null)
    -cn_idx : index of column containing sequencing centre (default = null)
    -pl_idx : index of column containing platform (default = null)
    -pi_idx : index of column containing median insert size (deafult = null)
    -pm_idx : index of column containing platform model (default = null)
    -dt_idx : index of column containing run date (default = null)
resource filename arguments:
    -ix|--index_base : basename of index files in resdir (default = genome.fa)
    -fa|--fasta : name of FASTA file in resdir (default = genome.fa) 

"

while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
        -i|--sampleinfo)
            sample_info=$2
            shift
            ;;
        -w|--workdir)
            workdir=$2
            shift
            ;;
        -f|--fqdir)
            fqdir=$2
            shift
            ;;
        -r|--resdir)
            resdir=$2
            shift
            ;;
        -l|--logdir)
            logdir=$2
            shift
            ;;
        -d|--delim)
            delim=$2
            shift
            ;;
        -g|--genome_size)
            genome_size=$2
            shift
            ;;
        --extract)
            extract=$2
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
        -cs_idx)
            cs_idx=$2
            shift
            ;;
        -lb_idx) 
            lb_idx=$2
            shift
            ;;
        -pu_idx)
            pu_idx=$2
            shift
            ;;
        -pl_idx)
            pl_idx=$2
            shift
            ;;
        -cn_idx) 
            cn_idx=$2
            shift
            ;;
        -pi_idx)
            pi_idx=$2
            shift
            ;;
        -pm_idx)
            pm_idx=$2
            shift
            ;;
        -dt_idx)
            dt_idx=$2
            shift
            ;;
        -ix|--index_base)
            index=$2
            shift
            ;;
        -fa|--fasta)
            fasta=$2
            shift
            ;;
        --force)
            force=$2
            shift
            ;;
        *) 
            printf "\nError: Unrecognised argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
        ;;
    esac
shift
done

# check required arg
if [[ -z "${sample_info:-}" ]]; then
	printf "\nERROR: Sample info file not provided.\n"
	echo "$help_message"; exit 1
fi
if [[ ! -r "$workdir/$sample_info" ]]; then
	printf "\nERROR: Input sample file doesnt exist/isnt readable: %s/%s\n" $workdir $sample_info
	echo "$help_message"; exit 1
fi
if [[ ! -d "$workdir/$resdir" ]]; then
	printf "\nERROR: Resources directory doesnt exist: %s/%s\n" $workdir $resdir
	echo "$help_message"; exit 1
fi
if [[ !	-d "$workdir/$fqdir" ]]; then
    printf "\nERROR: FASTQ directory doesnt exist: %s/%s\n" $workdir $fqdir
    echo "$help_message"; exit 1
fi
if [[ !	-e "$workdir/$resdir/$index.bwt" ]]; then
    printf "\nERROR: Index file doesnt exist/isnt readable: %s/%s/%s\n" $workdir $resdir $index
    echo "$help_message"; exit 1
fi
if [[ ! -r "$workdir/$resdir/$fasta" ]]; then
    printf "\nERROR: FASTA file doesnt exist/isnt readable: %s/%s/%s\n" $workdir $resdir $fasta
    echo "$help_message"; exit 1
fi

mkdir -p $workdir/logs

# function to extract fields from sample info file
function extract_info {
	awk -F $delim -v v1=$sm_idx -v v2=$sample_name -v v3=$1 \
	    '$v1 == v2 {print $v3}' $workdir/$sample_info
}

# function to list dependencies
function get_dependencies {
    job_name=$(echo $sample_name | cut -c1-15)
    depend_list=$(qstat -wu $USER | grep $job_name | cut -d. -f1 | tr '\n' ',')
    depend_list=${depend_list//,/,afterok:}
    depend_list=${depend_list%,*}
    echo "afterok:$depend_list"
}

# get unique sample names
samples=$(tail -n +2 "$sample_info" | cut -d "$delim" -f "$sm_idx" | sort | uniq)
printf "\nSample list:\n" > $workdir/$logdir/$log_file
printf "\t%s\n" ${samples[@]} >> $workdir/$logdir/$log_file

# loop through sample names
for sample_name in ${samples[@]}; do

    echo "Processing $sample_name"
    printf "\nProcessing %s:\n" $sample_name >> $workdir/$logdir/$log_file

    # get all FASTQ files for current sample -> array
    fq1_array=($(extract_info $fq1_idx))

    # do same for optional header info
    if [[ -n ${fq2_idx:-} ]]; then fq2_array=($(extract_info $fq2_idx)); fi
    if [[ -n ${id_idx:-} ]]; then id_array=($(extract_info $id_idx)); fi
    if [[ -n ${cn_idx:-} ]]; then cn_array=($(extract_info $cn_idx)); fi    
    if [[ -n ${dt_idx:-} ]]; then dt_array=($(extract_info $dt_idx)); fi
    if [[ -n ${pu_idx:-} ]]; then pu_array=($(extract_info $pu_idx)); fi
    if [[ -n ${lb_idx:-} ]]; then lb_array=($(extract_info $lb_idx)); fi
    if [[ -n ${pl_idx:-} ]]; then pl_array=($(extract_info $pl_idx)); fi
    if [[ -n ${pi_idx:-} ]]; then pi_array=($(extract_info $pi_idx)); fi
    if [[ -n ${pm_idx:-} ]]; then pm_array=($(extract_info $pm_idx)); fi

    printf "\tNumber of runs/FASTQ detected: %s\n" ${#fq1_array[@]} >> $logdir/$log_file

    # loop through each FASTQ1
    unset merge
    unset depend
    for idx in ${!fq1_array[@]}; do

        # get basename
        fq_name=$(basename ${fq1_array[$idx]})
        fq_name=${fq_name%%.*}

        # skip if exists
        if [[ -e bam/$fq_name.bam ]]; then
            merge="${merge:-},bam/$fq_name.bam"
            continue
        fi

        # if flagged, extract info from read names
        if [[ $extract_info = 'yes' ]]; then
            readname=$(gzip -dc $fqdir/${fq1_array[$idx]} | head -n 1)
            fcid=$(echo $readname | cut -d":" -f3) # flow cell id
            fcln=$(echo $readname | cut -d":" -f4) # flow cell lane
            bc=$(echo $readname | cut -d":" -f10) # barcode
            rgid=$fcid.$fcln # read group ID
            pu=$fcid.$fcln.$bc # platform unit
        fi

        # assemble command for alignment script with any provided optional args
        align=("bash $scriptdir/align_fastq_bwa.sh" 
               "-fq1 $fqdir/${fq1_array[$idx]} -n $fq_name --sm $sample_name"
               "-F $resdir/$fasta -I $resdir/$index -o bam -l logs -q qc")

        if [[ -n ${fq2_array[$idx]:-} ]]; then align+=("-fq2 $fqdir/${fq2_array[$idx]}"); fi 
        if [[ -n ${lb_array[$idx]:-} ]]; then align+=("--lb ${lb_array[$idx]}"); fi
       	if [[ -n ${cn_array[$idx]:-} ]]; then align+=("--cn ${cn_array[$idx]}"); fi
       	if [[ -n ${pl_array[$idx]:-} ]]; then align+=("--pl ${pl_array[$idx]}"); fi
       	if [[ -n ${pi_array[$idx]:-} ]]; then align+=("--pi ${pi_array[$idx]}"); fi
       	if [[ -n ${pm_array[$idx]:-} ]]; then align+=("--pm ${pm_array[$idx]}"); fi
        if [[ -n ${dt_array[$idx]:-} ]]; then align+=("--dt ${dt_array[$idx]}"); fi
        if [[ -n ${pu_array[$idx]:-} ]]; then 
            align+=("--pu ${pu_array[$idx]}")
        elif [[ -n ${pu:-} ]]; then
            align+=("--pu $pu")
        fi
        if [[ -n ${id_array[$idx]:-} ]]; then
            align+=("--id ${id_array[$idx]}")
        elif [[ -n ${rgid:-} ]]; then
            align+=("--id $rgid")
        fi
	
	    # run alignment script
        printf "\tAligning FASTQ: %s with command:\n" ${fq_name} >> ${logdir}/${log_file}
        printf "\t\t%s %s \\ \n" ${align[@]} >> ${logdir}/${log_file}
	    ${align[@]} > tmp.log
       
        # check return value and add jobid/filename to merge input lists
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: alignment failed for %s. STDERR output:\n" $sample_name
            cat tmp.log; rm tmp.log; exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            depend="${depend:-},afterok:${jobid}"
            merge="${merge:-},bam/${fq_name}.bam"
        fi
    done

    if [[ -e bam/$sample_name.bam ]]; then
        :
    else
        # remove leading commas
        merge=${merge#*,}
        if [[ -n ${depend:-} ]]; then depend="--depend ${depend#*,} --check no"; fi

        # generate merge command
        merge_cmd=("bash $scriptdir/merge_fixtags_picard.sh"
                   "-i $merge -F $resdir/$fasta -o bam"
                   "-n $sample_name -l logs -q qc/metrics ${depend:-}")
    
        # run merge
        printf "\tMerging runs with command:\n" >> $logdir/$log_file
        printf "\t\t%s %s \\ \n" ${merge_cmd[@]} >> $logdir/$log_file
        ${merge_cmd[@]} > tmp.log

        # check return/jobids
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: merge job failed for %s. STDERR output:\n" $sample_name
            cat tmp.log; rm tmp.log; exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            depend="--depend afterok:$jobid"
            macs_depend="${macs_depend:-},afterok:$jobid"
        fi
    fi

    if [[ -e tracks/$sample_name.bedGraph ]]; then
        :
    else
        # generate tracks
        tracks_cmd=("bash $scriptdir/generate_track_deeptools.sh -b bam/$sample_name.bam"
                    "-o tracks -F $resdir/$fasta -n $sample_name --logdir logs ${depend:-}")
        printf "\tGenerating track with command:\n" >> $logdir/$log_file
        printf "\t\t%s %s \\ \n" ${tracks_cmd[@]} >> $logdir/$log_file
        ${tracks_cmd[@]} > tmp.log
    fi
done

# check if dependencies
if [[ -n ${macs_depend:-} ]]; then
    macs_depend=${macs_depend#*,}
    macs_depend="--depend ${macs_depend} --check no"
fi

# run peak calling if control samples index provided
if [[ -n ${cs_idx:-} ]]; then
    printf "\nCalling peaks with MACS2:\n" >> $logdir/$log_file
    for sample_name in ${samples[@]}; do

        # list control sample names
        cs_array=($(extract_info $cs_idx))

        # check control field isnt empty
        if [[ -z "${cs_array:-}" ]]; then
            printf "\tWARNING: No control sample for sample %s" $sample_name >> $logdir/$log_file
            printf ", skipping peak calling\n" >> $logdir/$log_file
            continue
        fi

        # get unique control samples
        cs_array=($(echo ${cs_array[@]} | tr ' ' '\n' | sort -u | tr '\n' ' '))
        
        # check only one control sample provided
        if [[ ${#cs_array[@]} -ne 1 ]]; then
            printf "\tWARNING: multiple control samples for %s:" $sample_name >> $logdir/$logfile
            printf "\t\t%s\n" ${cs_array[@]} >> $logdir/$logfile
            printf "\tWARNING: skipping peak calling" >> $logdir/$logfile
            continue
        else
            cs_name=${cs_array[@]}
        fi 

        # compile command
        call_peaks=("bash $scriptdir/call_peaks_macs.sh -t bam/$sample_name.bam"
                    "-c bam/$cs_name.bam -g $genome_size -n $sample_name"
                    "-e 200 -F $resdir/$fasta -q qc -o peaks ${macs_depend:-}")

        # run macs script
        printf "\tCalling peaks for sample %s with command:\n" $sample_name >> $logdir/$log_file
        printf "\t\t%s %s \\ \n" ${call_peaks[@]} >> $logdir/$log_file
        ${call_peaks[@]} > tmp.log
    done
else
    printf "\nNo control sample index provided - skipping peak calling\n" >> $logdir/$log_file
fi
