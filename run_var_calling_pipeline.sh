#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arg
workdir=$PWD
fqdir=fastq
bamdir=bam
resdir=resources
logdir=logs
scriptdir=scripts
delim=","
fasta=Homo_sapiens_assembly38.fasta
index=hg38bundle
dbsnp=dbsnp_144.hg38.vcf.gz
log_file=$(basename $0 .sh)_$(date "+%Y-%m-%d").log

# sample info index defaults
fq1_idx=1
fq2_idx=2
sm_idx=3

# help message
help_message="
usage:
    bash $(basename "$0") -s <SAMPLE_INFO>
purpose:
    # parse a sample info file and run variant calling
    # SNPs and indels are called with GATKs HaplotypeCaller/Mutect2
    # SV are called with Delly2
required arguments:
    -i|--info : a delimited text file containing sample information (see below)
optional arguments:
    -w|--workdir : working directory (default = pwd)
    -f|--fqdir : directory within workdir containing FASTQ files (default = fastq)
    -r|--resdir : directory within workdir containing resource files (default = resources)
    -l|--logdir : directory within workdir to output log files (default = logs)
    -s|--scriptdir : directory within workdir containing pipeline scripts (default = scripts)
    -d|--delim : delimiter used in sample info file (default = ',')
additional info:
    # sample info file must minimally contain FASTQ filenames and sample names
    # additional columns may be included as listed below:
        > Filename of FASTQ R1 [required]
        > Filename of FASTQ R2 [required - leave blank column for SE]
        > Sample name [required]
        > Control sample name [optional]
        > Read group ID [optional]
        > Platform unit [optional]
        > Library name [optional]
        > Sequencing centre [optional]
        > Platform [default = ILLUMINA]
        > Platform model [optional]
        > Run date [optional]
        > Predicted median insert size [optional] 
    # to include any of the above - provide column indexes arguments as follows:
        --f1_idx : index of column containing FASTQ R1 filenames (default = 1)
        --f2_idx : index of column containing FASTQ R2 filenames (default = 2)
        --sm_idx : index of column containing sample names (default = 3)
        --cs_idx : index of column containing sample name of control (default = null) 
        --id_idx : index of column containing read group ID (default = null)
        --pu_idx : index of column containing platform unit (default = null)
        --lb_idx : index of column containing library name (default = null)
        --cn_idx : index of column containing sequencing centre (default = null)
        --pl_idx : index of column containing platform (default = null)
        --pi_idx : index of column containing median insert size (deafult = null)
        --pm_idx : index of column containing platform model (default = null)
        --dt_idx : index of column containing run date (default = null)
    # NB: If control sample index is provided (and field is non-empty) both germline and
          somatic variant calling will be run, otherwise ONLY germline var are called
    # the following resource files must be present in resdir:
        > whole genome sequence [FASTA and .dict files]
        > BWA index files
        > dbSNP VCF file [VCF and index]
    # the default expected resource filenames can be altered with following arg:
        -ix|--index_base : basename of index files in resdir (default = genome.fa)
        -fa|--fasta : name of fasta file (default = genome.fa)

"

# parse arguments
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -i|--info)
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
        --f1_idx)
            fq1_idx=$2
            shift
            ;;
        --f2_idx)
            fq2_idx=$2
            shift
            ;;
        --sm_idx)
            sm_idx=$2
            shift
            ;;
        --cs_idx)
            cs_idx=$2
            shift
            ;;
        --lb_idx)
            lb_idx=$2
            shift
            ;;
        --pu_idx)
            pu_idx=$2
            shift
            ;;
        --pl_idx)
            pl_idx=$2
            shift
            ;;
        --cn_idx)
            cn_idx=$2
            shift
            ;;
        --pi_idx)
            pi_idx=$2
            shift
            ;;
        --pm_idx)
            pm_idx=$2
            shift
            ;;
        --dt_idx)
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
        -db|--dbsnp)
            dbsnp=$2
            shift
            ;;
        *)
            echo "Error: Illegal argument"
            echo "$help_message"; exit 1
        ;;
    esac
shift
done

# check required arg
if [[ -z "${sample_info:-}" ]]; then
    printf "\nERROR: Sample info file not provided.\n"
    echo "$help_message"; exit 2
fi
if [[ ! -r "$workdir/$sample_info" ]]; then
    printf "\nERROR: Input sample file doesnt exist/isnt readable: %s/%s\n" $workdir $sample_info
    echo "$help_message"; exit 2
fi
if [[ ! -d "$workdir/$resdir" ]]; then
    printf "\nERROR: Resources directory doesnt exist: %s/%s\n" $workdir $resdir
    echo "$help_message"; exit 2
fi
if [[ ! -d "$workdir/$fqdir" ]]; then
    printf "\nERROR: FASTQ directory doesnt exist: %s/%s\n" $workdir $fqdir
    echo "$help_message"; exit 2
fi
if [[ ! -r "$workdir/$resdir/$index.sa" ]]; then
    printf "\nERROR: Index file doesnt exist/isnt readable: %s/%s/%s\n" $workdir $resdir $index
    echo "$help_message"; exit 2
fi
if [[ ! -r "$workdir/$resdir/$fasta" ]]; then
    printf "\nERROR: FASTA file doesnt exist/isnt readable: %s/%s/%s\n" $workdir $resdir $fasta
    echo "$help_message"; exit 2
fi

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

# for each sample, align, merge runs, mark dup and recal
for sample_name in ${samples[@]}; do
    :
done

# get list of standard chromosomes from fasta file
chr_list=($(cat $workdir/$resdir/$fasta | grep -Eo "^>chr[0-9XYM]+\b" | grep -Eo "chr[0-9MYX]+"))

# call germline variants with haplotype caller
for sample_name in ${samples[@]}; do
    
    #declare var to store args
    v_arg=''
    depend=''

    # call each chr seperately
    for chr in ${chr_list[@]}; do
        
        # assemble haplotype caller command
        hc_call="bash $workdir/$scriptdir/call_var_hc.sh -b $bamdir/$sample_name.recal.bam"
        hc_call="${hc_call} -f $resdir/$fasta --gvcf yes --intervals $chr"
        hc_call="${hc_call} --dbsnp $resdir/$dbsnp --check no --name $sample_name.$chr"
        
        # run haplotype caller
        printf "\tCalling germline variants for sample %s, %s with command:\n" \
                $sample_name $chr >> $workdir/$logdir/$log_file
        printf "\t\t%s %s \\ \n" $hc_call >> $workdir/$logdir/$log_file
        $hc_call > tmp.log
        
        # check return value
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: haplotype caller script failed for sample: %s\n" $sample_name
            echo "stderr: "
            cat tmp.log
            exit 1
        else
            jobid=$(cat tmp.log | grep -E "^JOBID:" | grep -Eo "[0-9]+")
            depend=${depend}",afterok:${jobid}"
            v_arg=${v_arg}",variants/${sample_name}.${chr}.hc.g.vcf"
            rm tmp.log
        fi
    done

    # remove leading commas
    depend=${depend#*,}
    v_arg=${v_arg#*,}
    
    # combine chr var
    cat_call="bash $workdir/$scriptdir/cat_vcf_gatk.sh -v $v_arg --depend $depend"
    cat_call="${cat_call} -f $resdir/$fasta --check no --name $sample_name"
    printf "\tConcatenating chromosome variants for sample %s with command:\n" \
           $sample_name >> $workdir/$logdir/$log_file
    printf "\t\t%s %s \\ \n" $cat_call >> $workdir/$logdir/$log_file
    $cat_call > tmp.log

    exit

done



 

