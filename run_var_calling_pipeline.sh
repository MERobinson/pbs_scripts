#!/usr/bin/env bash
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
            printf "\nERROR: Illegal argument: %s %s\n" $1 $2
            echo "$help_message"; exit 1
        ;;
    esac
hift
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

# get unique sample names
samples=($(tail -n +2 "$sample_info" | cut -d "$delim" -f "$sm_idx" | \
           grep -E "[0-9]+C" | sort | uniq))
printf "\nSample list:\n" > $workdir/$logdir/$log_file
printf "\t%s\n" ${samples[@]} >> $workdir/$logdir/$log_file

# for each sample, align, merge runs, mark dup and recal
for sample_name in ${samples[@]}; do
    :
done

# get list of standard chromosomes from fasta file
chr_list=($(cat $workdir/$resdir/$fasta | grep -Eo "^>chr[0-9XYM]+\b" | grep -Eo "chr[0-9MYX]+"))

######################################################
########### Germline SNV/indel calling ###############
######################################################

echo "Calling germline variants"
ggvcf_g_arg=''
ggvcf_depend=''
for sample_name in ${samples[@]}; do
    
    echo "Processing sample: $sample_name"
    cat_v_arg=''
    cat_depend=''
    for chr in ${chr_list[@]}; do
        
        # assemble haplotype caller command
        hc_call="bash $workdir/$scriptdir/call_var_hc.sh --bam $bamdir/$sample_name.recal.bam"
        hc_call="${hc_call} --outdir variants --logdir $logdir"
        hc_call="${hc_call} --fasta $resdir/$fasta --gvcf yes --chr $chr"
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
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            cat_v_arg=${cat_v_arg}",variants/${sample_name}.${chr}.hc.g.vcf"
            cat_depend=${cat_depend}",afterok:${jobid}"
            rm tmp.log
        fi
    done
    cat_v_arg=${cat_v_arg#*,}
    cat_depend=${cat_depend#*,}
   
    
    # combine chr var
    ## edit output/logdir dir for cat script and update this ##
    cat_call="bash $workdir/$scriptdir/cat_vcf_gatk.sh -v $cat_v_arg --depend $cat_depend"
    cat_call="${cat_call} -f $resdir/$fasta --check no --name $sample_name"
    printf "\tConcatenating chromosome variants for sample %s with command:\n" \
           $sample_name >> $workdir/$logdir/$log_file
    printf "\t\t%s %s \\ \n" $cat_call >> $workdir/$logdir/$log_file
    $cat_call > tmp.log

    # check return value / record jobid & filename
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: cat vcf script failed for sample: %s\n" $sample_name
        echo "stderr: "
        cat tmp.log
        exit 1
    else
        jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
        ggvcf_g_arg="${ggvcf_g_arg},variants/${sample_name}.hc.g.vcf"
        ggvcf_depend="${ggvcf_depend},afterok:${jobid}"
        rm tmp.log
    fi
done
ggvcf_g_arg=${ggvcf_g_arg#*,}
ggvcf_depend=${ggvcf_depend#*,}

## joint genotyping
cat_depend=''; cat_v_arg=''
for chr in ${chr_list[@]}; do

    # assemble genotypeGVCFs command
    ggvcf_call="bash $workdir/$scriptdir/genotype_gvcf_gatk.sh"
    ggvcf_call="${ggvcf_call} --gvcf $ggvcf_g_arg --fasta $resdir/$fasta --chr $chr"
    ggvcf_call="${ggvcf_call} --outdir variants --logdir $logdir --depend $ggvcf_depend"
    ggvcf_call="${ggvcf_call} --check no --name $sample_name.$chr"

    # run haplotype caller
    printf "\tPerforming joint genotyping for sample %s, %s with command:\n" \
            $sample_name $chr >> $workdir/$logdir/$log_file
    printf "\t\t%s %s \\ \n" $ggvcf_call >> $workdir/$logdir/$log_file
    $ggvcf_call > tmp.log

    # check return value / record jobid & filename
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: GenotypeGVCF script failed for sample: %s\n" $sample_name
        echo "stderr: "
        cat tmp.log
        exit 1
    else
        jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
        cat_depend=${cat_depend}",afterok:${jobid}"
        cat_v_arg=${cat_v_arg}",variants/${sample_name}.${chr}.vcf"
        rm tmp.log
    fi
done
cat_depend=${cat_depend#*,}
cat_v_arg=${cat_v_arg#*,}

# combine chr var
## edit output/logdir dir for cat script and update this ##
cat_call="bash $workdir/$scriptdir/cat_vcf_gatk.sh -v $cat_v_arg --depend $cat_depend"
cat_call="${cat_call} --fasta $resdir/$fasta --check no --name $sample_name"
printf "\tConcatenating chromosome variants for sample %s with command:\n" \
        $sample_name >> $workdir/$logdir/$log_file
printf "\t\t%s %s \\ \n" $cat_call >> $workdir/$logdir/$log_file
$cat_call


######################################################
############## Germline SV calling ###################
######################################################

echo "Calling germline structural variants"
cat_arg=''
cat_depend=''
for var in INS DEL INV BND DUP; do

    # discover var
    echo "Processing sample: $sample_name"
    merge_arg=''
    merge_depend=''
    for sample_name in ${samples[@]}; do
        # assemble haplotype caller command
        delly_call="bash $workdir/$scriptdir/call_sv_delly.sh --test $bamdir/$sample_name.recal.bam"
        delly_call="${delly_call} -v $var --outdir variants --logdir $logdir --name $sample_name"
        delly_call="${delly_call} --fasta $resdir/$fasta --excl $resdir/$exclude --check no"

        # run delly to discover sv
        printf "\tCalling %s variants for sample %s with command:\n" \
                $var $sample_name >> $workdir/$logdir/$log_file
        printf "\t\t%s %s \\ \n" $delly_call >> $workdir/$logdir/$log_file
        $delly_call > tmp.log

        # check return value
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: call_sv_delly failed for sample: %s\n" $sample_name
            echo "stderr: "
            cat tmp.log
            exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            merge_arg=${merge_arg}",variants/${sample_name}.$var.delly.bcf"
            merge_depend=${merge_depend}",afterok:${jobid}"
            rm tmp.log
        fi
    done
    merge_arg=${merge_arg#*,}
    merge_depend=${merge_depend#*,}

    # merge var accross samples
    merge_call="bash $workdir/$scriptdir/merge_sv_delly.sh --var $var --depend $merge_depend"
    merge_call="${cat_call} --sv_list $merge_arg --check no --name $sample_name"
    merge_call="${merge_call} --outdir variants --logdir $logdir" 
    merge_call="${merge_call} --delly_arg '-m 500 -b 500 -r 0.5'"
    printf "\tMerging %s variants with command:\n" \
           $var >> $workdir/$logdir/$log_file
    printf "\t\t%s %s \\ \n" $merge_call >> $workdir/$logdir/$log_file
    $merge_call > tmp.log

    # check return value
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: call_sv_delly failed for sample: %s\n" $sample_name
        echo "stderr: "
        cat tmp.log; rm tmp.log; exit 1
    else
        recall_depend="afterok:${jobid}"
        rm tmp.log
    fi

    # joint genotyping of sv
    merge_arg=''
    merge_depend=''
    for sample_name in ${samples[@]}; do
        # assemble delly call command
        delly_call="bash $workdir/$scriptdir/call_sv_delly.sh --test $bamdir/$sample_name.recal.bam"
        delly_call="${delly_call} --var $var --outdir variants --logdir $logdir --name $sample_name"
        delly_call="${delly_call} --fasta $resdir/$fasta --excl $resdir/$exclude --check no"
        delly_call="${delly_call} --depend $recall_depend"
        delly_call="${delly_call} --vcf variants/$sample_name.$var.delly.merged.bcf"

        # run delly to discover sv
        printf "\tRe-calling %s variants for sample %s with command:\n" \
                $var $sample_name >> $workdir/$logdir/$log_file
        printf "\t\t%s %s \\ \n" $delly_call >> $workdir/$logdir/$log_file
        $delly_call > tmp.log

        # check return value
        if [[ $? -ne 0 ]]; then
            printf "\nERROR: call_sv_delly failed for sample: %s\n" $sample_name
            echo "stderr: "
            cat tmp.log; rm tmp.log; exit 1
        else
            jobid=$(cat tmp.log | grep -Eo "^JOBID: [0-9]+.cx" | grep -Eo "[0-9]+")
            merge_arg=${merge_arg}",variants/${sample_name}.$var.delly.bcf"
            merge_depend=${merge_depend}",afterok:${jobid}"
            rm tmp.log
        fi
    done
    merge_arg=${merge_arg#*,}
    merge_depend=${merge_depend#*,}

    # merge re-called var accross samples
    merge_call="bash $workdir/$scriptdir/merge_sv_delly.sh --var $var --depend $merge_depend"
    merge_call="${cat_call} --sv_list $merge_arg --check no --name $sample_name"
    merge_call="${merge_call} --outdir variants --logdir $logdir"
    merge_call="${merge_call} --delly_arg '-m 500 -b 500 -r 0.5'"
    printf "\tMerging %s variants with command:\n" \
           $var >> $workdir/$logdir/$log_file
    printf "\t\t%s %s \\ \n" $merge_call >> $workdir/$logdir/$log_file
    $merge_call > tmp.log

    # check return value
    if [[ $? -ne 0 ]]; then
        printf "\nERROR: call_sv_delly failed for sample: %s\n" $sample_name
        echo "stderr: "
        cat tmp.log; rm tmp.log; exit 1
    else
        filter_depend="afterok:${jobid}"
        rm tmp.log
    fi

    # filter final call set
    filter_call="bash $workdir/$scriptdor/filter_sv_delly.sh --vartype $var --depend $filter_depend"
    filter_call="${filter_call} --bcf $vardir/$sample_name.$var.delly.merged.bcf --check no"
    filter_call="${filter_call} --outdir $vardir --logdir logs --name $sample_name" 

done
