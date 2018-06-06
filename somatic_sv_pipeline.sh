#!/bin/bash
# usage: bash $(basename "$0") [filter|merge|call|cat|combine] [DEL|INS|INV|BND|DUP]

step=$1
vartype=$2

# temp - needs proper wrapper writing
fasta=resources/Homo_sapiens_assembly38.fasta
exclude=resources/exclude.tsv
chr_list=$(cat $fasta | grep -Eo ">chr[0-9XYM]+\b" | grep -Eo "chr[0-9XYM]+")

# filter variants from original Delly2 call
if [[ $step = "filter1" ]]; then
    for bcf in tmp_files/*$vartype.c1.delly.bcf; do
        name=$(basename "$bcf"); name=${name%%.*}
        control=${name%%T*}; control=${control}C
        echo "TUMOR: $name CONTROL: $control"
        printf "$name\ttumor\n$control\tcontrol\n" > tmp_files/$name.sample_info
        if [[ -e tmp_files/$name.$vartype.c1.filt.bcf ]]; then
            :
        else
            bash scripts/filter_sv_delly.sh --bcf $bcf --vartype $vartype \
                --outdir tmp_files --name $name.$vartype.c1.filt --logdir logs \
                --sample_info tmp_files/$name.sample_info
        fi
    done
elif [[ $step = "merge1" ]]; then
    sv_list=$(printf "tmp_files/%sT.$vartype.c1.filt.bcf," $(seq 1 10)); sv_list=${sv_list%,*}
    if [[ -e tmp_files/$vartype.merge.bcf ]]; then
        :
    else
        bash scripts/merge_sv_delly.sh --sv_list $sv_list --vartype $vartype \
            --outdir tmp_files --name $vartype.c1.merge --logdir logs \
            --min 250 --max 100000000 --overlap 0.6 --offset 500
    fi
elif [[ $step = "call2" ]]; then
    for bam in bam/*.recal.bam; do
        name=$(basename "$bam"); name=${name%%.*}
        if [[ -e tmp_files/$name.$vartype.c2.bcf ]]; then
            :
        else
            bash scripts/call_sv_delly.sh --bam $bam --vartype $vartype --fasta $fasta \
                --outdir tmp_files --logdir logs --exclude $exclude \
                --name $name.$vartype.c2 \
                --vcf tmp_files/$vartype.c1.merge.bcf
        fi  
    done
elif [[ $step = "merge2" ]]; then
    bcf_list=$(printf "tmp_files/%sC.$vartype.c2.bcf," $(seq 1 10))
    bcf_list=${bcf_list}$(printf "tmp_files/%sT.$vartype.c2.bcf," $(seq 1 10))
    bcf_list=${bcf_list%,*}
    bash scripts/merge_bcf_bcftools.sh --bcf_list $bcf_list --name somatic_sv.$vartype.c2.merge \
        --outdir tmp_files --logdir logs
elif [[ $step = "filter2" ]]; then
    printf "%sC\tcontrol\n" {1..10} > $vartype.c2.sample_info
    printf "%sT\ttumor\n" {1..10} >> $vartype.c2.sample_info
    bash scripts/filter_sv_delly.sh --bcf tmp_files/somatic_sv.$vartype.c2.merge.bcf \
        --vartype $vartype --outdir variants \
        --name somatic_sv.$vartype --logdir logs \
        --sample_info $vartype.c2.sample_info
fi
