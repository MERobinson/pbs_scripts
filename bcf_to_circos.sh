#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

# default arguments
workdir=$PWD
vardir=variants
outdir=circos

# help message
help_message="
usage:
    bash $(basename "$0") [] -b <BCF>
purpose
    # parse Delly2 BCF output files to format suitable for input to circos plots.
required arguments:
    -b|--bcf : input bcf filename to parse
optional arguments:
    -w|--workdir : working directory - all other paths should be relative to workdir (default = pwd)
    -v|--vardir : variants directory containing BCF file (default = variants)
    -o|--outdir : output directory to write circos txt files (default = circos)
    -n|--name : filename prefix output files (default = extracted from bcf file)
    
"

# parse arguments
while [[ $# -gt 1 ]]; do
    key=$1
    case $key in
        -b|--bcf)
            bcf=$2
            shift
            ;;
        -w|--workdir)
            workdir=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            shift
            ;;
        -n|--name)
            name=$2
            shift
            ;;
        *)
            echo "ERROR: Illegal argument supplied"
            echo "$help_message"; exit 1
            ;;
    esac
shift
done

# check required arguments
if [[ ! -d "$workdir/$vardir" ]]; then
    printf "\nERROR: variant directory does not exist/isnt writable: %s/%s\n" $workdir $vardir 
    echo "$help_message"; exit 2
fi
if [[ -z "${bcf:-}" ]]; then
    printf "\nERROR: No --bcf argument provided.\n"
    echo "$help_message"; exit 2
elif [[ ! -r "$workdir/$vardir/$bcf" ]]; then
    printf "\nERROR: Input BCF file does not exist/isnt readable: %s/%s/%s\n" $workdir $vardir $bcf
    echo "$help_message"; exit 2
fi

# get name if not provided
if [[ -z ${name:-} ]]; then
    name=${bcf%.*}
fi 

# create outdir if needed
mkdir -p $workdir/$outdir

# load modules
module load bcftools

# output bcf info to temp file
bcftools query -f '%CHROM\t%POS\t%POS\t%INFO/CHR2\t%END\t%END\t%INFO/SVTYPE\n' "$workdir/$vardir/$bcf" > tmp_file

# loop through each SV type and print appropriately 
grep "BND" tmp_file | sed 's/chr/hs/g' | cut -f 1-6 > $outdir/${name}.BND.circos.txt
declare -a av_arr=("DEL" "DUP" "INS" "INV")
for sv in ${av_arr[@]}; do
    if grep -q $sv tmp_file; then
        grep "$sv" tmp_file | sed 's/chr/hs/g' | cut -f 1,3,6 > $outdir/${name}.${sv}.circos.txt
    else
        > $outdir/${name}.${sv}.circos.txt
    fi
done

# delete temp
rm tmp_file
