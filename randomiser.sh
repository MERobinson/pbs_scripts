#!/bin/bash

if [[ -z $1 ]] || [[ -z $2 ]]; then
    echo "ERROR: input and output dir required"
    printf "\nUsage: \nbash %s <INDIR> <OUTDIR>\n\n" $0
    exit 1
else
    indir=$1
    outdir=$2
fi

if [[ -d $outdir ]]; then
    echo "ERROR: output dir already exists"
    exit 1
else
    mkdir -p $outdir
    > $outdir/imgcodes.txt
fi

for f in $indir/*; do
    code=$(echo $f | md5sum | cut -f1 -d' ').tiff
    cp $f $outdir/$code
    printf "%s\t%s\n" $f $outdir/$code > $outdir/imgcodes.txt
done
