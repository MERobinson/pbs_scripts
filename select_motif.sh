#!/bin/bash

# usage: $0 <motifDB> <motif,list> 
# takes motif DB and a comma sep list of motifs
# outputs filtered DB to stdout

IFS=',' read -r -a motif_array <<< $2

contains () {
    local e match="$1"
    shift
    for e; do [[ "$e" == "$match" ]] && return 0; done
    return 1
}

keep=F
while read line; do
    if [[ ${line:0:1} == ">" ]]; then
        header=($line)
        if contains "${header[1]}" "${motif_array[@]}"; then
            echo "$line"
            keep=T
        else
            keep=F
        fi
    elif [[ $keep = T ]]; then
        echo "$line"
    fi
done < $1

exit 0
