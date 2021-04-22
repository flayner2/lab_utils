#!/bin/sh

# Provide two arguments when running:
input=$1 # The input file with the ids
outdir=$2 # The directory where to put the fastq files

ids=$(cat $input)

for id in $ids 
do 
        prefetch $id
        fastq-dump --split-files --outdir $outdir $id
        printf "Reads for $id:\n" && grep -c ">" ${outdir}/${id}*
        rm -f ${outdir}/${id}*
done
