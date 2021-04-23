#!/bin/sh

# Provide two arguments when running:
input=$1 # The input file with the ids
ref=$2 # The reference genome
fastq_dir=$3 # The directory where to put the fastq files
outdir=$4 # The directory where to put the sam file

ids=$(cat $input)

for id in $ids 
do
        prefetch $id
        fastq-dump --split-files --fastq_dir $fastq_dir $id
        raw_reads_1=$(grep -c -h ">" ${fastq_dir}/${id}_1*)
        raw_reads_2=$(grep -c -h ">" ${fastq_dir}/${id}_2*)
        printf "# of raw reads:\n1: $raw_reads_1\n2: $raw_reads_2"
        bwa mem $ref ${fastq_dir}/${id}_1* ${fastq_dir}/${id}_2* > ${outdir}/{$id}_bwa_out.sam
        #rm -f ${fastq_dir}/${id}*
        #rm -f $SRA_CACHE/*
done
