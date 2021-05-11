#!/bin/bash

dbs_list=$(ls /home/maycon/Documents/LAB/eusociality/local_data/genes/raw_seqs)
genomes_path="/mnt/2A2A43C42A438C2F/lab_data/genomes/2021-04-19/new/to_annotate"
blast_res_path="/mnt/2A2A43C42A438C2F/lab_data/genes_blast/2021-05-07"
found_genes_path="/mnt/2A2A43C42A438C2F/lab_data/genes_results/2021-05-10"

#for genome_file in `ls $genomes_path/*.fna`; do
        #for db in $dbs_list; do
                #blastn -query $genome_file -db $db -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen gaps" -out ${genome_file}_${db}.blastn -num_threads 4 -evalue 1e-5
        #done
#done

#mv ${genomes_path}/*.blastn $blast_res_path

for blast_result in `ls ${blast_res_path}/*.blastn`; do
        #awk -F '\t' 'OFS="\t" {print $0,($4-$14)/$13}' $blast_result > ${blast_result}.newcol
        LC_ALL="C" sort -t$'\t' -k 15,15 -r -k 3,3 -r -k 11,11 -sg ${blast_result}.newcol > ${blast_result}.sorted
        #rm ${blast_result} ${blast_result}.newcol
        #mv ${blast_result}.sorted ${blast_result}
done

#python get_genes_from_blast_result.py $blast_res_path $genomes_path $found_genes_path
