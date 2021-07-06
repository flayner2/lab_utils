#!/bin/bash

queries_list=$(ls /home/maycon/Documents/LAB/eusociality/local_data/genes/raw_seqs/2021-07-01/partial/*.fas)
subjects_list=$(ls /mnt/2A2A43C42A438C2F/lab_data/longest_cds/2021-07-06/originals/*.fasta)
blast_res_path="/mnt/2A2A43C42A438C2F/lab_data/genes_blast/2021-07-06"

for query_path in $queries_list; do
        query=${query_path##*/}
        for subject_path in $subjects_list; do
                subject=${subject_path##*/}
                blastn -query $query_path -subject $subject_path -outfmt 6 -out $blast_res_path/${subject}_blast_${query}_results.ts
        done
done
