#!/bin/bash

queries_list=$(ls /mnt/2A2A43C42A438C2F/lab_data/genes_results/2021-07-08/originals_complete_cds/*.fas)
subjects_list=$(ls /mnt/2A2A43C42A438C2F/lab_data/longest_cds/2021-07-06/new_annotated/*.longest)
blast_res_path="/mnt/2A2A43C42A438C2F/lab_data/genes_blast/2021-07-08"

for query_path in $queries_list; do
        query=${query_path##*/}
        for subject_path in $subjects_list; do
                subject=${subject_path##*/}
                blastn -query $query_path -subject $subject_path -outfmt 6 -out $blast_res_path/${subject}_blast_${query}_results.ts
        done
done
