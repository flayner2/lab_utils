#!/bin/sh

# Provide the following arguments when running:
input=$1 # The input file with the ids
ref=$2 # The reference genome
fastq_dir=$3 # The directory where to put the fastq files
outdir=$4 # The directory where to put the result files
threads=${5:-1} # The number of threads to run certain programs
mapq_qual=${6:-30}

ids=$(cat $input)

# Check if the indexed database files don't already exist
# before creating them with `bwa index`
if [ ! -f "${ref}.*" ]; then 
        bwa index $ref
fi

for id in $ids 
do
        # Download the SRA dataset, dump the fastq files and clear the prefetch cache
        #prefetch $id
        #fasterq-dump $id --split-files --outdir $fastq_dir -e $threads
        #rm -f $SRA_CACHE/* # Export this variable with the path to your prefetch downloads folder
        
        # Calculate some stats to print out later
        raw_reads_1=$(grep -c -h "${id}" ${fastq_dir}/${id}_1*)
        raw_reads_2=$(grep -c -h "${id}" ${fastq_dir}/${id}_2*)
        temp=$(tail -n +2 ${fastq_dir}/${id}_1* | head -n 1 | wc -m)
        length=$(($temp - 1))
        raw_nucleotides=$(($length * $raw_reads_1))

        
        # Run bwa to map the reads to the reference and remove the original fastq files
        bwa mem -t $threads -T $mapq_qual $ref ${fastq_dir}/${id}_1* ${fastq_dir}/${id}_2* > ${outdir}/${id}_bwa_out.sam
        #rm -f ${fastq_dir}/${id}*

        ## Sort the SAM file, converting it to BAM and removing the unsorted SAM
        samtools sort -n -o ${outdir}/${id}_bwa_out.sorted.bam -@ $threads ${outdir}/${id}_bwa_out.sam
        #rm -f ${outdir}/${id}_bwa_out.sam

        ## Use samtools view to filter the sorted BAM file, removing non-propper pairs. Remove the unfiltered BAM file
        samtools view -F 0x04 -f 0x2 -q $mapq_qual -o ${outdir}/${id}_bwa_out.sorted.filtered.bam -@ $threads -b ${outdir}/${id}_bwa_out.sorted.bam
        rm -f ${outdir}/${id}_bwa_out.sorted.bam
        
        # Convert BAM to fastq, removing the BAM files
        # This conversion only outputs paired reads. Singletons and non-paired reads are discarded
        samtools fastq -1 ${outdir}/${id}_enriched_1.fq -2 ${outdir}/${id}_enriched_2.fq -0 /dev/null -s /dev/null -n ${outdir}/${id}_bwa_out.sorted.filtered.bam
        rm -f ${outdir}/${id}_bwa_out.*

        # Calculate the needed stats and print the report out
        enriched_reads_1=$(grep -c -h "${id}" ${outdir}/${id}_enriched_1*)
        enriched_reads_2=$(grep -c -h "${id}" ${outdir}/${id}_enriched_2*)
        enriched_nucleotides=$(($length * $enriched_reads_1))
        ratio_1=$(echo "scale=6 ; $enriched_reads_1 / $raw_reads_1" | bc)
        ratio_2=$(echo "scale=6 ; $enriched_reads_2 / $raw_reads_2" | bc)
        bases_ratio=$(echo "scale=6 ; $enriched_nucleotides / $raw_nucleotides" | bc)
        printf "# of raw reads:\n1: ${raw_reads_1}\n2: ${raw_reads_2}\n" >> ${outdir}/${id}_report.txt
        printf "# of raw bases: ${raw_nucleotides}\n" >> ${outdir}/${id}_report.txt
        printf "# of enriched reads:\n1: ${enriched_reads_1}\n2: ${enriched_reads_2}\n" >> ${outdir}/${id}_report.txt
        printf "# of enriched bases: ${enriched_nucleotides}\n" >> ${outdir}/${id}_report.txt
        printf "Proportion of enriched reads on the original SRA dataset:\n1: ${ratio_1}\n2: ${ratio_2}\n" >> ${outdir}/${id}_report.txt
        printf "Proportion of enriched nucleotides on the original SRA dataset: ${bases_ratio}\n" >> ${outdir}/${id}_report.txt

        # Prepare the environment and get some data to run MitoZ
        cd $outdir
                
        # Run MitoZ in quick mode
        docker run -v $PWD:$PWD -w $PWD --rm guanliangmeng/mitoz:2.3 /app/release_MitoZ_v2.3/MitoZ.py all --genetic_code 5 --clade Arthropoda --outprefix ${id}_mit --thread_number $threads --fastq1 ${id}_enriched_1* --fastq2 ${id}_enriched_2* --fastq_read_length $length --run_mode 2 --requiring_taxa 'Arthropoda' 
        #sudo rm -rf tmp/

        #python3 create_mitoz_multikmer_files.py ${outdir}/${id}_mit.result/summary.txt $outdir
        #missing=$(cat missing_genes.txt)

        ## Run MitoZ in multi-kmer mode
        #if test -f "${id}_mit.result/*besthit.sim.filtered.high_abundance_10.0X.reformat.sorted"; then 
                #docker run -v $PWD:$PWD -w $PWD --rm guanliangmeng/mitoz:2.3 /app/release_MitoZ_v2.3/MitoZ.py all2 --genetic_code 5 --clade Arthropoda --outprefix ${id}_mit_2 --thread_number $threads --fastq1 ${id}_enriched_1* --fastq2 ${id}_enriched_2* --fastq_read_length $length --run_mode 3 --requiring_taxa 'Arthropoda' --quick_mode_seq_file ${id}_mit.result/work71.mitogenome.fa --quick_mode_fa_genes_file found_genes.txt --missing_PCGs $missing --quick_mode_score_file ${id}_mit.result/work71.hmmtblout.besthit.sim.filtered.high_abundance_10.0X.reformat.sorted --quick_mode_prior_seq_file ${id}_mit.result/work71.hmmtblout.besthit.sim.filtered.fa
        #else
                #docker run -v $PWD:$PWD -w $PWD --rm guanliangmeng/mitoz:2.3 /app/release_MitoZ_v2.3/MitoZ.py all2 --genetic_code 5 --clade Arthropoda --outprefix ${id}_mit_2 --thread_number $threads --fastq1 ${id}_enriched_1* --fastq2 ${id}_enriched_2* --fastq_read_length $length --run_mode 3 --requiring_taxa 'Arthropoda' --quick_mode_seq_file ${id}_mit.result/work71.mitogenome.fa --quick_mode_fa_genes_file found_genes.txt --missing_PCGs $missing
        #fi
done
