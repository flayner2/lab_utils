#!/bin/sh

# Provide the following arguments when running:
input=$1 # The input file with the ids
ref=$2 # The reference genome
fastq_dir=$3 # The directory where to put the fastq files
outdir=$4 # The directory where to put the result files
threads=${5:-1} # The number of threads to run certain programs

ids=$(cat $input)

# Check if the indexed database files don't already exist
# before creating them with `bwa index`
if [ ! -f "${ref}.*" ]; then 
        bwa index $ref
fi

for id in $ids 
do
        # Download the SRA dataset, dump the fastq files and clear the prefetch cache
        prefetch $id
        fasterq-dump $id --split-files --outdir $fastq_dir -e $threads
        rm -f $SRA_CACHE/* # Export this variable with the path to your prefetch downloads folder
        
        # Calculate some stats to print out later
        raw_reads_1=$(grep -c -h "${id}" ${fastq_dir}/${id}_1*)
        raw_reads_2=$(grep -c -h "${id}" ${fastq_dir}/${id}_2*)
        
        # Run bwa to map the reads to the reference and remove the original fastq files
        bwa mem -t $threads $ref ${fastq_dir}/${id}_1* ${fastq_dir}/${id}_2* > ${outdir}/${id}_bwa_out.sam
        rm -f ${fastq_dir}/${id}*

        # Use samtools to convert SAM to BAM and remove the SAM file
        samtools view -b -o ${outdir}/${id}_bwa_out.bam -@ $threads ${outdir}/${id}_bwa_out.sam
        rm -f ${outdir}/${id}_bwa_out.sam
        
        # Sort and index the BAM files, removing the unsorted BAM
        samtools sort -n -o ${outdir}/${id}_bwa_out.sorted.bam -@ $threads ${outdir}/${id}_bwa_out.bam
        rm -f ${outdir}/${id}_bwa_out.bam
        samtools index -b -@ $threads ${outdir}/${id}_bwa_out.sorted.bam
        
        # Convert BAM to fastq, removing the BAM and BAI files
        # This conversion only outputs paired reads. Singletons and non-paired reads are discarded
        samtools fastq -1 ${outdir}/${id}_enriched_1.fq -2 ${outdir}/${id}_enriched_2.fq -0 /dev/null -s /dev/null -n ${outdir}/${id}_bwa_out.bam.bai
        rm -f ${outdir}/${id}_bwa_out.*

        # Calculate the needed stats and print the report out
        enriched_reads_1=$(grep -c -h "${id}" ${outdir}/${id}_enriched_1*)
        enriched_reads_2=$(grep -c -h "${id}" ${outdir}/${id}_enriched_2*)
        ratio_1=$(echo "scale=4 ; $enriched_reads_1 / $raw_reads_1" | bc)
        ratio_2=$(echo "scale=4 ; $enriched_reads_2 / $raw_reads_2" | bc)
        printf "# of raw reads:\n1: ${raw_reads_1}\n2: ${raw_reads_2}\n"
        printf "# of enriched reads:\n1: ${enriched_reads_1}\n2: ${enriched_reads_2}\n"
        printf "Proportion of mitochondrial reads on the original SRA dataset:\n1: ${ratio_1}\n2: ${ratio_2}\n"
done
