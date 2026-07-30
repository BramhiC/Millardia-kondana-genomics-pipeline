#!/bin/bash

#====================================================
# Millardia kondana whole-genome sequencing pipeline
# Sequence processing
#====================================================

# Populations
# K  = Sinhagad
# R  = Rajgad
# T  = Torna
# RR = Raireshwar
# MEL = M. meltada


#====================================================
# Software installation
#====================================================

conda create -n mkondana_wgs -y
conda activate mkondana_wgs

conda install -c bioconda fastqc multiqc fastp bwa-mem2 samtools bcftools gatk4 plink -y


#====================================================
# Create project structure
#====================================================

mkdir -p raw_fastq trimmed fastqc_raw fastqc_trimmed
mkdir -p multiqc_raw multiqc_trimmed reference
mkdir -p alignment bam variants stats metadata


#====================================================
# Rename downloaded files if needed
#====================================================

for file in *\\?download=1; do
    newname=$(echo "$file" | sed "s/?download=1//")
    mv "$file" "$newname"
done


#====================================================
# Count sequencing reads
#====================================================

for file in raw_fastq/*.fastq.gz; do
  read_count=$(( $(zcat "$file" | wc -l) / 4 ))
  echo "$(basename "$file"): $read_count reads"
done


#====================================================
# Reads shorter than 150 bp
#====================================================

for file in raw_fastq/*.fastq.gz; do
  reads=$(zcat "$file" | awk 'NR % 4 == 2 {print length($0)}' | awk '$1<150' | wc -l)
  echo "$(basename "$file"): $reads reads shorter than 150bp"
done


#====================================================
# Reads longer than 150 bp
#====================================================

for file in raw_fastq/*.fastq.gz; do
  reads=$(zcat "$file" | awk 'NR % 4 == 2 {print length($0)}' | awk '$1>150' | wc -l)
  echo "$(basename "$file"): $reads reads longer than 150bp"
done


#====================================================
# Reads not equal to 150 bp
#====================================================

for file in raw_fastq/*.fastq.gz; do
  count=$(zcat "$file" | awk 'NR % 4 == 2 {if (length($0) != 150) print length($0)}' | wc -l)
  echo "$(basename "$file"): $count reads not equal to 150bp"
done


#====================================================
# Check adapter contamination
#====================================================

for file in raw_fastq/*.fastq.gz; do
  count=$(zcat "$file" | awk 'NR % 4 == 2' | grep -c 'CTGTCTCTTATACACATCT')
  echo "$(basename "$file"): $count adapter-contaminated reads"
done


#====================================================
# FastQC on raw reads
#====================================================

fastqc raw_fastq/*.fastq.gz -o fastqc_raw -t 16

multiqc fastqc_raw -o multiqc_raw


#====================================================
# Adapter and quality trimming
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  fastp \
    -i raw_fastq/${sample}_R1.fastq.gz \
    -I raw_fastq/${sample}_R2.fastq.gz \
    -o trimmed/${sample}_R1.trim.fastq.gz \
    -O trimmed/${sample}_R2.trim.fastq.gz \
    -h trimmed/${sample}.html \
    -j trimmed/${sample}.json \
    --detect_adapter_for_pe \
    --thread 16
done


#====================================================
# FastQC after trimming
#====================================================

fastqc trimmed/*.fastq.gz -o fastqc_trimmed -t 16

multiqc fastqc_trimmed -o multiqc_trimmed


#====================================================
# Prepare reference genome
#====================================================

# Place Rrattus.fa in reference/

bwa-mem2 index reference/Rrattus.fa

samtools faidx reference/Rrattus.fa

gatk CreateSequenceDictionary \
-R reference/Rrattus.fa \
-O reference/Rrattus.dict


#====================================================
# Read mapping
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  bwa-mem2 mem \
    -t 16 \
    reference/Rrattus.fa \
    trimmed/${sample}_R1.trim.fastq.gz \
    trimmed/${sample}_R2.trim.fastq.gz \
    > alignment/${sample}.sam
done


#====================================================
# Convert SAM to sorted BAM
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  samtools sort -@ 16 \
    -o bam/${sample}.sorted.bam \
    alignment/${sample}.sam

  samtools index bam/${sample}.sorted.bam
done


#====================================================
# Mark PCR duplicates
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  gatk MarkDuplicates \
    -I bam/${sample}.sorted.bam \
    -O bam/${sample}.dedup.bam \
    -M stats/${sample}.dup_metrics.txt

  samtools index bam/${sample}.dedup.bam
done


#====================================================
# Alignment statistics
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  samtools flagstat bam/${sample}.dedup.bam > stats/${sample}.flagstat

  samtools stats bam/${sample}.dedup.bam > stats/${sample}.stats
done


#====================================================
# Mean coverage
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  samtools depth -a bam/${sample}.dedup.bam | \
  awk '{sum+=$3} END {print sum/NR}' \
  > stats/${sample}.mean_coverage.txt
done


#====================================================
# Create BAM list
#====================================================

ls bam/*.dedup.bam > bamlist.txt


#====================================================
# Joint variant calling
#====================================================

bcftools mpileup \
-f reference/Rrattus.fa \
-b bamlist.txt \
-Ou \
-q 20 \
-Q 20 \
-a AD,DP | \
bcftools call \
-m -v -Oz \
-o variants/variants.raw.vcf.gz

bcftools index variants/variants.raw.vcf.gz


#====================================================
# Retain SNPs only
#====================================================

bcftools view -v snps \
variants/variants.raw.vcf.gz \
-Oz -o variants/variants.snps.vcf.gz


#====================================================
# QUAL > 20
#====================================================

bcftools filter -i 'QUAL>20' \
variants/variants.snps.vcf.gz \
-Oz -o variants/variants.qual20.vcf.gz


#====================================================
# Remove rare variants (MAC > 3)
#====================================================

bcftools view -i 'MAC>3' \
variants/variants.qual20.vcf.gz \
-Oz -o variants/variants.mac3.vcf.gz


#====================================================
# Remove SNPs with >20% missing data
#====================================================

bcftools +fill-tags variants/variants.mac3.vcf.gz -- -t F_MISSING | \
bcftools view -i 'F_MISSING<0.2' \
-Oz -o variants/variants.filtered.vcf.gz


bcftools index variants/variants.filtered.vcf.gz


#====================================================
# Final files
#====================================================

# Clean BAMs:
# bam/*.dedup.bam

# Filtered SNPs:
# variants/variants.filtered.vcf.gz
