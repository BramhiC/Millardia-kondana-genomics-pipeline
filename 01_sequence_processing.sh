#!/bin/bash

#====================================================
# Millardia kondana whole-genome population genomics
# Reference-based mapping using Rattus rattus genome
#====================================================


#----------------------------------------------------
# Configuration
# Replace placeholder values before running
#----------------------------------------------------

REF=reference/Rrattus.fa

THREADS=<NUM_THREADS>

MIN_MAPQ=<MAPQ_THRESHOLD>
MIN_BASEQ=<BASEQ_THRESHOLD>
MIN_QUAL=<VARIANT_QUAL_THRESHOLD>
MIN_MAC=<MINOR_ALLELE_COUNT>
MAX_MISSING=<MISSINGNESS_THRESHOLD>


#----------------------------------------------------
# Populations
#----------------------------------------------------

# K   = Sinhagad
# R   = Rajgad
# T   = Torna
# RR  = Raireshwar
# MEL = Millardia meltada


#====================================================
# Software installation
#====================================================

conda create -n mkondana_wgs -y
conda activate mkondana_wgs

conda install -c bioconda \
fastqc multiqc fastp bwa-mem2 samtools bcftools gatk4 \
plink qualimap -y


#====================================================
# Create project structure
#====================================================

mkdir -p raw_fastq
mkdir -p trimmed
mkdir -p fastqc_raw
mkdir -p fastqc_trimmed
mkdir -p multiqc_raw
mkdir -p multiqc_trimmed
mkdir -p reference
mkdir -p alignment
mkdir -p bam
mkdir -p variants
mkdir -p stats
mkdir -p qualimap
mkdir -p metadata


#====================================================
# Quality control of raw reads
#====================================================

fastqc raw_fastq/*.fastq.gz -o fastqc_raw -t ${THREADS}

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
    --thread ${THREADS}
done


#====================================================
# Quality control after trimming
#====================================================

fastqc trimmed/*.fastq.gz -o fastqc_trimmed -t ${THREADS}

multiqc fastqc_trimmed -o multiqc_trimmed


#====================================================
# Prepare reference genome (Rattus rattus)
#====================================================

bwa-mem2 index ${REF}

samtools faidx ${REF}

gatk CreateSequenceDictionary \
-R ${REF} \
-O reference/Rrattus.dict


#====================================================
# Read mapping
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  bwa-mem2 mem \
    -t ${THREADS} \
    ${REF} \
    trimmed/${sample}_R1.trim.fastq.gz \
    trimmed/${sample}_R2.trim.fastq.gz \
    > alignment/${sample}.sam
done


#====================================================
# Convert SAM to sorted BAM
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  samtools sort -@ ${THREADS} \
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
# Mean sequencing coverage
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  samtools depth -a bam/${sample}.dedup.bam | \
  awk '{sum+=$3} END {print sum/NR}' \
  > stats/${sample}.mean_coverage.txt
done


#====================================================
# Qualimap alignment assessment
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  qualimap bamqc \
    -bam bam/${sample}.dedup.bam \
    -outdir qualimap/${sample} \
    -nt ${THREADS}
done


#====================================================
# Create BAM list for downstream analyses
#====================================================

ls bam/*.dedup.bam > bamlist.txt


#====================================================
# Joint variant calling
#====================================================

bcftools mpileup \
-f ${REF} \
-b bamlist.txt \
-Ou \
-q ${MIN_MAPQ} \
-Q ${MIN_BASEQ} \
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
# Filter by variant quality
#====================================================

bcftools filter -i "QUAL>${MIN_QUAL}" \
variants/variants.snps.vcf.gz \
-Oz -o variants/variants.qual.vcf.gz


#====================================================
# Remove rare variants
#====================================================

bcftools view -i "MAC>${MIN_MAC}" \
variants/variants.qual.vcf.gz \
-Oz -o variants/variants.mac.vcf.gz


#====================================================
# Remove SNPs with excessive missing data
#====================================================

bcftools +fill-tags variants/variants.mac.vcf.gz -- -t F_MISSING | \
bcftools view -i "F_MISSING<${MAX_MISSING}" \
-Oz -o variants/variants.filtered.vcf.gz

bcftools index variants/variants.filtered.vcf.gz


#====================================================
# Prepare autosomal reference for PSMC / MSMC
# Remove mitochondrial and sex-chromosome scaffolds
#====================================================

grep -v -E "chrX|chrY|MT|M|mitochond" ${REF}.fai | cut -f1 > autosomes.list

samtools faidx ${REF} $(cat autosomes.list) > reference/Rrattus.autosomes.fa


#====================================================
# Notes for demographic analyses
#====================================================

# - PSMC should be run on high-quality individuals with adequate coverage.
# - Minimum depth for PSMC consensus generation should be ≥10×.
# - Exact depth thresholds should be chosen after inspecting coverage distribution.
# - Sensitivity analysis can be performed using alternative -t and -p values.
# - If multiple high-quality phased genomes are available, MSMC can be evaluated
#   as an alternative or complementary approach to PSMC.


#====================================================
# Final outputs
#====================================================

# Clean BAM files:
# bam/*.dedup.bam

# Final filtered SNPs:
# variants/variants.filtered.vcf.gz

# QC reports:
# multiqc_raw/multiqc_report.html
# multiqc_trimmed/multiqc_report.html

# Alignment QC:
# qualimap/*/genome_results.txt
