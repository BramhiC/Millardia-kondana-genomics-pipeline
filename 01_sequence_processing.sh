Yes. Your first page should also be made more compact and lab-style.

I removed:

* tutorial-style explanations,

* “expected outputs” prose,

* long descriptive comments,

* AI-sounding sentences.

I kept only:

* section headers,

* essential notes,

* commands.

Here is the cleaned `01_sequence_processing.sh` in the same style as the revised Page 2.

Bash

```
#!/bin/bash

#====================================================
# Millardia kondana whole-genome population genomics
# Reference-based mapping using Rattus rattus
#====================================================

REF=reference/Rrattus.fa

THREADS=<NUM_THREADS>

MIN_MAPQ=<MAPQ_THRESHOLD>
MIN_BASEQ=<BASEQ_THRESHOLD>
MIN_QUAL=<VARIANT_QUAL_THRESHOLD>
MIN_MAC=<MINOR_ALLELE_COUNT>
MAX_MISSING=<MISSINGNESS_THRESHOLD>


#====================================================
# Software
#====================================================

conda create -n mkondana_wgs -y
conda activate mkondana_wgs

conda install -c bioconda \
fastqc multiqc fastp bwa-mem2 samtools bcftools gatk4 \
plink qualimap -y


#====================================================
# Directories
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
# FastQC: raw reads
#====================================================

fastqc raw_fastq/*.fastq.gz -o fastqc_raw -t ${THREADS}

multiqc fastqc_raw -o multiqc_raw


#====================================================
# fastp
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
# FastQC: trimmed reads
#====================================================

fastqc trimmed/*.fastq.gz -o fastqc_trimmed -t ${THREADS}

multiqc fastqc_trimmed -o multiqc_trimmed


#====================================================
# Reference
#====================================================

bwa-mem2 index ${REF}

samtools faidx ${REF}

gatk CreateSequenceDictionary \
-R ${REF} \
-O reference/Rrattus.dict


#====================================================
# BWA-MEM2
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
# SAM -> sorted BAM
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  samtools sort -@ ${THREADS} \
    -o bam/${sample}.sorted.bam \
    alignment/${sample}.sam

  samtools index bam/${sample}.sorted.bam
done


#====================================================
# Mark duplicates
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
# Qualimap
#====================================================

for sample in K1 K2 K3 R1 R2 T1 T2 RR1 RR2 MEL1
do
  qualimap bamqc \
    -bam bam/${sample}.dedup.bam \
    -outdir qualimap/${sample} \
    -nt ${THREADS}
done


#====================================================
# BAM list
#====================================================

ls bam/*.dedup.bam > bamlist.txt


#====================================================
# Variant calling
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
# SNPs only
#====================================================

bcftools view -v snps \
variants/variants.raw.vcf.gz \
-Oz -o variants/variants.snps.vcf.gz


#====================================================
# QUAL filter
#====================================================

bcftools filter -i "QUAL>${MIN_QUAL}" \
variants/variants.snps.vcf.gz \
-Oz -o variants/variants.qual.vcf.gz


#====================================================
# MAC filter
#====================================================

bcftools view -i "MAC>${MIN_MAC}" \
variants/variants.qual.vcf.gz \
-Oz -o variants/variants.mac.vcf.gz


#====================================================
# Missingness filter
#====================================================

bcftools +fill-tags variants/variants.mac.vcf.gz -- -t F_MISSING | \
bcftools view -i "F_MISSING<${MAX_MISSING}" \
-Oz -o variants/variants.filtered.vcf.gz

bcftools index variants/variants.filtered.vcf.gz


#====================================================
# Autosomal reference for MSMC / PSMC
#====================================================

grep -v -E "chrX|chrY|MT|M|mitochond" ${REF}.fai | cut -f1 > autosomes.list

samtools faidx ${REF} $(cat autosomes.list) > reference/Rrattus.autosomes.fa
```
