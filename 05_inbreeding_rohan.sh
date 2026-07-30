#!/bin/bash

#====================================================
# Millardia kondana whole-genome sequencing pipeline
# Inbreeding and runs of homozygosity (ROH)
#====================================================

# ROHan is used to estimate:
# 1. Genome-wide inbreeding coefficient (FROH)
# 2. Runs of homozygosity (ROH)
# 3. Local heterozygosity

# ROHan works directly from BAM files and a reference genome.


#====================================================
# Software installation
#====================================================

conda activate mkondana_wgs

conda install -c bioconda rohan -y


#====================================================
# Create directories
#====================================================

mkdir -p rohan
mkdir -p rohan/K
mkdir -p rohan/R
mkdir -p rohan/T
mkdir -p rohan/RR
mkdir -p rohan/MEL


#====================================================
# Check reference index
#====================================================

samtools faidx reference/Rrattus.fa


#====================================================
# Run ROHan for a single individual
#====================================================

rohan \
--bam bam/K3.dedup.bam \
--ref reference/Rrattus.fa \
--out rohan/K/K3 \
--threads 16


# Output files:
# K3.global.txt
# K3.local.txt
# K3.roh.bed


#====================================================
# Run ROHan for all individuals
#====================================================

for sample in K1 K2 K3 K4 K5 K6 K7 K8 K9 K10
do
  rohan \
  --bam bam/${sample}.dedup.bam \
  --ref reference/Rrattus.fa \
  --out rohan/K/${sample} \
  --threads 16
done


for sample in R1 R2 R3 R4 R5 R6 R7 R8 R9 R10
do
  rohan \
  --bam bam/${sample}.dedup.bam \
  --ref reference/Rrattus.fa \
  --out rohan/R/${sample} \
  --threads 16
done


for sample in T1 T2 T3 T4 T5 T6 T7 T8 T9 T10
do
  rohan \
  --bam bam/${sample}.dedup.bam \
  --ref reference/Rrattus.fa \
  --out rohan/T/${sample} \
  --threads 16
done


for sample in RR1 RR2 RR3 RR4 RR5 RR6 RR7 RR8 RR9 RR10
do
  rohan \
  --bam bam/${sample}.dedup.bam \
  --ref reference/Rrattus.fa \
  --out rohan/RR/${sample} \
  --threads 16
done


for sample in MEL1 MEL2 MEL3 MEL4 MEL5 MEL6 MEL7 MEL8 MEL9 MEL10
do
  rohan \
  --bam bam/${sample}.dedup.bam \
  --ref reference/Rrattus.fa \
  --out rohan/MEL/${sample} \
  --threads 16
done


#====================================================
# Extract genome-wide inbreeding coefficient
#====================================================

grep 'Froh' rohan/K/*.global.txt > rohan/K/Froh_summary.txt

grep 'Froh' rohan/R/*.global.txt > rohan/R/Froh_summary.txt

grep 'Froh' rohan/T/*.global.txt > rohan/T/Froh_summary.txt

grep 'Froh' rohan/RR/*.global.txt > rohan/RR/Froh_summary.txt

grep 'Froh' rohan/MEL/*.global.txt > rohan/MEL/Froh_summary.txt


#====================================================
# Extract ROH lengths
#====================================================

awk '{sum += $3-$2} END {print sum}' rohan/K/K3.roh.bed

# Total ROH length for K3


#====================================================
# Calculate ROH statistics for all individuals
#====================================================

for file in rohan/K/*.roh.bed
do
  sample=$(basename $file .roh.bed)
  total=$(awk '{sum += $3-$2} END {print sum}' $file)
  echo -e "$sample\t$total"
done > rohan/K/K.roh_lengths.txt


# Repeat for other populations if needed.


#====================================================
# Local heterozygosity
#====================================================

head rohan/K/K3.local.txt


# Columns include:
# chromosome
# start
# end
# heterozygosity
# ROH/non-ROH state


#====================================================
# Final outputs
#====================================================

# Genome-wide inbreeding:
# rohan/*/Froh_summary.txt

# ROH coordinates:
# rohan/*/*.roh.bed

# Local heterozygosity:
# rohan/*/*.local.txt

# Total ROH length:
# rohan/*/*.roh_lengths.txt
