#!/bin/bash

#====================================================
# Millardia kondana whole-genome sequencing pipeline
# Historical and contemporary effective population size
#====================================================

# Historical Ne  : PSMC
# Current Ne     : GONE


#====================================================
# Software installation
#====================================================

conda activate mkondana_wgs

conda install -c bioconda psmc seqtk -y


#====================================================
# Create directories
#====================================================

mkdir -p psmc
mkdir -p gone


#====================================================
# Select high-coverage samples
#====================================================

# Choose the highest-coverage individuals from each population
# based on stats/*.mean_coverage.txt

ls stats/*.mean_coverage.txt


# Example selected samples:
# K3   K18
# R7   R21
# T5   T19
# RR4  RR15


#====================================================
# Generate consensus FASTQ for PSMC
#====================================================

bcftools mpileup \
-f reference/Rrattus.fa \
-q 20 -Q 20 \
bam/K3.dedup.bam | \
bcftools call -c | \
vcfutils.pl vcf2fq -d 5 -D 100 > psmc/K3.fq


# -d 5   : minimum coverage
# -D 100 : maximum coverage


#====================================================
# Convert FASTQ to PSMC format
#====================================================

fq2psmcfa -q20 psmc/K3.fq > psmc/K3.psmcfa


#====================================================
# Run PSMC
#====================================================

psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
-o psmc/K3.psmc \
psmc/K3.psmcfa


#====================================================
# Bootstrap analysis
#====================================================

splitfa psmc/K3.psmcfa > psmc/K3.split.psmcfa

for i in {1..100}
do
  seqtk sample -s$i psmc/K3.split.psmcfa 100 > psmc/bootstrap.$i.psmcfa

  psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" \
  -o psmc/bootstrap.$i.psmc \
  psmc/bootstrap.$i.psmcfa
done


#====================================================
# Plot historical Ne
#====================================================

psmc_plot.pl -u 2.5e-8 -g 0.5 psmc/K3 psmc/K3.psmc


# -u : mutation rate per site per generation
# -g : generation time (years)


# Repeat the same workflow for:
# K18 R7 R21 T5 T19 RR4 RR15


#====================================================
# Contemporary effective population size (GONE)
#====================================================

# Create population sample lists

grep "Sinhagad" metadata/sample_metadata.tsv | cut -f1 > K.list
grep "Rajgad" metadata/sample_metadata.tsv | cut -f1 > R.list
grep "Torna" metadata/sample_metadata.tsv | cut -f1 > T.list
grep "Raireshwar" metadata/sample_metadata.tsv | cut -f1 > RR.list


#====================================================
# Extract population-specific VCFs
#====================================================

bcftools view -S K.list variants/variants.filtered.vcf.gz \
-Oz -o variants/K.vcf.gz

bcftools view -S R.list variants/variants.filtered.vcf.gz \
-Oz -o variants/R.vcf.gz

bcftools view -S T.list variants/variants.filtered.vcf.gz \
-Oz -o variants/T.vcf.gz

bcftools view -S RR.list variants/variants.filtered.vcf.gz \
-Oz -o variants/RR.vcf.gz


#====================================================
# Convert VCF to PLINK format
#====================================================

plink --vcf variants/K.vcf.gz --make-bed --out gone/K
plink --vcf variants/R.vcf.gz --make-bed --out gone/R
plink --vcf variants/T.vcf.gz --make-bed --out gone/T
plink --vcf variants/RR.vcf.gz --make-bed --out gone/RR


# --make-bed creates:
# .bed : binary genotypes
# .bim : SNP information
# .fam : sample information


#====================================================
# Run GONE
#====================================================

cd gone

bash script_GONE.sh K
bash script_GONE.sh R
bash script_GONE.sh T
bash script_GONE.sh RR

cd ..


#====================================================
# Final outputs
#====================================================

# Historical Ne:
# psmc/*.psmc
# psmc/*.pdf
# psmc/*.png

# Contemporary Ne:
# gone/K*
# gone/R*
# gone/T*
# gone/RR*
