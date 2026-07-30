#!/bin/bash

#====================================================
# Effective population size analyses
# MSMC (primary) | PSMC (optional) | GONE (recent Ne)
#====================================================

REF=reference/Rrattus.autosomes.fa
VCF=variants/variants.filtered.vcf.gz


#====================================================
# Coverage distribution
#====================================================

for sample in K1 K2 K3 K4
do
  samtools depth -a bam/${sample}.dedup.bam | \
  awk '{print $3}' > stats/${sample}.depth.txt
done


#====================================================
# Consensus generation
#====================================================

for sample in K1 K2 K3 K4
do
  samtools mpileup -C50 -uf ${REF} bam/${sample}.dedup.bam | \
  bcftools call -c | \
  vcfutils.pl vcf2fq -d 10 -D 60 > ${sample}.fq
done


#====================================================
# MSMC input
#====================================================

for sample in K1 K2 K3 K4
do
  generate_multihetsep.py ${sample}.fq > ${sample}.multihetsep.txt
done


#====================================================
# MSMC
#====================================================

msmc2 \
-I 0,1,2,3 \
-o K_MSMC \
K1.multihetsep.txt \
K2.multihetsep.txt \
K3.multihetsep.txt \
K4.multihetsep.txt


#====================================================
# Sensitivity analysis
#====================================================

msmc2 -t 8 -I 0,1,2,3 -o K_MSMC_t8 \
K1.multihetsep.txt K2.multihetsep.txt K3.multihetsep.txt K4.multihetsep.txt

msmc2 -t 16 -I 0,1,2,3 -o K_MSMC_t16 \
K1.multihetsep.txt K2.multihetsep.txt K3.multihetsep.txt K4.multihetsep.txt


#====================================================
# Optional PSMC
#====================================================

fq2psmcfa -q20 K1.fq > K1.psmcfa

psmc -N25 -t15 -r5 -p "4+25*2+4+6" \
-o K1.psmc K1.psmcfa


#====================================================
# GONE: Sinhagad
#====================================================

bcftools view -S metadata/sinhagad.samples ${VCF} -Oz -o K.filtered.vcf.gz

plink --vcf K.filtered.vcf.gz \
      --make-bed \
      --out K_population

GONE K_population


#====================================================
# GONE: Rajgad
#====================================================

bcftools view -S metadata/rajgad.samples ${VCF} -Oz -o R.filtered.vcf.gz

plink --vcf R.filtered.vcf.gz \
      --make-bed \
      --out R_population

GONE R_population


#====================================================
# GONE: Torna
#====================================================

bcftools view -S metadata/torna.samples ${VCF} -Oz -o T.filtered.vcf.gz

plink --vcf T.filtered.vcf.gz \
      --make-bed \
      --out T_population

GONE T_population


#====================================================
# GONE: Raireshwar
#====================================================

bcftools view -S metadata/raireshwar.samples ${VCF} -Oz -o RR.filtered.vcf.gz

plink --vcf RR.filtered.vcf.gz \
      --make-bed \
      --out RR_population

GONE RR_population
