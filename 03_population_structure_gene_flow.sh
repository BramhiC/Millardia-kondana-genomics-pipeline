#!/bin/bash

#====================================================
# Millardia kondana whole-genome sequencing pipeline
# Population structure, FST and gene flow
#====================================================

# PCA        : PCAngsd
# Admixture  : NGSadmix
# FST        : ANGSD + realSFS
# Gene flow  : FEIMS


#====================================================
# Software installation
#====================================================

conda activate mkondana_wgs

conda install -c bioconda angsd pcangsd ngsadmix -y


#====================================================
# Create directories
#====================================================

mkdir -p angsd
mkdir -p pca
mkdir -p admixture
mkdir -p fst
mkdir -p feims


#====================================================
# Create genotype likelihoods for all samples
#====================================================

angsd \
-b bamlist.txt \
-ref reference/Rrattus.fa \
-GL 1 \
-doGlf 2 \
-doMajorMinor 1 \
-doMaf 1 \
-minMaf 0.05 \
-minMapQ 30 \
-minQ 20 \
-out angsd/mkondana


# -GL 1          : SAMtools genotype likelihood model
# -doGlf 2       : write Beagle file
# -doMajorMinor  : infer major/minor alleles
# -doMaf         : estimate allele frequencies


# Output:
# angsd/mkondana.beagle.gz


#====================================================
# Principal Component Analysis (PCA)
#====================================================

pcangsd.py \
-beagle angsd/mkondana.beagle.gz \
-o pca/mkondana_pca \
-threads 16


# Output:
# pca/mkondana_pca.cov
# pca/mkondana_pca.eigenvec
# pca/mkondana_pca.eigenval


#====================================================
# Admixture analysis
#====================================================

for K in 1 2 3 4 5 6
do
  NGSadmix \
    -likes angsd/mkondana.beagle.gz \
    -K $K \
    -P 16 \
    -o admixture/K${K}
done


# Choose the K with the best likelihood.


#====================================================
# Create population BAM lists
#====================================================

grep 'Sinhagad' metadata/sample_metadata.tsv | cut -f1 | sed 's#^#bam/#; s#$#.dedup.bam#' > fst/K.bamlist

grep 'Rajgad' metadata/sample_metadata.tsv | cut -f1 | sed 's#^#bam/#; s#$#.dedup.bam#' > fst/R.bamlist

grep 'Torna' metadata/sample_metadata.tsv | cut -f1 | sed 's#^#bam/#; s#$#.dedup.bam#' > fst/T.bamlist

grep 'Raireshwar' metadata/sample_metadata.tsv | cut -f1 | sed 's#^#bam/#; s#$#.dedup.bam#' > fst/RR.bamlist


#====================================================
# Estimate allele frequency likelihoods
#====================================================

angsd -b fst/K.bamlist  -ref reference/Rrattus.fa -GL 1 -doSaf 1 -anc reference/Rrattus.fa -out fst/K

angsd -b fst/R.bamlist  -ref reference/Rrattus.fa -GL 1 -doSaf 1 -anc reference/Rrattus.fa -out fst/R

angsd -b fst/T.bamlist  -ref reference/Rrattus.fa -GL 1 -doSaf 1 -anc reference/Rrattus.fa -out fst/T

angsd -b fst/RR.bamlist -ref reference/Rrattus.fa -GL 1 -doSaf 1 -anc reference/Rrattus.fa -out fst/RR


# .saf files contain site allele frequency likelihoods.


#====================================================
# Pairwise FST: Sinhagad vs Rajgad
#====================================================

realSFS fst/K.saf.idx fst/R.saf.idx > fst/K_R.sfs

realSFS fst index fst/K.saf.idx fst/R.saf.idx \
-sfs fst/K_R.sfs \
-fstout fst/K_R

realSFS fst stats fst/K_R.fst.idx > fst/K_R.txt


#====================================================
# Pairwise FST: Sinhagad vs Torna
#====================================================

realSFS fst/K.saf.idx fst/T.saf.idx > fst/K_T.sfs

realSFS fst index fst/K.saf.idx fst/T.saf.idx \
-sfs fst/K_T.sfs \
-fstout fst/K_T

realSFS fst stats fst/K_T.fst.idx > fst/K_T.txt


#====================================================
# Pairwise FST: Sinhagad vs Raireshwar
#====================================================

realSFS fst/K.saf.idx fst/RR.saf.idx > fst/K_RR.sfs

realSFS fst index fst/K.saf.idx fst/RR.saf.idx \
-sfs fst/K_RR.sfs \
-fstout fst/K_RR

realSFS fst stats fst/K_RR.fst.idx > fst/K_RR.txt


#====================================================
# Pairwise FST: Rajgad vs Torna
#====================================================

realSFS fst/R.saf.idx fst/T.saf.idx > fst/R_T.sfs

realSFS fst index fst/R.saf.idx fst/T.saf.idx \
-sfs fst/R_T.sfs \
-fstout fst/R_T

realSFS fst stats fst/R_T.fst.idx > fst/R_T.txt


#====================================================
# Pairwise FST: Rajgad vs Raireshwar
#====================================================

realSFS fst/R.saf.idx fst/RR.saf.idx > fst/R_RR.sfs

realSFS fst index fst/R.saf.idx fst/RR.saf.idx \
-sfs fst/R_RR.sfs \
-fstout fst/R_RR

realSFS fst stats fst/R_RR.fst.idx > fst/R_RR.txt


#====================================================
# Pairwise FST: Torna vs Raireshwar
#====================================================

realSFS fst/T.saf.idx fst/RR.saf.idx > fst/T_RR.sfs

realSFS fst index fst/T.saf.idx fst/RR.saf.idx \
-sfs fst/T_RR.sfs \
-fstout fst/T_RR

realSFS fst stats fst/T_RR.fst.idx > fst/T_RR.txt


# Weighted FST values are written to:
# fst/*.txt


#====================================================
# Prepare population map for FEIMS
#====================================================

cut -f1,2 metadata/sample_metadata.tsv > feims/population_map.txt


#====================================================
# Gene flow analysis (FEIMS)
#====================================================

feims \
--vcf variants/variants.filtered.vcf.gz \
--popmap feims/population_map.txt \
--out feims/mkondana


# Output:
# migration edges / gene-flow network


#====================================================
# Final outputs
#====================================================

# PCA:
# pca/

# Admixture:
# admixture/

# FST:
# fst/*.txt

# Gene flow:
# feims/
