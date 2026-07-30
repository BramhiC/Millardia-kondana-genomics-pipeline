#!/bin/bash

#====================================================
# Millardia kondana whole-genome sequencing pipeline
# Genetic diversity
#====================================================

# Diversity statistics are estimated using ANGSD.
# This includes:
# 1. Site Frequency Spectrum (SFS)
# 2. Nucleotide diversity (pi)
# 3. Watterson's theta
# 4. Tajima's D

# M. meltada is analyzed separately for comparison.


#====================================================
# Software installation
#====================================================

conda activate mkondana_wgs

conda install -c bioconda angsd -y


#====================================================
# Create directories
#====================================================

mkdir -p diversity/K
mkdir -p diversity/R
mkdir -p diversity/T
mkdir -p diversity/RR
mkdir -p diversity/MEL


#====================================================
# Create BAM lists
#====================================================

grep 'Sinhagad' metadata/sample_metadata.tsv | cut -f1 | sed 's#^#bam/#; s#$#.dedup.bam#' > diversity/K/K.bamlist

grep 'Rajgad' metadata/sample_metadata.tsv | cut -f1 | sed 's#^#bam/#; s#$#.dedup.bam#' > diversity/R/R.bamlist

grep 'Torna' metadata/sample_metadata.tsv | cut -f1 | sed 's#^#bam/#; s#$#.dedup.bam#' > diversity/T/T.bamlist

grep 'Raireshwar' metadata/sample_metadata.tsv | cut -f1 | sed 's#^#bam/#; s#$#.dedup.bam#' > diversity/RR/RR.bamlist

grep 'Meltada' metadata/sample_metadata.tsv | cut -f1 | sed 's#^#bam/#; s#$#.dedup.bam#' > diversity/MEL/MEL.bamlist


#====================================================
# Estimate SAF for each population
#====================================================

angsd -b diversity/K/K.bamlist   -ref reference/Rrattus.fa -anc reference/Rrattus.fa -GL 1 -doSaf 1 -minMapQ 30 -minQ 20 -out diversity/K/K

angsd -b diversity/R/R.bamlist   -ref reference/Rrattus.fa -anc reference/Rrattus.fa -GL 1 -doSaf 1 -minMapQ 30 -minQ 20 -out diversity/R/R

angsd -b diversity/T/T.bamlist   -ref reference/Rrattus.fa -anc reference/Rrattus.fa -GL 1 -doSaf 1 -minMapQ 30 -minQ 20 -out diversity/T/T

angsd -b diversity/RR/RR.bamlist -ref reference/Rrattus.fa -anc reference/Rrattus.fa -GL 1 -doSaf 1 -minMapQ 30 -minQ 20 -out diversity/RR/RR

angsd -b diversity/MEL/MEL.bamlist -ref reference/Rrattus.fa -anc reference/Rrattus.fa -GL 1 -doSaf 1 -minMapQ 30 -minQ 20 -out diversity/MEL/MEL


#====================================================
# Estimate folded Site Frequency Spectrum (SFS)
#====================================================

realSFS diversity/K/K.saf.idx   > diversity/K/K.sfs
realSFS diversity/R/R.saf.idx   > diversity/R/R.sfs
realSFS diversity/T/T.saf.idx   > diversity/T/T.sfs
realSFS diversity/RR/RR.saf.idx > diversity/RR/RR.sfs
realSFS diversity/MEL/MEL.saf.idx > diversity/MEL/MEL.sfs


#====================================================
# Calculate theta statistics
#====================================================

angsd -b diversity/K/K.bamlist \
-ref reference/Rrattus.fa \
-anc reference/Rrattus.fa \
-GL 1 -doSaf 1 -pest diversity/K/K.sfs \
-doThetas 1 -doCounts 1 -out diversity/K/K

angsd -b diversity/R/R.bamlist \
-ref reference/Rrattus.fa \
-anc reference/Rrattus.fa \
-GL 1 -doSaf 1 -pest diversity/R/R.sfs \
-doThetas 1 -doCounts 1 -out diversity/R/R

angsd -b diversity/T/T.bamlist \
-ref reference/Rrattus.fa \
-anc reference/Rrattus.fa \
-GL 1 -doSaf 1 -pest diversity/T/T.sfs \
-doThetas 1 -doCounts 1 -out diversity/T/T

angsd -b diversity/RR/RR.bamlist \
-ref reference/Rrattus.fa \
-anc reference/Rrattus.fa \
-GL 1 -doSaf 1 -pest diversity/RR/RR.sfs \
-doThetas 1 -doCounts 1 -out diversity/RR/RR

angsd -b diversity/MEL/MEL.bamlist \
-ref reference/Rrattus.fa \
-anc reference/Rrattus.fa \
-GL 1 -doSaf 1 -pest diversity/MEL/MEL.sfs \
-doThetas 1 -doCounts 1 -out diversity/MEL/MEL


#====================================================
# Genome-wide diversity statistics
#====================================================

thetaStat do_stat diversity/K/K.thetas.idx > diversity/K/K.theta_stats.txt

thetaStat do_stat diversity/R/R.thetas.idx > diversity/R/R.theta_stats.txt

thetaStat do_stat diversity/T/T.thetas.idx > diversity/T/T.theta_stats.txt

thetaStat do_stat diversity/RR/RR.thetas.idx > diversity/RR/RR.theta_stats.txt

thetaStat do_stat diversity/MEL/MEL.thetas.idx > diversity/MEL/MEL.theta_stats.txt


# Output includes:
# Tajima's D
# Pairwise nucleotide diversity (pi)
# Watterson's theta
# Number of segregating sites


#====================================================
# Sliding-window diversity
#====================================================

thetaStat do_stat diversity/K/K.thetas.idx -win 50000 -step 10000 > diversity/K/K.windowed.txt

thetaStat do_stat diversity/R/R.thetas.idx -win 50000 -step 10000 > diversity/R/R.windowed.txt

thetaStat do_stat diversity/T/T.thetas.idx -win 50000 -step 10000 > diversity/T/T.windowed.txt

thetaStat do_stat diversity/RR/RR.thetas.idx -win 50000 -step 10000 > diversity/RR/RR.windowed.txt

thetaStat do_stat diversity/MEL/MEL.thetas.idx -win 50000 -step 10000 > diversity/MEL/MEL.windowed.txt


# -win  50000 : 50 kb window
# -step 10000 : 10 kb step


#====================================================
# Final outputs
#====================================================

# Genome-wide diversity:
# diversity/*/*.theta_stats.txt

# Sliding-window diversity:
# diversity/*/*.windowed.txt

# Site frequency spectra:
# diversity/*/*.sfs
