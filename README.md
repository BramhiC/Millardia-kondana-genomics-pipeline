# Millardia-kondana-genomics-pipeline
Population genetics pipeline for M.kondana

This repository contains the analysis workflow to use for whole-genome population genomics of the endangered Kondana soft-furred rat (*Millardia kondana*) across four isolated high-elevation plateau populations in the northern Western Ghats, India.

---
## Study populations

- **K** — Sinhagad
- **R** — Rajgad
- **T** — Torna
- **RR** — Raireshwar
- **MEL** — *Millardia meltada* (comparative population)

- ## Analysis objectives

1. **Current and historical effective population size**
2. **Gene flow and genetic differentiation between populations**
3. **Genetic diversity**
4. **Inbreeding coefficient and runs of homozygosity**

## Software summary

| Step | Software |
|---|---|
| Quality control | FastQC, MultiQC, Qualimap |
| Adapter and quality trimming | fastp |
| Read mapping | BWA-MEM2 |
| BAM sorting and indexing | SAMtools |
| Duplicate marking | GATK4 MarkDuplicates |
| Variant calling and filtering | bcftools |
| Population structure | ANGSD, PCAngsd, NGSadmix |
| Genetic differentiation | realSFS, FEIMS |
| Genetic diversity | ANGSD, thetaStat |
| Historical effective population size | MSMC, PSMC(optional) |
| Recent effective population size | PLINK, GONE |
| Inbreeding and ROH | ROHan |
| Visualization | R (ggplot2, ComplexHeatmap) |
| Optional QC | IGV |

## Script descriptions

### 01_sequence_processing.sh
Quality control, adapter trimming, mapping to the *Rattus rattus* reference genome, BAM processing, joint SNP calling, and variant filtering.

### 02_effective_population_size.sh
Historical effective population size using MSMC and recent effective population size using GONE.

### 03_population_structure_gene_flow.sh
Population structure using PCAngsd and NGSadmix, pairwise FST estimation using ANGSD/realSFS, and gene-flow analysis using FEIMS.

### 04_genetic_diversity.sh
Estimation of nucleotide diversity (π), Watterson’s θ, Tajima’s D, and sliding-window diversity statistics using ANGSD.

### 05_inbreeding_rohan.sh
Estimation of genome-wide inbreeding coefficient (FROH), runs of homozygosity (ROH), and local heterozygosity using ROHan.
