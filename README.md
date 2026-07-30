# Millardia-kondana-genomics-pipeline
Population genetics pipeline for M.kondana
# Millardia kondana whole-genome population genomics pipeline

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

## Pipeline overview

```text
Raw Illumina FASTQ files
        │
        ▼
Quality control
(FastQC, MultiQC)
        │
        ▼
Adapter and quality trimming
(fastp)
        │
        ▼
Read mapping to Rattus rattus reference
(BWA-MEM2)
        │
        ▼
BAM processing
(SAMtools sort/index, GATK MarkDuplicates)
        │
        ▼
Coverage and mapping statistics
(SAMtools)
        │
        ▼
Joint SNP calling
(bcftools)
        │
        ▼
Variant filtering
(SNPs, QUAL, MAC, missingness)
        │
        ▼
Final filtered VCF + deduplicated BAMs
        │
        ├───────────────┬────────────────┬─────────────────┬──────────────────┐
        ▼               ▼                ▼                 ▼                  ▼
Historical Ne       Recent Ne        Population        Genetic diversity    Inbreeding
(PSMC)              (GONE)           structure & FST   (ANGSD)              (ROHan)
                                        │
                                        ▼
                                   Gene flow
                                     (FEIMS)
```
## Software summary

| Category | Software |
|---|---|
| Quality control | FastQC, MultiQC, fastp |
| Read mapping | BWA-MEM2 |
| BAM processing | SAMtools, GATK4 |
| Variant calling | bcftools |
| Population structure | ANGSD, PCAngsd, NGSadmix |
| Genetic differentiation | realSFS, FEIMS |
| Genetic diversity | ANGSD, thetaStat |
| Historical Ne | PSMC |
| Recent Ne | GONE |
| Inbreeding and ROH | ROHan |
| Visualization | R (ggplot2, ComplexHeatmap) |
| Optional QC | Qualimap, IGV |

## Script descriptions

### 01_sequence_processing.sh
Quality control, adapter trimming, mapping to the *Rattus rattus* reference genome, BAM processing, joint SNP calling, and variant filtering.

### 02_effective_population_size.sh
Historical effective population size using PSMC and recent effective population size using GONE.

### 03_population_structure_gene_flow.sh
Population structure using PCAngsd and NGSadmix, pairwise FST estimation using ANGSD/realSFS, and gene-flow analysis using FEIMS.

### 04_genetic_diversity.sh
Estimation of nucleotide diversity (π), Watterson’s θ, Tajima’s D, and sliding-window diversity statistics using ANGSD.

### 05_inbreeding_rohan.sh
Estimation of genome-wide inbreeding coefficient (FROH), runs of homozygosity (ROH), and local heterozygosity using ROHan.
