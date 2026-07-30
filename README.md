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

Raw Illumina FASTQ
        |
        v
FastQC + MultiQC
(Initial quality assessment)
        |
        v
fastp
(Adapter and quality trimming)
        |
        v
BWA-MEM2
(Map to Rattus rattus reference)
        |
        v
SAMtools sort + index
        |
        v
GATK MarkDuplicates
        |
        v
SAMtools flagstat + depth
(Coverage and mapping statistics)
        |
        v
bcftools mpileup + call
(Joint SNP calling)
        |
        v
SNP filtering
(QUAL, MAC, missingness)
        |
        v
+--------------------------------------------------+
|         Final BAMs + Filtered VCF                |
+--------------------------------------------------+
           /                              \
          /                                \
         v                                  v

BAM branch (genotype likelihoods)     VCF branch (hard genotypes)
---------------------------------------------------------------

ANGSD                                PLINK
   |                                    |
   +--> PCAngsd (PCA)                   +--> GONE (recent Ne)
   |
   +--> NGSadmix (admixture)
   |
   +--> realSFS (pairwise FST)
   |
   +--> thetaStat (π, θW, Tajima's D)
   |
   +--> PSMC (historical Ne)
   |
   +--> ROHan (FROH, ROH)

Filtered VCF
   |
   +--> FEIMS (gene flow)


Final outputs
-------------
- PCA plots
- Admixture plots
- Pairwise FST matrix
- Gene-flow network
- Nucleotide diversity (π)
- Watterson's θ
- Tajima's D
- Historical Ne trajectories
- Recent Ne estimates
- FROH values
- ROH coordinates

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
