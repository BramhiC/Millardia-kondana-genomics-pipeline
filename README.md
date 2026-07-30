# Millardia-kondana-genomics-pipeline
Population genetics pipeline for M.kondana
# Millardia kondana whole-genome population genomics pipeline

This repository contains the analysis workflow used for whole-genome population genomics of the endangered Kondana soft-furred rat (*Millardia kondana*) across four isolated high-elevation plateau populations in the northern Western Ghats, India.

---

## Study system

- **K** — Sinhagad
- **R** — Rajgad
- **T** — Torna
- **RR** — Raireshwar
- **MEL** — *Millardia meltada* (comparative population)

Illumina paired-end whole-genome sequencing data were generated for all individuals. Because no reference genome is currently available for *M. kondana*, reads were mapped to the **Rattus rattus** reference genome.

---

## Analysis objectives

1. **Current and historical effective population size**
2. **Gene flow and genetic differentiation between populations**
3. **Genetic diversity**
4. **Inbreeding coefficient and runs of homozygosity**

---

# Pipeline overview

```text
Raw FASTQ
   ↓
fastp
   ↓
FastQC / MultiQC
   ↓
BWA-MEM2
   ↓
SAMtools sort/index
   ↓
GATK MarkDuplicates
   ↓
bcftools variant calling
   ↓
SNP filtering (QUAL, MAC, missingness)
   ↓
Final filtered VCF
```

---

# Repository structure

```text
.
├── README.md
├── 01_sequence_processing.sh
├── 02_effective_population_size.sh
├── 03_population_structure_gene_flow.sh
├── 04_genetic_diversity.sh
├── 05_inbreeding_rohan.sh
├── reference/
├── metadata/
├── scripts/
└── results/
```

---

# Step 1 — Sequence processing

**File:** `01_sequence_processing.sh`

### Main software

- fastp
- FastQC
- MultiQC
- BWA-MEM2
- SAMtools
- GATK4
- bcftools

### Output

- Deduplicated BAM files
- Filtered SNP VCF

### Filtering

- SNPs only
- QUAL > 20
- MAC > 3
- Missingness < 20%

---

# Step 2 — Effective population size

**File:** `02_effective_population_size.sh`

### Historical Ne

- PSMC
- seqtk
- bcftools

### Contemporary Ne

- GONE
- PLINK

### Output

- PSMC demographic trajectories
- Recent Ne estimates for K, R, T and RR

---

# Step 3 — Population structure and gene flow

**File:** `03_population_structure_gene_flow.sh`

### Population structure

- ANGSD
- PCAngsd
- NGSadmix

### Genetic differentiation

- ANGSD
- realSFS

### Gene flow

- FEIMS

### Output

- PCA
- Admixture plots
- Pairwise FST
- Gene-flow network

---

# Step 4 — Genetic diversity

**File:** `04_genetic_diversity.sh`

### Software

- ANGSD
- realSFS
- thetaStat

### Statistics

- Nucleotide diversity (π)
- Watterson’s θ
- Tajima’s D
- Sliding-window diversity

### Comparative analysis

- *M. kondana* vs *M. meltada*

---

# Step 5 — Inbreeding and ROH

**File:** `05_inbreeding_rohan.sh`

### Software

- ROHan

### Statistics

- FROH
- ROH coordinates
- Local heterozygosity

### Output

- Genome-wide inbreeding coefficients
- ROH BED files
- Window-based heterozygosity

---

# Data flow

```text
BAM files
 ├── ANGSD → PCAngsd / NGSadmix / FST / Diversity
 ├── ROHan → Inbreeding / ROH
 └── bcftools → Filtered VCF → GONE / FEIMS
```

---

# Notes

- The workflow is designed for **moderate-coverage Illumina whole-genome data** from a **non-model species**.
- Genotype-likelihood approaches (ANGSD, ROHan) are preferred over hard-genotype methods where appropriate.
- Exact FEIMS command syntax may vary depending on the installed version.

---

# Citation

If this pipeline is used in publications or derivative analyses, please cite the relevant software packages and the associated thesis or manuscript describing the *Millardia kondana* dataset.
