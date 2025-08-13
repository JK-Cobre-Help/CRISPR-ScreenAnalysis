![GitHub Release](https://img.shields.io/github/v/release/KirklandLab/CRISPR-ScreenAnalysis)
![GitHub Release Date](https://img.shields.io/github/release-date/KirklandLab/CRISPR-ScreenAnalysis)
![GitHub repo size](https://img.shields.io/github/repo-size/KirklandLab/CRISPR-ScreenAnalysis)
![GitHub last commit](https://img.shields.io/github/last-commit/KirklandLab/CRISPR-ScreenAnalysis)
![GitHub contributors](https://img.shields.io/github/contributors/KirklandLab/CRISPR-ScreenAnalysis)
![GitHub Downloads (all assets, all releases)](https://img.shields.io/github/downloads/KirklandLab/CRISPR-ScreenAnalysis/total)
![GitHub commits since latest release](https://img.shields.io/github/commits-since/KirklandLab/CRISPR-ScreenAnalysis/latest)
[![DOI](https://zenodo.org/badge/1029909681.svg)](https://doi.org/10.5281/zenodo.16809383)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# CRISPR-ScreenAnalysis

![CRISPR-Screen](/images/CRISPR-Screen.png)
*(Image generated with DALL-E. OpenAI, 2025: Scientific data visualization: Genome Wide CRISPR Screens in bioinformatics)*

---

## 1) Project Description

**CRISPR-ScreenAnalysis** is an automated Snakemake workflow for pooled CRISPR screen analysis using MAGeCK, supporting both single-end and paired-end data. The pipeline performs FASTQ quality control, sgRNA counting with a user-provided library, automated or custom design matrix generation, and gene-level β score estimation via MAGeCK MLE. This workflow is designed for reproducibility and scalability in Slurm-managed HPC environments.

### Key Features

+ **Quality Control**
  + Runs **FastQC** on raw FASTQs and summarizes results with **MultiQC**
  + Works with both single-end and paired-end reads

+ **Optional FASTQ Staging**
  + Copies input FASTQ files from slower archive storage to a writable, fast-access location before analysis
  + Controlled by `stage_fastqs` in `config.yml`

+ **MAGeCK Count**
  + Reads directly from raw or staged (copied) FASTQs
  + Supports both paired-end and single-end modes
  + Assigns reads to sgRNAs from a user-provided library
  + Outputs both `.count.txt` and `.countsummary.txt` for all samples

+ **Automatic Design Matrix Generation**
  + Builds a valid MAGeCK MLE design matrix from `samples.csv`  
  + Includes:
    + Baseline column
    + One-hot encoded factors for experimental conditions
  + Option to bypass auto-build with a custom design matrix file

+ **MAGeCK MLE Analysis**
  + Estimates gene-level β scores using experimental factors defined in the design matrix
  + Accepts optional `control-sgRNA` file for normalization
  + Produces gene summary tables for ranking

+ **Reproducible Configuration**
  + All analysis parameters, file paths, and environment module versions stored in `config.yml`

---

## 2) Intended Use Case

This pipeline is intended for **researchers performing genome-wide or targeted CRISPR screens** who want:

+ Automated **QC**, **sgRNA counting**, and **statistical modeling**
+ Consistent output for downstream gene ranking and pathway analysis
+ A workflow that runs efficiently on **Slurm HPC systems**
+ Easy modification of experimental designs without altering code

---

## 3) Dependencies & Configuration

Tool versions and paths are defined in `config/config.yml`.  
Example key fields:

+ **samples_csv**: Path to `samples.csv` file  
+ **stage_fastqs**: `true`/`false` to enable staging  
+ **mageck_paired**: Set to `true` for paired-end data  
+ **sgRNA_library**: Path to sgRNA library file  
+ **counts_prefix**: Output name for count files  
+ **mle_enabled**: Enable or disable MAGeCK MLE step  
+ **mle_prefix**: Output name for MLE results  
+ **mle_control_sgrna**: Optional control sgRNA file  
+ **custom_design_matrix**: Path to user-provided design matrix file

**Example:**
```yaml
samples_csv: "config/samples.csv"
stage_fastqs: true
mageck_paired: false
sgRNA_library: "resources/calabreseA_library.txt"
counts_prefix: "sample1"
mle_enabled: true
mle_prefix: "results/mle/CalabreseA"
mle_control_sgrna: ""
custom_design_matrix: ""
```

---

## 3) Tools and Modules

---

## 5) Example Data

---

## 6) Explanation of `samples.csv`
**Required Columns**:  
+ sample – unique sample ID
+ fastq1 – path to R1 FASTQ
+ fastq2 – path to R2 FASTQ (or blank if single-end)
+ include_mle – whether to include this sample in MLE design matrix (true/false)
+ factor – experimental condition name (none if baseline)

**Example**:
```csv
sample,fastq1,fastq2,include_mle,factor
N1_D0,/path/N1_D0_R1.fastq.gz,/path/N1_D0_R2.fastq.gz,true,none
N2_D0,/path/N2_D0_R1.fastq.gz,/path/N2_D0_R2.fastq.gz,true,none
N3_D0,/path/N3_D0_R1.fastq.gz,/path/N3_D0_R2.fastq.gz,true,none
N1_D14_DMSO,/path/N1_D14_DMSO_R1.fastq.gz,/path/N1_D14_DMSO_R2.fastq.gz,true,D14_DMSO
N2_D14_DMSO,/path/N2_D14_DMSO_R1.fastq.gz,/path/N2_D14_DMSO_R2.fastq.gz,true,D14_DMSO
N3_D14_DMSO,/path/N3_D14_DMSO_R1.fastq.gz,/path/N3_D14_DMSO_R2.fastq.gz,true,D14_DMSO
N1_D14_PAC,/path/N1_D14_PAC_R1.fastq.gz,/path/N1_D14_PAC_R2.fastq.gz,true,D14_PAC
N2_D14_PAC,/path/N2_D14_PAC_R1.fastq.gz,/path/N2_D14_PAC_R2.fastq.gz,true,D14_PAC
N3_D14_PAC,/path/N3_D14_PAC_R1.fastq.gz,/path/N3_D14_PAC_R2.fastq.gz,true,D14_PAC
```

---

## 7) Output Overview
|  Category          | 	Output Location                         |
|--------------------|------------------------------------------|
| FastQC Reports     | `results/qc/fastqc/`                     |
| MultiQC Report     | `results/qc/multiqc/multiqc_report.html` |
| MAGeCK Counts      | `results/counts/*.count.txt`             |
| MAGeCK Count Summ. | `results/counts/*.countsummary.txt`      |
| Design Matrix      | `results/mle/design_matrix.txt`          |
| MLE Results        | `results/mle/*.gene_summary.txt`         |

---

## 8) Example Output Plots

---

## 9) Instructions to run on Slurm managed HPC  

9A. Download version controlled repository
```
git clone https://github.com/KirklandLab/CRISPR-ScreenAnalysis.git
```
9B. Load modules
```
module purge
module load slurm python/3.10 pandas/2.2.3 numpy/1.22.3 matplotlib/3.7.1
```
9C. Modify samples and config file
```
vim config/samples.csv
vim config/config.yml
```
9D. Dry Run
```
snakemake -npr
```

9E.  Make a DAG diagram
```
snakemake --dag | dot -Tpdf > dag.pdf
```

9F. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 999 --use-envmodules --rerun-incomplete --latency-wait 300 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {threads} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output}'"
```

---

## 10) Citation

+ **Boyd, K.A.** (2025). *CRISPR-ScreenAnalysis: A reproducible Snakemake workflow for pooled CRISPR screening data analysis*. *Zenodo*. https://doi.org/10.5281/zenodo.16809383  
+ [![DOI](https://zenodo.org/badge/1029909681.svg)](https://doi.org/10.5281/zenodo.16809383)

---

## 11) Authorship & Contributions

+ **Kevin A. Boyd** – Designed and implemented the Snakemake workflow for a Slurm-managed HPC environment, modularized the pipeline structure, implemented all processing steps, integrated peak consensus method, designed plots, and created the documentation.  
+ **Jacob Kirkland** – Principal Investigator; provided experimental data and validation of activation screen logic.
+ **Christopher L. Sansam** – Principal Investigator; provided additional experimental data and validation of knockout screen logic.  

This workflow was developed as part of a COBRE-funded collaborative effort. While the pipeline was built specifically for use within the Kirkland Lab, it is broadly applicable to CRISPR-Screen data analysis in other research settings.  

---

## 12) License

This project is licensed under the **Apache 2.0**. See the [LICENSE](LICENSE) file for details.  

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
