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

**CRISPR-ScreenAnalysis** is a fully automated Snakemake workflow for analyzing pooled CRISPR screens using **MAGeCK** that supports both single-end and paired-end sequencing data. The pipeline performs quality control on raw FASTQ files, counts sgRNA reads using a user-provided library, generates an experimental design matrix (automatically or from a custom file), and estimates gene-level β scores via MAGeCK MLE. This workflow is designed for reproducibility, scalability, and allows for easy user defined modifications of configuration settings in a Slurm-managed HPC environment.

### Key Features

+ **Quality Control**
  + Runs **FastQC** on raw FASTQs and compiles a summary report with **MultiQC**
  + Works with both single-end and paired-end sequencing

+ **Optional FASTQ Staging**
  + Copies input FASTQ files from slower archive storage to a writable, fast-access location before processing
  + Controlled by `stage_fastqs` option in `config.yml`

+ **MAGeCK Count**
  + Processes FASTQs directly or from staged copies
  + Supports both paired-end and single-end modes
  + Assigns reads to sgRNAs from a user-provided library file
  + Outputs both `.count.txt` and `.countsummary.txt` for all samples

+ **Automatic Design Matrix Generation**
  + Builds a valid **MAGeCK MLE** design matrix from `samples.csv`  
  + Includes:
    + Baseline column
    + One-hot encoded factors for experimental conditions
  + Option to bypass auto-build with a custom design matrix file

+ **MAGeCK MLE Analysis**
  + Estimates gene-level β scores using experimental factors defined in the design matrix
  + Accepts optional `control-sgRNA` file for normalization
  + Produces gene summary tables and ranking outputs

+ **Comprehensive Plotting & Visualization**
  + Integrates **MAGeCKFlute** for downstream visualization of MLE results
  + Generates publication-ready plots including:
  + β score histograms and normalization comparisons
  + β vs β scatterplots (with and without top-gene labels)
  + Volcano plots (effect size vs significance)
  + QC plots for Gini index, zero-count sgRNAs, and mapping rates
  + Produces selection tables of significantly enriched or depleted genes, filtered by configurable FDR threshold in the `config.yml`

+ **Reproducible Configuration**
  + All file paths, parameters, and tool versions are controlled in a single `config.yml`
  + Easily adapted for both CRISPR activation and knockout screens

---

## 2) Intended Use Case

This workflow is intended for researchers performing **genome-wide** or **targeted CRISPR screens** who need an end-to-end, reproducible solution for data processing, statistical modeling, and visualization. It is suited for use in Slurm-managed HPC environments where scalability and efficiency are essential.

It is ideal for researchers who want:  

+ Automated **QC**, **sgRNA counting**, **design matrix generation**, and **MAGeCK MLE statistical modeling**  
+ Consistent, publication-ready output for **gene ranking**, **pathway enrichment**, and **screen performance assessment**  
+ Integrated visualizations (via MAGeCKFlute) to quickly interpret screen results, including **β score distributions**, **volcano plots**, **scatter plots**, **QC summaries**, and **pathway analysis**  
+ A workflow that can efficiently handle **single-end** or **paired-end** sequencing from small targeted libraries to large genome-wide screens 
+ Easy modification of experimental designs and plotting parameters without altering code  

---

## 3) Dependencies & Configuration

All tool versions, paths, configurations, and parameters are defined in `config/config.yml`.  

### Key Configuration Fields
| Field                  | Description                                                  | Example                                         |
|------------------------|--------------------------------------------------------------|-------------------------------------------------|
| `samples_csv`          | Path to `samples.csv` with all input samples and metadata    | `config/samples.csv`                            |
| `stage_fastqs`         | Copy FASTQs to local storage before processing               | `true`                                          |
| `mageck_paired`        | Set to `true` for paired-end sequencing                      | `false`                                         |
| `sgRNA_library`        | Path to sgRNA reference library                              | `resources/calabreseA_library.txt`              |
| `counts_prefix`        | Output name for count files                                  | `sample1`                                       |
| `mle_enabled`          | Enable or disable MAGeCK MLE step                            | `true`                                          |
| `mle_prefix`           | Output name for MLE results                                  | `CalabreseA`                                    |
| `mle_control_sgrna`    | Path to optional control sgRNA file (blank for none)         | `resources/calabreseA_nonTargeting.txt` or `""` |
| `custom_design_matrix` | Path to user-supplied design matrix (skip auto-build)        | `resources/example_design_matrix.txt` or `""`   |
| `proj_name`            | Project name prefix for MAGeCKFlute outputs                  | `CalabreseA`                                    |
| `organism`             | Organism code for MAGeCKFlute enrichment analysis            | `hsa`                                           |
| `norm_method`          | Normalization method for MAGeCKFlute                         | `none`, `cc`, or `loess`                        |
| `fdr_threshold`        | FDR cutoff for significant genes in MAGeCKFlute plots/tables | `0.10`                                          |

Example `config.yml`:
```yaml
samples_csv: "config/samples.csv"
stage_fastqs: true
mageck_paired: false
sgRNA_library: "resources/calabreseA_library.txt"
counts_prefix: "sample1"
mle_enabled: true
mle_prefix: "CalabreseA"
mle_control_sgrna: ""
custom_design_matrix: ""
proj_name: "CalabreseA"
organism: "hsa"
norm_method: "none"
fdr_threshold: 0.10
```

---

## 4) Tools and Modules

The pipeline requires both **command-line tools** and **R packages** to run all steps, including plotting.  

### Core Command-Line Tools
| Tool / Module           | Purpose                                        | Example Version |
|-------------------------|------------------------------------------------|-----------------|
| `fastqc`                | FASTQ quality control                          | 0.12.1          |
| `multiqc`               | Summarizes QC results                          | 1.21            |
| `mageck`                | sgRNA counting & MLE analysis                  | 0.5.9.2         |
| `apptainer`             | Containerized execution (optional)             | 1.3.6           |
| `R` + `bioconductor`    | Required for plotting and MAGeCKFlute          | 4.4.1           |

### Required R/Bioconductor Packages
These must be installed in the R environment specified in your `config.yml` for plotting and MAGeCKFlute integration:  

+ **MAGeCKFlute**: Downstream visualization of MAGeCK results
+ **cowplot**: Plot arrangement
+ **ggplot2**: General plotting
+ **ggrepel**: Non-overlapping gene labels in scatter plots
+ **rlang**: Tidy evaluation helpers
+ **dplyr**: Data manipulation
+ **tibble**: Data frame handling

When using environment modules, ensure the R module loads with Bioconductor support (e.g., `bioconductor/3.19 gsl`).

Example `config.yml`:
```yaml
fastqc: "fastqc/0.12.1"
multiqc: "multiqc/1.21"
mageck: "mageck/0.5.9.2"
apptainer: "apptainer/1.3.6"
R: "R/4.4.1-mkl"
bioconductor: "bioconductor/3.19 gsl"
```

---

## 5) Example Data

A minimal test dataset can be placed in a `resources/` folder (not included currently). Update `samples.csv` to point to these FASTQs for a quick test run. Once confirmed, replace with your personal CRISPR-Screen data.

---

## 6) Explanation of `samples.csv`

The `samples.csv` file defines all input samples and their associated metadata.  
**Correct formatting of this file is essential** as any errors in column names, order, or values will cause the workflow to fail.  

### Required Columns
| Column         | Description |
|----------------|-------------|
| `sample`       | Unique short name for the sample (no spaces, special characters, or duplicates). Used in file naming and output tracking. |
| `fastq1`       | Full path to the R1 FASTQ file. Must point to an existing file readable by the compute environment. |
| `fastq2`       | Full path to the R2 FASTQ file (paired-end) or left blank if single-end. |
| `include_mle`  | Boolean flag (`true` or `false`) indicating whether this sample should be included in the MAGeCK MLE analysis and design matrix. |
| `factor`       | Experimental condition label for the sample. Use `none` for baseline/control samples; other values should be short condition names (no spaces). |

### Rules for Valid Entries
+ **Unique `sample` values** – duplicate names will overwrite outputs and break downstream steps. Do not use the exact same name as factor.    
+ **Case-sensitive paths** – ensure `fastq1` and `fastq2` match the exact case of filenames on disk.  
+ **Boolean formatting** – `include_mle` must be one of: `true`, `false`, `1`, `0`, `yes`, `no` (case-insensitive).  
+ **Factor naming** – avoid spaces; use underscores or short identifiers (e.g., `D4_DMSO`, `D14_PAC`, etc.). Do not use the exact same name as sample name.  
+ **Baseline handling** – only one factor should be set to `none`; all other conditions will be encoded relative to this baseline in the design matrix.  

**Important:** The `factor` value must not be identical to the `sample` name to avoid conflicts in downstream analysis.

### Example: Paired-End Design
**Note**: Even if paired-end file paths are provided the workflow can still be run as single-end.  

```csv
sample,fastq1,fastq2,include_mle,factor
N1_D0,/path/to/N1_D0_S19_R1_001.fastq.gz,/path/to/N1_D0_S19_R2_001.fastq.gz,true,none
N2_D0,/path/to/N2_D0_S23_R1_001.fastq.gz,/path/to/N2_D0_S23_R2_001.fastq.gz,true,none
N3_D0,/path/to/N3_D0_S27_R1_001.fastq.gz,/path/to/N3_D0_S27_R2_001.fastq.gz,true,none
N1_D14_DMSO,/path/to/N1_D14_DMSO_S21_R1_001.fastq.gz,/path/to/N1_D14_DMSO_S21_R2_001.fastq.gz,true,D14_DMSO
N2_D14_DMSO,/path/to/N2_D14_DMSO_S25_R1_001.fastq.gz,/path/to/N2_D14_DMSO_S25_R2_001.fastq.gz,true,D14_DMSO
N3_D14_DMSO,/path/to/N3_D14_DMSO_S29_R1_001.fastq.gz,/path/to/N3_D14_DMSO_S29_R2_001.fastq.gz,true,D14_DMSO
N1_D14_PAC,/path/to/N1_D14_PAC_S22_R1_001.fastq.gz,/path/to/N1_D14_PAC_S22_R2_001.fastq.gz,true,D14_PAC
N2_D14_PAC,/path/to/N2_D14_PAC_S26_R1_001.fastq.gz,/path/to/N2_D14_PAC_S26_R2_001.fastq.gz,true,D14_PAC
N3_D14_PAC,/path/to/N3_D14_PAC_S30_R1_001.fastq.gz,/path/to/N3_D14_PAC_S30_R2_001.fastq.gz,true,D14_PAC
```

### How This File is Used in the Workflow
+ **QC Steps** – `fastq1`/`fastq2` are passed directly to FastQC and MultiQC.  
+ **MAGeCK Count** – `sample` names become the `--sample-label` list; FASTQ paths are passed to `--fastq` / `--fastq-2`.  
+ **Design Matrix Building** – Only rows with `include_mle = true` are included. The `factor` column determines grouping for β score estimation.  
+ **Plotting** – `factor` labels are used to determine control vs treatment in MAGeCKFlute visualizations. The first `factor` provided will be set as the control (x-axis) while the last `factor` provided will be the treatment (y-axis) of the β score plots.  

### Troubleshooting Common Errors
| Problem | Cause | Fix |
|---------|-------|-----|
| Workflow fails with `KeyError` on a sample name | Duplicate `sample` values in CSV | Ensure each sample name is unique |
| `No such file or directory` error during FastQC or MAGeCK count | `fastq1`/`fastq2` paths are incorrect or have wrong capitalization | Verify file paths and match exact case on disk |
| Design matrix is missing samples | `include_mle` set to `false` or mistyped | Ensure intended samples have `true` in `include_mle` |
| MAGeCK MLE produces only one factor column | All `factor` values set to `none` or baseline | Assign correct condition labels for treatments |
| Plotting script fails to find FDR columns | Incorrect `factor` names or mismatch between design matrix and MAGeCK output | Check `factor` column values match experimental conditions exactly |

---

## 7) Output Overview  

| Category | Output Location | Description |
|----------|-----------------|-------------|
| **FastQC Reports** | `results/qc/fastqc/` | Per-sample HTML and ZIP files from FastQC (`*_R1_fastqc.html/.zip`, `*_R2_fastqc.html/.zip`) |
| **MultiQC Report** | `results/qc/multiqc/multiqc_report.html` | Aggregated QC report summarizing all FastQC results |
| **MAGeCK Counts** | `results/counts/*.count.txt` | Main sgRNA count table from `mageck count` |
| **MAGeCK Count Summary** | `results/counts/*.countsummary.txt` | Summary statistics from `mageck count` |
| **Design Matrix** | `results/mle/design_matrix.txt` *(or custom path)* | Auto-generated or user-supplied MAGeCK MLE design matrix |
| **MLE-Dependent Outputs** | *-* | **All outputs below are only generated if** `mle_enabled: true` in `config.yml` |
| **MLE Results** | `results/mle/*.gene_summary.txt` | Gene-level β scores and FDRs from `mageck mle` |
| **Selection Table** | `results/plots/selection_table.norm_<method>.tsv` | Table of positive/negative selected genes filtered by FDR threshold |
| **Significant Hits Table** | `results/plots/<proj_name>_sig_hits_FDR_<thr>.tsv` | Final ranked list of significant genes based on β score and FDR |
| **QC Plots** | `results/plots/qc_gini.png`, `qc_zero_counts.png`, `qc_maprates.png` | Plots showing screen quality metrics (Gini index, zero counts, mapping rates) |
| **β Score Plots** | `results/plots/beta_hist.png`, `beta_norm.png`, `beta_scatter.png`, `beta_scatter_labeled.png` | Histograms, normalization comparisons, scatter plots (with and without labels) |
| **Volcano Plot** | `results/plots/volcano_diff_vs_neglog10FDR.png`, `results/plots/custom_volcano_diff_vs_neglog10FDR.png` | Effect size vs statistical significance for treatments |
| **MAGeCKFlute Directory** | `results/plots/MAGeCKFlute_<proj_name>/` | MAGeCKFlute output including enrichment plots and pathway analysis results |

---

## 8) Example Output Plots

If `mle_enabled: true`, the workflow generates a variety of plots for visualizing CRISPR screen results and quality metrics.  
All plots are saved in `results/plots/` unless otherwise specified.

| Plot | Example Filename | Description |
|------|------------------|-------------|
| **QC Plot: Mapping Rates** | `qc_maprates.png` | Mapping rate for each sample from the count summary. |
| **QC Plot: Zero Count sgRNAs** | `qc_zero_counts.png` | Bar plot showing the number of sgRNAs with zero reads per sample. |
| **QC Plot: Gini Index** | `qc_gini.png` | Bar plot showing Gini index for each sample (evenness of sgRNA distribution). |
| **β Score Histograms** | `beta_hist.png` | Side-by-side histograms of β score distributions for control and treatment samples. |
| **β Score Normalization Comparison** | `beta_norm.png` | Density plots comparing raw β scores with cell-cycle and loess normalization methods. |
| **β vs β Scatter Plot** | `beta_scatter.png` | Scatter plot comparing β scores for control (x-axis) vs treatment (y-axis), colored by FDR. |
| **β vs β Scatter Plot (Labeled)** | `beta_scatter_labeled.png` | Same as above, but with top-ranked genes labeled using `ggrepel`. |
| **Volcano Plot** | `volcano_diff_vs_neglog10FDR.png` | Effect size (β score difference) vs -log10(FDR) for treatment samples. |
| **Custom Volcano Plot** | `custom_volcano_diff_vs_neglog10FDR.png` | Effect size (β score difference) vs -log10(FDR) for treatment samples with 20 lowest treatment FDR labeled. |

**MAGeCKFlute Directory**  
+ **Note**: The workflow also generates a directory:  
  + `results/plots/MAGeCKFlute_<proj_name>/`  
  + This folder contains additional pathway enrichment plots, QC metrics, and summary visualizations created by **MAGeCKFlute**. These are not shown here in full due to the large number of outputs, but they provide further insight into screen performance and biological interpretation.

Below are example plots generated by this pipeline.  

| 1. **QC: Mapping Rates**                                                              | 2. **QC: Zero Count sgRNAs**                                                               |
| :-----------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------: |
| <img src="/images/qc_maprates.png" width="300">                                       | <img src="/images/qc_zero_counts.png" width="300">                                         |
| *Per-sample mapping rates from count summary*                                         | *Number of sgRNAs with zero reads per sample*                                              |

| 3. **QC: Gini Index**                                                                 |
| :-----------------------------------------------------------------------------------: |
| <img src="/images/qc_gini.png" width="300">                                           |
| *Evenness of sgRNA read distribution per sample*                                      |

| 4. **β Score Histograms**                                                             | 5. **β Score Normalization Comparison**                                                    |
| :-----------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------: |
| <img src="/images/beta_hist.png" width="300">                                         | <img src="/images/beta_norm.png" width="300">                                              |
| *Side-by-side histograms for control and treatment beta scores*                       | *Comparison of raw, cell-cycle, and loess normalized beta scores*                          |

| 6. **β vs β Scatter Plot**                                                            | 7. **β vs β Scatter Plot (Labeled)**                                                       |
| :-----------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------: |
| <img src="/images/beta_scatter.png" width="300">                                      | <img src="/images/beta_scatter_labeled.png" width="300">                                   |
| *Control (x-axis) vs treatment (y-axis) Beta Scores, colored by FDR*                  | *Control vs Treatment Beta Scores with lowest 20 treatment FDR genes labeled*              |

| 8. **Volcano Plot**                                                                   | 9. **Custom Volcano Plot**                                                                 |
| :-----------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------: |
| <img src="/images/volcano_diff_vs_neglog10FDR.png" width="300">                       | <img src="/images/custom_volcano_diff_vs_neglog10FDR.png" width="300">                     |
| *β difference (treat − ctrl) vs −log10(FDR) with 20 top hits labeled*                 | *β difference (treat − ctrl) vs −log10(FDR) with lowest 20 treatment FDR genes labeled*    |

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

+ **Kevin A. Boyd** – Designed and implemented the Snakemake workflow for a Slurm-managed HPC environment, modularized the pipeline structure, implemented all processing steps, designed custom volcano/beta score plots, and created the documentation.  
+ **Jacob Kirkland** – Principal Investigator; provided experimental data and validation of activation screen logic.
+ **Christopher L. Sansam** – Principal Investigator; experimental data and validation of knockout screen logic.
+ **Dean Dawson** - Principal Investigator; provided experimental data used in example plots.  

This workflow was developed as part of a COBRE-funded collaborative effort. While the pipeline was built specifically for use within the Kirkland Lab, it is broadly applicable to CRISPR-Screen data analysis in other research settings.  

---

## 12) License

This project is licensed under the **Apache 2.0**. See the [LICENSE](LICENSE) file for details.  

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
