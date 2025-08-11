# CRISPR-ScreenAnalysis




## 6) Explanation of `samples.csv`

**IMPORTANT**: Always check the config/samples.csv with your sample name, FASTQs locations, include_mle option, and factors defined before running.

This pipeline builds a MAGeCK‐MLE design matrix **directly from `config/samples.csv`** (unless you supply a custom matrix).  
The CSV is simple but flexible, and it produces a design matrix that matches how MAGeCK models log fold changes.

### Required columns in `samples.csv`

| column        | type    | meaning |
|---------------|---------|---------|
| `sample`      | string  | Unique sample label. This must match the sample labels used in `mageck count` (i.e., row names in the count matrix). |
| `fastq1`      | path    | R1 FASTQ path. |
| `fastq2`      | path    | R2 FASTQ path (paired) or a placeholder/empty if single-end (the pipeline still expects the column). |
| `include_mle` | bool    | Whether this row participates in the **MLE** model. Rows with `false` are ignored by the design builder. |
| `factor`      | string  | The treatment/condition label for this sample. Use the literal sentinel **`none`** (or `NA`, empty) to mark **baseline** samples. |

+ **Note**:  
  + **`include_mle`** is **case‐insensitive** — any of the following values will be interpreted as **True**:  
    + `true`, `TRUE`, `True`, `1`, `t`, `T`, `yes`, `Yes`, `y`, `Y`  
    + All other values will be interpreted as **False**.  
  + **`factor`** values are normalized — the following will all be converted to **`none`**:  
    + `none`, `NONE`, `None`, `NA`, `na`, `NaN`, `null`, or an empty cell.  
    + All other non‐`none` values are treated as factor names for one‐hot encoding in the design matrix.  

### MAGeCK MLE Model  

MAGeCK MLE fits gene‐level effects using a **linear design** per sample:
+ effect(sample) = β_baseline
+ β_factor1 * I(sample ∈ factor1)
+ β_factor2 * I(sample ∈ factor2)

+ **Baseline**: all samples contribute a baseline column (all **1**s).
+ **Factors**: each distinct non-`none` value in the `factor` column becomes **one binary column** (one-hot encoding).
+ There are **no interaction terms**; it’s baseline + one-hot factors.

### Why a single `factor` column?  
+ Keeps the CSV compact and readable.  
+ Still supports many designs — just use unique factor labels (e.g., `D14_DMSO`, `D14_PAC`, `Day6_Auxin`).  
+ For more complex designs (interactions, multi-factor models), you can prebuild your own matrix.  

### How the pipeline builds the design matrix

1. **Read & filter**: Keep only rows with `include_mle = true`.
2. **Order rows**: Row order matches the `sample` column (must match counts file).
3. **Normalize factors**: Replace `NA`, empty, `none` → `"none"`.
4. **Collect factors**: Unique factor values (excluding `none`) in first appearance order.
5. **Emit matrix**:  
   + Column 1: `group` (sample name)  
   + Column 2: `baseline` (all `1`)  
   + One column per factor; a row gets `1` if its `factor` matches that column, else `0`.

Matrix is saved as:
+ results/mle/design_matrix.txt

**Unless `custom_design_matrix` is set in `config.yml`.**

### Example1 — Endpoint vs baseline

**CSV**

```csv
sample,fastq1,fastq2,include_mle,factor
N1_D0,...,TRUE,none
N2_D0,...,TRUE,none
N3_D0,...,TRUE,none
N1_D14_DMSO,...,TRUE,D14_DMSO
N2_D14_DMSO,...,TRUE,D14_DMSO
N3_D14_DMSO,...,TRUE,D14_DMSO
N1_D14_PAC,...,TRUE,D14_PAC
N2_D14_PAC,...,TRUE,D14_PAC
N3_D14_PAC,...,TRUE,D14_PAC
```

group         baseline  D14_DMSO  D14_PAC
N1_D0                1         0        0
N2_D0                1         0        0
N3_D0                1         0        0
N1_D14_DMSO          1         1        0
N2_D14_DMSO          1         1        0
N3_D14_DMSO          1         1        0
N1_D14_PAC           1         0        1
N2_D14_PAC           1         0        1
N3_D14_PAC           1         0        1


### Best practices
+ Sample names: must match mageck count output exactly.
+ Baseline: Use none for controls; multiple baselines are fine.
+ Factor names: No spaces; underscores recommended.
+ Balance: Keep replicate counts even across factors.
+ Validation: Always check results/mle/design_matrix.txt before interpreting results.

### MLE Output Interpretation
+ The gene summary includes columns for each factor (plus baseline).
+ β_baseline: shared drift/essentiality effect.
+ β_factorX: treatment effect relative to baseline.

**For endpoint treatment effects, baseline vs. treatment factors is the cleanest approach.**

### Checklist before running
+ Unique sample names
+ include_mle correctly set
+ Baselines labeled none
+ Paths valid
+ Inspect generated matrix in results/mle/design_matrix.txt
 
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
```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

9F. Run on HPC with config.yml options
```
sbatch --wrap="snakemake -j 999 --use-envmodules --rerun-incomplete --latency-wait 300 --cluster-config config/cluster_config.yml --cluster 'sbatch -A {cluster.account} -p {cluster.partition} --cpus-per-task {cluster.cpus-per-task}  -t {cluster.time} --mem {cluster.mem} --output {cluster.output} --job-name {cluster.name}'"
```
