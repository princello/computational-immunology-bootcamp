# Final Practice: Capstone Project

**Estimated time: 2-3 weeks**

Congratulations. You have made it to the final stage of the Computational Immunology
Bootcamp. Over the past four weeks, you installed a bioinformatics environment from
scratch, learned the fundamentals of immunology, ran a full single-cell RNA-seq pipeline,
analyzed TCR and BCR repertoires, and integrated multi-omics data. Those are real skills,
and you should be proud of the work it took to get here.

This capstone is where you bring everything together into a project that looks like real
research. The output will be three things: clean, well-commented Jupyter notebooks;
publication-quality figures; and a short written report with the same structure scientists
use in their papers. When you are done, you will have a portfolio piece that demonstrates
genuine computational biology ability.

---

## Table of Contents

1. [Overview and Track Selection](#1-overview-and-track-selection)
2. [Project Setup](#2-project-setup)
3. [Track 1: RNA + ATAC Integration](#3-track-1-rna--atac-integration)
   - [Step 1: Choose and Download Your Data](#step-1-choose-and-download-your-data)
   - [Step 2: RNA Quality Control and Processing](#step-2-rna-quality-control-and-processing)
   - [Step 3: ATAC Quality Control and Processing](#step-3-atac-quality-control-and-processing)
   - [Step 4: Integration and Co-embedding](#step-4-integration-and-co-embedding)
   - [Step 5: Cell Type Annotation](#step-5-cell-type-annotation)
   - [Step 6: Regulatory Inference (Challenge)](#step-6-regulatory-inference-challenge)
   - [Step 7: Biological Summary](#step-7-biological-summary)
4. [Track 2: TCR Specificity Prediction](#4-track-2-tcr-specificity-prediction)
   - [Step 1: Understand the Problem](#step-1-understand-the-problem)
   - [Step 2: Get the Data](#step-2-get-the-data)
   - [Step 3: Feature Engineering](#step-3-feature-engineering)
   - [Step 4: Build a Baseline Model](#step-4-build-a-baseline-model)
   - [Step 5: Improve with Better Models](#step-5-improve-with-better-models)
   - [Step 6: Evaluate and Visualize](#step-6-evaluate-and-visualize)
   - [Step 7: Proper Data Splitting (Critical)](#step-7-proper-data-splitting-critical)
   - [Step 8: Deep Learning Approach (Challenge)](#step-8-deep-learning-approach-challenge)
5. [Report Template](#5-report-template)
6. [Minimum Figure Set](#6-minimum-figure-set)
7. [Self-Check and Completion Criteria](#7-self-check-and-completion-criteria)

---

## 1. Overview and Track Selection

You will choose **one** of two tracks for your capstone project. Each track draws on
different parts of what you have learned, and each emphasizes different skills.

| | Track 1: RNA + ATAC Integration | Track 2: TCR Specificity Prediction |
|---|---|---|
| **What you build** | Joint analysis of gene expression and chromatin accessibility | Machine learning model that predicts TCR-epitope binding |
| **Skills emphasized** | Biology, data integration, visualization | Coding, machine learning, model evaluation |
| **Difficulty** | More biology-heavy; moderate coding | More coding/ML-heavy; moderate biology |
| **Draws on** | Week 2 (scRNA-seq pipeline) and Week 4 (multi-omics) | Week 3 (TCR/BCR analysis) and general Python |
| **Output** | Integrated UMAP, cell type annotations, regulatory insights | Trained models, ROC/PR curves, fair evaluation |

**How to choose:** If you found Week 4 (multi-omics integration) the most interesting part
of the bootcamp and enjoy interpreting biological results, choose Track 1. If you found
Week 3 (repertoire analysis) fascinating and want to learn machine learning, choose
Track 2. If you are ambitious and have the time, you may complete both tracks for extra
credit.

Read through both track descriptions before committing. You do not need to decide
immediately.

---

## 2. Project Setup

**[REQUIRED]**

Before starting either track, set up a clean project directory. Good project organization
is a hallmark of professional research. Every file should have a logical home.

### Create the Project Structure

Open your terminal and run:

```bash
mkdir -p ~/bootcamp/final_project/{data/{raw,processed},notebooks,src,figs,reports}
cd ~/bootcamp/final_project
```

This creates the following structure:

```
~/bootcamp/final_project/
  data/
    raw/              # Original downloaded data (never modify these files)
    processed/        # Cleaned, filtered, transformed data
  notebooks/
    track1_rna_atac_integration.ipynb    (if doing Track 1)
    track2_tcr_specificity.ipynb         (if doing Track 2)
  src/
    utils.py          # Helper functions you write and reuse
  figs/               # All saved figures go here
  reports/
    final_report.md   # Your written report
  environment.yml     # Conda environment specification (for reproducibility)
  README.md           # Brief description of your project
```

### Create the Conda Environment File

Create a file called `environment.yml` in your project root. This records exactly which
software versions you used, so anyone (including your future self) can reproduce your work.

For **Track 1**, your `environment.yml` should contain:

```yaml
name: capstone
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.10
  - scanpy>=1.9
  - muon>=0.1.5
  - scvi-tools>=0.20
  - episcanpy>=0.4
  - matplotlib>=3.7
  - seaborn>=0.12
  - pandas>=1.5
  - numpy>=1.23
  - jupyter
  - pip
  - pip:
    - celltypist
```

For **Track 2**, your `environment.yml` should contain:

```yaml
name: capstone
channels:
  - conda-forge
dependencies:
  - python=3.10
  - scanpy>=1.9
  - scikit-learn>=1.2
  - matplotlib>=3.7
  - seaborn>=0.12
  - pandas>=1.5
  - numpy>=1.23
  - jupyter
  - pip
  - pip:
    - xgboost
```

Create the environment:

```bash
conda env create -f environment.yml
conda activate capstone
```

**If you get errors:** Try creating a simpler environment first and adding packages one at
a time:

```bash
conda create -n capstone python=3.10 scanpy matplotlib seaborn pandas numpy jupyter
conda activate capstone
pip install scvi-tools muon episcanpy celltypist   # Track 1
# or
pip install scikit-learn xgboost                   # Track 2
```

### Create Your Notebook

Launch Jupyter and create a new notebook in the `notebooks/` folder:

```bash
cd ~/bootcamp/final_project
jupyter notebook
```

Navigate to the `notebooks/` folder and create a new Python 3 notebook. Name it
`track1_rna_atac_integration.ipynb` or `track2_tcr_specificity.ipynb` depending on your
chosen track.

In the first cell of your notebook, add a title and description:

```python
# Final Capstone Project
# Track 1: RNA + ATAC Integration  (or Track 2: TCR Specificity Prediction)
# Author: [Your Name]
# Date: [Today's Date]
#
# This notebook contains my capstone analysis for the
# Computational Immunology Bootcamp.
```

`[CHECKPOINT 0]` Your project directory exists, your conda environment is active, and you
have an empty notebook ready to fill. Verify by running `ls ~/bootcamp/final_project/` in
the terminal -- you should see `data`, `notebooks`, `src`, `figs`, `reports`.

---

## 3. Track 1: RNA + ATAC Integration

**Goal:** Reproduce the core idea of joint scRNA-seq and scATAC-seq analysis using a
public dataset. You will process both RNA (gene expression) and ATAC (chromatin
accessibility) data from the same cells, integrate them into a single representation,
annotate cell types using both modalities, and explore regulatory relationships.

This track builds directly on the multi-omics concepts from Week 4 and the scRNA-seq
pipeline from Week 2. The key new skill you will learn is handling ATAC-seq data, which
measures which parts of the genome are physically "open" and accessible in each cell. Open
chromatin near a gene suggests that gene might be actively regulated.

---

### Step 1: Choose and Download Your Data

**[REQUIRED]**

You need a dataset where RNA and ATAC were measured from the **same cells**. This is called
a "multiome" dataset. The best publicly available option is from 10x Genomics.

#### Option A: 10x Genomics Multiome PBMC 10k (Recommended)

This dataset contains approximately 10,000 peripheral blood mononuclear cells (PBMCs) with
both RNA and ATAC data from each cell. You already know PBMCs from Week 1 -- they include
T cells, B cells, NK cells, and monocytes.

**How to download:**

1. Go to: https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0
2. You may need to create a free 10x Genomics account
3. Download the following files:
   - **Filtered feature barcode matrix (HDF5)** -- this is the main data file containing
     both RNA and ATAC counts
   - **ATAC Per fragment information file (TSV.GZ)** -- this contains the raw ATAC
     fragment locations
   - **ATAC Per fragment information index (TSV.GZ index)** -- the index file for the
     fragments
4. Save all downloaded files into `~/bootcamp/final_project/data/raw/`

**What are these files?**

- The **.h5 file** (HDF5 format) is a compressed binary file containing the cells-by-features
  matrix. "Features" includes both gene expression counts (RNA) and peak counts (ATAC). This
  is the same type of file you used in Week 2, but now with two modalities inside.
- The **fragments.tsv.gz file** contains the raw location of every ATAC fragment in the
  genome. Each line tells you: this barcode had a fragment of open chromatin between
  position X and position Y on chromosome Z. You need this for some ATAC quality control
  steps.
- The **.tbi file** is an index that lets software quickly look up fragments for a
  specific genomic region without reading the entire file.

```python
# After downloading, verify the files exist
import os

data_dir = os.path.expanduser("~/bootcamp/final_project/data/raw")
expected_files = [
    "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
    "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
    "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi"
]

for f in expected_files:
    path = os.path.join(data_dir, f)
    if os.path.exists(path):
        size_mb = os.path.getsize(path) / 1e6
        print(f"  Found: {f} ({size_mb:.1f} MB)")
    else:
        print(f"  MISSING: {f}")
```

**Expected output:** Three files found. The HDF5 file will be around 60-80 MB. The
fragments file will be much larger (1-3 GB).

#### Option B: Use muon to Load from AnnotatedData

If the 10x download is not working or the files are too large, you can try loading a
smaller preprocessed version using the `muon` package:

```python
import muon as mu

# muon provides some example datasets
# Check the muon documentation for available datasets:
# https://muon.readthedocs.io/

# If a multiome PBMC dataset is available:
mdata = mu.read("path/to/dataset")
```

#### Option C: Fallback -- Use Separate RNA and ATAC Datasets

If you cannot get a multiome dataset (same cells for both modalities), you can still do
this project using separate RNA and ATAC PBMC datasets. The analysis changes slightly
because you will need to "align" cells across datasets instead of using paired
measurements.

For this fallback:
- Use the pbmc3k dataset from scanpy for RNA: `adata_rna = sc.datasets.pbmc3k()`
- Download a PBMC ATAC dataset from 10x Genomics (search for "PBMC ATAC" on their
  datasets page)

The integration steps in Step 4 will use label transfer instead of direct pairing.

`[CHECKPOINT 1]` Your data files are downloaded and saved in `data/raw/`. You know what
each file contains. If using the multiome dataset, you should have at least the HDF5 file.

---

### Step 2: RNA Quality Control and Processing

**[REQUIRED]**

This step should feel familiar from Week 2. You will load the RNA portion of the data,
compute quality metrics, visualize them, and filter out low-quality cells.

```python
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor='white')
```

#### Load the RNA Data

```python
# Load the filtered feature barcode matrix
# This file contains BOTH RNA and ATAC features
# We will separate them in this step

adata = sc.read_10x_h5(
    "data/raw/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
    gex_only=False  # Load all features, not just gene expression
)

print(f"Full dataset: {adata.n_obs} cells x {adata.n_vars} features")
print(f"Feature types: {adata.var['feature_types'].value_counts().to_dict()}")
```

**Expected output:** You should see something like:
```
Full dataset: ~10000 cells x ~150000 features
Feature types: {'Gene Expression': ~36000, 'Peaks': ~110000}
```

The dataset contains two types of features: Gene Expression (RNA) and Peaks (ATAC). We
need to separate them.

```python
# Separate RNA and ATAC features
# RNA features have feature_types == 'Gene Expression'
# ATAC features have feature_types == 'Peaks'

adata_rna = adata[:, adata.var['feature_types'] == 'Gene Expression'].copy()
adata_atac = adata[:, adata.var['feature_types'] == 'Peaks'].copy()

print(f"RNA: {adata_rna.n_obs} cells x {adata_rna.n_vars} genes")
print(f"ATAC: {adata_atac.n_obs} cells x {adata_atac.n_vars} peaks")

# Make gene names unique (some gene names appear more than once)
adata_rna.var_names_make_unique()
```

#### RNA Quality Control

```python
# Calculate QC metrics -- same approach as Week 2
# Mitochondrial genes start with 'MT-' in human data
adata_rna.var['mt'] = adata_rna.var_names.str.startswith('MT-')

# Ribosomal genes start with 'RPS' or 'RPL'
adata_rna.var['ribo'] = adata_rna.var_names.str.startswith(('RPS', 'RPL'))

sc.pp.calculate_qc_metrics(
    adata_rna,
    qc_vars=['mt', 'ribo'],
    inplace=True
)

# Print summary statistics
print("RNA QC Summary:")
print(f"  Genes per cell: {adata_rna.obs['n_genes_by_counts'].median():.0f} median "
      f"(range: {adata_rna.obs['n_genes_by_counts'].min()}-"
      f"{adata_rna.obs['n_genes_by_counts'].max()})")
print(f"  Counts per cell: {adata_rna.obs['total_counts'].median():.0f} median")
print(f"  Mitochondrial %: {adata_rna.obs['pct_counts_mt'].median():.1f}% median")
```

#### Visualize RNA QC Metrics

```python
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Panel 1: Genes per cell
axes[0].hist(adata_rna.obs['n_genes_by_counts'], bins=50, color='steelblue',
             edgecolor='white')
axes[0].set_title('Genes per Cell')
axes[0].set_xlabel('Number of genes detected')
axes[0].set_ylabel('Number of cells')
# Draw threshold lines in red
axes[0].axvline(x=500, color='red', linestyle='--', linewidth=1, label='min=500')
axes[0].axvline(x=5000, color='red', linestyle='--', linewidth=1, label='max=5000')
axes[0].legend(fontsize=8)

# Panel 2: Total counts per cell
axes[1].hist(adata_rna.obs['total_counts'], bins=50, color='coral',
             edgecolor='white')
axes[1].set_title('Total Counts per Cell')
axes[1].set_xlabel('Total UMI counts')
axes[1].set_ylabel('Number of cells')

# Panel 3: Mitochondrial percentage
axes[2].hist(adata_rna.obs['pct_counts_mt'], bins=50, color='mediumseagreen',
             edgecolor='white')
axes[2].set_title('% Mitochondrial Reads')
axes[2].set_xlabel('Mitochondrial %')
axes[2].set_ylabel('Number of cells')
axes[2].axvline(x=20, color='red', linestyle='--', linewidth=1, label='max=20%')
axes[2].legend(fontsize=8)

plt.tight_layout()
plt.savefig('figs/rna_qc.png', dpi=150, bbox_inches='tight')
plt.show()

print("Figure saved to figs/rna_qc.png")
```

**Expected output:** Three histograms showing the distribution of QC metrics. The genes
per cell histogram should have a main peak and possibly a small tail of low-quality cells
on the left. The mitochondrial percentage should be mostly low (under 10%) with a tail
extending higher.

#### Apply RNA Filters

```python
# Filter cells based on QC thresholds
# These thresholds are guidelines -- adjust if your distributions suggest
# different cutoffs
n_before = adata_rna.n_obs

adata_rna = adata_rna[
    (adata_rna.obs['n_genes_by_counts'] > 500) &
    (adata_rna.obs['n_genes_by_counts'] < 5000) &
    (adata_rna.obs['pct_counts_mt'] < 20)
].copy()

n_after = adata_rna.n_obs
print(f"Cells before RNA QC: {n_before}")
print(f"Cells after RNA QC:  {n_after}")
print(f"Removed: {n_before - n_after} cells ({100*(n_before-n_after)/n_before:.1f}%)")
```

**Expected output:** You should retain the majority of cells (typically 80-95%). If you
lose more than 30%, revisit your thresholds -- they may be too strict for this dataset.

#### Standard RNA Processing

```python
# Normalize, log-transform, find variable genes, scale, PCA
# This is the standard pipeline from Week 2

# Filter genes expressed in very few cells
sc.pp.filter_genes(adata_rna, min_cells=10)

# Store raw counts for later use (some tools need them)
adata_rna.layers['counts'] = adata_rna.X.copy()

# Normalize to 10,000 counts per cell
sc.pp.normalize_total(adata_rna, target_sum=1e4)

# Log transform
sc.pp.log1p(adata_rna)

# Find highly variable genes
sc.pp.highly_variable_genes(adata_rna, n_top_genes=3000)
print(f"Highly variable genes: {adata_rna.var['highly_variable'].sum()}")

# Scale the data (zero mean, unit variance)
sc.pp.scale(adata_rna, max_value=10)

# PCA
sc.tl.pca(adata_rna, n_comps=50, use_highly_variable=True)

# Quick check: how much variance do the first PCs explain?
print(f"Variance explained by PC1: {adata_rna.uns['pca']['variance_ratio'][0]*100:.1f}%")
print(f"Variance explained by PC1-10: {sum(adata_rna.uns['pca']['variance_ratio'][:10])*100:.1f}%")
```

#### RNA UMAP (Preliminary)

```python
# Compute neighbors and UMAP for the RNA data alone
sc.pp.neighbors(adata_rna, n_pcs=30)
sc.tl.umap(adata_rna)

# Cluster
sc.tl.leiden(adata_rna, resolution=0.8)

# Plot
sc.pl.umap(adata_rna, color='leiden', title='RNA-only UMAP (Leiden clusters)')
plt.savefig('figs/rna_umap_preliminary.png', dpi=150, bbox_inches='tight')
plt.show()
```

**Expected output:** A UMAP plot with colored clusters. You should see several distinct
groups corresponding to different immune cell types. This is your RNA-only baseline -- we
will compare it to the integrated result later.

`[CHECKPOINT 2]` RNA data is loaded, quality-filtered, normalized, and reduced to a UMAP.
You have a preliminary clustering. Save your notebook.

---

### Step 3: ATAC Quality Control and Processing

**[REQUIRED]**

ATAC-seq data requires different quality control metrics and a different processing
pipeline compared to RNA. This is because ATAC measures a fundamentally different thing:
instead of counting mRNA molecules, it counts where the genome is physically accessible
(open chromatin). The data is binary-like -- a peak is either accessible or not in a given
cell.

#### ATAC-Specific QC Concepts

Before writing code, understand these three key ATAC QC metrics:

**TSS Enrichment (Transcription Start Site enrichment):** In a good ATAC experiment, you
expect to see more open chromatin near the places where genes begin (transcription start
sites). This is because active genes need to have their promoter regions open. TSS
enrichment compares the amount of signal near TSSs versus background. A value of 4 or
higher is typical for good-quality cells. If this number is low (below 2), the cell's data
may not be reliable.

**Nucleosome Signal:** DNA in the nucleus wraps around proteins called histones in units
called nucleosomes. ATAC fragments come in characteristic sizes: short fragments (~100-150
bp) come from nucleosome-free regions (NFR, the truly open parts), and longer fragments
(~200 bp) come from regions where the enzyme cut around a single nucleosome. The
nucleosome signal is the ratio of mononucleosomal to nucleosome-free fragments. Lower is
better -- a high nucleosome signal means less of the data comes from truly open chromatin.

**FRiP (Fraction of Reads in Peaks):** What fraction of all ATAC reads fall within called
peaks (regions of significant enrichment)? Higher is better, because it means the
experiment captured mostly signal rather than noise. A FRiP above 0.2 (20%) is acceptable;
above 0.4 is good.

#### ATAC Processing

```python
# ATAC processing is different from RNA processing
# The key steps are:
# 1. Filter peaks present in very few cells
# 2. Apply TF-IDF normalization
# 3. Perform LSI dimensionality reduction

# First, let's filter the ATAC data to only keep cells that also passed RNA QC
# Since this is a multiome dataset, the cell barcodes match between RNA and ATAC
shared_cells = adata_rna.obs_names.intersection(adata_atac.obs_names)
print(f"Cells in RNA after QC: {adata_rna.n_obs}")
print(f"Cells in ATAC (before QC): {adata_atac.n_obs}")
print(f"Shared cells: {len(shared_cells)}")

adata_atac = adata_atac[shared_cells].copy()
print(f"ATAC cells after matching to RNA: {adata_atac.n_obs}")
```

#### Filter ATAC Peaks

```python
# Remove peaks that are accessible in very few cells
# These are likely noise
min_cells_for_peak = 20
peak_counts = np.array((adata_atac.X > 0).sum(axis=0)).flatten()

print(f"Total peaks before filtering: {adata_atac.n_vars}")
adata_atac = adata_atac[:, peak_counts >= min_cells_for_peak].copy()
print(f"Total peaks after filtering (>={min_cells_for_peak} cells): {adata_atac.n_vars}")
```

**Expected output:** You should retain most peaks, removing perhaps 10-20% of very rare
peaks.

#### TF-IDF Normalization

```python
# TF-IDF stands for Term Frequency - Inverse Document Frequency
# It comes from text mining, where:
#   - "Term Frequency" = how often a word appears in a document
#   - "Inverse Document Frequency" = words that appear in many documents are less informative
#
# For ATAC-seq:
#   - "Term" = a peak (genomic region)
#   - "Document" = a cell
#   - A peak that is open in many cells is less informative about cell identity
#   - A peak that is open in only a few cells is highly informative
#
# TF-IDF normalization accounts for both the depth of each cell
# and the informativeness of each peak

from sklearn.preprocessing import normalize
from scipy.sparse import diags, issparse
import scipy.sparse as sp

def tfidf_normalize(adata_atac):
    """
    Apply TF-IDF normalization to ATAC peak count matrix.
    This is the standard normalization for binary/count ATAC data.
    """
    # Binarize the matrix (any count > 0 becomes 1)
    X = adata_atac.X.copy()
    if issparse(X):
        X.data[:] = 1  # Set all nonzero values to 1
    else:
        X = (X > 0).astype(float)

    # Term Frequency: normalize each cell to sum to 1
    tf = normalize(X, norm='l1', axis=1)

    # Inverse Document Frequency: log(N / n_j) where N = total cells,
    # n_j = number of cells with peak j open
    n_cells = X.shape[0]
    peak_frequencies = np.array(X.sum(axis=0)).flatten()
    # Add 1 to avoid division by zero
    idf = np.log1p(n_cells / (peak_frequencies + 1))

    # TF-IDF = TF * IDF
    if issparse(tf):
        tfidf = tf.multiply(idf)
    else:
        tfidf = tf * idf

    return tfidf

adata_atac.X = tfidf_normalize(adata_atac)
print("TF-IDF normalization complete.")
```

#### LSI Dimensionality Reduction

```python
# LSI stands for Latent Semantic Indexing
# It is the ATAC equivalent of PCA for RNA
# Technically, LSI is just SVD (Singular Value Decomposition) applied to
# the TF-IDF-normalized matrix
#
# We skip the first component because it usually captures sequencing depth
# (technical variation) rather than biological variation

from sklearn.decomposition import TruncatedSVD

# Compute LSI with 50 components
n_components = 50
svd = TruncatedSVD(n_components=n_components, random_state=42)
lsi = svd.fit_transform(adata_atac.X)

# The first component often correlates with sequencing depth -- check this
depth = np.array(adata_atac.X.sum(axis=1)).flatten()
corr_with_depth = np.corrcoef(lsi[:, 0], depth)[0, 1]
print(f"Correlation of LSI component 1 with sequencing depth: {corr_with_depth:.3f}")

if abs(corr_with_depth) > 0.5:
    print("  -> High correlation. Dropping component 1 (standard practice).")
    lsi_filtered = lsi[:, 1:]
else:
    print("  -> Low correlation. Keeping all components.")
    lsi_filtered = lsi

# Store in the AnnData object
adata_atac.obsm['X_lsi'] = lsi_filtered

# Check variance explained
print(f"\nVariance explained by first 10 LSI components: "
      f"{svd.explained_variance_ratio_[1:11].sum()*100:.1f}%")
```

#### ATAC UMAP (Preliminary)

```python
# Compute UMAP for the ATAC data alone
sc.pp.neighbors(adata_atac, use_rep='X_lsi', n_neighbors=30)
sc.tl.umap(adata_atac)

# We don't have labels yet, so just show the structure
sc.tl.leiden(adata_atac, resolution=0.5)
sc.pl.umap(adata_atac, color='leiden', title='ATAC-only UMAP (Leiden clusters)')
plt.savefig('figs/atac_umap_preliminary.png', dpi=150, bbox_inches='tight')
plt.show()
```

**Expected output:** A UMAP plot from ATAC data alone. The structure may look different
from the RNA UMAP -- that is expected and interesting. ATAC captures regulatory state,
which can differ from transcriptional state.

`[CHECKPOINT 3]` ATAC data is filtered, TF-IDF normalized, and reduced via LSI. You have
UMAP plots for both RNA and ATAC modalities separately. Save your notebook.

---

### Step 4: Integration and Co-embedding

**[REQUIRED]**

Now comes the core of this capstone: combining RNA and ATAC into a single integrated
representation. Because this is a multiome dataset (RNA and ATAC from the same cells), you
can directly combine the two latent spaces. You do not need to infer which cells correspond
across modalities -- they are already paired by barcode.

There are two main approaches. We will implement the simpler one (concatenation + WNN-like
weighting) and describe the more advanced one (MultiVI) as an option.

#### Approach A: Concatenated Latent Spaces

The simplest approach: concatenate the RNA PCA embeddings and ATAC LSI embeddings for each
cell, then compute a shared UMAP.

```python
# Ensure both objects have the same cells in the same order
shared_cells = adata_rna.obs_names.intersection(adata_atac.obs_names)
adata_rna_shared = adata_rna[shared_cells].copy()
adata_atac_shared = adata_atac[shared_cells].copy()

print(f"Shared cells for integration: {len(shared_cells)}")

# Get the RNA PCA embeddings (first 30 components)
rna_embedding = adata_rna_shared.obsm['X_pca'][:, :30]

# Get the ATAC LSI embeddings (first 30 components, or all if fewer)
n_lsi = min(30, adata_atac_shared.obsm['X_lsi'].shape[1])
atac_embedding = adata_atac_shared.obsm['X_lsi'][:, :n_lsi]

print(f"RNA embedding shape:  {rna_embedding.shape}")
print(f"ATAC embedding shape: {atac_embedding.shape}")

# Concatenate the two embeddings
# Normalize each embedding to have similar scale
from sklearn.preprocessing import StandardScaler

scaler_rna = StandardScaler()
scaler_atac = StandardScaler()

rna_scaled = scaler_rna.fit_transform(rna_embedding)
atac_scaled = scaler_atac.fit_transform(atac_embedding)

combined_embedding = np.concatenate([rna_scaled, atac_scaled], axis=1)
print(f"Combined embedding shape: {combined_embedding.shape}")

# Create a new AnnData object for the integrated analysis
adata_integrated = adata_rna_shared.copy()
adata_integrated.obsm['X_integrated'] = combined_embedding
```

```python
# Compute neighbors and UMAP on the integrated embedding
sc.pp.neighbors(adata_integrated, use_rep='X_integrated', n_neighbors=30)
sc.tl.umap(adata_integrated)
sc.tl.leiden(adata_integrated, resolution=0.8)

# Plot the integrated UMAP
sc.pl.umap(adata_integrated, color='leiden',
           title='Integrated RNA+ATAC UMAP (Leiden clusters)')
plt.savefig('figs/integrated_umap.png', dpi=150, bbox_inches='tight')
plt.show()
```

**Expected output:** A UMAP that captures information from both modalities. The clusters
should be well-separated, and you may see finer structure than in either single-modality
UMAP.

#### Approach B: MultiVI (Advanced, Optional)

If you want a more principled integration, scvi-tools provides MultiVI, a deep generative
model designed specifically for multiome data. It learns a shared latent space that
accounts for the different distributions of RNA and ATAC data.

```python
# --- OPTIONAL: Advanced integration with MultiVI ---
# Only attempt this if Approach A is working and you want to explore further

import scvi

# MultiVI requires raw counts for RNA and binary peaks for ATAC
# Prepare the data according to MultiVI documentation

# Create a combined AnnData with both RNA and ATAC features
# See: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/multimodal/MultiVI_tutorial.html

# Set up the model
# scvi.model.MULTIVI.setup_anndata(
#     adata_combined,
#     batch_key=None,
# )

# Train
# model = scvi.model.MULTIVI(adata_combined)
# model.train(max_epochs=100)

# Get the latent representation
# adata_combined.obsm['X_multivi'] = model.get_latent_representation()

# Then use this representation for UMAP and clustering
# sc.pp.neighbors(adata_combined, use_rep='X_multivi')
# sc.tl.umap(adata_combined)
```

#### Compare All Three UMAPs

```python
# Create a side-by-side comparison of RNA-only, ATAC-only, and integrated UMAPs

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel 1: RNA-only
sc.pl.umap(adata_rna_shared, color='leiden', ax=axes[0], show=False,
           title='RNA-only')

# Panel 2: ATAC-only
# Transfer RNA leiden labels to ATAC for comparison
adata_atac_shared.obs['rna_leiden'] = adata_rna_shared.obs['leiden'].values
sc.pl.umap(adata_atac_shared, color='rna_leiden', ax=axes[1], show=False,
           title='ATAC-only (colored by RNA clusters)')

# Panel 3: Integrated
sc.pl.umap(adata_integrated, color='leiden', ax=axes[2], show=False,
           title='Integrated RNA+ATAC')

plt.tight_layout()
plt.savefig('figs/umap_comparison.png', dpi=150, bbox_inches='tight')
plt.show()
```

**Expected output:** Three UMAPs side by side. Look for how the structure changes when
you add ATAC information. Do clusters become more distinct? Do any clusters split or merge?

`[CHECKPOINT 4]` Integration is complete. You have a combined RNA+ATAC embedding and a
comparison of single-modality versus integrated UMAPs. Save your notebook.

---

### Step 5: Cell Type Annotation

**[REQUIRED]**

Now assign biological cell type labels to your clusters. You will use a combination of
known marker genes and automated annotation.

#### Marker-Based Annotation

```python
# Define canonical immune cell markers
# You learned many of these in Weeks 1-2
markers = {
    'CD4 T cells':  ['CD3D', 'CD3E', 'CD4', 'IL7R'],
    'CD8 T cells':  ['CD3D', 'CD3E', 'CD8A', 'CD8B'],
    'NK cells':     ['NKG7', 'GNLY', 'KLRD1', 'NCAM1'],
    'B cells':      ['CD79A', 'CD79B', 'MS4A1', 'CD19'],
    'Monocytes':    ['CD14', 'LYZ', 'S100A8', 'S100A9'],
    'DCs':          ['FCER1A', 'CST3', 'CLEC10A'],
    'Platelets':    ['PPBP', 'PF4'],
}

# Dot plot of marker genes across clusters
sc.pl.dotplot(
    adata_integrated,
    var_names=markers,
    groupby='leiden',
    standard_scale='var',
    title='Marker Gene Expression by Cluster'
)
plt.savefig('figs/marker_dotplot.png', dpi=150, bbox_inches='tight')
plt.show()
```

**Expected output:** A dot plot where each row is a cluster and each column is a marker
gene. The size of the dot indicates the fraction of cells expressing the gene, and the
color indicates the expression level. You should see clear patterns -- for example, CD3D
should be high in T cell clusters and low everywhere else.

```python
# Based on the dot plot, manually assign cell type labels
# YOU MUST ADJUST THIS based on YOUR data -- look at the dot plot!
# The cluster numbers may be different in your analysis

cluster_to_celltype = {
    '0': 'CD4 T cells',
    '1': 'CD8 T cells',
    '2': 'Monocytes',
    '3': 'B cells',
    '4': 'NK cells',
    '5': 'CD4 T cells',
    '6': 'Monocytes',
    '7': 'DCs',
    '8': 'Platelets',
    # Add more clusters as needed based on YOUR dot plot
}

# Apply the labels
adata_integrated.obs['cell_type'] = (
    adata_integrated.obs['leiden']
    .map(cluster_to_celltype)
    .fillna('Unknown')
)

# Print cell type counts
print("Cell type counts:")
print(adata_integrated.obs['cell_type'].value_counts())
```

#### Automated Annotation with CellTypist (Optional)

```python
# CellTypist uses a pre-trained model to automatically label cells
# This is a good sanity check for your manual annotations

import celltypist
from celltypist import models

# Download a pre-trained immune cell model
models.download_models(model='Immune_All_Low.pkl')

# Run prediction
# CellTypist needs log-normalized data
predictions = celltypist.annotate(
    adata_integrated,
    model='Immune_All_Low.pkl',
    majority_voting=True
)

# Add predictions to the AnnData
adata_integrated.obs['celltypist_label'] = predictions.predicted_labels['majority_voting']

# Compare manual and automated labels
print("\nCellTypist predictions:")
print(adata_integrated.obs['celltypist_label'].value_counts().head(10))
```

#### Annotated UMAP

```python
# Create the final annotated UMAP
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

sc.pl.umap(adata_integrated, color='cell_type', ax=axes[0], show=False,
           title='Manual Annotation', frameon=True)

sc.pl.umap(adata_integrated, color='celltypist_label', ax=axes[1], show=False,
           title='CellTypist Annotation', frameon=True)

plt.tight_layout()
plt.savefig('figs/cell_type_umap.png', dpi=150, bbox_inches='tight')
plt.show()
```

**Expected output:** Two UMAPs showing your manual annotations and the automated
annotations side by side. They should broadly agree. If there are major disagreements,
investigate -- look at marker genes in the disagreeing clusters.

`[CHECKPOINT 5]` Cells are annotated with cell type labels. You have verified your
annotations using marker genes and (optionally) an automated tool. Save your notebook.

---

### Step 6: Regulatory Inference (Challenge)

**[CHALLENGE]**

This section goes beyond standard analysis into the realm of regulatory genomics. The idea
is to use the ATAC data to infer which regulatory elements (enhancers, promoters) control
which genes. This is where the real power of joint RNA+ATAC analysis shines.

You should attempt this section only after completing Steps 1-5 and writing your report for
those steps.

#### Gene Activity Scores

Gene activity is a way to estimate gene expression from ATAC data alone. The logic: if the
chromatin near a gene's promoter is open, the gene is more likely to be transcribed.

```python
# Gene activity: sum ATAC signal in peaks near each gene's promoter and body
#
# For each gene:
#   1. Find its genomic location (chromosome, start, end)
#   2. Find all ATAC peaks within a window around the gene (e.g., 2 kb upstream)
#   3. Sum the accessibility of those peaks for each cell
#
# This gives a "gene activity" score per cell per gene

# In practice, tools like episcanpy or ArchR compute gene activity scores
# Here is a conceptual implementation:

# Parse peak names to get genomic coordinates
# Peak names look like: chr1:9777-10668

def parse_peak_name(peak_name):
    """Parse a peak name like 'chr1:9777-10668' into (chrom, start, end)."""
    try:
        chrom, coords = peak_name.split(':')
        start, end = coords.split('-')
        return chrom, int(start), int(end)
    except (ValueError, AttributeError):
        return None, None, None

# Example: Check what peak names look like in your data
print("Sample peak names:")
for p in adata_atac_shared.var_names[:5]:
    print(f"  {p}")
    chrom, start, end = parse_peak_name(p)
    if chrom:
        print(f"    -> {chrom}: {start:,} - {end:,}")
```

```python
# A simplified gene activity calculation for a few key genes
# In practice, you would use a gene annotation file (GTF) to get gene coordinates
# and compute this genome-wide

# For now, let's visualize which peaks are accessible in each cell type
# and how that relates to gene expression

# Plot a heatmap of key marker gene expression across cell types
marker_genes = ['CD3D', 'CD8A', 'CD79A', 'CD14', 'NKG7', 'FCER1A']

sc.pl.heatmap(
    adata_integrated,
    var_names=marker_genes,
    groupby='cell_type',
    standard_scale='var',
    cmap='viridis',
    show_gene_labels=True
)
plt.savefig('figs/gene_activity_heatmap.png', dpi=150, bbox_inches='tight')
plt.show()
```

#### TF Motif Enrichment (Conceptual)

```python
# Transcription Factor (TF) motif analysis asks:
# "Which transcription factors have their binding motifs enriched
#  in the open chromatin regions of each cell type?"
#
# The idea:
# 1. For each cell type, identify the peaks that are uniquely accessible
#    (open in that cell type but not others)
# 2. Search those peak sequences for known TF binding motifs
#    (from databases like JASPAR)
# 3. TFs whose motifs are enriched are likely active regulators
#    in that cell type
#
# Tools for this:
# - chromVAR (R package): computes per-cell TF deviation scores
# - pycistopic (Python): topic modeling for ATAC, includes motif analysis
# - SCENIC+ (Python): comprehensive regulatory network inference
#
# This is computationally intensive and requires additional reference files
# (genome sequence, motif databases). If you want to attempt it:
#
# 1. Install pycistopic or use episcanpy's motif functions
# 2. Download the JASPAR motif database
# 3. Run motif enrichment on cell-type-specific peaks
#
# Expected result: T cell peaks enriched for TCF/LEF motifs,
# B cell peaks enriched for PAX5/EBF1 motifs, monocyte peaks
# enriched for SPI1 (PU.1) and CEBPA motifs.

print("TF motif analysis is an advanced topic.")
print("If you complete it, include the results in your report as a bonus figure.")
```

#### Peak-to-Gene Links

```python
# Peak-to-gene links ask:
# "For each ATAC peak, is there a nearby gene whose expression
#  correlates with the peak's accessibility?"
#
# The approach:
# 1. For each peak, find all genes within a window (e.g., 250 kb)
# 2. Across cells, correlate peak accessibility with gene expression
# 3. Statistically significant positive correlations suggest that
#    the peak is a regulatory element (enhancer/promoter) for that gene
#
# This requires paired RNA+ATAC data from the same cells -- exactly
# what we have in the multiome dataset!
#
# Tools: ArchR (R) or custom correlation analysis
# The compute is straightforward but involves many peak-gene pairs

# Conceptual code:
# for each peak:
#     for each gene within 250kb:
#         correlation = pearsonr(peak_accessibility, gene_expression)
#         if correlation.pvalue < 0.01 and correlation.statistic > 0:
#             store as significant link

print("Peak-to-gene link analysis is computationally intensive.")
print("Consider running it on a subset of peaks (e.g., top variable peaks)")
print("and genes (e.g., highly variable genes) to make it feasible.")
```

---

### Step 7: Biological Summary

**[REQUIRED]**

Write up what you found. This becomes the core of your final report. Answer these
questions in your notebook (as markdown cells) and then transfer the answers to your
report:

1. **What cell types are present in the PBMC sample?** List them with approximate
   percentages. Do the proportions make biological sense for healthy human blood?

2. **How do the RNA-only and ATAC-only UMAPs compare?** Are the same clusters visible in
   both? Are there cell types that are better separated by one modality than the other?

3. **What does integration add?** Compare the integrated UMAP to the single-modality
   UMAPs. Did integration resolve any ambiguities or reveal finer structure?

4. **Do the manual and automated annotations agree?** If there are disagreements, which do
   you trust more and why?

5. **(If you did Step 6) What regulatory insights did you find?** Which transcription
   factors are active in which cell types? Any interesting peak-to-gene links?

#### Required Figures for Track 1

Create a final figure panel for your report:

```python
# Create a summary figure panel
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Row 1: QC and processing
# Panel 1a: RNA QC violin plots
sc.pl.violin(adata_integrated, ['n_genes_by_counts'], ax=axes[0, 0], show=False)
axes[0, 0].set_title('RNA: Genes per Cell')

# Panel 1b: Integrated UMAP by cluster
sc.pl.umap(adata_integrated, color='leiden', ax=axes[0, 1], show=False,
           title='Integrated Clusters')

# Panel 1c: Integrated UMAP by cell type
sc.pl.umap(adata_integrated, color='cell_type', ax=axes[0, 2], show=False,
           title='Cell Type Annotation')

# Row 2: Biological insights
# Panel 2a: Marker gene dot plot (save separately -- too complex for subplot)
axes[1, 0].text(0.5, 0.5, 'See marker_dotplot.png', ha='center', va='center',
                fontsize=12, style='italic', transform=axes[1, 0].transAxes)
axes[1, 0].set_title('Marker Genes')

# Panel 2b: Cell type proportions
celltype_counts = adata_integrated.obs['cell_type'].value_counts()
axes[1, 1].barh(celltype_counts.index, celltype_counts.values, color='steelblue')
axes[1, 1].set_xlabel('Number of Cells')
axes[1, 1].set_title('Cell Type Composition')

# Panel 2c: Specific gene on UMAP
sc.pl.umap(adata_integrated, color='CD3D', ax=axes[1, 2], show=False,
           title='CD3D Expression (T cell marker)')

plt.tight_layout()
plt.savefig('figs/summary_panel.png', dpi=200, bbox_inches='tight')
plt.show()

print("Summary panel saved to figs/summary_panel.png")
```

#### Save Processed Data

```python
# Save the final annotated dataset for future use
adata_integrated.write('data/processed/integrated_annotated.h5ad')
print("Integrated dataset saved to data/processed/integrated_annotated.h5ad")
```

`[CHECKPOINT 6]` Your Track 1 analysis is complete. You have QC panels, integrated UMAPs,
cell type annotations, and a written biological summary. Proceed to the Report Template
section (Section 5) to write your final report.

---

## 4. Track 2: TCR Specificity Prediction

**Goal:** Build a machine learning model that predicts whether a T cell receptor (TCR)
binds a specific epitope (a short peptide fragment from a pathogen). This is one of the
grand challenges in computational immunology, and you will tackle a simplified version
of it.

This track draws heavily on the TCR biology you learned in Weeks 1 and 3, and introduces
machine learning concepts that are fundamental to modern computational biology.

---

### Step 1: Understand the Problem

**[REQUIRED]**

Before writing any code, make sure you understand what you are trying to predict and why
it matters.

#### The Problem

Every T cell carries a unique T Cell Receptor (TCR) on its surface. The TCR recognizes
specific peptide fragments (epitopes) presented by MHC molecules on other cells. When a
TCR binds its matching epitope, the T cell activates and mounts an immune response. You
learned about this in Week 1 (antigen presentation) and Week 3 (repertoire analysis).

The question is: **Given a TCR sequence and an epitope sequence, can we predict whether
they will bind?**

- **Input:** A TCR CDR3 sequence (the most variable part of the receptor, typically 10-20
  amino acids long) and an epitope sequence (typically 8-11 amino acids for MHC class I
  epitopes)
- **Output:** A binary prediction -- "yes, this TCR binds this epitope" or "no, it does
  not"

#### Why This Matters

If we could accurately predict TCR-epitope binding from sequence alone, we could:

- **Design better vaccines:** predict which T cells will respond to a vaccine candidate
- **Develop immunotherapies:** identify TCRs that recognize tumor-specific epitopes for
  cancer treatment
- **Diagnose infections:** determine what a patient's T cells are responding to without
  expensive wet-lab assays
- **Understand autoimmunity:** find TCRs that mistakenly target the body's own proteins

#### The Analogy

Think of it as a lock-and-key problem. The TCR is the lock, and the epitope is the key.
You are trying to predict whether a given key fits a given lock, based only on knowing the
shape of the key (epitope sequence) and the shape of the lock (TCR CDR3 sequence). The
catch is that these "shapes" are encoded as strings of amino acid letters, and the
relationship between sequence and binding is complex and nonlinear.

#### What Makes This Hard

This is considered an unsolved problem in immunology, and for good reason:

1. **Sparse data:** We only know the binding partners for a tiny fraction of all possible
   TCR-epitope pairs
2. **Imbalanced data:** For every TCR that binds a given epitope, there are millions that
   do not
3. **Generalization:** A model that memorizes known pairs is useless. We need it to predict
   for new, unseen TCRs and epitopes
4. **Sequence flexibility:** Two TCRs with very different sequences can bind the same
   epitope, and vice versa

You will encounter all of these challenges in this project.

`[CHECKPOINT 1]` You can explain the TCR specificity prediction problem in your own words:
what the input is, what the output is, why it matters, and what makes it hard.

---

### Step 2: Get the Data

**[REQUIRED]**

You will use data from **VDJdb**, a curated database of TCR sequences with known epitope
specificities. Each entry in VDJdb tells us: "This TCR CDR3 sequence was experimentally
shown to bind this epitope."

#### Download the Data

```python
import pandas as pd
import numpy as np
import os
import zipfile
import urllib.request

# --- Option A: Download from VDJdb ---
# VDJdb maintains a regularly updated database on GitLab

url = "https://gitlab.com/vdjdb/vdjdb-db/-/raw/master/latest-version.zip"
zip_path = "data/raw/vdjdb_latest.zip"
extract_dir = "data/raw/vdjdb"

# Download if not already present
if not os.path.exists(zip_path):
    print("Downloading VDJdb database...")
    urllib.request.urlretrieve(url, zip_path)
    print(f"Downloaded to {zip_path}")
else:
    print(f"Already downloaded: {zip_path}")

# Extract
os.makedirs(extract_dir, exist_ok=True)
with zipfile.ZipFile(zip_path, 'r') as z:
    z.extractall(extract_dir)
    print(f"Extracted to {extract_dir}")
    print(f"Files: {z.namelist()}")
```

```python
# Load the data
# VDJdb has a tab-separated file with many columns
# Look for a file ending in .tsv in the extracted directory

tsv_files = [f for f in os.listdir(extract_dir) if f.endswith('.tsv')]
print(f"TSV files found: {tsv_files}")

# Load the main database file (usually the largest .tsv file)
# Adjust the filename if needed based on what you see above
df = pd.read_csv(os.path.join(extract_dir, tsv_files[0]), sep='\t')
print(f"\nDataset shape: {df.shape}")
print(f"Columns: {list(df.columns)}")
```

**If the download fails:** You can manually download VDJdb from https://vdjdb.cdr3.net/
by clicking "Browse" and then the download button. Save the file in `data/raw/` and adjust
the path in the code above.

#### Explore the Data

```python
# Look at the first few rows
print(df.head())
```

```python
# Key columns we care about:
# - cdr3: the CDR3 amino acid sequence of the TCR
# - antigen.epitope: the epitope (peptide) the TCR binds
# - gene: TRA (alpha chain) or TRB (beta chain)
# - species: human, mouse, etc.
# - vdjdb.score: confidence score (0-3, higher is better)

# Let's filter to get a clean dataset
# 1. Human TCRs only
# 2. Beta chain only (TRB -- this is the most commonly sequenced chain)
# 3. Confidence score >= 1 (at least some experimental evidence)
# 4. Remove entries with missing CDR3 or epitope

df_filtered = df[
    (df['species'] == 'HomoSapiens') &
    (df['gene'] == 'TRB') &
    (df['vdjdb.score'] >= 1) &
    (df['cdr3'].notna()) &
    (df['antigen.epitope'].notna())
].copy()

# Keep only the columns we need
df_filtered = df_filtered[['cdr3', 'antigen.epitope', 'vdjdb.score',
                            'antigen.gene', 'antigen.species']].copy()
df_filtered.columns = ['cdr3', 'epitope', 'score', 'antigen_gene', 'antigen_species']

print(f"Filtered dataset: {len(df_filtered)} entries")
print(f"Unique TCR CDR3 sequences: {df_filtered['cdr3'].nunique()}")
print(f"Unique epitopes: {df_filtered['epitope'].nunique()}")
```

```python
# What are the most common epitopes?
print("Top 15 epitopes (by number of known TCR binders):")
epitope_counts = df_filtered['epitope'].value_counts()
for i, (epitope, count) in enumerate(epitope_counts.head(15).items()):
    source = df_filtered[df_filtered['epitope'] == epitope]['antigen_species'].iloc[0]
    print(f"  {i+1:2d}. {epitope:15s} (n={count:5d}) -- from {source}")
```

**Expected output:** You should see epitopes from common pathogens like CMV (NLVPMVATV),
EBV (GILGFVFTL from Influenza is also common), SARS-CoV-2, and others. The most studied
epitopes will have thousands of known TCR binders.

```python
# For a tractable ML problem, let's focus on epitopes with enough data
# We need at least ~50 known binders to learn meaningful patterns
min_binders = 50
frequent_epitopes = epitope_counts[epitope_counts >= min_binders].index.tolist()

df_ml = df_filtered[df_filtered['epitope'].isin(frequent_epitopes)].copy()
print(f"\nEpitopes with >= {min_binders} binders: {len(frequent_epitopes)}")
print(f"Entries after filtering: {len(df_ml)}")
```

#### Create Negative Examples

The database only contains positive pairs (TCR X binds epitope Y). To train a classifier,
we also need negative examples (TCR X does NOT bind epitope Y). We create these by
mismatching -- pairing a TCR with a random epitope it was not observed to bind.

```python
# Create negative pairs by mismatching TCRs and epitopes
# For each positive pair (TCR_i, epitope_j), create a negative pair
# by replacing epitope_j with a random different epitope

np.random.seed(42)  # For reproducibility

positive_pairs = df_ml[['cdr3', 'epitope']].drop_duplicates()
positive_pairs['label'] = 1
print(f"Unique positive pairs: {len(positive_pairs)}")

# Create a set of positive pairs for quick lookup
positive_set = set(zip(positive_pairs['cdr3'], positive_pairs['epitope']))

# Generate negatives
all_epitopes = positive_pairs['epitope'].unique()
negative_rows = []

for _, row in positive_pairs.iterrows():
    tcr = row['cdr3']
    true_epitope = row['epitope']

    # Pick a random epitope that this TCR is NOT known to bind
    attempts = 0
    while attempts < 100:
        random_epitope = np.random.choice(all_epitopes)
        if (tcr, random_epitope) not in positive_set:
            negative_rows.append({'cdr3': tcr, 'epitope': random_epitope, 'label': 0})
            break
        attempts += 1

negative_pairs = pd.DataFrame(negative_rows)
print(f"Negative pairs generated: {len(negative_pairs)}")

# Combine positive and negative pairs
df_dataset = pd.concat([positive_pairs, negative_pairs], ignore_index=True)
df_dataset = df_dataset.sample(frac=1, random_state=42).reset_index(drop=True)  # Shuffle

print(f"\nFinal dataset: {len(df_dataset)} pairs")
print(f"  Positive: {(df_dataset['label']==1).sum()} ({100*(df_dataset['label']==1).mean():.1f}%)")
print(f"  Negative: {(df_dataset['label']==0).sum()} ({100*(df_dataset['label']==0).mean():.1f}%)")
```

**Expected output:** A balanced dataset with roughly equal positive and negative pairs. The
total size will depend on VDJdb's current version but should be in the thousands to tens of
thousands range.

```python
# Save the processed dataset
df_dataset.to_csv('data/processed/tcr_epitope_dataset.csv', index=False)
print("Dataset saved to data/processed/tcr_epitope_dataset.csv")
```

`[CHECKPOINT 2]` Data is downloaded, filtered, and processed into a clean dataset with
positive and negative TCR-epitope pairs. Save your notebook.

---

### Step 3: Feature Engineering

**[REQUIRED]**

Machine learning models cannot read amino acid sequences directly. You need to convert
the sequences into numerical features (vectors of numbers) that the model can work with.
This process is called **feature engineering**, and it is often the most important part of
a machine learning project.

We will implement three approaches, each representing a different level of sophistication.
You should try all three and compare their performance.

#### Approach 1: K-mer Features

The simplest approach: treat the amino acid sequence like text and count subsequences
(k-mers) of length k. This is the same idea as the "bag of words" model in natural
language processing.

```python
from sklearn.feature_extraction.text import CountVectorizer

def seq_to_kmers(seq, k=3):
    """
    Convert an amino acid sequence into a string of k-mers (overlapping subsequences).

    Example: seq_to_kmers("CASSF", k=3) -> "CAS ASS SSF"

    This is like breaking a sentence into overlapping groups of k letters.
    The machine learning model will then count how often each k-mer appears.
    """
    return ' '.join([seq[i:i+k] for i in range(len(seq) - k + 1)])

# Test it
print("Example k-mer conversion:")
print(f"  Input:  CASSLGQAYEQYF")
print(f"  k=3:    {seq_to_kmers('CASSLGQAYEQYF', k=3)}")
print(f"  k=2:    {seq_to_kmers('CASSLGQAYEQYF', k=2)}")
```

```python
# Create k-mer features for both TCR and epitope, then combine them
# We concatenate the TCR k-mers and epitope k-mers so the model can learn
# patterns in both sequences

df_dataset['tcr_kmers'] = df_dataset['cdr3'].apply(lambda x: seq_to_kmers(str(x), k=3))
df_dataset['epi_kmers'] = df_dataset['epitope'].apply(lambda x: seq_to_kmers(str(x), k=3))
df_dataset['combined_kmers'] = df_dataset['tcr_kmers'] + ' ||| ' + df_dataset['epi_kmers']

# Use CountVectorizer to convert k-mer strings into feature vectors
vec_kmer = CountVectorizer(analyzer='word')
X_kmer = vec_kmer.fit_transform(df_dataset['combined_kmers'])

print(f"K-mer feature matrix shape: {X_kmer.shape}")
print(f"  ({X_kmer.shape[0]} samples x {X_kmer.shape[1]} features)")
print(f"  Each feature is a unique 3-mer (or the separator)")
```

**Expected output:** A sparse matrix with thousands of rows (one per TCR-epitope pair) and
hundreds to thousands of columns (one per unique k-mer).

#### Approach 2: Amino Acid Properties

Instead of treating amino acids as abstract letters, encode their physical and chemical
properties. Each amino acid has measurable characteristics like hydrophobicity (does it
like water?), charge (positive, negative, or neutral), and molecular weight.

```python
# Physical/chemical properties of the 20 standard amino acids
# [hydrophobicity (Kyte-Doolittle scale), charge at pH 7, molecular weight in Da]
aa_properties = {
    'A': [ 1.8,  0,  89],  # Alanine: small, nonpolar
    'R': [-4.5,  1, 174],  # Arginine: large, positive charge
    'N': [-3.5,  0, 132],  # Asparagine: polar, uncharged
    'D': [-3.5, -1, 133],  # Aspartic acid: negative charge
    'C': [ 2.5,  0, 121],  # Cysteine: can form disulfide bonds
    'E': [-3.5, -1, 147],  # Glutamic acid: negative charge
    'Q': [-3.5,  0, 146],  # Glutamine: polar, uncharged
    'G': [-0.4,  0,  75],  # Glycine: smallest amino acid
    'H': [-3.2,  0.5, 155],# Histidine: can be positive at low pH
    'I': [ 4.5,  0, 131],  # Isoleucine: nonpolar, branched
    'L': [ 3.8,  0, 131],  # Leucine: nonpolar, branched
    'K': [-3.9,  1, 146],  # Lysine: positive charge
    'M': [ 1.9,  0, 149],  # Methionine: nonpolar, contains sulfur
    'F': [ 2.8,  0, 165],  # Phenylalanine: large, aromatic
    'P': [-1.6,  0, 115],  # Proline: rigid, causes kinks
    'S': [-0.8,  0, 105],  # Serine: small, polar
    'T': [-0.7,  0, 119],  # Threonine: polar, branched
    'W': [-0.9,  0, 204],  # Tryptophan: largest amino acid, aromatic
    'Y': [-1.3,  0, 181],  # Tyrosine: aromatic, polar
    'V': [ 4.2,  0, 117],  # Valine: nonpolar, branched
}

def encode_sequence_properties(seq, max_len=20):
    """
    Encode an amino acid sequence as a vector of physical/chemical properties.

    Each amino acid is represented by 3 numbers (hydrophobicity, charge, weight).
    The sequence is padded or truncated to max_len amino acids.
    Result: a vector of length max_len * 3.
    """
    encoded = []
    for aa in str(seq)[:max_len]:
        if aa in aa_properties:
            encoded.extend(aa_properties[aa])
        else:
            encoded.extend([0, 0, 0])  # Unknown amino acids get zeros

    # Pad to fixed length if the sequence is shorter than max_len
    while len(encoded) < max_len * 3:
        encoded.append(0)

    return encoded[:max_len * 3]  # Truncate if somehow too long

# Test it
print("Encoding example:")
seq = "CASSLGQAYEQYF"
enc = encode_sequence_properties(seq, max_len=20)
print(f"  Input sequence: {seq} ({len(seq)} amino acids)")
print(f"  Encoded vector length: {len(enc)}")
print(f"  First 9 values (first 3 amino acids):")
print(f"    C: {enc[0:3]}  (hydrophobic, neutral, 121 Da)")
print(f"    A: {enc[3:6]}  (slightly hydrophobic, neutral, 89 Da)")
print(f"    S: {enc[6:9]}  (hydrophilic, neutral, 105 Da)")
```

```python
# Create property-based features for the full dataset
# Encode TCR CDR3 (max 20 amino acids) and epitope (max 15 amino acids) separately
# then concatenate

tcr_features = np.array([encode_sequence_properties(s, max_len=20)
                         for s in df_dataset['cdr3']])
epi_features = np.array([encode_sequence_properties(s, max_len=15)
                         for s in df_dataset['epitope']])

X_props = np.concatenate([tcr_features, epi_features], axis=1)

print(f"Property feature matrix shape: {X_props.shape}")
print(f"  TCR features: {tcr_features.shape[1]} (20 positions x 3 properties)")
print(f"  Epitope features: {epi_features.shape[1]} (15 positions x 3 properties)")
print(f"  Combined: {X_props.shape[1]}")
```

#### Approach 3: One-Hot Encoding

Represent each amino acid as a binary vector of length 20 (one position per possible amino
acid). This preserves the identity of each amino acid without imposing assumptions about
properties.

```python
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # 20 standard amino acids

def one_hot_encode(seq, max_len=20):
    """
    One-hot encode an amino acid sequence.

    Each amino acid becomes a binary vector of length 20 (one per possible amino acid).
    For example, Alanine (A) = [1, 0, 0, 0, ..., 0]
                 Cysteine (C) = [0, 1, 0, 0, ..., 0]

    Result: a flattened vector of length max_len * 20.
    """
    encoding = np.zeros((max_len, len(amino_acids)))
    for i, aa in enumerate(str(seq)[:max_len]):
        if aa in amino_acids:
            encoding[i, amino_acids.index(aa)] = 1
    return encoding.flatten()

# Test it
print("One-hot encoding example:")
seq = "CAS"
enc = one_hot_encode(seq, max_len=5)
print(f"  Input: {seq}")
print(f"  Vector length: {len(enc)}  (5 positions x 20 amino acids)")
# Show the encoding for position 0 (C)
enc_2d = enc.reshape(5, 20)
print(f"  Position 0 (C): {enc_2d[0].astype(int)}")
print(f"  Position 1 (A): {enc_2d[1].astype(int)}")
print(f"  Position 2 (S): {enc_2d[2].astype(int)}")
```

```python
# Create one-hot features for the full dataset
tcr_onehot = np.array([one_hot_encode(s, max_len=20) for s in df_dataset['cdr3']])
epi_onehot = np.array([one_hot_encode(s, max_len=15) for s in df_dataset['epitope']])
X_onehot = np.concatenate([tcr_onehot, epi_onehot], axis=1)

print(f"One-hot feature matrix shape: {X_onehot.shape}")
print(f"  TCR: {tcr_onehot.shape[1]} (20 positions x 20 amino acids)")
print(f"  Epitope: {epi_onehot.shape[1]} (15 positions x 20 amino acids)")
print(f"  Combined: {X_onehot.shape[1]}")
```

#### Summary of Feature Approaches

```python
print("Feature Engineering Summary:")
print(f"  Approach 1 (K-mers):     {X_kmer.shape[1]:6d} features")
print(f"  Approach 2 (Properties): {X_props.shape[1]:6d} features")
print(f"  Approach 3 (One-hot):    {X_onehot.shape[1]:6d} features")
print()
print("Each approach captures different information:")
print("  K-mers:     local sequence patterns (which 3-letter motifs are present)")
print("  Properties:  physical/chemical character (is the sequence hydrophobic? charged?)")
print("  One-hot:     exact amino acid identity at each position")
```

`[CHECKPOINT 3]` You have created three different feature representations of your
TCR-epitope pairs. Save your notebook.

---

### Step 4: Build a Baseline Model

**[REQUIRED]**

Start with the simplest reasonable model: logistic regression. This gives you a baseline
to compare against. In machine learning, you always want to know: "How much better is my
fancy model than a simple one?"

```python
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
import time

# Get labels
y = df_dataset['label'].values

# We will test each feature set with logistic regression
feature_sets = {
    'K-mers': X_kmer,
    'AA Properties': X_props,
    'One-hot': X_onehot,
}

# 5-fold stratified cross-validation
# "Stratified" means each fold has the same ratio of positives to negatives
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

results = {}

for feat_name, X in feature_sets.items():
    print(f"\n--- Logistic Regression with {feat_name} features ---")
    aucs, aps = [], []
    start_time = time.time()

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        # Handle both sparse and dense matrices
        if hasattr(X, 'toarray'):
            X_train, X_test = X[train_idx], X[test_idx]
        else:
            X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Train logistic regression
        model = LogisticRegression(max_iter=300, random_state=42, solver='lbfgs')
        model.fit(X_train, y_train)

        # Predict probabilities
        y_pred_proba = model.predict_proba(X_test)[:, 1]

        # Evaluate
        auc = roc_auc_score(y_test, y_pred_proba)
        ap = average_precision_score(y_test, y_pred_proba)
        aucs.append(auc)
        aps.append(ap)
        print(f"  Fold {fold+1}: AUROC={auc:.3f}, AP={ap:.3f}")

    elapsed = time.time() - start_time
    print(f"  Mean AUROC: {np.mean(aucs):.3f} +/- {np.std(aucs):.3f}")
    print(f"  Mean AP:    {np.mean(aps):.3f} +/- {np.std(aps):.3f}")
    print(f"  Time: {elapsed:.1f}s")

    results[feat_name] = {
        'model': 'Logistic Regression',
        'auroc_mean': np.mean(aucs),
        'auroc_std': np.std(aucs),
        'ap_mean': np.mean(aps),
        'ap_std': np.std(aps),
    }
```

**Expected output:** AUROC values for logistic regression. For random guessing, AUROC
would be 0.5. You should see values significantly above 0.5 (probably 0.7-0.9 depending
on the feature set). K-mer features often perform surprisingly well on this task.

#### Interpret the Metrics

A brief explanation of the two metrics you are using:

- **AUROC (Area Under the ROC Curve):** Measures how well the model separates positives
  from negatives. 0.5 = random guessing, 1.0 = perfect. A value of 0.85 means that a
  randomly chosen positive example will be scored higher than a randomly chosen negative
  example 85% of the time.

- **Average Precision (AP):** Similar to AUROC but more sensitive to performance on the
  positive class. This matters when positives are rare (which they often are in real TCR
  prediction problems). Higher is better.

`[CHECKPOINT 4]` You have a baseline model trained and evaluated. You know which feature
set works best with logistic regression. Save your notebook.

---

### Step 5: Improve with Better Models

**[REQUIRED]**

Now try more powerful models. Tree-based models (Random Forest, Gradient Boosting) often
perform well on tabular data like this.

```python
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

# We will use the best feature set from Step 4
# (Or try all of them -- your choice)
# For this example, we'll use the k-mer features
X_best = X_kmer  # Change this based on your Step 4 results
best_feat_name = 'K-mers'

models_to_try = {
    'Logistic Regression': LogisticRegression(max_iter=300, random_state=42),
    'Random Forest': RandomForestClassifier(
        n_estimators=200,     # Number of trees in the forest
        max_depth=None,       # Trees can grow as deep as needed
        min_samples_leaf=5,   # Prevent overfitting on very small groups
        random_state=42,
        n_jobs=-1             # Use all CPU cores
    ),
    'Gradient Boosting': GradientBoostingClassifier(
        n_estimators=200,     # Number of boosting stages
        learning_rate=0.1,    # How much each tree contributes
        max_depth=5,          # Shallow trees (boosting prefers shallow trees)
        random_state=42
    ),
}

# Try to import XGBoost (a faster, more powerful version of Gradient Boosting)
try:
    from xgboost import XGBClassifier
    models_to_try['XGBoost'] = XGBClassifier(
        n_estimators=200,
        learning_rate=0.1,
        max_depth=5,
        random_state=42,
        eval_metric='logloss',
        use_label_encoder=False,
    )
    print("XGBoost is available and will be tested.")
except ImportError:
    print("XGBoost not installed. Skipping. Install with: pip install xgboost")

print(f"\nModels to evaluate: {list(models_to_try.keys())}")
print(f"Feature set: {best_feat_name} ({X_best.shape[1]} features)")
```

```python
# Evaluate all models with cross-validation
all_results = {}

for model_name, model in models_to_try.items():
    print(f"\n--- {model_name} ---")
    aucs, aps = [], []
    start_time = time.time()

    for fold, (train_idx, test_idx) in enumerate(skf.split(X_best, y)):
        if hasattr(X_best, 'toarray'):
            X_train, X_test = X_best[train_idx], X_best[test_idx]
        else:
            X_train, X_test = X_best[train_idx], X_best[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        model_copy = model.__class__(**model.get_params())
        model_copy.fit(X_train, y_train)

        y_pred_proba = model_copy.predict_proba(X_test)[:, 1]
        auc = roc_auc_score(y_test, y_pred_proba)
        ap = average_precision_score(y_test, y_pred_proba)
        aucs.append(auc)
        aps.append(ap)
        print(f"  Fold {fold+1}: AUROC={auc:.3f}, AP={ap:.3f}")

    elapsed = time.time() - start_time
    mean_auc = np.mean(aucs)
    mean_ap = np.mean(aps)
    print(f"  Mean AUROC: {mean_auc:.3f} +/- {np.std(aucs):.3f}")
    print(f"  Mean AP:    {mean_ap:.3f} +/- {np.std(aps):.3f}")
    print(f"  Time: {elapsed:.1f}s")

    all_results[model_name] = {
        'auroc_mean': mean_auc,
        'auroc_std': np.std(aucs),
        'ap_mean': mean_ap,
        'ap_std': np.std(aps),
        'time': elapsed
    }
```

```python
# Create a comparison table
print("\n" + "="*70)
print("MODEL COMPARISON TABLE")
print("="*70)
print(f"{'Model':<25} {'AUROC':>12} {'AP':>12} {'Time':>8}")
print("-"*70)

for model_name, res in sorted(all_results.items(),
                               key=lambda x: x[1]['auroc_mean'], reverse=True):
    print(f"{model_name:<25} "
          f"{res['auroc_mean']:.3f}+/-{res['auroc_std']:.3f} "
          f"{res['ap_mean']:.3f}+/-{res['ap_std']:.3f} "
          f"{res['time']:>6.1f}s")

print("="*70)
best_model = max(all_results, key=lambda x: all_results[x]['auroc_mean'])
print(f"\nBest model: {best_model} "
      f"(AUROC={all_results[best_model]['auroc_mean']:.3f})")
```

**Expected output:** A table comparing all models. Tree-based models (Random Forest,
Gradient Boosting, XGBoost) usually outperform logistic regression on this task, sometimes
by a substantial margin. XGBoost or Gradient Boosting typically gives the best results.

`[CHECKPOINT 5]` You have compared multiple models and identified the best performer. Save
your notebook.

---

### Step 6: Evaluate and Visualize

**[REQUIRED]**

Now create publication-quality evaluation plots for your best model. These are the figures
that will go in your report.

```python
from sklearn.metrics import (roc_curve, precision_recall_curve,
                             confusion_matrix, ConfusionMatrixDisplay)

# Train the best model on the full training set and evaluate on a held-out set
# Use the last CV fold for final visualization

# Retrain on the last fold's train/test split
last_train_idx, last_test_idx = list(skf.split(X_best, y))[-1]

if hasattr(X_best, 'toarray'):
    X_train = X_best[last_train_idx]
    X_test = X_best[last_test_idx]
else:
    X_train = X_best[last_train_idx]
    X_test = X_best[last_test_idx]

y_train = y[last_train_idx]
y_test = y[last_test_idx]

# Retrain the best model
# (Replace with your best model class and parameters)
final_model = RandomForestClassifier(n_estimators=200, min_samples_leaf=5,
                                     random_state=42, n_jobs=-1)
final_model.fit(X_train, y_train)

y_pred_proba = final_model.predict_proba(X_test)[:, 1]
y_pred = (y_pred_proba >= 0.5).astype(int)

final_auc = roc_auc_score(y_test, y_pred_proba)
final_ap = average_precision_score(y_test, y_pred_proba)
print(f"Final model on held-out fold:")
print(f"  AUROC: {final_auc:.3f}")
print(f"  AP:    {final_ap:.3f}")
```

#### ROC and Precision-Recall Curves

```python
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel 1: ROC Curve
fpr, tpr, thresholds_roc = roc_curve(y_test, y_pred_proba)
axes[0].plot(fpr, tpr, color='steelblue', lw=2,
             label=f'Model (AUC={final_auc:.3f})')
axes[0].plot([0, 1], [0, 1], 'k--', lw=1, label='Random (AUC=0.500)')
axes[0].set_xlabel('False Positive Rate', fontsize=12)
axes[0].set_ylabel('True Positive Rate', fontsize=12)
axes[0].set_title('ROC Curve', fontsize=14)
axes[0].legend(fontsize=10)
axes[0].set_xlim([0, 1])
axes[0].set_ylim([0, 1.02])

# Panel 2: Precision-Recall Curve
precision, recall, thresholds_pr = precision_recall_curve(y_test, y_pred_proba)
axes[1].plot(recall, precision, color='coral', lw=2,
             label=f'Model (AP={final_ap:.3f})')
# Baseline: fraction of positives
baseline = y_test.mean()
axes[1].axhline(y=baseline, color='k', linestyle='--', lw=1,
                label=f'Random (AP={baseline:.3f})')
axes[1].set_xlabel('Recall', fontsize=12)
axes[1].set_ylabel('Precision', fontsize=12)
axes[1].set_title('Precision-Recall Curve', fontsize=14)
axes[1].legend(fontsize=10)
axes[1].set_xlim([0, 1])
axes[1].set_ylim([0, 1.02])

plt.tight_layout()
plt.savefig('figs/model_evaluation.png', dpi=150, bbox_inches='tight')
plt.show()

print("Figure saved to figs/model_evaluation.png")
```

**Expected output:** Two curves. The ROC curve should bow toward the upper-left corner
(farther from the diagonal = better). The Precision-Recall curve should stay high as recall
increases (farther from the bottom-right = better).

#### Confusion Matrix

```python
fig, ax = plt.subplots(figsize=(6, 5))

cm = confusion_matrix(y_test, y_pred)
disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                               display_labels=['Non-binder', 'Binder'])
disp.plot(ax=ax, cmap='Blues', values_format='d')
ax.set_title('Confusion Matrix', fontsize=14)

plt.tight_layout()
plt.savefig('figs/confusion_matrix.png', dpi=150, bbox_inches='tight')
plt.show()

# Print confusion matrix interpretation
tn, fp, fn, tp = cm.ravel()
print(f"\nConfusion Matrix Interpretation:")
print(f"  True Negatives (correctly predicted non-binder):  {tn}")
print(f"  False Positives (predicted binder, actually not): {fp}")
print(f"  False Negatives (predicted non-binder, actually is): {fn}")
print(f"  True Positives (correctly predicted binder):      {tp}")
print(f"\n  Accuracy: {(tp+tn)/(tp+tn+fp+fn):.3f}")
print(f"  Precision: {tp/(tp+fp):.3f}  (of predicted binders, how many are real)")
print(f"  Recall: {tp/(tp+fn):.3f}  (of real binders, how many did we find)")
```

#### Feature Importance

```python
# What did the model learn? Which features are most important?
# This only works for tree-based models

if hasattr(final_model, 'feature_importances_'):
    importances = final_model.feature_importances_

    # Get feature names
    if hasattr(X_best, 'toarray'):
        # For k-mer features from CountVectorizer
        feature_names = vec_kmer.get_feature_names_out()
    else:
        # For property or one-hot features, create generic names
        feature_names = [f'feature_{i}' for i in range(X_best.shape[1])]

    # Sort by importance and show top 20
    top_n = 20
    top_indices = np.argsort(importances)[-top_n:][::-1]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(range(top_n), importances[top_indices], color='steelblue')
    ax.set_yticks(range(top_n))
    ax.set_yticklabels([feature_names[i] for i in top_indices])
    ax.set_xlabel('Feature Importance', fontsize=12)
    ax.set_title(f'Top {top_n} Most Important Features', fontsize=14)
    ax.invert_yaxis()  # Most important at top

    plt.tight_layout()
    plt.savefig('figs/feature_importance.png', dpi=150, bbox_inches='tight')
    plt.show()

    print("\nTop 10 features:")
    for i, idx in enumerate(top_indices[:10]):
        print(f"  {i+1:2d}. {feature_names[idx]:15s} (importance: {importances[idx]:.4f})")
else:
    print("Feature importance not available for this model type.")
```

**Expected output:** A bar chart showing which features the model relies on most. For
k-mer features, look for biologically meaningful patterns -- certain 3-mers may correspond
to known TCR binding motifs.

#### Model Comparison Visualization

```python
# Create a visual comparison of all models tested

model_names = list(all_results.keys())
aurocs = [all_results[m]['auroc_mean'] for m in model_names]
auroc_stds = [all_results[m]['auroc_std'] for m in model_names]
aps_vals = [all_results[m]['ap_mean'] for m in model_names]
ap_stds = [all_results[m]['ap_std'] for m in model_names]

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# AUROC comparison
x = np.arange(len(model_names))
axes[0].bar(x, aurocs, yerr=auroc_stds, color='steelblue', capsize=5,
            edgecolor='white')
axes[0].set_xticks(x)
axes[0].set_xticklabels(model_names, rotation=30, ha='right')
axes[0].set_ylabel('AUROC', fontsize=12)
axes[0].set_title('Model Comparison: AUROC', fontsize=14)
axes[0].axhline(y=0.5, color='red', linestyle='--', linewidth=1, label='Random')
axes[0].set_ylim([0.4, 1.0])
axes[0].legend()

# AP comparison
axes[1].bar(x, aps_vals, yerr=ap_stds, color='coral', capsize=5,
            edgecolor='white')
axes[1].set_xticks(x)
axes[1].set_xticklabels(model_names, rotation=30, ha='right')
axes[1].set_ylabel('Average Precision', fontsize=12)
axes[1].set_title('Model Comparison: Average Precision', fontsize=14)
axes[1].set_ylim([0.3, 1.0])

plt.tight_layout()
plt.savefig('figs/model_comparison.png', dpi=150, bbox_inches='tight')
plt.show()

print("Figure saved to figs/model_comparison.png")
```

`[CHECKPOINT 6]` You have ROC curves, PR curves, a confusion matrix, feature importance
plots, and a model comparison chart. Save your notebook.

---

### Step 7: Proper Data Splitting (Critical)

**[REQUIRED]**

This is the most important conceptual step in the entire project. Everything you did above
has a hidden flaw that you need to understand and fix.

#### The Problem: Data Leakage

In Steps 4-6, you used `StratifiedKFold` to split data into train and test sets. This
split was **random** -- any TCR-epitope pair could end up in either the train or test set.

Here is the problem: if the same epitope (or even the same TCR) appears in both the
training set and the test set, the model can "cheat." Instead of learning general rules
about how TCRs recognize epitopes, it memorizes which TCRs go with which epitopes.

**Example of leakage:**

```
Training set:  TCR_1 + NLVPMVATV -> binds (positive)
               TCR_2 + NLVPMVATV -> binds (positive)
               TCR_3 + GILGFVFTL -> binds (positive)

Test set:      TCR_4 + NLVPMVATV -> binds (positive)    <-- LEAKAGE!
               TCR_5 + GILGFVFTL -> binds (positive)    <-- LEAKAGE!
```

The model has already seen the epitopes NLVPMVATV and GILGFVFTL during training. It can
score highly on the test set just by recognizing "oh, this is an NLVPMVATV-binding TCR,
and I saw lots of those in training." This is not real generalization.

In the real world, we want our model to work on **new epitopes** it has never seen before
(e.g., from a newly emerged virus). The random split gives an inflated performance
estimate.

#### The Fix: Split by Epitope

Use `GroupKFold` to ensure that all examples of a given epitope go into either the training
set or the test set, never both. This forces the model to generalize to new epitopes.

```python
from sklearn.model_selection import GroupKFold

# Groups = epitope identity
# All pairs involving epitope X will be in the same fold
groups = df_dataset['epitope'].values

gkf = GroupKFold(n_splits=5)

print("="*60)
print("FAIR EVALUATION: Split by Epitope")
print("="*60)

fair_results = {}

for model_name, model in models_to_try.items():
    print(f"\n--- {model_name} (epitope-split) ---")
    aucs, aps = [], []

    for fold, (train_idx, test_idx) in enumerate(gkf.split(X_best, y, groups)):
        if hasattr(X_best, 'toarray'):
            X_train, X_test = X_best[train_idx], X_best[test_idx]
        else:
            X_train, X_test = X_best[train_idx], X_best[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Show which epitopes are in train vs test for the first fold
        if fold == 0:
            train_epitopes = set(df_dataset.iloc[train_idx]['epitope'])
            test_epitopes = set(df_dataset.iloc[test_idx]['epitope'])
            overlap = train_epitopes.intersection(test_epitopes)
            print(f"  Fold 1: {len(train_epitopes)} train epitopes, "
                  f"{len(test_epitopes)} test epitopes, "
                  f"{len(overlap)} overlap (should be 0)")

        model_copy = model.__class__(**model.get_params())
        model_copy.fit(X_train, y_train)

        y_pred_proba = model_copy.predict_proba(X_test)[:, 1]
        auc = roc_auc_score(y_test, y_pred_proba)
        ap = average_precision_score(y_test, y_pred_proba)
        aucs.append(auc)
        aps.append(ap)
        print(f"  Fold {fold+1}: AUROC={auc:.3f}, AP={ap:.3f}")

    mean_auc = np.mean(aucs)
    mean_ap = np.mean(aps)
    print(f"  Mean AUROC: {mean_auc:.3f} +/- {np.std(aucs):.3f}")
    print(f"  Mean AP:    {mean_ap:.3f} +/- {np.std(aps):.3f}")

    fair_results[model_name] = {
        'auroc_mean': mean_auc,
        'auroc_std': np.std(aucs),
        'ap_mean': mean_ap,
        'ap_std': np.std(aps),
    }
```

```python
# Compare random split vs epitope split
print("\n" + "="*70)
print("COMPARISON: Random Split vs Epitope Split")
print("="*70)
print(f"{'Model':<25} {'Random AUROC':>14} {'Epitope AUROC':>14} {'Drop':>8}")
print("-"*70)

for model_name in all_results:
    if model_name in fair_results:
        random_auc = all_results[model_name]['auroc_mean']
        fair_auc = fair_results[model_name]['auroc_mean']
        drop = random_auc - fair_auc
        print(f"{model_name:<25} "
              f"{random_auc:>11.3f}    "
              f"{fair_auc:>11.3f}    "
              f"{drop:>+6.3f}")

print("="*70)
print("\nNotice: the 'fair' (epitope-split) AUROC is LOWER than the random split.")
print("This is expected and honest. The random split was overestimating performance")
print("because the model could memorize epitope-specific patterns from training.")
print("\nThe epitope-split score reflects real-world generalization:")
print("how well the model predicts binding for UNSEEN epitopes.")
```

**Expected output:** The epitope-split AUROC will be noticeably lower (often 0.05-0.15
lower) than the random-split AUROC. This drop is the data leakage you just eliminated. The
epitope-split score is the honest estimate of how well your model would work in practice.

#### Why This Matters Beyond This Project

This data leakage problem is not unique to TCR prediction. It appears in many areas of
computational biology:

- Drug-target interaction prediction (split by target, not randomly)
- Protein function prediction (split by protein family)
- Patient outcome prediction (split by hospital or time period)

Understanding data leakage and knowing how to prevent it is one of the most important
skills in machine learning for biology. Many published papers have been criticized for
inflated results due to improper data splitting.

`[CHECKPOINT 7]` You understand data leakage, have implemented a proper epitope-split
evaluation, and can explain why the fair score is lower. Save your notebook.

---

### Step 8: Deep Learning Approach (Challenge)

**[CHALLENGE]**

This section introduces a simple neural network for TCR-epitope prediction. Only attempt
this if you are comfortable with everything in Steps 1-7 and have time to spare.

Deep learning can potentially learn complex, nonlinear patterns in sequence data that
simpler models miss. However, it requires more code, more data, and more tuning.

```python
# This section requires PyTorch
# Install if needed: pip install torch

try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import Dataset, DataLoader
    print(f"PyTorch version: {torch.__version__}")
    print(f"GPU available: {torch.cuda.is_available()}")
except ImportError:
    print("PyTorch is not installed.")
    print("Install with: pip install torch")
    print("Then restart your notebook and re-run this cell.")
```

#### Define the Dataset

```python
class TCREpitopeDataset(Dataset):
    """
    A PyTorch dataset that holds TCR-epitope pairs.
    Each sequence is converted to integer indices (one per amino acid).
    """
    def __init__(self, cdr3_seqs, epitope_seqs, labels, max_tcr_len=20, max_epi_len=15):
        self.aa_to_idx = {aa: i+1 for i, aa in enumerate('ACDEFGHIKLMNPQRSTVWY')}
        # Index 0 is reserved for padding
        self.max_tcr_len = max_tcr_len
        self.max_epi_len = max_epi_len
        self.cdr3_seqs = cdr3_seqs
        self.epitope_seqs = epitope_seqs
        self.labels = labels

    def encode(self, seq, max_len):
        """Convert amino acid sequence to integer indices."""
        indices = [self.aa_to_idx.get(aa, 0) for aa in str(seq)[:max_len]]
        # Pad with zeros
        indices += [0] * (max_len - len(indices))
        return torch.tensor(indices, dtype=torch.long)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        tcr = self.encode(self.cdr3_seqs[idx], self.max_tcr_len)
        epi = self.encode(self.epitope_seqs[idx], self.max_epi_len)
        label = torch.tensor(self.labels[idx], dtype=torch.float32)
        return tcr, epi, label
```

#### Define the Model

```python
class SimpleTCRModel(nn.Module):
    """
    A simple 1D CNN model for TCR-epitope binding prediction.

    Architecture:
    1. Embedding layer: converts amino acid indices to dense vectors
    2. 1D Convolution: captures local sequence patterns (like k-mers but learned)
    3. Adaptive max pooling: get a fixed-size representation regardless of length
    4. Fully connected layers: combine TCR and epitope representations
    5. Sigmoid output: probability of binding
    """
    def __init__(self, vocab_size=21, embed_dim=32, hidden_dim=64):
        super().__init__()
        # Embedding: each amino acid becomes a vector of length embed_dim
        self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=0)

        # 1D Convolution: slides a filter across the sequence
        # kernel_size=3 means it looks at 3 amino acids at a time (like 3-mers!)
        self.conv1 = nn.Conv1d(embed_dim, hidden_dim, kernel_size=3, padding=1)
        self.conv2 = nn.Conv1d(hidden_dim, hidden_dim, kernel_size=3, padding=1)

        # Max pooling: take the strongest signal from anywhere in the sequence
        self.pool = nn.AdaptiveMaxPool1d(1)

        # Fully connected layers to combine TCR and epitope
        self.fc1 = nn.Linear(hidden_dim * 2, 64)  # *2 because TCR + epitope
        self.fc2 = nn.Linear(64, 1)
        self.dropout = nn.Dropout(0.3)

    def encode_sequence(self, x):
        """Encode a sequence through embedding -> conv -> pool."""
        x = self.embedding(x)          # (batch, seq_len, embed_dim)
        x = x.permute(0, 2, 1)         # (batch, embed_dim, seq_len) for Conv1d
        x = torch.relu(self.conv1(x))  # (batch, hidden_dim, seq_len)
        x = torch.relu(self.conv2(x))  # (batch, hidden_dim, seq_len)
        x = self.pool(x).squeeze(-1)   # (batch, hidden_dim)
        return x

    def forward(self, tcr, epitope):
        tcr_feat = self.encode_sequence(tcr)
        epi_feat = self.encode_sequence(epitope)
        combined = torch.cat([tcr_feat, epi_feat], dim=1)  # (batch, hidden_dim*2)
        x = self.dropout(torch.relu(self.fc1(combined)))
        x = torch.sigmoid(self.fc2(x))
        return x.squeeze(-1)
```

#### Training Loop

```python
def train_model(model, train_loader, val_loader, n_epochs=30, lr=0.001):
    """Train the model and track performance on validation set."""
    optimizer = optim.Adam(model.parameters(), lr=lr)
    criterion = nn.BCELoss()  # Binary Cross-Entropy for binary classification
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)

    history = {'train_loss': [], 'val_auc': []}

    for epoch in range(n_epochs):
        # Training
        model.train()
        epoch_loss = 0
        for tcr, epi, labels in train_loader:
            tcr, epi, labels = tcr.to(device), epi.to(device), labels.to(device)
            optimizer.zero_grad()
            outputs = model(tcr, epi)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        avg_loss = epoch_loss / len(train_loader)
        history['train_loss'].append(avg_loss)

        # Validation
        model.eval()
        all_preds, all_labels = [], []
        with torch.no_grad():
            for tcr, epi, labels in val_loader:
                tcr, epi, labels = tcr.to(device), epi.to(device), labels.to(device)
                outputs = model(tcr, epi)
                all_preds.extend(outputs.cpu().numpy())
                all_labels.extend(labels.cpu().numpy())

        val_auc = roc_auc_score(all_labels, all_preds)
        history['val_auc'].append(val_auc)

        if (epoch + 1) % 5 == 0:
            print(f"  Epoch {epoch+1}/{n_epochs}: "
                  f"loss={avg_loss:.4f}, val_AUROC={val_auc:.3f}")

    return history

# Prepare data for PyTorch
# Use the last fold split from the epitope-split evaluation
last_train_idx, last_test_idx = list(gkf.split(X_best, y, groups))[-1]

train_dataset = TCREpitopeDataset(
    df_dataset.iloc[last_train_idx]['cdr3'].values,
    df_dataset.iloc[last_train_idx]['epitope'].values,
    y[last_train_idx]
)

test_dataset = TCREpitopeDataset(
    df_dataset.iloc[last_test_idx]['cdr3'].values,
    df_dataset.iloc[last_test_idx]['epitope'].values,
    y[last_test_idx]
)

train_loader = DataLoader(train_dataset, batch_size=128, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=128, shuffle=False)

# Create and train the model
model = SimpleTCRModel(vocab_size=21, embed_dim=32, hidden_dim=64)
print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")
print(f"Training samples: {len(train_dataset)}")
print(f"Test samples: {len(test_dataset)}")
print()

history = train_model(model, train_loader, test_loader, n_epochs=30, lr=0.001)
```

```python
# Plot training history
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].plot(history['train_loss'], color='steelblue')
axes[0].set_xlabel('Epoch')
axes[0].set_ylabel('Training Loss')
axes[0].set_title('Training Loss')

axes[1].plot(history['val_auc'], color='coral')
axes[1].set_xlabel('Epoch')
axes[1].set_ylabel('Validation AUROC')
axes[1].set_title('Validation AUROC')
axes[1].axhline(y=fair_results.get(best_model, {}).get('auroc_mean', 0.5),
                color='gray', linestyle='--', label='Best sklearn model')
axes[1].legend()

plt.tight_layout()
plt.savefig('figs/deep_learning_training.png', dpi=150, bbox_inches='tight')
plt.show()

print(f"Final deep learning AUROC: {history['val_auc'][-1]:.3f}")
print(f"Best sklearn AUROC (epitope-split): "
      f"{fair_results.get(best_model, {}).get('auroc_mean', 'N/A')}")
```

**Expected output:** Training loss should decrease over epochs. Validation AUROC should
increase and then plateau (or even decrease if the model starts overfitting). Compare the
final deep learning AUROC with your best sklearn model from Step 7.

Note: With this small dataset and simple architecture, the deep learning model may not
outperform Gradient Boosting. That is normal. Deep learning shines with large datasets and
complex architectures, which are beyond the scope of this bootcamp.

---

## 5. Report Template

**[REQUIRED for both tracks]**

Create a file called `reports/final_report.md` in your project directory. Use the template
below as a starting point. Replace the placeholder text with your own content.

Your report should be at least 2 pages when printed (roughly 800-1200 words, plus
figures). Write in complete sentences. Do not just paste code output -- explain what it
means.

```markdown
# [Your Title Here]
# e.g., "Joint scRNA-seq and scATAC-seq Analysis of Human PBMCs"
# or "Predicting TCR-Epitope Binding Using Machine Learning"

**Author:** [Your Name]
**Date:** [Date]
**Bootcamp:** Computational Immunology Capstone Project

---

## 1. Background and Objectives

[2-3 paragraphs explaining:]
[- What biological question are you addressing?]
[- Why does this question matter for immunology or medicine?]
[- What is your specific objective in this project?]
[- Briefly mention what approaches you will use.]

## 2. Data

[Describe the dataset(s) you used:]
[- Where did you download the data from? Include URLs.]
[- How many cells/samples/entries?]
[- What organism and tissue/source?]
[- What filtering or quality control did you apply? What were the thresholds?]
[- How many cells/entries remained after filtering?]

## 3. Methods

[Describe each step of your analysis in order:]
[- For Track 1: RNA processing, ATAC processing, integration method, cell type annotation]
[- For Track 2: feature engineering, model types, cross-validation strategy, data splitting]
[- Include specific parameter choices (e.g., "I used 30 principal components",]
[  "I trained a Random Forest with 200 trees")]
[- Explain WHY you made key choices, not just what you did]

## 4. Results

[Present your main findings. For each figure:]
[- Reference the figure by name (e.g., "Figure 1 shows..."")]
[- Describe what the figure shows]
[- State the key finding or observation]

[For Track 1:]
[- How many cells passed QC?]
[- How many clusters did you find?]
[- What cell types are present and in what proportions?]
[- How does the integrated UMAP compare to single-modality UMAPs?]

[For Track 2:]
[- What is the baseline performance (logistic regression)?]
[- Which model performed best?]
[- How much did performance drop with the proper epitope-split evaluation?]
[- What features were most important?]

## 5. Discussion

[2-3 paragraphs reflecting on:]
[- What did you learn from this analysis?]
[- What were the main limitations?]
[  - Track 1: small sample size, single donor, limited cell types, etc.]
[  - Track 2: data leakage, limited epitope diversity, simple features, etc.]
[- How do your results compare to published work?]
[- What surprised you?]

## 6. Reproducibility

[- Conda environment: list the environment.yml file or key package versions]
[- Random seeds: what seed did you use?]
[- Data version: when did you download the data? What version?]
[- Hardware: what computer did you run this on? How long did the analysis take?]

## 7. Next Steps

[What would you do with more time? Some ideas:]
[- Track 1: add more donors, use SCENIC+ for regulatory networks, validate with ChIP-seq]
[- Track 2: try protein language model embeddings (ESM), use attention-based architectures,]
[  train on more epitopes, add structural features]

---

## Figures

[List all figures with captions:]

**Figure 1.** [Description of figure 1]

**Figure 2.** [Description of figure 2]

[... etc.]
```

---

## 6. Minimum Figure Set

These are the figures you must produce. Save all of them in the `figs/` directory with
descriptive filenames.

### Track 1: RNA + ATAC Integration

| Figure | Filename | Description |
|--------|----------|-------------|
| 1 | `rna_qc.png` | Three-panel histogram of RNA QC metrics (genes per cell, total counts, % mitochondrial) |
| 2 | `rna_umap_preliminary.png` | RNA-only UMAP colored by Leiden cluster |
| 3 | `atac_umap_preliminary.png` | ATAC-only UMAP colored by Leiden cluster |
| 4 | `integrated_umap.png` | Integrated RNA+ATAC UMAP colored by Leiden cluster |
| 5 | `umap_comparison.png` | Side-by-side comparison of RNA-only, ATAC-only, and integrated UMAPs |
| 6 | `cell_type_umap.png` | Integrated UMAP colored by cell type annotation |
| 7 | `marker_dotplot.png` | Dot plot of marker genes across clusters |
| 8 | `summary_panel.png` | Multi-panel summary figure for the report |
| (optional) | `gene_activity_heatmap.png` | Gene activity heatmap for key immune markers |
| (optional) | `tf_motif_enrichment.png` | TF motif enrichment across cell types |
| (optional) | `peak_gene_links.png` | Peak-to-gene correlation examples |

### Track 2: TCR Specificity Prediction

| Figure | Filename | Description |
|--------|----------|-------------|
| 1 | `epitope_distribution.png` | Bar chart of top epitopes by number of known TCR binders |
| 2 | `model_evaluation.png` | Two-panel figure: ROC curve and Precision-Recall curve |
| 3 | `confusion_matrix.png` | Confusion matrix for the best model |
| 4 | `model_comparison.png` | Bar chart comparing AUROC and AP across models |
| 5 | `feature_importance.png` | Top 20 most important features |
| 6 | `fair_vs_random.png` | Comparison of random-split vs epitope-split performance |
| (optional) | `deep_learning_training.png` | Training loss and validation AUROC curves |
| (optional) | `per_epitope_performance.png` | AUROC broken down by epitope |
| (optional) | `calibration_plot.png` | Calibration curve (predicted probability vs actual fraction) |

To create the optional epitope distribution figure (Track 2) and the fair-vs-random
comparison figure:

```python
# Epitope distribution figure
fig, ax = plt.subplots(figsize=(10, 6))
top_epitopes = df_filtered['epitope'].value_counts().head(15)
ax.barh(range(len(top_epitopes)), top_epitopes.values, color='steelblue')
ax.set_yticks(range(len(top_epitopes)))
ax.set_yticklabels(top_epitopes.index)
ax.set_xlabel('Number of Known TCR Binders')
ax.set_title('Most Studied Epitopes in VDJdb')
ax.invert_yaxis()
plt.tight_layout()
plt.savefig('figs/epitope_distribution.png', dpi=150, bbox_inches='tight')
plt.show()
```

```python
# Fair vs Random comparison figure
fig, ax = plt.subplots(figsize=(10, 5))

model_names = [m for m in all_results if m in fair_results]
x = np.arange(len(model_names))
width = 0.35

random_aurocs = [all_results[m]['auroc_mean'] for m in model_names]
fair_aurocs = [fair_results[m]['auroc_mean'] for m in model_names]

bars1 = ax.bar(x - width/2, random_aurocs, width, label='Random Split',
               color='steelblue', edgecolor='white')
bars2 = ax.bar(x + width/2, fair_aurocs, width, label='Epitope Split (fair)',
               color='coral', edgecolor='white')

ax.set_ylabel('AUROC', fontsize=12)
ax.set_title('Effect of Proper Data Splitting on Performance', fontsize=14)
ax.set_xticks(x)
ax.set_xticklabels(model_names, rotation=30, ha='right')
ax.legend()
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=1)
ax.set_ylim([0.4, 1.0])

plt.tight_layout()
plt.savefig('figs/fair_vs_random.png', dpi=150, bbox_inches='tight')
plt.show()
```

---

## 7. Self-Check and Completion Criteria

Go through this checklist before you consider your capstone complete.

### Technical Checklist

- [ ] My notebook runs from top to bottom without errors (Kernel -> Restart & Run All)
- [ ] All figures are saved as PNG files in the `figs/` directory
- [ ] My processed data is saved in the `data/processed/` directory
- [ ] My `environment.yml` accurately reflects the packages I used
- [ ] I set random seeds wherever randomness is involved (for reproducibility)

### Report Checklist

- [ ] My report has all 7 sections filled in
- [ ] My report is at least 2 pages with all sections completed
- [ ] Every figure is referenced and described in the report
- [ ] I explained my parameter choices (not just listed them)
- [ ] I discussed limitations honestly

### Understanding Checklist

- [ ] I can explain every step of my analysis in my own words
- [ ] I understand why each QC filter was applied (Track 1)
- [ ] I understand the difference between TF-IDF and log-normalization (Track 1)
- [ ] I understand what data leakage is and why epitope-split evaluation matters (Track 2)
- [ ] I can explain what AUROC and Average Precision measure (Track 2)
- [ ] I know what the main limitations of my analysis are

### Final File Inventory

When you are done, your project directory should contain:

```
~/bootcamp/final_project/
  data/
    raw/                    # Downloaded data files
    processed/              # Saved processed data (.h5ad or .csv)
  notebooks/
    track1_rna_atac_integration.ipynb   (and/or)
    track2_tcr_specificity.ipynb
  src/
    utils.py                # Any helper functions you created
  figs/
    rna_qc.png              # (Track 1) or epitope_distribution.png (Track 2)
    ...                     # All other figures from the minimum figure set
  reports/
    final_report.md         # Your written report
  environment.yml           # Conda environment specification
  README.md                 # Brief project description
```

---

## Final Words

Congratulations. You have completed the Computational Immunology Bootcamp. You went from
zero programming experience to analyzing real single-cell data, building machine learning
models, and writing a scientific report. These are the same skills used by researchers at
universities and biotech companies around the world. Whether you continue in computational
biology or take another path, you now have a powerful foundation in data science and
biological thinking. Well done.
