# Week 2: Single-Cell RNA-seq Analysis Pipeline

**Time needed:** 8-10 hours
**Goal:** Go from raw count data to a publication-quality UMAP showing different immune cell types in human blood.

Last week you learned about the immune system and loaded your first dataset. This week, you will do what real computational biologists do every day: take messy raw data, clean it up, find structure in it, and identify the immune cell types hiding inside thousands of individual cells. By the end, you will have a labeled map of 2,700 human immune cells -- the exact kind of figure you see in published research papers.

---

## Table of Contents

1. [What You Will Build This Week](#1-what-you-will-build-this-week)
2. [Understanding the Data Format](#2-understanding-the-data-format)
3. [Step 1: Quality Control (QC)](#3-step-1-quality-control-qc) [REQUIRED]
4. [Step 2: Normalization](#4-step-2-normalization) [REQUIRED]
5. [Step 3: Feature Selection (Highly Variable Genes)](#5-step-3-feature-selection-highly-variable-genes) [REQUIRED]
6. [Step 4: Dimensionality Reduction](#6-step-4-dimensionality-reduction) [REQUIRED]
7. [Step 5: Clustering](#7-step-5-clustering) [REQUIRED]
8. [Step 6: Cell Type Annotation](#8-step-6-cell-type-annotation) [REQUIRED]
9. [Step 7: Publication-Quality Figures](#9-step-7-publication-quality-figures) [REQUIRED]
10. [Automated Annotation with CellTypist](#10-automated-annotation-with-celltypist) [CHALLENGE]
11. [Understanding scVI](#11-understanding-scvi) [DEEP DIVE]
12. [Self-Check](#12-self-check)
13. [Key Resources](#13-key-resources)
14. [What's Coming Next](#14-whats-coming-next)

---

## 1. What You Will Build This Week

By the end of this tutorial, you will have produced a colored UMAP plot -- a two-dimensional map where each dot is one immune cell, and the colors tell you which type of immune cell it is. This is one of the most common figures in single-cell biology papers published in journals like Nature, Cell, and Science.

Here is the pipeline you will follow, from start to finish:

```
Raw Data --> QC --> Normalization --> Feature Selection --> Dimensionality Reduction --> Clustering --> Annotation --> Visualization
```

Each step builds on the previous one. Here is what each step does in one sentence:

| Step | What It Does |
|------|--------------|
| **Raw Data** | Load the unprocessed count matrix (how many times each gene was detected in each cell). |
| **QC (Quality Control)** | Remove low-quality cells and genes that would add noise to the analysis. |
| **Normalization** | Adjust the data so cells with different total RNA amounts can be fairly compared. |
| **Feature Selection** | Keep only the ~2,000 most informative genes out of 20,000+ to reduce noise. |
| **Dimensionality Reduction** | Compress the data from ~2,000 dimensions into 2 dimensions you can visualize. |
| **Clustering** | Group similar cells together automatically, like sorting coins by denomination. |
| **Annotation** | Assign biological labels (CD4 T cell, B cell, etc.) to each cluster using known marker genes. |
| **Visualization** | Create publication-quality plots showing your results. |

Do not worry if this looks like a lot. You will work through each step one at a time, and every piece of code is explained line by line.

---

## 2. Understanding the Data Format

Before you write any code, you need to understand how single-cell data is stored. Scanpy (the Python library you will use) stores everything in an object called **AnnData** (short for "Annotated Data"). Think of it as a smart spreadsheet that holds not just the data, but also information about the cells, the genes, and any analysis results.

Here are the key parts of an AnnData object:

| Part | What It Stores | Analogy |
|------|----------------|---------|
| `adata.X` | The count matrix (cells x genes) | The actual spreadsheet of numbers |
| `adata.obs` | Metadata about each cell (one row per cell) | A sidebar with notes about each row |
| `adata.var` | Metadata about each gene (one row per gene) | A header with notes about each column |
| `adata.obsm` | Embeddings like PCA and UMAP coordinates | Extra sheets with compressed versions of the data |
| `adata.uns` | Unstructured info (colors, settings, etc.) | A miscellaneous drawer |

Here is how these pieces fit together visually:

```
              genes (adata.var)
              |
              v
    cells --> [  count matrix  ] = adata.X
    (adata.obs)

    Example:
                   Gene_A  Gene_B  Gene_C  ...
    Cell_001  [     3       0       12     ... ]
    Cell_002  [     0       7        1     ... ]
    Cell_003  [     5       2        0     ... ]
    ...
```

Each number in the matrix tells you: "In this cell, this gene was detected this many times." A zero means the gene was not detected in that cell (which is very common -- most of the matrix is zeros, making it "sparse").

You will see these parts of AnnData referenced throughout the tutorial. Come back to this section if you forget what `adata.obs` or `adata.var` means.

---

## 3. Step 1: Quality Control (QC) [REQUIRED]

### Getting Started

Open your terminal, activate your environment, and start Jupyter:

```bash
conda activate immunology
cd ~/bootcamp
jupyter notebook
```

Create a new notebook and name it `week2_scrna_pipeline.ipynb`. This is where all your code for this week will go.

### Loading the Data

In the first cell of your notebook, type and run the following:

```python
# Import the libraries we need
import scanpy as sc        # Main single-cell analysis library
import numpy as np         # Numerical computing (arrays, math)
import matplotlib.pyplot as plt  # Plotting library

# Set scanpy to show progress messages (verbosity=2 means "informative")
sc.settings.verbosity = 2

# Set figure parameters: resolution and background color
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Load the PBMC 3k dataset (peripheral blood mononuclear cells from a healthy donor)
# This is a famous dataset from 10x Genomics with ~2,700 cells
adata = sc.datasets.pbmc3k()

# Print a summary of the data
print(adata)
```

**Expected output:** You should see something like:

```
AnnData object with n_obs x n_vars = 2700 x 32738
    var: 'gene_ids', 'feature_types'
```

This tells you there are 2,700 cells (n_obs) and 32,738 genes (n_vars). That is a matrix with over 88 million entries -- and you are about to analyze it.

> **Common mistake:** If you see an error about the dataset not being found, make sure you have an internet connection. Scanpy will download it the first time you run this command.

### Why Do We Need Quality Control?

Not every cell in your data is good. Some cells were damaged during the experiment, some "cells" are actually empty droplets that captured stray RNA, and some are doublets (two cells stuck together). If you include these in your analysis, they will distort your results -- like trying to calculate class averages when some answer sheets are blank or have two students' answers mixed together.

We check three main things:

1. **Number of genes detected per cell** -- A real cell should express hundreds of genes. If a "cell" has very few genes detected (say, fewer than 200), it is probably an empty droplet. If it has an unusually high number (say, more than 2,500), it might be a doublet.

2. **Total counts per cell** -- The total number of RNA molecules detected. Very low total counts suggest a dead or empty cell.

3. **Percentage of mitochondrial genes** -- This is the key biological insight. When a cell is dying, its outer membrane breaks open and the contents of the cytoplasm leak out. But mitochondria (the energy factories inside cells) have their own membrane, so their RNA stays behind. A dying cell therefore has a disproportionately high percentage of mitochondrial RNA. In human data, mitochondrial genes start with "MT-".

### Computing QC Metrics

Add a new cell and run:

```python
# Mark mitochondrial genes
# Human mitochondrial gene names start with "MT-" (e.g., MT-CO1, MT-ND1)
adata.var['mt'] = adata.var_names.str.startswith('MT-')

# Calculate QC metrics and add them to adata.obs
# qc_vars=['mt'] tells scanpy to compute mitochondrial percentage
# percent_top=None skips computing the percentage of counts in top N genes
# log1p=False means do not log-transform the metrics
# inplace=True means add the results directly to adata
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=['mt'],
    percent_top=None,
    log1p=False,
    inplace=True
)

# Now adata.obs has three new columns:
# n_genes_by_counts = number of genes detected in each cell
# total_counts = total number of RNA molecules detected in each cell
# pct_counts_mt = percentage of counts that come from mitochondrial genes
print(adata.obs[['n_genes_by_counts', 'total_counts', 'pct_counts_mt']].describe())
```

**Expected output:** A table showing count, mean, std, min, 25%, 50%, 75%, and max for each of the three QC metrics. You should see that most cells have a few hundred to a couple thousand genes detected, total counts in the low thousands, and mitochondrial percentages mostly under 5%.

### Visualizing QC Metrics

Before filtering, you should always look at the distributions. This helps you decide where to set your thresholds. Add a new cell and run:

```python
# Create a figure with three side-by-side histograms
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Plot 1: Distribution of genes detected per cell
axes[0].hist(adata.obs['n_genes_by_counts'], bins=50, color='steelblue', edgecolor='white')
axes[0].set_xlabel('Genes per Cell')
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('Genes Detected per Cell')
# Red dashed lines mark our filtering thresholds
axes[0].axvline(x=200, color='red', linestyle='--', label='Min threshold')
axes[0].axvline(x=2500, color='red', linestyle='--', label='Max threshold')
axes[0].legend()

# Plot 2: Distribution of total RNA counts per cell
axes[1].hist(adata.obs['total_counts'], bins=50, color='coral', edgecolor='white')
axes[1].set_xlabel('Total Counts')
axes[1].set_ylabel('Number of Cells')
axes[1].set_title('Total RNA Counts per Cell')

# Plot 3: Distribution of mitochondrial gene percentage
axes[2].hist(adata.obs['pct_counts_mt'], bins=50, color='mediumseagreen', edgecolor='white')
axes[2].set_xlabel('% Mitochondrial')
axes[2].set_ylabel('Number of Cells')
axes[2].set_title('Mitochondrial Gene %')
# Red dashed line marks our filtering threshold
axes[2].axvline(x=5, color='red', linestyle='--', label='Threshold')
axes[2].legend()

# Adjust spacing so labels do not overlap
plt.tight_layout()

# Save the figure as a PNG file
plt.savefig('qc_plots.png', dpi=150)

# Display the figure in the notebook
plt.show()
```

**Expected output:** Three histograms side by side. The genes-per-cell plot should show a rough bell shape with a peak around 500-1000. The total counts plot should look similar. The mitochondrial percentage plot should show most cells clustered near 0-5%, with a long tail of cells at higher percentages (those are the dying cells).

The red dashed lines show where you will cut. Cells outside these thresholds get removed.

### Applying the Filters

Now apply the filters:

```python
# Record how many cells we start with
print(f"Cells before filtering: {adata.n_obs}")

# Filter 1: Remove cells with fewer than 200 genes detected
# These are likely empty droplets
sc.pp.filter_cells(adata, min_genes=200)

# Filter 2: Remove genes detected in fewer than 3 cells
# These genes are too rare to be informative
sc.pp.filter_genes(adata, min_cells=3)

# Filter 3: Remove cells with more than 5% mitochondrial counts
# These are likely dying cells
adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

# Report how many cells survived
print(f"Cells after filtering: {adata.n_obs}")
```

**Expected output:**

```
Cells before filtering: 2700
Cells after filtering: ~2600 (the exact number may vary slightly)
```

You should lose around 100 cells. That is normal -- those were the low-quality cells that would have added noise to your analysis.

**Why `.copy()`?** When you subset an AnnData object with bracket notation (`adata[...]`), scanpy creates a "view" (a reference to the original data), not an independent copy. Calling `.copy()` ensures you get a fully independent object, which prevents confusing warnings later.

> **[CHECKPOINT 1]** You have loaded the raw PBMC 3k data, computed QC metrics, visualized the distributions, and filtered out low-quality cells. You should have a saved file called `qc_plots.png` and around 2,600 cells remaining. If your cell count after filtering is between 2,400 and 2,700, you are on track.

---

## 4. Step 2: Normalization [REQUIRED]

### Why Normalize?

Imagine you are comparing two classrooms to see which one likes math more. Classroom A has 30 students and Classroom B has 10 students. If you count total hands raised for "I like math," Classroom A will almost always have more -- not because they like math more, but because there are more students. To compare fairly, you would use percentages instead: "What fraction of students in each classroom raised their hand?"

The same problem exists in single-cell data. Different cells capture different amounts of total RNA during the experiment. A cell with 10,000 total counts and 100 counts of gene X is not necessarily expressing more of gene X than a cell with 2,000 total counts and 50 counts of gene X. In fact, the second cell has gene X as a higher fraction of its total (50/2000 = 2.5% vs 100/10000 = 1%). Normalization fixes this by putting every cell on the same scale.

### Performing Normalization

```python
# Normalize each cell so that its total counts sum to 10,000
# This is like converting raw numbers to "per 10,000" (similar to percentages)
# After this step, gene expression values are comparable across cells
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform the data: replace each value x with log(x + 1)
# Why? Gene expression data is very skewed -- a few genes have enormous counts
# while most have small counts. Log-transform compresses the high values and
# spreads out the low values, making the data easier to analyze
# The "+1" (that is the "1p" in log1p) prevents log(0), which is undefined
sc.pp.log1p(adata)

# Save the current state as "raw" for later use
# Some plotting functions (like dotplot) need access to the normalized-but-not-subset data
# Storing it in adata.raw means we can always access the full gene set later
adata.raw = adata

print("Normalization complete.")
print(f"Data shape: {adata.shape}")
```

**Expected output:**

```
Normalization complete.
Data shape: (~2600, ~13000)
```

The number of genes is now lower than the original 32,738 because `filter_genes(min_cells=3)` in the QC step removed very rare genes.

> **[CHECKPOINT 2]** Your data is now normalized and log-transformed. The values in `adata.X` are no longer raw counts -- they are normalized log-counts that are comparable across cells.

---

## 5. Step 3: Feature Selection (Highly Variable Genes) [REQUIRED]

### Why Select Features?

Your data currently has around 13,000 genes. But most of those genes are expressed at roughly the same level in every cell -- they are "housekeeping" genes that all cells need to survive. These genes tell you nothing about what makes a T cell different from a B cell. They just add noise.

The genes that matter are the ones that **vary** between cells. A gene that is high in T cells but low in B cells is informative. A gene that is medium in every cell is not. By selecting the ~2,000 most variable genes (called **highly variable genes** or HVGs), you keep the signal and throw away the noise.

### Finding Highly Variable Genes

```python
# Identify highly variable genes
# min_mean=0.0125: ignore genes with very low average expression (too rare)
# max_mean=3: ignore genes with very high average expression (housekeeping)
# min_disp=0.5: keep genes with high dispersion (lots of variation between cells)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# How many genes were selected?
print(f"Number of highly variable genes: {adata.var.highly_variable.sum()}")

# Visualize the selection
# Black dots = selected (highly variable), grey dots = not selected
sc.pl.highly_variable_genes(adata)
```

**Expected output:** You should see that around 1,500-2,000 genes are marked as highly variable. The plot will show a scatter of genes, with the selected ones highlighted. The selected genes have both reasonable mean expression and high dispersion (variation).

```python
# Keep only the highly variable genes for downstream analysis
# We are reducing from ~13,000 genes to ~1,500-2,000
adata = adata[:, adata.var.highly_variable].copy()

print(f"Data shape after filtering: {adata.shape}")
```

**Expected output:**

```
Data shape after filtering: (~2600, ~1800)
```

Your matrix just got much smaller and more informative.

> **[CHECKPOINT 3]** You have selected the most informative genes. Your data now has around 1,500-2,000 genes instead of 13,000+. The plot of highly variable genes should show a clear separation between selected and non-selected genes.

---

## 6. Step 4: Dimensionality Reduction [REQUIRED]

### What Is Dimensionality Reduction?

Even after selecting ~2,000 genes, you cannot visualize 2,000 dimensions on a screen. You need to compress the data into fewer dimensions while preserving the important structure. This happens in two stages:

**Stage 1: PCA (Principal Component Analysis)**

Imagine you have data about students: height, weight, arm length, shoe size, and leg length. These five measurements are not all independent -- tall students tend to also have long arms, big shoes, and long legs. PCA figures out that most of the variation can be explained by just two "directions": something like "overall body size" and "body proportions." These compressed directions are called **principal components (PCs)**.

In our case, PCA compresses ~2,000 genes into ~50 principal components. Each PC captures a pattern of genes that vary together across cells.

**Stage 2: UMAP (Uniform Manifold Approximation and Projection)**

Fifty dimensions is still too many to visualize. UMAP takes those 50 PCs and squishes them down to 2 dimensions that you can plot on a screen. The key property of UMAP is that it keeps similar cells close together: if two cells have similar gene expression in 50D space, they will be near each other in the 2D UMAP.

Think of UMAP like making a map. Tokyo and Osaka are close on a map because they are close in real life. New York and Tokyo are far apart on a map because they are far apart in real life. UMAP does the same thing, but the "distance" between cells is measured by gene expression similarity, not geography.

> **Important caveat:** Distances within a UMAP cluster are meaningful (nearby points are similar), but distances between clusters can be misleading. Two clusters that look far apart on UMAP might actually be quite similar. UMAP preserves local structure well but can distort global structure.

### Running PCA

```python
# Scale the data: for each gene, subtract the mean and divide by the standard deviation
# This ensures that every gene has mean=0 and variance=1
# Without scaling, genes with naturally higher expression would dominate the PCA
# max_value=10 clips extreme values to prevent outliers from having too much influence
sc.pp.scale(adata, max_value=10)

# Run PCA: compress ~2,000 genes into 50 principal components
# svd_solver='arpack' is an efficient algorithm for computing PCA
# n_comps=50 means compute 50 principal components
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# Plot the variance explained by each PC (called a "scree plot" or "elbow plot")
# This tells you how much information each PC captures
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
```

**Expected output:** A plot showing variance explained by each PC. The first few PCs explain the most variance, and the curve drops off sharply (forming an "elbow"). After about 30-40 PCs, additional PCs explain very little. This tells you that 40 PCs capture most of the meaningful variation -- the rest is noise.

**How to read the elbow plot:** Look for where the curve flattens out. The PCs before the "elbow" carry real biological signal. The PCs after the elbow are mostly noise. You do not need to find the exact elbow -- using 30-40 PCs is a safe default for most datasets.

### Building the Neighborhood Graph and Running UMAP

```python
# Build a neighborhood graph
# This computes which cells are most similar to each other based on their PCA coordinates
# n_neighbors=10: for each cell, find the 10 most similar cells
# n_pcs=40: use the first 40 principal components for distance calculations
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Run UMAP: compress the neighborhood graph into 2 dimensions for visualization
sc.tl.umap(adata)

# Plot the UMAP -- no colors yet, just showing the structure
sc.pl.umap(adata, color=None, title='UMAP of PBMC 3k')
```

**Expected output:** A scatter plot where each dot is one cell. You should see several distinct clusters (groups of dots close together) separated by gaps. These clusters represent different cell types, but we do not know which is which yet -- that comes in the next steps.

If you see one big blob with no structure, something went wrong in an earlier step. Go back and check your QC filtering and normalization.

> **Common mistake:** If your UMAP looks very different each time you run it, that is normal. UMAP has randomness built in, so the exact layout changes, but the clusters should always be the same groups of cells. If you want reproducible results, add `random_state=42` to the `sc.tl.umap()` call.

> **[CHECKPOINT 4]** You should see a UMAP with several distinct clusters of points. The cells have organized themselves into groups based on gene expression similarity. You have gone from ~2,000 dimensions to 2 dimensions while preserving the structure. This is one of the most satisfying moments in single-cell analysis.

---

## 7. Step 5: Clustering [REQUIRED]

### What Is Clustering?

Clustering means grouping similar cells together automatically. Imagine you have a pile of coins from different countries: US quarters, Japanese yen, British pounds, and euros. Without knowing anything about the coins, you could still sort them into piles based on size, color, and weight. That is what clustering does with cells -- it groups them based on gene expression patterns, without knowing what the cell types are.

### The Leiden Algorithm

Scanpy uses the **Leiden algorithm** for clustering. Here is the intuition:

1. You already built a neighborhood graph (in the previous step) that connects each cell to its most similar neighbors.
2. The Leiden algorithm finds "communities" in this graph -- groups of cells that are densely connected to each other but sparsely connected to cells in other groups.
3. It is like finding friend groups in a social network: people within a friend group know each other well, but they have fewer connections to people in other groups.

### The Resolution Parameter

The `resolution` parameter controls how many clusters you get:

- **Higher resolution (e.g., 1.0 or 2.0):** More clusters, finer-grained grouping. Like sorting coins by country AND denomination.
- **Lower resolution (e.g., 0.3 or 0.5):** Fewer clusters, coarser grouping. Like sorting coins by just metal type (gold, silver, copper).

There is no single "correct" resolution. It depends on what you want to learn. For the PBMC 3k dataset, `resolution=0.5` gives a reasonable number of clusters that correspond to major immune cell types.

### Running the Clustering

```python
# Run Leiden clustering
# resolution=0.5: moderate resolution, good starting point for PBMC data
# flavor='igraph': use the igraph implementation (faster for small datasets)
# n_iterations=2: number of optimization iterations (2 is usually enough)
sc.tl.leiden(
    adata,
    resolution=0.5,
    flavor='igraph',
    n_iterations=2
)

# Plot the UMAP, now colored by cluster assignment
# Each cluster gets a different color and a number label
sc.pl.umap(adata, color='leiden', title='Clusters (Leiden)', legend_loc='on data')

# How many clusters did we get?
print(f"Number of clusters: {adata.obs['leiden'].nunique()}")

# How many cells are in each cluster?
print(adata.obs['leiden'].value_counts())
```

**Expected output:** A UMAP plot with each cluster colored differently and labeled with a number (0, 1, 2, ...). You should see roughly 7-9 clusters, depending on random variation. Each cluster corresponds to a group of cells with similar gene expression. The value counts will show you how many cells are in each cluster -- some clusters will be large (hundreds of cells) and others small (tens of cells).

> **Common mistake:** If you get an error about `leidenalg` not being installed, run `pip install leidenalg` in your terminal (with the immunology environment activated) and restart your notebook kernel.

> **[CHECKPOINT 5]** Your UMAP should now show colored, numbered clusters. Each cluster is a group of cells that the algorithm thinks are similar to each other. You do not know what these clusters ARE yet (are they T cells? B cells?), but you can see that the algorithm has found structure in the data.

---

## 8. Step 6: Cell Type Annotation [REQUIRED]

### From Clusters to Cell Types

Clustering tells you WHICH cells group together, but not WHAT they are. That is like sorting your coin pile into groups but not knowing which group is quarters and which is dimes. To identify cell types, you need to look at **marker genes** -- specific genes that are known to be expressed by specific cell types.

This is where your immunology knowledge from Week 1 comes in. Each immune cell type has a characteristic set of genes that it expresses. By checking which marker genes are active in each cluster, you can figure out what cell type that cluster represents.

### Marker Gene Reference Table

Here are the marker genes for the major immune cell types you will find in PBMC (peripheral blood mononuclear cells). Save this table -- you will refer to it repeatedly.

| Cell Type | Key Markers | What to Look For |
|-----------|------------|-----------------|
| CD4 T cells | IL7R, CD4 | IL7R is the most reliable; CD4 is also expressed by monocytes at low levels |
| CD8 T cells | CD8A, CD8B | Both should be high; distinguish from NK cells by lack of GNLY |
| B cells | MS4A1 (CD20), CD79A | MS4A1 is the classic B cell marker |
| NK cells | GNLY, NKG7 | Natural killer cells; high cytotoxic gene expression |
| CD14+ Monocytes | CD14, LYZ | The most common monocyte subtype in blood |
| FCGR3A+ Monocytes | FCGR3A, MS4A7 | Also called CD16+ monocytes; a less common subtype |
| Dendritic cells | FCER1A, CST3 | Antigen-presenting cells; relatively rare in blood |
| Platelets | PPBP | Megakaryocyte-derived cell fragments; very small cluster |

### Visualizing Marker Gene Expression

The **dot plot** is one of the most informative ways to check which markers are active in which cluster. The size of each dot shows what fraction of cells in that cluster express the gene, and the color shows the average expression level.

```python
# Define the marker genes we want to check
marker_genes = [
    'IL7R',    # CD4 T cells
    'CD8A',    # CD8 T cells
    'MS4A1',   # B cells (also known as CD20)
    'CD79A',   # B cells
    'GNLY',    # NK cells
    'NKG7',    # NK cells
    'CD14',    # CD14+ Monocytes
    'LYZ',     # CD14+ Monocytes (lysozyme)
    'FCGR3A',  # FCGR3A+ Monocytes (also known as CD16)
    'MS4A7',   # FCGR3A+ Monocytes
    'FCER1A',  # Dendritic cells
    'CST3',    # Dendritic cells
    'PPBP'     # Platelets
]

# Dot plot: rows are clusters, columns are marker genes
# Dot size = fraction of cells expressing the gene
# Dot color = average expression level
sc.pl.dotplot(adata, marker_genes, groupby='leiden', standard_scale='var')
```

**Expected output:** A grid of dots. For each cluster (row), you should see large, dark dots for the marker genes of that cluster's cell type, and small or absent dots for other markers. For example, the cluster that is B cells should have large dark dots at MS4A1 and CD79A, but small/light dots at CD14 and CD8A.

Now look at the markers directly on the UMAP to see where each marker is expressed spatially:

```python
# Feature plots: overlay gene expression on the UMAP
# Yellow/bright = high expression, dark/purple = low expression
sc.pl.umap(
    adata,
    color=['CD14', 'MS4A1', 'CD8A', 'IL7R', 'GNLY', 'FCER1A'],
    ncols=3,       # Arrange plots in 3 columns
    cmap='viridis'  # Color map: purple (low) to yellow (high)
)
```

**Expected output:** Six UMAP plots, each colored by the expression of one gene. You should see that each gene "lights up" a different part of the UMAP:
- CD14 should light up one or two clusters (monocytes)
- MS4A1 should light up a different cluster (B cells)
- CD8A should light up another cluster (CD8 T cells)
- IL7R should light up the largest cluster(s) (CD4 T cells)
- GNLY should light up a cluster near the T cells (NK cells)
- FCER1A should light up a small cluster (dendritic cells)

### Assigning Cell Type Labels

Now comes the part where YOU make a scientific judgment. Look at your dot plot and feature plots, and decide which cluster corresponds to which cell type. Here is the process:

1. Look at each cluster number in the dot plot
2. See which marker genes have the biggest dots in that cluster
3. Match those markers to the table above
4. Assign the cell type label

```python
# Based on marker gene expression, assign cell type labels to each cluster
# IMPORTANT: Look at YOUR dot plot and feature plots to decide!
# The cluster numbers below are an example -- your numbers may differ
# because clustering has some randomness

cluster_to_celltype = {
    '0': 'CD4 T cells',
    '1': 'CD14 Monocytes',
    '2': 'CD4 T cells',
    '3': 'B cells',
    '4': 'CD8 T cells',
    '5': 'NK cells',
    '6': 'FCGR3A Monocytes',
    '7': 'Dendritic cells',
}

# TIP: Your cluster numbers will very likely not match these exactly!
# You MUST look at your own dot plot and adjust the mapping above.
# For example, if YOUR cluster 0 has high MS4A1 and CD79A,
# then change '0' to 'B cells' in the dictionary above.
#
# If you have more clusters than listed here (e.g., cluster 8),
# check its markers and add it to the dictionary.
# If a cluster does not match any known type well,
# you can label it 'Unknown' for now.

# Apply the mapping: create a new column 'cell_type' in adata.obs
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)

# Check if any cells were not assigned (NaN means a cluster was missing from your dictionary)
unassigned = adata.obs['cell_type'].isna().sum()
if unassigned > 0:
    print(f"WARNING: {unassigned} cells have no assigned cell type.")
    print("Check if all cluster numbers are in your cluster_to_celltype dictionary.")
    # Fill any missing with 'Unknown'
    adata.obs['cell_type'] = adata.obs['cell_type'].fillna('Unknown')

# Plot the final annotated UMAP
sc.pl.umap(
    adata,
    color='cell_type',
    title='Cell Types in PBMC 3k',
    legend_loc='on data',      # Put labels directly on the clusters
    legend_fontsize=8,         # Keep the labels readable but not too large
    frameon=False              # Remove the box around the plot for a cleaner look
)

# Save this figure -- it is your main result for this week
plt.savefig('pbmc3k_celltype_umap.png', dpi=150, bbox_inches='tight')
plt.show()
```

**Expected output:** A beautiful UMAP plot where each cluster is labeled with an immune cell type name. You should see:
- One or two large clusters of CD4 T cells (the most common immune cell in blood)
- A cluster of CD14 Monocytes
- A cluster of B cells
- A cluster of CD8 T cells
- A cluster of NK cells
- Smaller clusters for FCGR3A Monocytes, Dendritic cells, and possibly Platelets

If a cluster does not clearly match any cell type, it is okay to label it "Unknown" and come back to it later. Real researchers face this situation regularly.

> **[CHECKPOINT 6]** You just identified immune cell types from real human blood data. You took 2,700 unlabeled cells, filtered them, normalized them, found structure in the data, clustered them, and identified what kind of immune cell each one is -- using nothing but gene expression patterns and your knowledge of marker genes. This is exactly what researchers do in published papers. Your saved file `pbmc3k_celltype_umap.png` is a publication-style figure.

---

## 9. Step 7: Publication-Quality Figures [REQUIRED]

A single UMAP is a good start, but research papers typically include multiple complementary visualizations. In this section, you will create violin plots and a heatmap to strengthen your results.

### Violin Plots

Violin plots show the distribution of gene expression for each cell type. They are useful for showing that your marker genes are genuinely specific to the expected cell types.

```python
# Violin plots of key marker genes, grouped by cell type
# Each "violin" shows the distribution of expression values for that gene in that cell type
sc.pl.violin(
    adata,
    ['CD14', 'LYZ', 'MS4A1', 'CD79A', 'CD8A', 'IL7R', 'GNLY'],
    groupby='cell_type',
    rotation=45  # Rotate x-axis labels so they do not overlap
)
```

**Expected output:** A panel of violin plots. For each gene, you should see that one or two cell types have high expression (a tall, wide violin shape) while the others have low expression (a flat line near zero). For example, CD14 and LYZ should be high in CD14 Monocytes and low in everything else. MS4A1 and CD79A should be high in B cells only.

### Heatmap of Top Marker Genes

The heatmap shows the top differentially expressed genes for each cell type. Scanpy can automatically find these genes using a statistical test.

```python
# Find the top marker genes for each cell type using the Wilcoxon rank-sum test
# This test asks: "Is this gene significantly more expressed in this cell type
# compared to all other cell types?"
sc.tl.rank_genes_groups(adata, 'cell_type', method='wilcoxon')

# Plot a heatmap of the top 5 genes per cell type
# Each row is a gene, each column is a cell (grouped by cell type)
# Yellow = high expression, purple = low expression
sc.pl.rank_genes_groups_heatmap(
    adata,
    n_genes=5,                # Show the top 5 genes per cell type
    groupby='cell_type',
    show_gene_labels=True     # Display gene names on the plot
)
```

**Expected output:** A heatmap with clear blocks of yellow along the diagonal. Each cell type should have a set of 5 genes that are bright yellow (highly expressed) in that cell type but dark purple (low expression) in all others. This confirms that your cell type assignments are biologically meaningful.

> **Tip:** If the gene labels are hard to read, try adding `figsize=(12, 8)` or adjusting `show_gene_labels=True` in the plot call. You can also view the top genes as a table:

```python
# Print the top 5 marker genes for each cell type (as a table)
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
top_genes = {group: [result['names'][group][i] for i in range(5)] for group in groups}

for cell_type, genes in top_genes.items():
    print(f"{cell_type}: {', '.join(genes)}")
```

**Expected output:** A list showing the top 5 marker genes for each cell type. You should recognize some of them from the marker gene table earlier. You may also discover new markers you did not know about -- this is how researchers identify novel marker genes.

> **[CHECKPOINT 7]** You now have multiple visualizations supporting your cell type annotations: a labeled UMAP, violin plots, and a heatmap. These are the core figures you would include in a research paper's supplementary materials or main figures.

---

## 10. Automated Annotation with CellTypist [CHALLENGE]

Manual annotation (looking at marker genes and deciding cell types yourself) is the gold standard, but it requires expertise and is time-consuming. **CellTypist** is a tool that automates this process using a pre-trained machine learning model. It has been trained on millions of cells with known labels and can predict cell types automatically.

This section is optional. If you are short on time, skip to the Self-Check.

### Running CellTypist

```python
import celltypist
from celltypist import models

# Download a pre-trained model
# 'Immune_All_Low.pkl' is trained on immune cells at a low (broad) resolution
# It can identify major immune cell types like T cells, B cells, monocytes, etc.
models.download_models(model='Immune_All_Low.pkl')

# Run cell type prediction
# majority_voting=True uses a smoothing step: if most cells in a neighborhood
# are predicted as one type, isolated differently-labeled cells get corrected
# This reduces noisy predictions
predictions = celltypist.annotate(
    adata,
    model='Immune_All_Low.pkl',
    majority_voting=True
)

# Convert predictions to an AnnData object with the labels added
adata_ct = predictions.to_adata()

# Plot the CellTypist predictions on the UMAP
sc.pl.umap(
    adata_ct,
    color='majority_voting',
    title='CellTypist Predictions',
    legend_fontsize=6,
    frameon=False
)
```

**Expected output:** A UMAP colored by CellTypist's predicted cell types. The labels will be more granular than your manual annotations -- for example, CellTypist might distinguish "Classical monocytes" from "Non-classical monocytes" or "Naive CD4 T cells" from "Memory CD4 T cells."

### Comparing with Your Manual Annotation

Take a moment to compare the CellTypist predictions with your manual annotations. Do they agree on the major cell types? Where do they differ? CellTypist might identify subtypes that you grouped together (e.g., it might split your "CD4 T cells" into naive and memory subsets).

Neither approach is "right" or "wrong" -- they serve different purposes. Manual annotation gives you full control and understanding, while automated annotation is faster and can catch subtypes you might miss.

> **[CHECKPOINT 8 -- CHALLENGE]** You have run an automated cell type prediction tool and compared its results with your manual analysis. If you completed this section, you now understand both manual and automated approaches to cell type annotation.

---

## 11. Understanding scVI [DEEP DIVE]

This section is purely informational. You do not need to run any code.

**scVI** (single-cell Variational Inference) is a deep learning framework for single-cell data analysis, developed by the Yosef Lab at UC Berkeley. It uses a type of neural network called a **variational autoencoder** (VAE) to learn a compressed representation of single-cell data.

### Why Would You Use scVI Instead of the Standard Pipeline?

The pipeline you ran today (normalize, select HVGs, PCA, UMAP) works great for many datasets. But it has limitations:

- **Batch effects:** If you combine data from multiple experiments, technical differences between experiments can overshadow biological differences. scVI can correct for these batch effects during the analysis.
- **Count data modeling:** The standard pipeline log-transforms the data, which is a rough approximation. scVI directly models the count data using a negative binomial distribution, which is a better statistical fit.
- **Uncertainty:** scVI can tell you how confident it is in its results, which the standard pipeline cannot.
- **Scalability:** scVI handles millions of cells efficiently because neural networks are good at processing large datasets in batches.

### When Will You Use scVI?

In later weeks and in real research projects, you may encounter datasets with multiple batches or very large cell counts. That is when scVI becomes essential. For now, the standard pipeline you learned today is perfectly appropriate.

If you are curious, the scvi-tools documentation is excellent: [https://scvi-tools.org/](https://scvi-tools.org/)

---

## 12. Self-Check

### Review Questions

Answer these questions to make sure you understand what you did. Try to answer from memory before looking back at the tutorial.

1. **What does QC stand for, and why do we do it?**

   QC stands for Quality Control. We do it to remove low-quality cells (empty droplets, dying cells, doublets) that would add noise and distort our analysis.

2. **Why do we filter out cells with high mitochondrial gene percentage?**

   Dying cells have damaged outer membranes, so their cytoplasmic RNA leaks out. Mitochondrial RNA stays behind because mitochondria have their own membrane. A high mitochondrial percentage is therefore a signature of cell death.

3. **What is the purpose of normalization?**

   Different cells capture different amounts of total RNA during the experiment. Normalization puts all cells on the same scale so that gene expression levels are comparable across cells.

4. **What does UMAP show you?**

   UMAP is a 2D visualization where each dot is one cell. Cells that are close together on the UMAP have similar gene expression profiles. It lets you see the structure in your data at a glance.

5. **How many cell types did you identify?**

   You should have identified 7-8 major cell types: CD4 T cells, CD8 T cells, B cells, NK cells, CD14 Monocytes, FCGR3A Monocytes, Dendritic cells, and possibly Platelets.

6. **Can you name 3 marker genes and the cell types they mark?**

   Examples: IL7R marks CD4 T cells, MS4A1 (CD20) marks B cells, CD14 marks CD14+ Monocytes, GNLY marks NK cells, CD8A marks CD8 T cells.

### Output Verification

Check that you have produced these files:

- [ ] `qc_plots.png` -- Three histograms showing the distributions of genes per cell, total counts, and mitochondrial percentage
- [ ] `pbmc3k_celltype_umap.png` -- A UMAP colored by cell type labels

### Achievement

You just analyzed 2,700 human immune cells and identified their types -- this is exactly what researchers do in real published papers. The pipeline you learned today (QC, normalization, feature selection, dimensionality reduction, clustering, annotation) is the foundation of virtually every single-cell RNA-seq study. Whether the paper is in Nature, Cell, or Science, it almost certainly includes this pipeline or a close variant.

You should feel confident about what you accomplished. This is real computational biology.

---

## 13. Key Resources

These resources can help you go deeper or troubleshoot issues:

| Resource | What It Covers |
|----------|---------------|
| [Single-cell best practices (Chapters 1-21)](https://www.sc-best-practices.org/) | Comprehensive guide to scRNA-seq analysis with best practices and pitfalls |
| [Scanpy tutorials](https://scanpy-tutorials.readthedocs.io/) | Official Scanpy tutorials, including the PBMC 3k walkthrough |
| [CellTypist](https://www.celltypist.org/) | Automated cell type annotation tool and documentation |
| [Kaggle](https://www.kaggle.com/) | Practice datasets and competitions for data science skills |
| [scikit-learn documentation](https://scikit-learn.org/) | Machine learning library documentation (useful for understanding PCA, clustering concepts) |
| [scvi-tools](https://scvi-tools.org/) | Deep learning framework for single-cell analysis |

---

## 14. What's Coming Next

In **Week 3: TCR and BCR Analysis**, you will explore the other side of single-cell immune data: the receptors. Every T cell and B cell carries a unique receptor (TCR or BCR) that determines which pathogens it can recognize. By combining gene expression data (what you learned this week) with receptor sequence data, you can ask powerful questions like:

- Which T cells are expanding in response to an infection? (Clonal expansion)
- Do T cells that recognize the same antigen also have similar gene expression?
- How diverse is the immune repertoire in a patient?

You will use **scirpy** (a companion library to scanpy, specifically for immune receptor data) to analyze TCR and BCR sequences from single-cell experiments.

Everything you learned this week -- AnnData objects, UMAP, clustering, cell type annotation -- will carry over directly. Week 3 builds on top of this foundation.

---

*You have finished Week 2. Take a break, then come back for Week 3 when you are ready.*
