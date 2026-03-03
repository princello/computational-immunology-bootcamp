# Week 4: Multi-omics Integration

**Time needed:** 8-10 hours
**Goal:** Learn how different measurement technologies capture different aspects of a cell, and integrate multiple data types to get a more complete picture of immune biology.

By the end of this week, you will understand the multi-omics landscape, work with CITE-seq data (RNA + protein), explore chromatin accessibility concepts, and integrate multiple modalities using modern computational tools. This is one of the most in-demand skills in computational biology today.

---

## Table of Contents

1. [What You Will Learn This Week](#1-what-you-will-learn-this-week)
2. [The Multi-omics Landscape](#2-the-multi-omics-landscape)
3. [CITE-seq: Gene Expression + Protein in the Same Cell](#3-cite-seq-gene-expression--protein-in-the-same-cell)
4. [Hands-On: Analyzing a CITE-seq Dataset](#4-hands-on-analyzing-a-cite-seq-dataset) [REQUIRED]
5. [Understanding scATAC-seq (Conceptual)](#5-understanding-scatac-seq-conceptual) [REQUIRED]
6. [Multi-omics Integration Concepts](#6-multi-omics-integration-concepts) [REQUIRED]
7. [Hands-On: Working with MuData (Multi-modal Data)](#7-hands-on-working-with-mudata-multi-modal-data) [REQUIRED]
8. [Hands-On: Integration with totalVI](#8-hands-on-integration-with-totalvi) [REQUIRED]
9. [Comparing Modalities](#9-comparing-modalities) [REQUIRED]
10. [Regulatory Inference from ATAC Data](#10-regulatory-inference-from-atac-data) [CHALLENGE]
11. [The Chopp et al. Paper](#11-the-chopp-et-al-paper) [DEEP DIVE]
12. [Key Papers](#12-key-papers) [DEEP DIVE]
13. [Self-Check](#13-self-check)
14. [What's Coming Next](#14-whats-coming-next)

---

## 1. What You Will Learn This Week

Over the past three weeks, you have built a solid foundation: you can set up a computing environment, you understand immune cell types, you have analyzed single-cell RNA-seq data from scratch, and you have explored TCR/BCR repertoires. All of that work used **one type of measurement at a time**.

This week, we go further. In real biology, a single cell is doing many things simultaneously -- expressing genes, displaying proteins on its surface, opening and closing regions of its DNA, and carrying a unique immune receptor. No single technology captures the full picture. **Multi-omics** means measuring multiple of these layers at once and combining them computationally.

### The Analogy

Imagine you are trying to understand a person. If you could only read their emails, you would know what they are working on right now (that is like **RNA** -- the cell's current to-do list). But you would miss a lot:

- Their **social media profile** shows how they present themselves to the world (like **surface proteins** -- what the cell actually displays)
- Their **diary** reveals what they *could* think about, even if they are not thinking about it right now (like **chromatin accessibility** -- which genes are bookmarked as "available to use")
- Their **contacts list** shows who they interact with (like **immune receptor sequences** -- which antigens the cell can recognize)

To truly understand someone, you want all of these. The same is true for cells.

### By the End of This Week

You will be able to:

- Explain what different omics technologies measure and why combining them matters
- Load and process CITE-seq data (RNA + protein)
- Create and manipulate MuData objects for multi-modal data
- Integrate RNA and protein data using totalVI
- Compare gene expression and protein levels for the same markers
- Describe how chromatin accessibility (scATAC-seq) reveals regulatory potential
- Discuss integration strategies for different experimental designs

---

## 2. The Multi-omics Landscape

Before we write any code, let's understand the different technologies and what they each bring to the table.

### Technologies at a Glance

| Technology | What It Measures | Analogy | Typical # of Features |
|-----------|-----------------|---------|----------------------|
| scRNA-seq | Which genes are actively being transcribed | Reading the cell's to-do list | ~20,000 genes |
| CITE-seq | Surface protein levels (alongside RNA) | Seeing the cell's uniform and badges | ~100-300 proteins |
| scATAC-seq | Which DNA regions are physically accessible | Seeing which pages of the instruction manual are bookmarked | ~100,000-500,000 peaks |
| scTCR/BCR-seq | Immune receptor sequences | Reading the cell's ID card | 2-4 chains per cell |
| 10x Multiome | RNA + ATAC together from the same cell | Getting the to-do list AND the bookmarks from the same cell | ~20,000 genes + ~100,000 peaks |

### Why Integration Matters

Each data type answers a different question:

- **RNA** tells you what the cell is doing RIGHT NOW -- which genes it is actively transcribing into messenger RNA
- **Protein** confirms whether those instructions actually made it to the cell surface. This is the "ground truth" for surface markers because proteins are what other cells and the immune system actually see
- **ATAC** (chromatin accessibility) reveals what the cell COULD do. Open chromatin regions are like bookmarked pages -- the cell has made those genes available for rapid activation, even if they are not being used yet
- **Receptor sequences** (TCR/BCR) link each cell to its antigen specificity -- what that particular cell is designed to recognize

When you combine these layers, you get insights that no single measurement could provide. For example, you might find that a group of T cells all have similar RNA profiles, but their chromatin patterns reveal that some are "primed" to become one cell type while others are headed in a different direction. That kind of insight is invisible to RNA alone.

### A Concrete Example

Consider a developing T cell in the thymus that is deciding whether to become a CD4+ helper T cell or a CD8+ killer T cell:

- **RNA-seq alone** might show the cell expressing both CD4 and CD8 genes (it hasn't decided yet)
- **Protein measurement** might show CD4 protein on the surface but CD8 protein fading (the decision is being made)
- **ATAC-seq** might show that the chromatin around CD4-related genes is opening up while CD8-related regions are closing (the cell has committed, even if the RNA hasn't caught up)

This is exactly what researchers found in the Chopp et al. paper we will discuss later -- chromatin changes predicted cell fate before gene expression did.

---

## 3. CITE-seq: Gene Expression + Protein in the Same Cell

### What Is CITE-seq?

**CITE-seq** stands for **C**ellular **I**ndexing of **T**ranscriptomes and **E**pitopes by **seq**uencing. It was developed by Stoeckius et al. (Nature Methods, 2017) and is one of the most widely used multi-modal single-cell technologies.

### How It Works

1. Antibodies that bind to specific surface proteins are tagged with unique DNA barcodes (these are called **Antibody-Derived Tags**, or **ADTs**)
2. These tagged antibodies are mixed with cells
3. Each antibody binds to its target protein on the cell surface
4. The cells go through the standard 10x Genomics (or similar) droplet-based capture
5. Inside each droplet, both the cell's mRNA AND the antibody DNA barcodes are captured
6. After sequencing, you get two readouts per cell: gene expression counts AND protein abundance counts

### Why Not Just Use RNA?

This is an important question. If RNA eventually gets translated into protein, why bother measuring protein separately?

Several reasons:

1. **Timing mismatch.** A cell might have already stopped making CD4 mRNA but still have CD4 protein on its surface (proteins last longer than mRNA)
2. **Translation is not guaranteed.** Not all mRNA molecules get translated into protein. Post-transcriptional regulation can block translation
3. **Surface display.** Some proteins are made but stored inside the cell. Only surface proteins are visible to other cells and the immune system
4. **Sensitivity.** For many surface markers (like CD4, CD8, CD3), the protein measurement is cleaner and less noisy than the RNA measurement
5. **Ground truth for cell typing.** Immunologists have spent decades defining cell types by surface protein expression using flow cytometry. CITE-seq lets you use those same protein markers in a single-cell sequencing experiment

### Key Differences Between RNA and Protein Data

| Aspect | RNA (Gene Expression) | Protein (ADT) |
|--------|----------------------|----------------|
| Number of features | ~20,000 genes | ~100-300 proteins |
| Noise level | High (dropout, low counts) | Lower for measured proteins |
| Coverage | Genome-wide | Only proteins with antibodies |
| Normalization | Log-normalization or scran | CLR (Centered Log-Ratio) |
| What it tells you | What the cell is transcribing | What the cell is displaying |

---

## 4. Hands-On: Analyzing a CITE-seq Dataset [REQUIRED]

Time to write code. Open a new Jupyter notebook and name it `week4_multiomics.ipynb`.

### 4.1 Install Required Packages

Before starting, make sure you have the necessary packages. Run this in your terminal (not in the notebook):

```bash
# Make sure your environment is activated
conda activate immunology

# Install muon for multi-modal data handling
pip install muon

# scvi-tools should already be installed from Week 0
# If not, run: pip install scvi-tools
```

### 4.2 Set Up Your Notebook

```python
# Cell 1: Imports and settings
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=80, facecolor='white')
sc.settings.verbosity = 1  # Show some progress messages

print("Scanpy version:", sc.__version__)
print("AnnData version:", ad.__version__)
print("Setup complete.")
```

**Expected output:**
```
Scanpy version: 1.x.x
AnnData version: 0.x.x
Setup complete.
```

### 4.3 Load a CITE-seq Dataset

We will use the `scvi-tools` built-in PBMC CITE-seq dataset. This is a dataset of peripheral blood mononuclear cells (PBMCs) with both RNA and protein measurements.

> **Note on data loading:** Multi-omics datasets can be large and servers sometimes go down. We provide multiple fallback options below. If one approach does not work, do not worry -- try the next one.

```python
# Cell 2: Load the CITE-seq dataset
# Option A: Use scvi-tools built-in dataset (recommended)
import scvi

# This downloads a PBMC dataset with both RNA and protein (ADT) measurements
# It may take a minute to download the first time
adata = scvi.data.pbmc_seurat_v4_cite_seq(
    apply_filters=True,    # Pre-filtered for quality
    normalize=False         # We will handle normalization ourselves
)

print("Full dataset:")
print(adata)
print()
print("Shape:", adata.shape)
print("This means:", adata.shape[0], "cells and", adata.shape[1], "genes")
```

**Expected output (approximately):**
```
Full dataset:
AnnData object with n_obs x n_vars = XXXX x XXXX
    obs: 'celltype', ...
    obsm: 'protein_expression'
    ...
Shape: (XXXX, XXXX)
This means: XXXX cells and XXXX genes
```

> **If Option A does not work**, try this alternative:

```python
# Cell 2 (Alternative): If scvi.data.pbmc_seurat_v4_cite_seq() fails,
# try loading the dataset manually

# Option B: Download from a direct URL
# The scvi-tools team hosts datasets that are usually reliable
try:
    adata = scvi.data.pbmc_seurat_v4_cite_seq(apply_filters=True, normalize=False)
    print("Loaded scvi dataset successfully.")
except Exception as e:
    print(f"scvi dataset failed: {e}")
    print("Trying alternative approach...")

    # Option C: Use a smaller demonstration dataset
    # We will create a multi-modal dataset from the PBMC3k data you already know
    # This won't have real protein data, but it lets you learn the workflow

    adata_rna = sc.datasets.pbmc3k_processed()
    print("Loaded PBMC3k as fallback.")
    print("Shape:", adata_rna.shape)
    print()
    print("NOTE: This fallback dataset does not have real protein data.")
    print("We will simulate protein data so you can still learn the workflow.")
    print("The concepts and code patterns are the same with real data.")

    # Simulate protein data based on known marker genes
    # In real CITE-seq, these would come from antibody measurements
    protein_genes = ['CD3D', 'CD4', 'CD8A', 'CD8B', 'CD14', 'CD19',
                     'NCAM1', 'CD27', 'CD38', 'ITGAX']
    protein_names = ['CD3', 'CD4', 'CD8a', 'CD8b', 'CD14', 'CD19',
                     'CD56', 'CD27', 'CD38', 'CD11c']

    # Extract expression of marker genes and add noise to simulate protein
    protein_data = np.zeros((adata_rna.n_obs, len(protein_genes)))
    for i, gene in enumerate(protein_genes):
        if gene in adata_rna.var_names:
            # Get RNA expression and add realistic noise
            expr = adata_rna[:, gene].X.toarray().flatten() if hasattr(adata_rna[:, gene].X, 'toarray') else adata_rna[:, gene].X.flatten()
            protein_data[:, i] = np.maximum(0, expr + np.random.normal(0, 0.3, len(expr)))

    adata_rna.obsm['protein_expression'] = pd.DataFrame(
        protein_data,
        index=adata_rna.obs_names,
        columns=protein_names
    )
    adata = adata_rna
    print("Simulated protein data added to adata.obsm['protein_expression']")
```

**Common mistake:** If you get a `ModuleNotFoundError` for `scvi`, go back to your terminal and run `pip install scvi-tools`. Then restart your notebook kernel (Kernel menu > Restart).

### 4.4 Explore the Dataset

```python
# Cell 3: Explore what we have
print("=== Dataset Overview ===")
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes (RNA): {adata.n_vars}")
print()

# Check what is stored in the object
print("=== Observation metadata (per-cell info) ===")
print(adata.obs.columns.tolist())
print()

# Check for protein data
if 'protein_expression' in adata.obsm:
    protein_df = adata.obsm['protein_expression']
    if isinstance(protein_df, pd.DataFrame):
        print(f"=== Protein (ADT) Data ===")
        print(f"Number of proteins measured: {protein_df.shape[1]}")
        print(f"Protein names: {protein_df.columns.tolist()[:20]}")
        if len(protein_df.columns) > 20:
            print(f"  ... and {len(protein_df.columns) - 20} more")
    else:
        print(f"Protein data shape: {protein_df.shape}")
else:
    print("No protein data found in adata.obsm['protein_expression']")
    print("Available keys in obsm:", list(adata.obsm.keys()))
```

**Expected output:** You should see the number of cells, genes, and a list of measured proteins. CITE-seq panels typically include 100-300 proteins, though some panels are smaller.

```python
# Cell 4: Look at cell type annotations (if available)
if 'celltype' in adata.obs.columns:
    ct_col = 'celltype'
elif 'cell_type' in adata.obs.columns:
    ct_col = 'cell_type'
else:
    # Find the most likely cell type column
    ct_candidates = [c for c in adata.obs.columns if 'type' in c.lower() or 'cell' in c.lower()]
    ct_col = ct_candidates[0] if ct_candidates else None
    print("Available columns:", adata.obs.columns.tolist())

if ct_col:
    print(f"Cell types (column: '{ct_col}'):")
    print(adata.obs[ct_col].value_counts())
else:
    print("No cell type column found. We will annotate cells ourselves.")
```

**Expected output:** A table showing different immune cell types and how many cells belong to each type. You should see familiar names from Week 2 like T cells, B cells, monocytes, and NK cells.

### 4.5 Process the RNA Data

This should feel familiar from Week 2. We are doing the same standard pipeline.

```python
# Cell 5: RNA preprocessing
# Store raw counts for later use with scvi-tools
adata.layers['counts'] = adata.X.copy()

# Standard RNA preprocessing (same as Week 2)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3',
                             layer='counts')
print(f"Found {adata.var['highly_variable'].sum()} highly variable genes")

# PCA on RNA
sc.tl.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

print("RNA preprocessing complete.")
```

**Expected output:**
```
Found 3000 highly variable genes
RNA preprocessing complete.
```

```python
# Cell 6: Visualize RNA-only UMAP
# Use whatever cell type column we found, or cluster if no annotation exists
if ct_col:
    sc.pl.umap(adata, color=ct_col, title='RNA-only UMAP (by cell type)',
               frameon=False)
else:
    sc.tl.leiden(adata, resolution=0.8)
    sc.pl.umap(adata, color='leiden', title='RNA-only UMAP (by cluster)',
               frameon=False)
```

**Expected output:** A UMAP plot colored by cell type or cluster. You should see distinct groups of cells separating in the plot. This is your baseline -- what RNA alone can tell you.

### 4.6 Process the Protein (ADT) Data

Now we process the protein data. This is where things get different from Week 2.

```python
# Cell 7: Extract and examine protein data
protein_data = adata.obsm['protein_expression']

if isinstance(protein_data, pd.DataFrame):
    protein_matrix = protein_data.values
    protein_names = protein_data.columns.tolist()
elif isinstance(protein_data, np.ndarray):
    protein_matrix = protein_data
    protein_names = [f'Protein_{i}' for i in range(protein_data.shape[1])]
else:
    # Might be sparse
    protein_matrix = np.array(protein_data.todense()) if hasattr(protein_data, 'todense') else np.array(protein_data)
    protein_names = [f'Protein_{i}' for i in range(protein_matrix.shape[1])]

print(f"Protein data shape: {protein_matrix.shape}")
print(f"Number of proteins: {len(protein_names)}")
print(f"First 10 proteins: {protein_names[:10]}")
print()

# Basic statistics
print("=== Protein expression statistics ===")
print(f"Min value: {protein_matrix.min():.2f}")
print(f"Max value: {protein_matrix.max():.2f}")
print(f"Mean value: {protein_matrix.mean():.2f}")
```

**Expected output:** You will see the protein matrix dimensions and basic statistics. Notice how there are far fewer proteins (~100-300) than genes (~20,000).

### 4.7 Normalize Protein Data with CLR

Protein data requires different normalization than RNA data. The standard approach is **Centered Log-Ratio (CLR)** transformation.

```python
# Cell 8: CLR normalization for protein data
# CLR = Centered Log-Ratio transformation
# This is the standard normalization for ADT (protein) data in CITE-seq
# It accounts for the fact that total protein signal varies between cells

def clr_normalize(data):
    """
    Centered Log-Ratio normalization.
    For each cell: log(x + 1) - mean(log(x + 1)) across all proteins.
    This removes cell-to-cell differences in total protein capture.
    """
    # Add pseudocount to avoid log(0)
    data_pseudo = data + 1

    # Log transform
    log_data = np.log(data_pseudo)

    # Subtract the geometric mean (mean of log values) per cell
    geometric_mean = log_data.mean(axis=1, keepdims=True)
    clr_data = log_data - geometric_mean

    return clr_data

protein_clr = clr_normalize(protein_matrix)

print("CLR normalization complete.")
print(f"Before CLR - range: [{protein_matrix.min():.2f}, {protein_matrix.max():.2f}]")
print(f"After CLR  - range: [{protein_clr.min():.2f}, {protein_clr.max():.2f}]")
print()
print("Why CLR? Unlike RNA normalization (which scales to total counts),")
print("CLR handles the fact that protein measurements have different background")
print("levels and different antibody efficiencies.")
```

**Expected output:** You should see that after CLR, the data range is centered around zero, with some positive and some negative values. This is expected -- CLR transforms the data relative to the cell's overall protein signal.

### 4.8 Create a Protein-only UMAP

```python
# Cell 9: Create a protein-only embedding
# Make a separate AnnData object for the protein data
adata_protein = ad.AnnData(
    X=protein_clr,
    obs=adata.obs.copy()
)
adata_protein.var_names = protein_names

# PCA and UMAP on protein data
# Note: with only ~100-300 proteins, we use fewer PCs
n_pcs_protein = min(30, protein_clr.shape[1] - 1)
sc.tl.pca(adata_protein, n_comps=n_pcs_protein)
sc.pp.neighbors(adata_protein, n_neighbors=15, n_pcs=n_pcs_protein)
sc.tl.umap(adata_protein)

# Store protein UMAP coordinates back in the main object
adata.obsm['X_umap_protein'] = adata_protein.obsm['X_umap']

print(f"Protein UMAP computed using {n_pcs_protein} PCs.")
```

### 4.9 Compare RNA-only vs Protein-only UMAPs

This is one of the most instructive visualizations in multi-omics analysis.

```python
# Cell 10: Side-by-side comparison
color_col = ct_col if ct_col else 'leiden'

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# RNA-only UMAP
sc.pl.umap(adata, color=color_col, ax=axes[0], show=False,
           title='RNA-only UMAP', frameon=False)

# Protein-only UMAP
# We need to temporarily swap the UMAP coordinates
original_umap = adata.obsm['X_umap'].copy()
adata.obsm['X_umap'] = adata.obsm['X_umap_protein']

sc.pl.umap(adata, color=color_col, ax=axes[1], show=False,
           title='Protein-only UMAP', frameon=False)

# Restore original RNA UMAP
adata.obsm['X_umap'] = original_umap

plt.tight_layout()
plt.savefig('rna_vs_protein_umap.png', dpi=150, bbox_inches='tight')
plt.show()

print("Figure saved as 'rna_vs_protein_umap.png'")
print()
print("NOTICE: The two UMAPs look different! Some cell types separate")
print("better in protein space (especially T cell subtypes defined by")
print("surface markers like CD4 and CD8), while RNA space captures")
print("more transcriptional nuance.")
```

**Expected output:** Two side-by-side UMAP plots. You should notice:
- Some cell types cluster more tightly in the protein UMAP (e.g., CD4+ vs CD8+ T cells, because their defining markers are surface proteins)
- The RNA UMAP may show more gradual transitions between states (because it captures ~20,000 features vs ~100-300)
- Some clusters may appear in one plot but not the other

### 4.10 The WNN Concept

Before we move on, let's understand a key idea: **Weighted Nearest Neighbors (WNN)**.

WNN was introduced by Hao et al. (Cell, 2021) as part of Seurat v4. The core insight is simple but powerful:

- For some cells, RNA is more informative (e.g., a rare transcriptional state that is not captured by the limited protein panel)
- For other cells, protein is more informative (e.g., CD4 vs CD8 T cells are cleanly separated by protein but noisy in RNA)
- WNN learns a **per-cell weight** for each modality, then combines them

Think of it this way: if you are trying to identify someone and you have both a photo and a voice recording, WNN says "for this person, the photo is more helpful" and "for that person, the voice is more helpful," rather than treating both equally for everyone.

We will not implement WNN from scratch here (it is built into Seurat v5 in R), but you will use a similar idea through totalVI in Section 8.

> **Checkpoint 4.1:** You have loaded a CITE-seq dataset, processed both the RNA and protein components, created UMAPs for each, and can see how they differ. If your two UMAP plots look different from each other, that is correct -- it means the two modalities are capturing different information.

---

## 5. Understanding scATAC-seq (Conceptual) [REQUIRED]

We will not run a full scATAC-seq pipeline this week (that is part of the Final Practice), but understanding the concepts is essential for multi-omics integration.

### 5.1 What Is Chromatin Accessibility?

Every cell in your body contains the same DNA -- the complete instruction manual. But different cell types only read certain chapters. How does this work?

DNA is not floating freely in the nucleus. It is wrapped tightly around proteins called **histones**, like thread wound around spools. This compacted structure is called **chromatin**.

- **Open chromatin** = regions where the DNA is loosely wrapped and physically accessible to the cell's machinery. These regions CAN be read, meaning the genes there can potentially be turned on
- **Closed chromatin** = regions where the DNA is tightly compacted and inaccessible. These genes are effectively silenced

This is like a library where some books are on the shelves (open, accessible) and others are locked in storage (closed, inaccessible). Just because a book is on the shelf does not mean someone is reading it right now -- but it COULD be read. That distinction between "could be used" and "is being used" is what makes chromatin accessibility so informative alongside gene expression.

### 5.2 How scATAC-seq Works

**scATAC-seq** stands for **s**ingle-**c**ell **A**ssay for **T**ransposase-**A**ccessible **C**hromatin using **seq**uencing. Here is the simplified version:

1. An enzyme called **Tn5 transposase** is added to the cells
2. Tn5 can only insert itself into DNA that is accessible (open)
3. Where Tn5 inserts, it cuts the DNA and adds sequencing adapters
4. The resulting fragments are sequenced
5. By mapping where the fragments came from, you know which regions were open

### 5.3 Key Concepts in scATAC-seq

**Peaks:** When you pile up all the fragments across many cells, certain genomic regions have many overlapping fragments. These pile-ups are called "peaks" -- they mark consistently open chromatin regions. Peak calling is analogous to finding the most popular books in the library.

**Fragments:** The raw data from scATAC-seq. Each fragment represents a piece of accessible DNA from a single cell. A typical cell yields 5,000-50,000 fragments.

**Gene activity:** An estimated measure of gene expression derived from chromatin data. The idea is: if the chromatin around gene X is open, gene X is likely to be expressed. This is calculated by counting fragments that overlap the gene body and promoter region. Gene activity is noisier than actual RNA measurement but provides a bridge between ATAC and RNA data.

**TSS enrichment:** A quality control metric. TSS = Transcription Start Site (the beginning of a gene). In good-quality data, fragments should pile up near TSS regions because promoters are generally in open chromatin. A TSS enrichment score above 4-5 indicates good data quality.

**Nucleosome signal:** Another QC metric. DNA wraps around nucleosomes in ~147 base-pair units. In scATAC-seq, you see a characteristic pattern of fragment lengths: short fragments (~150 bp) that fit between nucleosomes, and longer fragments (~300 bp, ~450 bp) that span 1 or 2 nucleosomes. A clear pattern means good data.

**TF motifs:** Transcription factors (TFs) are proteins that bind to specific DNA sequences to turn genes on or off. Each TF recognizes a specific short sequence pattern called a **motif** (e.g., the GATA family of TFs binds to the sequence "GATA"). By searching for these motif patterns in open chromatin regions, you can infer which TFs are likely active in each cell. This is a powerful way to understand gene regulation.

### 5.4 Why scATAC-seq Matters for Immunology

In the immune system, cells frequently change state: naive T cells become activated, then become memory cells or effector cells. These transitions involve massive chromatin remodeling -- opening new genomic regions and closing old ones.

The remarkable finding is that chromatin changes often happen BEFORE the corresponding changes in gene expression. This means scATAC-seq can reveal a cell's future trajectory -- where it is headed -- before RNA-seq can. This is like seeing someone pack a suitcase before they announce their trip.

For T cell development specifically (which you will explore in the Final Practice), scATAC-seq has revealed:
- The chromatin landscape of CD4/CD8 lineage commitment changes before the mRNA profile does
- Specific transcription factor networks (like ThPOK for CD4, Runx3 for CD8) create distinct chromatin signatures
- Regulatory elements far from gene bodies control lineage-specific gene expression

> **Checkpoint 5.1:** You can explain in your own words: (a) what chromatin accessibility is, (b) how scATAC-seq measures it, and (c) why it matters alongside RNA-seq. Write a 2-3 sentence explanation and compare it with the summary above.

---

## 6. Multi-omics Integration Concepts [REQUIRED]

Now that you understand what different technologies measure, the question is: how do you combine them computationally? The answer depends on your experimental design.

### 6.1 Three Scenarios

#### Scenario 1: Same Cells, Multiple Measurements

**Examples:** CITE-seq (RNA + protein), 10x Multiome (RNA + ATAC)

This is the easiest scenario. Because RNA and protein (or RNA and ATAC) come from the SAME cell, you have a natural pairing. Cell #1234 has both an RNA profile and a protein profile. You can directly compare and combine them.

**Integration methods for this scenario:**
- WNN (Weighted Nearest Neighbors) -- weights each modality per cell
- totalVI -- jointly models RNA + protein
- MOFA+ -- finds shared and modality-specific factors

#### Scenario 2: Different Cells, Same Sample

**Example:** You ran scRNA-seq on one aliquot of your blood sample and scATAC-seq on another aliquot. The cells are different, but they come from the same biological sample.

This is harder because there is no one-to-one cell mapping. You need computational "bridges" to link cells across modalities.

**Integration methods for this scenario:**
- **Label transfer:** Annotate cell types in one modality, then predict labels in the other
- **Gene activity scores:** Convert ATAC data into estimated gene expression, then find nearest neighbors between the two
- **Anchor-based integration:** Find pairs of cells that are similar across modalities (Seurat v3/v5)
- **MultiVI:** Deep learning model that handles missing modalities (scvi-tools)

#### Scenario 3: Different Samples

**Example:** You have RNA-seq from Patient A and ATAC-seq from Patient B.

This is the hardest scenario because you need to account for both modality differences AND biological/batch differences between samples. You typically need batch correction on top of cross-modality integration.

### 6.2 Integration Methods Explained Simply

#### WNN (Weighted Nearest Neighbors)

**What it does:** For each cell, WNN determines how much to trust the RNA data vs the protein (or ATAC) data. It then builds a combined nearest-neighbor graph using per-cell weights.

**Analogy:** You are a teacher grading a student based on a written exam and an oral exam. For one student, the written exam is more informative (they are a great writer). For another, the oral exam is better (they explain things well verbally). WNN adjusts the weights per student.

**Tool:** Seurat v5 (R). There is no direct Python implementation, but muon has a similar approach.

#### totalVI

**What it does:** A deep generative model (variational autoencoder) that jointly models RNA counts and protein counts. It learns a shared low-dimensional representation (latent space) where both modalities are combined. It also handles batch effects.

**Analogy:** Imagine an artist who can create a portrait of a person by combining their photograph and their voice recording into a single representation. The artist learns which features from each source are most important.

**Tool:** scvi-tools (Python). This is what we will use in Section 8.

#### MultiVI

**What it does:** Extends totalVI to handle RNA + ATAC data, and can even work when some cells have only RNA and others have only ATAC (the "missing modality" problem from Scenario 2).

**Analogy:** Same artist, but now they can also work with DNA accessibility data, and they can fill in the blanks even if some people only provided a photo without a voice recording.

**Tool:** scvi-tools (Python).

#### Label Transfer

**What it does:** Train a cell type classifier on one dataset (the "reference") and predict cell type labels in another dataset (the "query"). This is the simplest form of integration -- you are not combining the data itself, just transferring annotations.

**Analogy:** You have a well-labeled photo album (reference). When someone shows you a new photo (query), you find the most similar photo in your album and use its label.

**Tool:** Seurat, CellTypist, scArches.

### 6.3 How to Choose an Integration Method

| Your Situation | Recommended Method |
|---------------|-------------------|
| Same cells, RNA + protein | totalVI or WNN |
| Same cells, RNA + ATAC | WNN, MOFA+, or directly combine latent spaces |
| Different cells, same sample | MultiVI, label transfer, or anchor-based integration |
| Different cells, different samples | MultiVI with batch correction, or scArches |

> **Checkpoint 6.1:** You can name the three integration scenarios and describe one method for each. You understand why same-cell multi-modal data is easier to integrate than data from separate experiments.

---

## 7. Hands-On: Working with MuData (Multi-modal Data) [REQUIRED]

### 7.1 What Is MuData?

In Weeks 2 and 3, you worked with **AnnData** objects -- the standard container for single-modal single-cell data. **MuData** (from the `muon` package) extends this to multiple modalities.

Think of it this way:
- **AnnData** = a single spreadsheet with cells as rows and features as columns
- **MuData** = a binder containing multiple spreadsheets (one per modality), all sharing the same rows (cells)

```
MuData
  |-- 'rna'     -> AnnData (cells x genes)
  |-- 'protein'  -> AnnData (cells x proteins)
  |-- 'atac'    -> AnnData (cells x peaks)  [if you have it]
```

### 7.2 Create a MuData Object

```python
# Cell 11: Import muon and create MuData
import muon as mu
from muon import MuData

# We already have our RNA data in adata and protein data in adata_protein
# Let's create a proper MuData object

# First, create a clean RNA AnnData
adata_rna = adata.copy()

# Create protein AnnData from the CLR-normalized data
adata_prot = ad.AnnData(
    X=protein_clr,
    obs=adata.obs[[ct_col]].copy() if ct_col else adata.obs[['leiden']].copy()
)
adata_prot.var_names = protein_names
adata_prot.obs_names = adata.obs_names.copy()

# Create the MuData object
mdata = MuData({
    'rna': adata_rna,
    'protein': adata_prot
})

print("=== MuData Object ===")
print(mdata)
print()
print(f"Modalities: {list(mdata.mod.keys())}")
print(f"RNA shape: {mdata['rna'].shape}")
print(f"Protein shape: {mdata['protein'].shape}")
print(f"Shared cells: {mdata.n_obs}")
```

**Expected output:**
```
=== MuData Object ===
MuData object with n_obs x n_vars = XXXX x XXXXX
  2 modalities
    rna:     XXXX x XXXXX
    protein: XXXX x XXX

Modalities: ['rna', 'protein']
RNA shape: (XXXX, XXXXX)
Protein shape: (XXXX, XXX)
Shared cells: XXXX
```

### 7.3 Explore the MuData Object

```python
# Cell 12: Navigate the MuData structure
# Access individual modalities
rna_data = mdata['rna']      # or mdata.mod['rna']
prot_data = mdata['protein']  # or mdata.mod['protein']

print("=== Accessing Modalities ===")
print(f"RNA genes (first 10): {rna_data.var_names[:10].tolist()}")
print(f"Proteins (first 10): {prot_data.var_names[:10].tolist()}")
print()

# MuData keeps track of which cells are in which modality
print("=== Shared Observations ===")
print(f"Cells in RNA: {rna_data.n_obs}")
print(f"Cells in Protein: {prot_data.n_obs}")
print(f"Cells in both: {mdata.n_obs}")
print()

# Access shared metadata
print("=== Shared Metadata ===")
print(mdata.obs.head())
```

### 7.4 Visualize Individual Modalities from MuData

```python
# Cell 13: Plot from MuData modalities
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot a few marker genes from RNA
if 'CD3D' in mdata['rna'].var_names:
    sc.pl.umap(mdata['rna'], color='CD3D', ax=axes[0], show=False,
               title='CD3D Gene Expression (RNA)')
elif 'CD3E' in mdata['rna'].var_names:
    sc.pl.umap(mdata['rna'], color='CD3E', ax=axes[0], show=False,
               title='CD3E Gene Expression (RNA)')
else:
    axes[0].set_title('(CD3 gene not found in RNA data)')

# Plot a protein if available
if 'CD3' in mdata['protein'].var_names:
    sc.pl.umap(mdata['protein'], color='CD3', ax=axes[1], show=False,
               title='CD3 Protein Level')
elif len(mdata['protein'].var_names) > 0:
    first_prot = mdata['protein'].var_names[0]
    sc.pl.umap(mdata['protein'], color=first_prot, ax=axes[1], show=False,
               title=f'{first_prot} Protein Level')

plt.tight_layout()
plt.savefig('mudata_modalities.png', dpi=150, bbox_inches='tight')
plt.show()
print("Figure saved as 'mudata_modalities.png'")
```

### 7.5 MuData Tips and Tricks

```python
# Cell 14: Useful MuData operations

# Subset cells (applies to all modalities at once)
# For example, keep only T cells:
if ct_col and ct_col in mdata.obs.columns:
    t_cell_mask = mdata.obs[ct_col].str.contains('T cell|CD4|CD8', case=False, na=False)
    if t_cell_mask.sum() > 0:
        mdata_tcells = mdata[t_cell_mask].copy()
        print(f"Subset to T cells: {mdata_tcells.n_obs} cells across {len(mdata_tcells.mod)} modalities")
    else:
        print("No T cells found with the current filter. This is OK -- the column names may differ.")
        print(f"Available cell types: {mdata.obs[ct_col].unique()[:10]}")

# Save MuData to disk (for later use)
# mdata.write("my_cite_seq_data.h5mu")
# Load it back: mdata = mu.read("my_cite_seq_data.h5mu")
print()
print("MuData objects can be saved as .h5mu files (like .h5ad for AnnData).")
print("Use mdata.write('filename.h5mu') to save and mu.read('filename.h5mu') to load.")
```

> **Checkpoint 7.1:** You have created a MuData object with RNA and protein modalities, accessed individual modalities, and understand the relationship between MuData and AnnData. Try accessing `mdata['rna'].var_names` and `mdata['protein'].var_names` to confirm you can navigate the structure.

---

## 8. Hands-On: Integration with totalVI [REQUIRED]

Now we get to the exciting part: using a deep learning model to integrate RNA and protein data into a single representation.

### 8.1 What Is totalVI?

**totalVI** (Total Variational Inference) is a generative model from the scvi-tools suite. It was developed by Gayoso et al. (Nature Methods, 2021) and is specifically designed for CITE-seq data.

What it does:
1. Takes raw RNA counts and protein counts as input
2. Uses a neural network (variational autoencoder) to learn a shared low-dimensional representation
3. This "latent space" captures information from BOTH modalities
4. It also handles batch effects if your data comes from multiple experiments

Why it works well:
- It models the data-generating process: RNA counts follow a negative binomial distribution, and protein counts follow a mixture of negative binomial distributions (to account for background signal)
- By modeling both together, it finds patterns that neither modality would reveal alone

### 8.2 Prepare Data for totalVI

```python
# Cell 15: Set up data for totalVI
import scvi

# totalVI needs:
# 1. Raw RNA counts (not normalized) in adata.X or a layer
# 2. Protein data stored in adata.obsm

# Let's prepare a fresh copy of our data
adata_total = adata.copy()

# Make sure we have raw counts
if 'counts' in adata_total.layers:
    adata_total.X = adata_total.layers['counts'].copy()
    print("Using raw counts from layers['counts']")
else:
    print("WARNING: No raw counts found. totalVI works best with raw counts.")
    print("Using current .X (which may be normalized).")

# Make sure protein data is in obsm
if 'protein_expression' not in adata_total.obsm:
    # Store it there
    adata_total.obsm['protein_expression'] = protein_matrix
    print("Stored protein data in obsm['protein_expression']")

# Filter to highly variable genes to speed up training
if 'highly_variable' in adata_total.var.columns:
    adata_total = adata_total[:, adata_total.var['highly_variable']].copy()

print(f"Data ready for totalVI: {adata_total.shape[0]} cells, {adata_total.shape[1]} genes")
print(f"Protein features: {adata_total.obsm['protein_expression'].shape[1]}")
```

### 8.3 Set Up and Train totalVI

```python
# Cell 16: Set up totalVI model
# Tell scvi-tools where to find the RNA and protein data
scvi.model.TOTALVI.setup_anndata(
    adata_total,
    protein_expression_obsm_key="protein_expression",
    layer="counts" if "counts" in adata_total.layers else None,
    # batch_key="batch"  # Uncomment if your data has batches
)

# Create the model
model = scvi.model.TOTALVI(
    adata_total,
    latent_distribution="normal",  # Use normal distribution for the latent space
    n_latent=20                    # 20-dimensional latent space
)

print("Model created:")
print(model)
```

**Expected output:** A summary of the totalVI model showing the number of input genes, proteins, and latent dimensions.

```python
# Cell 17: Train the model
# This is the computationally intensive step
# On a laptop, it may take 5-15 minutes depending on dataset size
# On a GPU, it would take 1-2 minutes

print("Training totalVI model...")
print("This may take 5-15 minutes. Watch the loss decrease over epochs.")
print()

model.train(
    max_epochs=200,       # Maximum number of training epochs
    early_stopping=True,  # Stop early if the model converges
    early_stopping_patience=10,  # Wait 10 epochs before stopping
    train_size=0.9        # Use 90% for training, 10% for validation
)

print("Training complete!")
```

**Expected output:** You will see training progress with the loss value decreasing over epochs. It should start around 2000-5000 and decrease to a lower value. If you see `early_stopping` triggered, that means the model converged before reaching 200 epochs -- that is normal and good.

**Common mistakes:**
- If training fails with a CUDA/GPU error, that is fine -- totalVI works on CPU too, it is just slower
- If you get a memory error, try reducing the data: `adata_total = adata_total[:5000, :]` to use only 5000 cells
- If training loss goes to NaN or infinity, try reducing the learning rate: `model.train(max_epochs=200, lr=1e-4)`

```python
# Cell 18: Check training history
# Plot the training loss to make sure the model converged
train_history = model.history

fig, ax = plt.subplots(1, 1, figsize=(8, 4))
ax.plot(train_history['elbo_train']['elbo_train'], label='Training loss')
if 'elbo_validation' in train_history:
    ax.plot(train_history['elbo_validation']['elbo_validation'], label='Validation loss')
ax.set_xlabel('Epoch')
ax.set_ylabel('ELBO Loss')
ax.set_title('totalVI Training History')
ax.legend()
plt.tight_layout()
plt.savefig('totalvi_training.png', dpi=150, bbox_inches='tight')
plt.show()
print("Figure saved as 'totalvi_training.png'")
print()
print("The loss should decrease and then flatten out. If it is still decreasing")
print("at the end, you may want to increase max_epochs.")
```

**Expected output:** A line plot showing the loss decreasing over epochs and eventually flattening. Both training and validation loss should follow similar trends. If the validation loss starts increasing while training loss continues decreasing, that indicates overfitting (but totalVI is fairly robust to this).

### 8.4 Extract the Integrated Representation

```python
# Cell 19: Get the latent representation
# This is the integrated embedding that combines RNA and protein information
latent = model.get_latent_representation()
adata_total.obsm["X_totalvi"] = latent

print(f"Latent representation shape: {latent.shape}")
print(f"This means: {latent.shape[0]} cells in a {latent.shape[1]}-dimensional space")
print()
print("This 20-dimensional representation captures information from BOTH")
print("RNA (~3000 genes) and protein (~100-300 markers) in a compact form.")
```

### 8.5 Compute the Integrated UMAP

```python
# Cell 20: Compute UMAP from the integrated representation
sc.pp.neighbors(adata_total, use_rep="X_totalvi", n_neighbors=15)
sc.tl.umap(adata_total)

# Also cluster the integrated data
sc.tl.leiden(adata_total, resolution=0.8)

print("Integrated UMAP computed.")
print(f"Number of clusters: {adata_total.obs['leiden'].nunique()}")
```

```python
# Cell 21: Visualize the integrated UMAP
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Color by cell type (if available)
if ct_col and ct_col in adata_total.obs.columns:
    sc.pl.umap(adata_total, color=ct_col, ax=axes[0], show=False,
               title='totalVI Integrated UMAP (cell type)', frameon=False)
else:
    sc.pl.umap(adata_total, color='leiden', ax=axes[0], show=False,
               title='totalVI Integrated UMAP (clusters)', frameon=False)

# Color by Leiden clusters from integrated data
sc.pl.umap(adata_total, color='leiden', ax=axes[1], show=False,
           title='totalVI Integrated UMAP (Leiden clusters)', frameon=False)

plt.tight_layout()
plt.savefig('totalvi_integrated_umap.png', dpi=150, bbox_inches='tight')
plt.show()
print("Figure saved as 'totalvi_integrated_umap.png'")
```

**Expected output:** A UMAP plot that should show well-separated cell populations. Because totalVI uses information from both RNA and protein, you should see:
- Cleaner separation of cell types compared to RNA-only or protein-only UMAPs
- T cell subtypes (CD4, CD8) should be well-resolved (protein helps here)
- Monocyte and B cell populations should also be distinct (RNA helps here)

### 8.6 What Just Happened?

Let's pause and understand what totalVI did:

1. **Input:** Raw RNA counts (~3000 genes) + raw protein counts (~100-300 proteins) for each cell
2. **Encoding:** A neural network compressed both data types into a shared 20-dimensional latent space
3. **Learning:** The model learned which RNA genes and which proteins are most informative for distinguishing cell types, and how to weight them
4. **Output:** A 20-dimensional vector per cell that captures the combined RNA + protein signal

This is fundamentally different from simply concatenating RNA and protein features, because:
- totalVI accounts for the different noise characteristics of each modality
- It handles the fact that RNA has ~3000 features while protein has ~100-300
- It models the data-generating process (count distributions), not just arbitrary numbers

> **Checkpoint 8.1:** You have trained a totalVI model and generated an integrated UMAP. The integrated UMAP should look as good as or better than either the RNA-only or protein-only UMAP. If the plot looks reasonable with distinct cell clusters, you are on track.

---

## 9. Comparing Modalities [REQUIRED]

One of the most valuable things you can do with multi-modal data is compare measurements across modalities. Where do they agree? Where do they disagree? The disagreements are often the most biologically interesting.

### 9.1 RNA vs Protein for Key Markers

```python
# Cell 22: Compare RNA and protein for CD4
# First, find the gene and protein names in your dataset

# Common marker pairs: (RNA gene name, Protein name)
marker_pairs = [
    ('CD4', 'CD4'),
    ('CD8A', 'CD8a'),
    ('CD14', 'CD14'),
    ('CD19', 'CD19'),
    ('NCAM1', 'CD56'),   # NCAM1 gene encodes CD56 protein
    ('CD3D', 'CD3'),
    ('CD3E', 'CD3'),
]

# Find which markers are available in both modalities
available_pairs = []
rna_var = set(adata_total.var_names)
prot_names_set = set(protein_names)

for rna_name, prot_name in marker_pairs:
    rna_found = rna_name in rna_var
    prot_found = prot_name in prot_names_set
    if rna_found and prot_found:
        available_pairs.append((rna_name, prot_name))
        print(f"  Found pair: {rna_name} (RNA) <-> {prot_name} (protein)")
    else:
        status = []
        if not rna_found:
            status.append(f"{rna_name} not in RNA")
        if not prot_found:
            status.append(f"{prot_name} not in protein")
        print(f"  Missing: {', '.join(status)}")

print(f"\nAvailable pairs for comparison: {len(available_pairs)}")
```

```python
# Cell 23: Side-by-side comparison plots
if len(available_pairs) > 0:
    n_pairs = min(len(available_pairs), 4)  # Show up to 4 pairs
    fig, axes = plt.subplots(n_pairs, 2, figsize=(12, 4 * n_pairs))

    if n_pairs == 1:
        axes = axes.reshape(1, 2)

    for i, (rna_name, prot_name) in enumerate(available_pairs[:n_pairs]):
        # RNA expression
        sc.pl.umap(adata_total, color=rna_name, ax=axes[i, 0], show=False,
                   title=f'{rna_name} Gene Expression (RNA)',
                   frameon=False, vmin=0)

        # Protein expression
        # We need to add the protein data as a column in obs temporarily
        prot_idx = protein_names.index(prot_name)
        adata_total.obs[f'protein_{prot_name}'] = protein_clr[:adata_total.n_obs, prot_idx]

        sc.pl.umap(adata_total, color=f'protein_{prot_name}', ax=axes[i, 1],
                   show=False, title=f'{prot_name} Surface Protein (ADT)',
                   frameon=False)

    plt.tight_layout()
    plt.savefig('rna_vs_protein_markers.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("Figure saved as 'rna_vs_protein_markers.png'")
else:
    print("No matching RNA-protein pairs found in this dataset.")
    print("This can happen with certain dataset configurations.")
    print("Try printing adata_total.var_names and protein_names to see what is available.")
```

**Expected output:** Side-by-side UMAP plots for each marker, with RNA on the left and protein on the right. You should notice:

- **Agreement:** For strong markers like CD3 (T cells) and CD14 (monocytes), RNA and protein should highlight the same clusters
- **Disagreement:** Some cells may have RNA but not protein (the gene was recently turned on but protein has not accumulated yet) or protein but not RNA (the gene was turned off but old protein is still on the surface)
- **Sharpness:** Protein measurements often look "sharper" -- cells are either positive or negative, with less in-between. RNA is noisier with more intermediate values

### 9.2 Scatter Plot: RNA vs Protein Correlation

```python
# Cell 24: Direct scatter comparison
if len(available_pairs) > 0:
    rna_name, prot_name = available_pairs[0]

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))

    # Get RNA expression for this gene
    if hasattr(adata_total[:, rna_name].X, 'toarray'):
        rna_expr = adata_total[:, rna_name].X.toarray().flatten()
    else:
        rna_expr = np.array(adata_total[:, rna_name].X).flatten()

    # Get protein expression
    prot_idx = protein_names.index(prot_name)
    prot_expr = protein_clr[:adata_total.n_obs, prot_idx]

    # Color by cell type if available
    if ct_col and ct_col in adata_total.obs.columns:
        categories = adata_total.obs[ct_col].astype('category')
        colors = plt.cm.tab20(np.linspace(0, 1, len(categories.cat.categories)))
        color_map = dict(zip(categories.cat.categories, colors))
        point_colors = [color_map[c] for c in categories]

        scatter = ax.scatter(rna_expr, prot_expr, c=point_colors,
                           alpha=0.3, s=5, rasterized=True)

        # Add legend
        handles = [plt.Line2D([0], [0], marker='o', color='w',
                  markerfacecolor=color_map[cat], markersize=8, label=cat)
                  for cat in list(categories.cat.categories)[:15]]
        ax.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left',
                 fontsize=7, frameon=False)
    else:
        ax.scatter(rna_expr, prot_expr, alpha=0.3, s=5, c='steelblue',
                  rasterized=True)

    ax.set_xlabel(f'{rna_name} RNA Expression', fontsize=12)
    ax.set_ylabel(f'{prot_name} Protein Level (CLR)', fontsize=12)
    ax.set_title(f'RNA vs Protein: {rna_name}/{prot_name}', fontsize=14)

    # Add correlation
    correlation = np.corrcoef(rna_expr, prot_expr)[0, 1]
    ax.text(0.05, 0.95, f'Pearson r = {correlation:.3f}',
            transform=ax.transAxes, fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('rna_protein_scatter.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("Figure saved as 'rna_protein_scatter.png'")
    print()
    print(f"Correlation between {rna_name} RNA and {prot_name} protein: r = {correlation:.3f}")
    print()
    if correlation > 0.5:
        print("Moderate to strong correlation. RNA and protein largely agree for this marker.")
    elif correlation > 0.2:
        print("Weak to moderate correlation. There is some agreement, but also meaningful")
        print("differences -- some cells with high RNA have low protein, and vice versa.")
    else:
        print("Weak correlation. RNA and protein tell quite different stories for this marker.")
        print("This is common and biologically meaningful -- it reflects post-transcriptional")
        print("regulation, protein turnover, and timing differences.")
```

**Expected output:** A scatter plot with RNA expression on the x-axis and protein level on the y-axis. The Pearson correlation will likely be moderate (r = 0.3-0.6) for most markers. This is normal and expected -- if RNA and protein agreed perfectly, there would be no point in measuring both.

### 9.3 Discussion: Why Do RNA and Protein Disagree?

After looking at your plots, you should understand these key points:

1. **Protein is the "ground truth" for surface markers.** When immunologists define cell types (CD4+ T cells, CD8+ T cells, etc.), they are referring to protein expression, not RNA. This is because proteins are what other cells actually interact with.

2. **RNA is a leading indicator.** When a cell starts to change state, RNA changes first (within minutes to hours), while protein changes lag behind (hours to days). This means RNA can detect cells in transition that protein has not caught up to yet.

3. **RNA covers more genes.** CITE-seq panels typically measure 100-300 proteins, while RNA captures ~20,000 genes. So RNA wins on breadth, while protein wins on accuracy for the markers it covers.

4. **Post-transcriptional regulation.** Cells actively regulate which mRNA molecules get translated into protein. MicroRNAs, RNA-binding proteins, and other mechanisms can block translation. This means RNA presence does not guarantee protein presence.

5. **Protein degradation.** Even after a cell stops making mRNA for a gene, the existing protein molecules may persist on the surface for hours or days. This creates a "memory" effect in the protein data.

These are not flaws in the data -- they are real biology. And the differences between RNA and protein are precisely why multi-omics is so valuable.

> **Checkpoint 9.1:** You have created RNA vs protein comparison plots and understand why the two modalities can disagree. You can explain at least three reasons why RNA and protein levels for the same gene might not match.

---

## 10. Regulatory Inference from ATAC Data [CHALLENGE]

This section is conceptual with code sketches. You will implement a full ATAC analysis pipeline in the Final Practice.

### 10.1 From Peaks to Gene Regulation

scATAC-seq gives you a map of open chromatin regions (peaks). But how do you go from "this DNA region is accessible" to "this gene is regulated"?

**Step 1: Link peaks to genes.** Peaks near a gene's promoter (the DNA region just before the gene starts) likely regulate that gene. Peaks far from any gene might be **enhancers** -- distant regulatory elements that can activate genes from thousands of base pairs away.

**Step 2: Estimate gene activity.** By counting all the ATAC fragments that fall within or near a gene, you can estimate how "active" that gene's locus is. This is called a **gene activity score**. It is a rough estimate, but it provides a bridge between ATAC and RNA data.

**Step 3: Find transcription factor motifs.** Open chromatin regions contain short DNA sequences (6-20 base pairs) where transcription factors can bind. By scanning open regions for known TF motif patterns, you can infer which TFs are likely active in each cell.

### 10.2 Code Sketch: ATAC Analysis Pipeline

The following code is a roadmap, not something to run this week. You will implement each step in the Final Practice.

```python
# ============================================================
# CONCEPTUAL CODE SKETCH -- ATAC Analysis Pipeline
# You will implement this in the Final Practice
# ============================================================

# Step 1: Load ATAC fragment data
# The raw input is a fragments.tsv.gz file from 10x Genomics
# Each row = one fragment (chromosome, start, end, cell barcode, count)

# import snapatac2 as snap  # or use ArchR in R
# adata_atac = snap.read_fragments("fragments.tsv.gz")

# Step 2: Quality control
# - TSS enrichment: should be > 4-5 for good quality
# - Number of fragments per cell: typically 5,000-50,000
# - Fraction of reads in peaks (FRiP): should be > 0.3
# - Nucleosome signal: characteristic banding pattern

# Step 3: Call peaks (find consistently open chromatin regions)
# - Use MACS2 or the built-in peak caller
# - Results in a "peaks x cells" matrix (like the "genes x cells" matrix in RNA)

# Step 4: Normalization
# - TF-IDF: Term Frequency - Inverse Document Frequency
#   (borrowed from text mining -- treats peaks like "words" and cells like "documents")
# - LSI: Latent Semantic Indexing (like PCA but for TF-IDF data)

# snap.tl.spectral(adata_atac)  # Compute LSI
# snap.tl.umap(adata_atac)      # Compute UMAP

# Step 5: Gene activity estimation
# - Count fragments overlapping each gene body + promoter
# - This creates a "pseudo-RNA" matrix from ATAC data
# - Useful for label transfer from RNA to ATAC

# gene_activity = snap.tl.gene_activity(adata_atac)

# Step 6: TF motif enrichment with chromVAR
# - For each cell, score the enrichment of each TF motif
#   in its accessible regions compared to the background
# - High score = that TF is likely active in this cell

# snap.tl.motif_enrichment(adata_atac, motifs="JASPAR2022")

# Step 7: Link peaks to genes (co-accessibility)
# - Find correlations between peak accessibility and nearby gene expression
# - Peaks that correlate with a gene are likely regulatory elements

# links = snap.tl.link_peaks_to_genes(adata_atac, adata_rna)
```

### 10.3 Why This Matters for Immunology

In immune cell development and activation, chromatin remodeling is not just a consequence of cell state changes -- it often DRIVES them. Key examples:

- **T cell development in the thymus:** The decision to become CD4+ or CD8+ involves opening lineage-specific enhancers and closing others. These chromatin changes happen before the corresponding gene expression changes, effectively "committing" the cell to its fate.

- **B cell activation:** When a naive B cell encounters its antigen, it undergoes rapid chromatin remodeling at antibody gene loci, enabling class-switch recombination (changing from IgM to IgG, IgA, or IgE).

- **Exhaustion in chronic infection:** T cells fighting prolonged infections (like chronic viral infections or cancer) gradually lose their effector functions. scATAC-seq has revealed that this "exhaustion" is written into the chromatin -- exhausted T cells have a distinct chromatin landscape that is difficult to reverse, explaining why some immunotherapies have limited success.

---

## 11. The Chopp et al. Paper [DEEP DIVE]

### 11.1 Paper Overview

**Chopp et al., "An integrated epigenomic and transcriptomic map of mouse and human alpha-beta T cell development," Immunity, 2020.**

This paper is central to the Final Practice, so let's understand what they did and found.

### 11.2 The Biological Question

T cells develop in the thymus, where they go through several stages. One critical decision point is **lineage commitment**: will a developing T cell become a CD4+ helper T cell or a CD8+ killer T cell? This decision depends on how the T cell receptor (TCR) interacts with MHC molecules during thymic selection (concepts from Week 1).

The question Chopp et al. asked: Can we understand this decision at the level of individual cells, using BOTH gene expression and chromatin accessibility?

### 11.3 What They Did

1. **Sorted developing T cells** from mouse and human thymuses at different developmental stages
2. **Performed bulk RNA-seq and ATAC-seq** on each stage (note: this was before the 10x Multiome technology existed, so they used sorted populations rather than single cells)
3. **Performed scRNA-seq** on unsorted thymocytes to get single-cell resolution
4. **Integrated** the data to map chromatin dynamics onto transcriptional trajectories

### 11.4 Key Findings

1. **Asymmetric lineage emergence:** CD4+ and CD8+ T cells do not emerge symmetrically from the double-positive (DP) stage. The transcriptional programs diverge at different rates.

2. **Chromatin precedes transcription:** Open chromatin regions associated with CD4-lineage genes appeared BEFORE those genes were highly expressed. Similarly for CD8-lineage genes. This means the cell "decides" its fate at the chromatin level before it fully commits at the transcriptional level.

3. **Transcription factor networks:** Specific TFs like ThPOK (encoded by *Zbtb7b*) for CD4 and Runx3 for CD8 create self-reinforcing chromatin loops. Once activated, they open more chromatin around their own target genes, which further activates them -- a positive feedback loop that makes the commitment irreversible.

4. **Conservation across species:** Many of the chromatin regulatory patterns were conserved between mouse and human, suggesting they are fundamental to T cell development.

### 11.5 Why This Paper Matters for Your Training

Chopp et al. demonstrated the power of combining transcriptomic and epigenomic data to understand a fundamental immunological process. In the Final Practice, you will apply similar methods (but using modern single-cell multi-modal data) to a related question.

---

## 12. Key Papers [DEEP DIVE]

Here are the foundational papers in multi-omics immunology. You do not need to read all of them, but knowing they exist and what they showed will help you put your work in context.

### CITE-seq and Multi-modal Methods

- **Stoeckius et al., Nature Methods, 2017** -- "Simultaneous epitope and transcriptome measurement in single cells." The original CITE-seq paper. Introduced the idea of DNA-barcoded antibodies for simultaneous RNA and protein measurement.

- **Hao et al., Cell, 2021** -- "Integrated analysis of multimodal single-cell data." Introduced Weighted Nearest Neighbors (WNN) in Seurat v4. Showed how to combine RNA, protein, and ATAC data in a principled way. This paper also released the PBMC CITE-seq reference dataset you may have used today.

- **Gayoso et al., Nature Methods, 2021** -- "Joint probabilistic modeling of single-cell multi-omic data with totalVI." The totalVI paper. Showed that probabilistic modeling of CITE-seq data outperforms simpler approaches.

### Multi-omics in Thymus and T Cell Development

- **Steier et al., Nature Immunology, 2023** -- "Single-cell multiomic analysis of thymocyte development reveals drivers of CD4/CD8 lineage commitment." Used CITE-seq to study thymocyte development, finding that CD4+ T cell differentiation is initiated first, with a second TCR signal determining CD8+ fate via the calcineurin-NFAT-GATA3 pathway.

- **Park et al., Science, 2020** -- "A cell atlas of human thymic development defines T cell repertoire formation." A comprehensive single-cell atlas of the human thymus across the lifespan, revealing new cell types, TCR recombination biases, and lineage commitment kinetics.

- **Chopp et al., Immunity, 2020** -- "An integrated epigenomic and transcriptomic map of mouse and human alpha-beta T cell development." Showed that chromatin changes precede and predict transcriptional changes during CD4/CD8 lineage commitment (discussed in detail in Section 11).

### Integration Methods

- **Stuart et al., Cell, 2019** -- "Comprehensive integration of single-cell data." Introduced anchor-based integration in Seurat v3 for combining datasets across modalities, batches, and technologies.

---

## 13. Self-Check

Before moving on, verify the following.

### Conceptual Questions

Answer these in your own words (write them in a markdown cell in your notebook):

1. **What does CITE-seq measure that scRNA-seq alone cannot?**
   Expected answer: CITE-seq measures surface protein levels alongside gene expression. scRNA-seq only measures mRNA, which is the instructions for making protein but not the protein itself.

2. **Why might RNA and protein levels disagree for the same gene?**
   Expected answer: Several reasons -- timing differences (RNA changes first, protein lags behind), post-transcriptional regulation (not all mRNA is translated), protein degradation rates differ from mRNA degradation rates, and proteins may be stored inside the cell rather than displayed on the surface.

3. **What is chromatin accessibility and why does it matter?**
   Expected answer: Chromatin accessibility refers to whether a region of DNA is physically open (unwound from histones) and available for the cell's machinery to read. It matters because it reveals the cell's *potential* -- what genes it could turn on -- rather than just what it is currently expressing. Chromatin changes often predict future cell state changes before RNA does.

4. **Name two integration methods and briefly describe how they work.**
   Expected answer: (a) totalVI -- a deep learning model that jointly models RNA and protein counts, learning a shared low-dimensional representation that captures both modalities. (b) WNN (Weighted Nearest Neighbors) -- assigns per-cell weights to each modality, so cells where protein is more informative weight protein more heavily, and vice versa. Other valid answers include MultiVI, label transfer, or anchor-based integration.

5. **What is the advantage of multi-modal data over single-modal data?**
   Expected answer: Different modalities capture different aspects of cell biology. RNA shows current transcription, protein confirms surface display, and ATAC reveals regulatory potential. By combining them, you get a more complete and accurate picture of each cell's state and trajectory. You can also resolve ambiguities that arise in single-modal analysis (e.g., two cell types that look similar in RNA may be clearly distinct in protein).

### Outputs to Verify

Make sure you have produced the following figures:

- [ ] RNA-only UMAP colored by cell type or cluster
- [ ] Protein-only UMAP colored by cell type or cluster
- [ ] Side-by-side comparison of RNA-only and protein-only UMAPs
- [ ] totalVI training history plot (loss decreasing over epochs)
- [ ] totalVI integrated UMAP colored by cell type or cluster
- [ ] RNA vs protein comparison for at least one marker (UMAP or scatter)

### Files You Should Have

Your `~/bootcamp/` folder should now contain:

```
~/bootcamp/
  week4_multiomics.ipynb         <-- This week's notebook
  rna_vs_protein_umap.png        <-- RNA vs protein UMAP comparison
  mudata_modalities.png          <-- MuData visualization
  totalvi_training.png           <-- Training history plot
  totalvi_integrated_umap.png    <-- Integrated UMAP
  rna_vs_protein_markers.png     <-- Marker comparison plots
  rna_protein_scatter.png        <-- RNA-protein scatter plot
  ... (files from previous weeks)
```

### Achievement Unlocked

You can now integrate multiple types of biological data from the same cells. You understand the multi-omics landscape, you have worked with CITE-seq data, you have used a deep learning model (totalVI) for integration, and you understand how chromatin accessibility adds another layer of information.

This is one of the most in-demand skills in modern computational biology. Researchers at top institutions are actively developing new methods for multi-omics integration, and the ability to work with these data types puts you at the frontier of the field.

---

## 14. What's Coming Next

You have now completed all four instructional weeks of the bootcamp:

| Week | What You Learned |
|------|-----------------|
| Week 0 | Set up your computing environment |
| Week 1 | Immunology foundations |
| Week 2 | Single-cell RNA-seq analysis (PBMC3k) |
| Week 3 | TCR/BCR repertoire analysis |
| **Week 4** | **Multi-omics integration (CITE-seq, scATAC concepts, totalVI)** |

Next up is the **[Final Practice](../final_practice/README.md)** -- a capstone project with two tracks where you apply everything you have learned:

- **Track 1: Recreate a "Chopp-style" result.** You will perform RNA + ATAC integration on a public dataset (like the 10x Multiome PBMC), including peak calling, gene activity estimation, TF motif analysis, and co-embedding. This puts into practice the scATAC-seq concepts from Section 5 and the integration methods from Section 6.

- **Track 2: TCR specificity prediction.** You will build machine learning models to predict whether a TCR binds a given epitope, using the repertoire analysis skills from Week 3 combined with feature engineering and model evaluation.

You can choose one track or attempt both. Either way, the Final Practice is your chance to produce a complete, publication-quality analysis from start to finish.

Take a break, review your notes, and come back ready to tie everything together.

---

## Common Problems and Solutions

| Problem | Solution |
|---------|----------|
| `ModuleNotFoundError: No module named 'muon'` | Run `pip install muon` in your terminal with the immunology environment activated |
| `ModuleNotFoundError: No module named 'scvi'` | Run `pip install scvi-tools` in your terminal |
| `scvi.data.pbmc_seurat_v4_cite_seq()` fails to download | Check your internet connection. If the server is down, use the fallback approach in Cell 2 (Alternative) |
| totalVI training is very slow | Reduce the number of cells: `adata_total = adata_total[:5000, :]`. On a laptop without GPU, expect 5-15 minutes |
| Training loss becomes NaN | Try `model.train(max_epochs=200, lr=1e-4)` to use a lower learning rate |
| Memory error during training | Close other applications, reduce data size, or try `model.train(batch_size=128)` for smaller batches |
| UMAP looks like a single blob | Check that your normalization worked correctly. Try different resolution values for Leiden clustering |
| Protein data not found in `obsm` | Print `list(adata.obsm.keys())` to see what keys are available. The protein data might be stored under a different key |
| CLR normalization gives unexpected values | Make sure you are applying CLR to the raw protein counts, not to already-normalized data |
| Mismatch in cell numbers between RNA and protein | Ensure both AnnData objects have the same `obs_names`. Use `adata_rna = adata_rna[common_cells]` to subset |

For more help, see the [Troubleshooting Guide](../TROUBLESHOOTING.md).
