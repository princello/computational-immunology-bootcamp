# Week 3: TCR and BCR Repertoire Analysis

**Time needed:** 8-10 hours
**Goal:** Analyze the immune receptor sequences of T cells and B cells to understand how the immune system recognizes specific threats.

By the end of this week, you will calculate diversity metrics, visualize clonal expansion, and map clonotypes onto gene expression clusters. These are the same analyses used by researchers developing cancer immunotherapies and next-generation vaccines.

**Prerequisites:** You should have completed Weeks 0-2. You can use Python and Scanpy, you understand the major immune cell types, and you can generate UMAPs from single-cell RNA-seq data. No prior knowledge of immune receptors is needed beyond knowing that T cells and B cells each have a unique receptor.

---

## Table of Contents

1. [What You Will Learn This Week](#1-what-you-will-learn-this-week)
2. [The Immune Receptor: Your Body's Lock-and-Key System](#2-the-immune-receptor-your-bodys-lock-and-key-system)
3. [VDJ Recombination: How Your Body Makes Billions of Unique Receptors](#3-vdj-recombination-how-your-body-makes-billions-of-unique-receptors)
4. [Clonotypes: Families of Identical Immune Cells](#4-clonotypes-families-of-identical-immune-cells)
5. [The Data: What TCR/BCR Sequencing Produces](#5-the-data-what-tcrbcr-sequencing-produces)
6. [Setting Up: Install and Import Tools](#6-setting-up-install-and-import-tools)
7. [Loading a Practice Dataset](#7-loading-a-practice-dataset)
8. [Exploring the TCR Data](#8-exploring-the-tcr-data)
9. [Clonotype Analysis](#9-clonotype-analysis)
10. [Diversity Analysis](#10-diversity-analysis)
11. [Mapping Clonotypes to Gene Expression](#11-mapping-clonotypes-to-gene-expression)
12. [Clonotype Overlap Between Groups](#12-challenge-clonotype-overlap-between-groups)
13. [Sequence Motif Analysis](#13-challenge-sequence-motif-analysis)
14. [Key Databases for TCR/BCR Research](#14-deep-dive-key-databases-for-tcrbcr-research)
15. [Self-Check](#15-self-check)
16. [Key Papers (Optional Reading)](#16-key-papers-optional-reading)
17. [What's Coming Next](#17-whats-coming-next)

---

## 1. What You Will Learn This Week

In Weeks 1 and 2, you learned to profile immune cells by their gene expression -- which genes are turned on in each cell. That told you *what kind* of cell it is (a T cell, a B cell, a monocyte). This week, you will go deeper. You will analyze the actual receptor sequences that make each T cell and B cell unique, answering a fundamentally different question: not just *what kind* of cell is this, but *what specific threat is this cell designed to fight?*

Here is what you will accomplish:

- **Understand immune receptor biology.** Learn what TCRs and BCRs are, how they are assembled from gene segments, and why each one is unique.
- **Load and explore TCR sequencing data.** Work with real paired single-cell RNA-seq and TCR-seq data from cancer patients.
- **Analyze clonal expansion.** Identify which T cell clones have multiplied, suggesting they recognized and responded to a specific antigen.
- **Calculate diversity metrics.** Quantify how varied the repertoire is using richness, Shannon entropy, and other measures.
- **Map clonotypes onto gene expression clusters.** Connect receptor identity to cell state -- see which cell types are expanding clonally.
- **Explore clonotype sharing and sequence patterns.** Investigate whether patients share clonotypes and what CDR3 sequences look like.

By the end, you will have produced diversity plots, clonal expansion visualizations, and clonotype-to-cluster maps -- the core analyses in any immune repertoire study.

---

## 2. The Immune Receptor: Your Body's Lock-and-Key System

### The Big Picture

In Week 1, you learned that your immune system has many cell types, each with a different job. Two of the most important are:

- **T cells**, which directly kill infected cells or coordinate immune responses
- **B cells**, which produce antibodies to neutralize pathogens

What makes T cells and B cells special among all immune cells is this: every single T cell carries a unique **T Cell Receptor (TCR)**, and every single B cell carries a unique **B Cell Receptor (BCR)**. These receptors sit on the cell surface and act like molecular antennae, scanning for foreign molecules called **antigens**.

Think of it as a lock-and-key system. Each antigen (from a virus, bacterium, or cancer cell) is a "lock." Each receptor is a "key." Your body manufactures billions of different keys, hoping that when a new lock shows up, at least one key will fit.

### Two Chains Make One Receptor

Each receptor is not a single molecule but a pair of two protein chains that work together:

- **TCR** = alpha chain + beta chain (written as TRA + TRB)
- **BCR** = heavy chain + light chain (written as IGH + IGL or IGK)

Both chains contribute to the receptor's ability to bind an antigen, but the most critical region is the **CDR3**.

### The CDR3 Region: The "Teeth" of the Key

The **Complementarity Determining Region 3**, or **CDR3**, is a short stretch of amino acids at the very tip of the receptor. It is the part that physically contacts the antigen. If the receptor is a key, the CDR3 is the teeth -- the precise shape that determines which lock it fits.

Key facts about CDR3:

- It is typically 10-20 amino acids long
- It is the most variable part of the receptor (no two people have the same CDR3 repertoire)
- It is the sequence most commonly used to identify and compare receptors in computational analyses
- When you see a "clonotype" defined in a dataset, it is almost always based on CDR3 identity

Here is an example of what a CDR3 amino acid sequence looks like:

```
CASSLAPGATNEKLFF
```

That string of letters represents a specific amino acid sequence. Each letter is one amino acid (C = cysteine, A = alanine, S = serine, and so on). This particular CDR3 belongs to a TCR beta chain, and it defines one unique receptor out of the billions your body can make.

### Why This Matters for Computational Analysis

When we sequence immune receptors from thousands of cells simultaneously, we get a table of CDR3 sequences. By analyzing these sequences, we can answer questions like:

- Has this person's immune system responded to a specific infection? (Look for expanded clonotypes.)
- How diverse is this patient's immune repertoire? (Calculate diversity metrics.)
- Do different patients respond to the same pathogen with similar receptors? (Look for shared clonotypes.)
- Which cell types are actively fighting? (Map clonotypes to gene expression clusters.)

These are exactly the analyses you will perform in the hands-on sections below.

---

## 3. VDJ Recombination: How Your Body Makes Billions of Unique Receptors

### The Problem

Your genome contains roughly 20,000 genes. But you need *billions* of unique receptors -- far more than the number of genes you have. How does your body pull this off?

### The Solution: Mix and Match Gene Segments

The answer is **VDJ recombination**, one of the most elegant solutions in all of biology. Instead of encoding each receptor as a single gene, your DNA encodes a set of interchangeable *segments* that get shuffled and combined randomly in each developing immune cell.

There are three types of segments:

| Segment | Name | Approximate Count | Role |
|---------|------|-------------------|------|
| **V** | Variable | ~50 options | The main body of the receptor's antigen-binding region |
| **D** | Diversity | ~25 options | An extra piece that adds diversity (heavy/beta chains only) |
| **J** | Joining | ~6 options | Connects the variable region to the rest of the receptor |

During development, each immune cell randomly picks one V, one D (for heavy/beta chains), and one J segment, then physically cuts and pastes its own DNA to join them together. This is an irreversible process -- once a cell has made its choice, it is locked in for life, and all its descendants will carry the same rearrangement.

### The Math of Combinatorial Diversity

Even with just the combinatorial math, the numbers are impressive:

```
TCR beta chain: ~50 V x 2 D x 13 J = ~1,300 combinations
TCR alpha chain: ~50 V x 50 J = ~2,500 combinations (no D segment)
Together: 1,300 x 2,500 = ~3.25 million TCR combinations
```

But wait -- that is only 3.25 million. We said we need billions. Where does the rest of the diversity come from?

### Junctional Diversity: The Secret Ingredient

When the V, D, and J segments are joined, the cell's molecular machinery does something remarkable: it randomly *deletes* a few nucleotides from the ends of each segment and *inserts* new random nucleotides at the junctions. This process is called **junctional diversity**, and it happens right at the CDR3 region.

This randomness is what makes the CDR3 so variable. It is also why the CDR3 is the "fingerprint" of each receptor: even two cells that chose the same V, D, and J segments will almost certainly have different CDR3 sequences because of these random insertions and deletions.

### The Analogy: A Combination Lock with Extra Randomness

Think of it like a combination lock with three dials (V, D, J). Each dial has many positions. That alone gives you thousands of combinations. But then, after you set the dials, someone randomly files down a few teeth on the lock and welds on a few new ones in unpredictable places. Now the number of possible locks is astronomical.

### The Final Numbers

| Receptor | Estimated Diversity |
|----------|-------------------|
| TCR (alpha + beta combined) | ~10^15 (one quadrillion) |
| BCR (heavy + light combined) | ~10^18 (one quintillion) |

For context, there are roughly 10^11 (100 billion) stars in the Milky Way. Your immune system can generate more unique receptors than there are stars in thousands of galaxies.

### What You Will See in the Data

When you load TCR sequencing data, each cell's receptor will be described by:

- **v_gene**: which V segment was chosen (e.g., TRBV5-1)
- **d_gene**: which D segment was chosen (e.g., TRBD1)
- **j_gene**: which J segment was chosen (e.g., TRBJ2-7)
- **cdr3**: the amino acid sequence of the CDR3 region
- **cdr3_nt**: the nucleotide sequence of the CDR3 region

Understanding VDJ recombination helps you interpret these columns. When two cells share the same V, D, J genes *and* the same CDR3 sequence, they almost certainly descended from the same original cell -- they are a **clonotype**.

---

## 4. Clonotypes: Families of Identical Immune Cells

### What Is a Clonotype?

A **clonotype** is a group of cells that all share the same immune receptor sequence. In practice, this usually means they share identical CDR3 sequences (and often the same V, D, J gene usage).

Why would multiple cells have the exact same receptor? Because of **clonal expansion**.

### Clonal Expansion: The Immune System's Response

Here is what happens when your body detects a threat:

1. Among your billions of T cells and B cells, each with a unique receptor, a few happen to have receptors that bind the invading antigen.
2. Those cells become activated -- they receive a signal that says "you found a match!"
3. The activated cells begin to divide rapidly, making many copies of themselves. Each copy carries the same receptor.
4. This army of identical cells then goes to work: T cells kill infected cells, B cells produce antibodies.

This process is called **clonal expansion**, and it is one of the most important concepts in immunology. When we see a clonotype with many cells in our sequencing data, it tells us that those cells recognized something and responded.

### What Clonal Expansion Looks Like in Data

| Clonotype ID | Number of Cells | Interpretation |
|-------------|----------------|----------------|
| clone_001 | 1 | Not expanded -- just one cell with this receptor |
| clone_002 | 1 | Not expanded |
| clone_003 | 47 | Highly expanded -- this receptor recognized something |
| clone_004 | 3 | Slightly expanded |
| clone_005 | 1 | Not expanded |

In a healthy person with no active infection, most clonotypes will appear only once or a few times. During an infection, cancer response, or autoimmune flare, you will see a few clonotypes balloon to large numbers.

### Clonal Diversity

**Clonal diversity** describes how varied the repertoire is overall:

- **High diversity** = many unique clonotypes, each appearing once or a few times. This is what you expect in a healthy, resting immune system. It means the body has a broad collection of receptors ready to respond to many different threats.
- **Low diversity** = a few clonotypes dominate the repertoire. This can indicate an active immune response (the expanded clones are crowding out the rare ones), a chronic infection, or immune exhaustion.

Measuring diversity is important in clinical settings. For example, cancer patients receiving checkpoint immunotherapy tend to have better outcomes when their TCR repertoire is more diverse, because it means their immune system has more "keys" available to fight the tumor.

---

## 5. The Data: What TCR/BCR Sequencing Produces

### How the Data Is Generated

Modern single-cell technologies, particularly those from 10x Genomics, can simultaneously capture two types of information from the same cell:

1. **Gene expression (scRNA-seq)**: The full transcriptome -- which genes are on and how active they are. This is what you analyzed in Week 2.
2. **Immune receptor sequences (scTCR-seq or scBCR-seq)**: The exact VDJ rearrangement and CDR3 sequence of the cell's receptor.

Because both measurements come from the same cell (linked by a shared cell barcode), you can do something powerful: overlay receptor information onto gene expression clusters. You can see not just that a cell is a CD8+ T cell, but which specific antigen receptor that CD8+ T cell carries and whether it has expanded.

### What the Output Table Looks Like

The receptor sequencing output is a table where each row is one chain from one cell. Here are the key columns:

| Column | Description | Example |
|--------|-------------|---------|
| `barcode` | Unique identifier for each cell (links to scRNA-seq data) | `AAACCTGCAGTTCATG-1` |
| `chain` | Which chain was sequenced | `TRA`, `TRB`, `IGH`, `IGL`, `IGK` |
| `v_gene` | The V segment used | `TRBV5-1` |
| `d_gene` | The D segment used (only for TRB/IGH) | `TRBD1` |
| `j_gene` | The J segment used | `TRBJ2-7` |
| `c_gene` | The constant region gene | `TRBC1` |
| `cdr3` | CDR3 amino acid sequence | `CASSLAPGATNEKLFF` |
| `cdr3_nt` | CDR3 nucleotide sequence | `TGTGCCAGCAGCTTG...` |
| `clonotype_id` | Group label for cells with identical receptors | `clonotype1` |

A few things to notice:

- Each cell typically has **two rows** (one for each chain: alpha + beta for TCR, heavy + light for BCR)
- Some cells may be missing one chain (technical dropout) or have extra chains (dual-receptor cells)
- The `barcode` column is what links the receptor data to the gene expression data

### Chain Naming Conventions

| Receptor | Chain 1 | Chain 2 | Notes |
|----------|---------|---------|-------|
| TCR | TRA (alpha) | TRB (beta) | Most T cells have one alpha and one beta |
| BCR | IGH (heavy) | IGL (lambda) or IGK (kappa) | B cells use either lambda or kappa as the light chain, never both |

You will see these abbreviations throughout the analysis. "VDJ" in the context of data columns usually refers to the heavy chain (IGH) or beta chain (TRB), because those are the chains that use all three segments (V, D, and J). The alpha chain (TRA) and light chains (IGL/IGK) only use V and J segments (no D segment).

---

## 6. Setting Up: Install and Import Tools [REQUIRED]

Time to get your hands on the keyboard. Open a Jupyter notebook (if you need a refresher on how, revisit [Week 0, Section 6](../week_0_setup/README.md#6-install-jupyter-notebook)).

### Install Scirpy

If you followed the Week 0 setup, Scirpy should already be installed. If not, run this in your terminal before opening Jupyter:

```bash
# Make sure your immunology environment is active
conda activate immunology

# Install scirpy (the main tool for immune receptor analysis)
pip install scirpy
```

### Import All Required Packages

In the first cell of your notebook, type and run the following:

```python
# ============================================================
# Week 3: TCR and BCR Repertoire Analysis
# ============================================================

# scanpy: the single-cell analysis framework you already know from Week 2
import scanpy as sc

# scirpy: a scverse package for immune receptor (TCR/BCR) analysis
# It integrates tightly with scanpy -- same data structures, same patterns
import scirpy as ir

# matplotlib: for creating plots and figures
import matplotlib.pyplot as plt

# pandas: for working with tabular data (DataFrames)
import pandas as pd

# suppress warning messages to keep the output clean
import warnings
warnings.filterwarnings('ignore')

# set the visual style for our figures
sc.settings.set_figure_params(dpi=80, facecolor='white')

# print version numbers to verify everything is installed correctly
print(f"Scanpy version: {sc.__version__}")
print(f"Scirpy version: {ir.__version__}")
```

**Expected output:**

```
Scanpy version: 1.10.x   (or similar)
Scirpy version: 0.18.x   (or similar)
```

The exact version numbers may differ, and that is fine. What matters is that both lines print without errors. If you get an `ImportError` or `ModuleNotFoundError`, go back to your terminal and install the missing package with `pip install scanpy` or `pip install scirpy`.

### A Note About Scirpy

Scirpy (pronounced "SKIR-py") is part of the **scverse** ecosystem -- the same family of tools that includes Scanpy, which you used in Week 2. It follows the same conventions:

- Data is stored in AnnData objects (`adata`)
- Preprocessing functions are in `ir.pp`
- Analysis/tool functions are in `ir.tl`
- Plotting functions are in `ir.pl`

If you are comfortable with Scanpy's `sc.pp`, `sc.tl`, and `sc.pl` pattern, Scirpy will feel familiar.

---

## 7. Loading a Practice Dataset [REQUIRED]

### The Dataset: Breast Cancer T Cells

Scirpy comes with a built-in tutorial dataset from **Wu et al., 2020** -- a study of T cells from breast cancer patients. This dataset contains both gene expression (RNA) and TCR sequence data from the same cells, which is exactly what we need.

The researchers sequenced T cells from multiple patients and tissue sites (tumor, adjacent normal tissue, blood), giving us a rich dataset to explore clonal expansion, diversity, and clonotype sharing.

### Load the Data

```python
# Load the built-in tutorial dataset from Wu et al., 2020
# This contains T cells from breast cancer patients with both
# gene expression (RNA) and TCR sequence data already integrated
adata = ir.datasets.wu2020()

# Print the basic structure of the data
print(adata)
print(f"\nNumber of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
```

**Expected output:**

```
AnnData object with n_obs x n_vars = XXXX x XXXX
    obs: 'cluster', ...
    var: ...
    obsm: 'X_umap', ...
    ...

Number of cells: XXXX
Number of genes: XXXX
```

The exact numbers will depend on the version of scirpy you have installed, but you should see an AnnData object with several thousand cells and several thousand genes.

### Understand the Data Structure

This AnnData object is the same type you worked with in Week 2, but it now has extra columns in `adata.obs` that store TCR information. Let's explore:

```python
# Look at the first few rows of the cell-level metadata
# This is where both the gene expression annotations (like cluster)
# and the TCR information are stored
print("Shape of metadata table:", adata.obs.shape)
print("\nAll columns in adata.obs:")
for i, col in enumerate(adata.obs.columns):
    print(f"  {i+1}. {col}")
```

**Expected output:**

A numbered list of column names. You will see familiar ones like `cluster` from gene expression analysis, and new ones that start with `IR_` -- those are the immune receptor columns added by scirpy.

```python
# Look at a few example cells and their metadata
# Selecting only the most interesting columns for readability
interesting_cols = [col for col in adata.obs.columns if col in [
    'cluster', 'patient', 'source'
] or col.startswith('IR_VDJ_1') or col.startswith('IR_VJ_1')]

# Show the first 5 cells
print(adata.obs[interesting_cols].head())
```

**Expected output:**

A table showing, for each cell, its cluster assignment, patient ID, tissue source, and the V/D/J gene segments and CDR3 sequences for its TCR chains. Some cells may have `NaN` values for the receptor columns, meaning no TCR was detected for that cell (this is normal -- the capture efficiency is not 100%).

### What the Column Prefixes Mean

Scirpy uses a specific naming convention for receptor data columns:

| Prefix | Meaning |
|--------|---------|
| `IR_VDJ_1_` | First VDJ chain (typically TRB for T cells, IGH for B cells) |
| `IR_VDJ_2_` | Second VDJ chain (if the cell has two beta/heavy chains) |
| `IR_VJ_1_` | First VJ chain (typically TRA for T cells, IGL/IGK for B cells) |
| `IR_VJ_2_` | Second VJ chain (if the cell has two alpha/light chains) |

Remember: "VDJ" chains use V, D, and J segments (beta and heavy chains), while "VJ" chains only use V and J segments (alpha and light chains). This naming convention helps you immediately know which chain type you are looking at.

> **[CHECKPOINT]** You have loaded the dataset successfully if `adata.n_obs` returns a number greater than zero and you can see `IR_` columns in `adata.obs`. Take a moment to look at the data. Can you find the CDR3 sequences? Can you tell which cells come from which patients? This is the raw material for everything that follows.

---

## 8. Exploring the TCR Data [REQUIRED]

Now that the data is loaded, let's explore the TCR information in detail.

### View the Receptor Columns

```python
# Scirpy stores TCR/BCR info in adata.obs with column names
# that start with "IR_". Let's see all of them.
print("Columns related to immune receptors:")
ir_cols = [c for c in adata.obs.columns if c.startswith('IR_')]
for col in ir_cols[:15]:
    # Show the column name and an example non-null value
    example = adata.obs[col].dropna().iloc[0] if adata.obs[col].notna().any() else "all NaN"
    print(f"  {col}: e.g., {example}")
```

**Expected output:**

A list of columns with example values like gene names (e.g., `TRBV5-1`), CDR3 sequences (e.g., `CASSLAPGATNEKLFF`), and chain types (e.g., `TRB`).

### Chain Pairing Quality

A critical quality check for TCR data is **chain pairing**: does each cell have both an alpha and a beta chain? In a perfect world, every T cell would have exactly one alpha chain (TRA) and one beta chain (TRB). In reality, some cells are missing one chain (technical dropout), some have extra chains (dual TCR cells), and some have no TCR detected at all.

```python
# Run the chain quality control analysis
# This examines each cell and categorizes its chain pairing status
ir.tl.chain_qc(adata)

# See the results
print("Chain pairing categories:")
print(adata.obs['chain_pairing'].value_counts())
```

**Expected output:**

```
chain_pairing
single pair       XXXX
orphan VDJ        XXXX
orphan VJ         XXXX
extra VDJ         XXXX
extra VJ          XXXX
two full chains   XXXX
no IR             XXXX
...
```

Here is what each category means:

| Category | Meaning | Is This Good? |
|----------|---------|---------------|
| `single pair` | One alpha + one beta chain detected | Ideal -- this is what we want |
| `orphan VDJ` | Only a beta chain detected, no alpha | Acceptable -- alpha chains have higher dropout |
| `orphan VJ` | Only an alpha chain detected, no beta | Acceptable but less common |
| `extra VDJ` | More than one beta chain detected | Unusual -- could be a doublet or dual-receptor cell |
| `extra VJ` | More than one alpha chain detected | Fairly common -- some T cells genuinely express two alpha chains |
| `two full chains` | Two complete alpha-beta pairs | Could be doublet (two cells in one droplet) |
| `no IR` | No receptor detected | Normal for non-T/B cells or technical dropout |
| `multichain` | Complex multi-chain pattern | Rare -- usually a technical artifact |

### Visualize Chain Pairing

```python
# Create a bar chart of chain pairing categories
fig, ax = plt.subplots(figsize=(8, 4))

# Count the cells in each category and plot them
pairing_counts = adata.obs['chain_pairing'].value_counts()
pairing_counts.plot(kind='barh', ax=ax, color='steelblue', edgecolor='white')

# Label the axes and title
ax.set_xlabel('Number of Cells')
ax.set_ylabel('Chain Pairing Category')
ax.set_title('TCR Chain Pairing Quality')

# Add count labels to the end of each bar for clarity
for i, (count, name) in enumerate(zip(pairing_counts.values, pairing_counts.index)):
    ax.text(count + 10, i, str(count), va='center', fontsize=9)

plt.tight_layout()
plt.savefig('chain_pairing.png', dpi=150)
plt.show()
print("Saved: chain_pairing.png")
```

**Expected output:**

A horizontal bar chart showing the distribution of chain pairing categories. The `single pair` bar should be the longest -- most cells should have a clean alpha+beta pair. You should see some `orphan` categories (where only one chain was captured) and possibly some cells with extra chains.

### What If Most Cells Have No IR?

If you see a large number of cells in the `no IR` category, do not worry. This can happen for two reasons:

1. The dataset includes non-T cells (monocytes, B cells, etc.) that do not have TCRs
2. Technical dropout -- some cells' TCR transcripts were not captured

For downstream analysis, scirpy will automatically focus on cells that have receptor information.

> **Common Mistake:** Do not confuse chain pairing quality with data quality overall. Having 20-30% of cells in the `orphan` categories is normal for 10x Genomics data. The alpha chain is particularly prone to dropout because of its genomic structure.

> **[CHECKPOINT]** You should have a `chain_pairing.png` file saved and be able to answer: What percentage of your cells have a clean "single pair" of alpha and beta chains? If the majority of cells with TCR data fall into `single pair` or `orphan` categories, your data is in good shape. Continue to the next section.

---

## 9. Clonotype Analysis [REQUIRED]

This is the core of immune repertoire analysis. You are about to identify clonotypes -- groups of cells that share identical receptors -- and measure how much each one has expanded.

### Step 1: Define Clonotypes

To group cells into clonotypes, we need to decide what "identical" means. The most common approach is to require exact identity of the CDR3 nucleotide sequence. Two cells are in the same clonotype if and only if their CDR3 nucleotide sequences match perfectly.

```python
# Step 1: Compute pairwise distances between receptor sequences
# metric='identity' means we require exact matches (distance = 0)
# sequence='nt' means we compare nucleotide sequences (more stringent than amino acid)
ir.pp.ir_dist(adata, metric='identity', sequence='nt')

# Step 2: Group cells into clonotypes based on those distances
# Cells with identical CDR3 nucleotide sequences get the same clone_id
ir.tl.define_clonotypes(adata)

# How many unique clonotypes did we find?
print(f"Number of unique clonotypes: {adata.obs['clone_id'].nunique()}")
print(f"Number of cells with clonotype info: {adata.obs['clone_id'].notna().sum()}")
```

**Expected output:**

```
Number of unique clonotypes: XXXX
Number of cells with clonotype info: XXXX
```

You should see a large number of unique clonotypes. Most will be "singletons" (appearing only once), but some will be expanded.

### Step 2: See the Most Expanded Clonotypes

```python
# Which clonotypes have the most cells?
# This tells us which receptors have undergone the most clonal expansion
clone_counts = adata.obs['clone_id'].value_counts()

print("Top 10 most expanded clonotypes:")
print(clone_counts.head(10))
print(f"\nTotal unique clonotypes: {len(clone_counts)}")
print(f"Singletons (appear only once): {(clone_counts == 1).sum()}")
print(f"Expanded (appear 2+ times): {(clone_counts >= 2).sum()}")
```

**Expected output:**

A table showing the top 10 clonotype IDs and their cell counts. You should see that the most expanded clonotype has significantly more cells than most others. The majority of clonotypes should be singletons.

Take a moment to think about what this means biologically. That top clonotype -- with dozens of cells -- represents a T cell that recognized some antigen (perhaps a tumor antigen in this breast cancer dataset) and rapidly divided to fight it. The singletons are T cells whose receptors have not yet encountered their matching antigen, or that responded only weakly.

### Step 3: Categorize Clonal Expansion

Scirpy can automatically categorize each cell's clonotype by its expansion level:

```python
# Categorize clonotypes by expansion level
# This adds a 'clonal_expansion' column to adata.obs
ir.tl.clonal_expansion(adata)

# See the categories and their counts
print("Clonal expansion categories:")
print(adata.obs['clonal_expansion'].value_counts())
```

**Expected output:**

```
clonal_expansion
1                  XXXX    (singletons -- not expanded)
2                  XXXX    (two cells with this clonotype)
3-10               XXXX    (moderately expanded)
>10                XXXX    (highly expanded)
...
```

The categories represent how many cells share each clonotype. Most cells will be singletons (clonal expansion = 1), reflecting the enormous diversity of the naive repertoire. The cells in the `>10` category are the most interesting -- they represent strong immune responses.

### Step 4: Visualize Clonal Expansion by Cell Cluster

Now the real power of paired data: we can see which *types* of T cells are expanding the most.

```python
# Plot clonal expansion broken down by cell cluster
# normalize=True shows proportions instead of raw counts
# clip_at=4 groups all clonotypes with 4+ cells together
ir.pl.clonal_expansion(
    adata,
    groupby='cluster',    # group by the cell type clusters from gene expression
    clip_at=4,            # group together clonotypes with 4 or more cells
    normalize=True        # show proportions so clusters of different sizes are comparable
)
plt.savefig('clonal_expansion.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved: clonal_expansion.png")
```

**Expected output:**

A stacked bar chart where each bar represents a cell cluster, and the colors represent clonal expansion categories. Clusters enriched for expanded clonotypes (more warm/dark colors) are the ones actively involved in immune responses.

Look for patterns:

- **CD8+ effector clusters** often show high clonal expansion (they are the "killers" actively fighting the tumor)
- **Naive T cell clusters** should have very little expansion (they have not yet encountered their antigen)
- **Regulatory T cell (Treg) clusters** may show moderate expansion

If you see these patterns, you are seeing real biology in the data. The computational analysis is revealing the immune response in action.

> **Common Mistake:** If you see an error about `clone_id` not existing, make sure you ran *both* `ir.pp.ir_dist()` and `ir.tl.define_clonotypes()` before this step. The order matters: distance computation must come before clonotype definition.

> **[CHECKPOINT]** You should have a `clonal_expansion.png` file showing clonal expansion by cluster. Ask yourself: which clusters have the most expanded clonotypes? Which have the least? Does this make biological sense given what you know about T cell subtypes from Week 1? Write down your observations -- you will revisit them in the self-check.

---

## 10. Diversity Analysis [REQUIRED]

### What Is Repertoire Diversity?

In ecology, scientists measure the "diversity" of species in an ecosystem. In immunology, we do the same thing, but instead of species, we measure the diversity of clonotypes in a repertoire. The concept is identical -- and in fact, we use the same mathematical formulas.

Diversity is important because it tells you about the *state* of the immune system:

- A diverse repertoire means the immune system has many different receptors ready to respond to many different threats
- A low-diversity repertoire means a few clonotypes dominate, which could indicate an ongoing immune response, chronic infection, or immune dysfunction

### Three Key Diversity Metrics

Let's walk through three commonly used metrics. We will use an analogy to make them intuitive.

Imagine you are standing in a forest and counting tree species.

#### 1. Richness

**What it measures:** The total number of unique clonotypes (species).

**Forest analogy:** You count 15 different tree species in the forest. The richness is 15.

**In immunology:** If a patient has 5,000 unique clonotypes, the richness is 5,000.

**Limitation:** Richness does not care about how evenly distributed the clonotypes are. A forest with 1,000 oaks and 1 maple has richness = 2, the same as a forest with 500 oaks and 500 maples.

#### 2. Shannon Entropy (Shannon Diversity Index)

**What it measures:** Both the number of unique clonotypes AND how evenly they are distributed.

**Forest analogy:**
- Forest A: 100 trees total, 50 oaks and 50 maples. Shannon entropy is high (two species, perfectly even).
- Forest B: 100 trees total, 99 oaks and 1 maple. Shannon entropy is low (two species, but extremely uneven -- practically a monoculture).

**In immunology:** A repertoire with 1,000 unique clonotypes, each appearing once, has high Shannon entropy. A repertoire where one clonotype has 990 cells and 10 others have one cell each has low Shannon entropy, even though it has 11 unique clonotypes.

**Why it matters:** Shannon entropy is the most commonly reported diversity metric in repertoire studies because it captures both richness and evenness in a single number.

#### 3. D50

**What it measures:** The minimum number of the most abundant clonotypes needed to account for 50% of all cells.

**In immunology:** If a patient's top 5 clonotypes (out of 1,000 total) make up 50% of all cells, then D50 = 5. This means the repertoire is highly skewed -- a few clonotypes dominate.

**Why it matters:** D50 is a quick way to assess how "top-heavy" a repertoire is. A low D50 means high clonal expansion; a high D50 means a more balanced repertoire.

### Calculate Diversity

```python
# Calculate Shannon diversity for each cluster
# This tells us which cell types have the most diverse repertoires
ir.tl.alpha_diversity(adata, groupby='cluster', metric='shannon')

# Print the diversity values
print("Shannon diversity by cluster:")
diversity_data = adata.obs.groupby('cluster')['alpha_diversity_shannon'].first()
print(diversity_data.sort_values(ascending=False))
```

**Expected output:**

A table showing Shannon diversity values for each cluster. Higher values mean more diverse repertoires. You should see that naive T cell clusters have higher diversity (many unique clonotypes, each appearing once) while effector/memory clusters have lower diversity (dominated by expanded clonotypes).

### Visualize Diversity

```python
# Create a violin plot of Shannon diversity by cluster
fig, ax = plt.subplots(figsize=(10, 5))
sc.pl.violin(
    adata,
    keys='alpha_diversity_shannon',    # the metric we just calculated
    groupby='cluster',                  # one violin per cluster
    ax=ax,
    rotation=45                         # rotate x-axis labels for readability
)
ax.set_title('TCR Repertoire Diversity by Cluster (Shannon Entropy)')
ax.set_ylabel('Shannon Diversity')
plt.tight_layout()
plt.savefig('diversity.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved: diversity.png")
```

**Expected output:**

A violin plot with one violin per cluster. The height of each violin represents the Shannon diversity. Look for the pattern: clusters with high diversity (tall violins, high values) are likely naive or resting T cells, while clusters with low diversity (short violins, low values) are likely effector cells in the middle of an immune response.

### Interpreting the Results

Here is a guide for interpreting diversity patterns:

| Pattern | What It Suggests |
|---------|-----------------|
| High diversity, low expansion | Naive or resting T cells -- broad repertoire, no dominant clones |
| Low diversity, high expansion | Active immune response -- a few clones are dominating |
| Low diversity across all clusters | Could indicate immune exhaustion, chronic infection, or aging |
| Different diversity between patients | Individual variation in immune history and genetics |

> **[CHECKPOINT]** You should have a `diversity.png` file and be able to explain what Shannon diversity measures. Here is a test: if Patient A has a Shannon diversity of 8.5 and Patient B has a Shannon diversity of 3.2, which patient has a more diverse TCR repertoire? (Answer: Patient A. Higher Shannon entropy means more diversity.) If you got that right, continue.

---

## 11. Mapping Clonotypes to Gene Expression [REQUIRED]

This section is where everything comes together. By overlaying clonotype information on the gene expression UMAP, you can see *which cell types* are expanding -- connecting receptor biology to cell identity.

### UMAP with Clonal Expansion Overlay

```python
# Create side-by-side UMAPs: one colored by cell cluster,
# one colored by clonal expansion status
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left panel: cell clusters (from gene expression analysis)
sc.pl.umap(
    adata,
    color='cluster',          # color by cell type cluster
    frameon=False,            # remove the box around the plot
    title='Cell Clusters',
    ax=axes[0],
    show=False                # do not display yet -- we will show both panels together
)

# Right panel: clonal expansion status (from TCR analysis)
sc.pl.umap(
    adata,
    color='clonal_expansion',    # color by expansion category
    frameon=False,
    title='Clonal Expansion',
    ax=axes[1],
    show=False
)

plt.tight_layout()
plt.savefig('umap_clonal.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved: umap_clonal.png")
```

**Expected output:**

Two UMAP plots side by side. The left one should look familiar -- colored by cell type, just like the UMAPs you made in Week 2. The right one is new: it shows the same cells, but colored by how expanded their clonotype is.

Look at the right panel carefully. You should see that clonal expansion is not evenly distributed across the UMAP. Some clusters are "hot" with expanded clonotypes (concentrated color), while others are almost entirely singletons.

### Cross-Tabulation: Which Cell Types Are Expanding?

Let's quantify what the UMAP shows:

```python
# Create a cross-tabulation table
# Rows: cell clusters
# Columns: clonal expansion categories
# Values: proportion of cells in each category
ct = pd.crosstab(
    adata.obs['cluster'],              # rows
    adata.obs['clonal_expansion'],     # columns
    normalize='index'                  # normalize by row (each row sums to 1.0)
)

# Display the table, rounded for readability
print("Proportion of cells in each expansion category, by cluster:")
print(ct.round(3))
```

**Expected output:**

A table where each row is a cell cluster and each column is an expansion category. The values show what fraction of each cluster's cells fall into each expansion category. For example:

```
clonal_expansion         1      2   3-10    >10
cluster
CD8_effector         0.320  0.150  0.280  0.250
CD8_naive            0.910  0.060  0.020  0.010
Treg                 0.650  0.120  0.150  0.080
...
```

(Your exact numbers and cluster names will differ.)

This table tells a clear story: CD8 effector cells are highly expanded (only ~32% are singletons), while CD8 naive cells are almost entirely singletons (~91%). This makes perfect biological sense:

- **Naive cells** have not yet encountered their antigen, so they have not expanded
- **Effector cells** have been activated by an antigen and have divided many times
- **Memory cells** often show moderate expansion from past immune responses
- **Regulatory T cells** can show variable expansion depending on the disease context

### Highlight the Top Clonotypes on UMAP

For an even more specific view, you can highlight where the most expanded clonotypes appear:

```python
# Find the single most expanded clonotype
top_clone = adata.obs['clone_id'].value_counts().index[0]
top_clone_size = adata.obs['clone_id'].value_counts().iloc[0]
print(f"Most expanded clonotype: {top_clone} ({top_clone_size} cells)")

# Mark cells belonging to this clonotype
adata.obs['is_top_clone'] = (adata.obs['clone_id'] == top_clone).astype(str)

# Plot
fig, ax = plt.subplots(figsize=(7, 5))
sc.pl.umap(
    adata,
    color='is_top_clone',
    frameon=False,
    title=f'Top Clonotype: {top_clone} ({top_clone_size} cells)',
    palette={'True': 'red', 'False': 'lightgrey'},
    ax=ax,
    show=False
)
plt.tight_layout()
plt.savefig('umap_top_clone.png', dpi=150, bbox_inches='tight')
plt.show()
print("Saved: umap_top_clone.png")
```

**Expected output:**

A UMAP where most cells are light grey, and the cells belonging to the most expanded clonotype are highlighted in red. Notice where the red dots cluster -- they should be concentrated in a specific region of the UMAP, typically within an effector or memory cluster. This tells you the cell-type identity of that highly expanded clone.

> **Common Mistake:** If your UMAP looks odd (all one color, or no colored dots), double-check that the `clone_id` column exists and has non-null values. You may need to re-run the clonotype definition step (Section 9) if something went wrong.

> **[CHECKPOINT]** You can now connect receptor sequences to cell identities. This is a powerful capability. You should have `umap_clonal.png` and `umap_top_clone.png` saved, and you should be able to answer: which cell clusters in your dataset show the most clonal expansion? Why does this make biological sense? If you can answer these questions, you have mastered the core of repertoire analysis. The remaining sections are optional enrichment.

---

## 12. [CHALLENGE] Clonotype Overlap Between Groups

An important question in immunology is whether different patients, tissues, or conditions share clonotypes. If two patients independently develop T cells with the same receptor, it suggests they are responding to the same antigen -- this is especially relevant in cancer and infectious disease research.

### Check Available Grouping Variables

```python
# See what grouping variables are available in our dataset
# Common ones: patient, tissue source, treatment group
for col in ['patient', 'source', 'sample', 'group']:
    if col in adata.obs.columns:
        print(f"\n{col}:")
        print(adata.obs[col].value_counts())
```

### Compute and Visualize Overlap

```python
# Calculate clonotype overlap between groups
# The Jaccard index measures the proportion of shared clonotypes:
# J = (shared clonotypes) / (total unique clonotypes in both groups)
# A value of 0 means no sharing; 1 means identical repertoires

# First, check which grouping variable to use
group_col = None
for col in ['patient', 'source', 'sample', 'group']:
    if col in adata.obs.columns:
        group_col = col
        break

if group_col is not None:
    print(f"Computing clonotype overlap by: {group_col}")

    # Compute overlap
    ir.tl.clonotype_overlap(adata, groupby=group_col, metric='jaccard')

    # Visualize as a heatmap
    ir.pl.clonotype_overlap(adata, groupby=group_col)
    plt.savefig('clonotype_overlap.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("Saved: clonotype_overlap.png")
else:
    print("No suitable grouping variable found in this dataset.")
    print("Available columns:", list(adata.obs.columns))
```

**Expected output:**

A heatmap showing the Jaccard overlap between each pair of groups. If using `patient` as the grouping variable, each row and column represents a patient. The diagonal will be 1.0 (each patient's repertoire is identical to itself). Off-diagonal values near 0 indicate little clonotype sharing between patients.

### Interpreting Clonotype Overlap

| Overlap Value | Interpretation |
|--------------|----------------|
| 0.0 | No shared clonotypes between the two groups |
| 0.0-0.05 | Very low overlap (typical between unrelated patients) |
| 0.05-0.20 | Moderate overlap (could indicate shared antigen targets) |
| > 0.20 | High overlap (unusual -- suggests convergent immune response or related samples) |

In most cancer studies, overlap between different patients is very low. This reflects the fact that each patient's tumor has unique mutations (neoantigens) that drive unique T cell responses. However, some shared antigens (like common viral epitopes) can produce shared clonotypes.

If your dataset has a `source` column (e.g., tumor vs. blood), you might see higher overlap between different tissue sites from the same patient than between different patients -- expanded clonotypes often circulate between tissues.

---

## 13. [CHALLENGE] Sequence Motif Analysis

Beyond just counting clonotypes, we can look at the actual amino acid sequences of CDR3 regions to understand patterns in how receptors are built.

### Extract and Explore CDR3 Sequences

```python
# Get the CDR3 amino acid sequences from the beta chain (VDJ_1)
# The beta chain CDR3 is the most commonly analyzed because it has
# the most diversity (V, D, and J segments + junctional diversity)
cdr3_col = 'IR_VDJ_1_junction_aa'
cdr3_seqs = adata.obs[cdr3_col].dropna()

print(f"Number of cells with CDR3 beta sequences: {len(cdr3_seqs)}")
print(f"Number of unique CDR3 beta sequences: {cdr3_seqs.nunique()}")
print(f"\nExample CDR3 sequences (first 10):")
for i, seq in enumerate(cdr3_seqs.unique()[:10]):
    print(f"  {i+1}. {seq} (length: {len(seq)})")
```

**Expected output:**

A list of 10 CDR3 amino acid sequences, each a short string of capital letters. Typical lengths are 10-20 characters. Notice how they all start with C (cysteine) and usually end with F (phenylalanine) or W (tryptophan) -- this is a conserved feature of CDR3 regions across almost all TCRs.

### CDR3 Length Distribution

The length of the CDR3 region affects which antigens a receptor can bind. Let's visualize the distribution:

```python
# Calculate the length of each CDR3 sequence
cdr3_lengths = cdr3_seqs.str.len()

# Create the histogram
plt.figure(figsize=(8, 4))
plt.hist(
    cdr3_lengths,
    bins=range(5, 30),        # bin edges from 5 to 29
    color='steelblue',
    edgecolor='white',
    align='left'              # align bars to the left edge of each bin
)
plt.xlabel('CDR3 Length (amino acids)')
plt.ylabel('Number of Cells')
plt.title('CDR3 Beta Chain Length Distribution')

# Add a vertical line at the median length
median_len = cdr3_lengths.median()
plt.axvline(median_len, color='red', linestyle='--', label=f'Median: {median_len:.0f} aa')
plt.legend()

plt.tight_layout()
plt.savefig('cdr3_length.png', dpi=150)
plt.show()
print(f"\nMedian CDR3 length: {median_len:.0f} amino acids")
print(f"Range: {cdr3_lengths.min()} to {cdr3_lengths.max()} amino acids")
print("Saved: cdr3_length.png")
```

**Expected output:**

A histogram showing a roughly bell-shaped distribution centered around 13-16 amino acids. The red dashed line marks the median. Most CDR3 sequences cluster in a narrow range because the basic structure of the CDR3 is constrained by the V and J gene segments on either side.

### V Gene Usage

Beyond CDR3 sequences, we can also look at which V gene segments are most commonly used:

```python
# Count how often each V gene is used (for the beta chain)
v_gene_col = 'IR_VDJ_1_v_call'
if v_gene_col in adata.obs.columns:
    v_gene_counts = adata.obs[v_gene_col].dropna().value_counts()

    # Plot the top 20 most common V genes
    plt.figure(figsize=(10, 5))
    v_gene_counts.head(20).plot(
        kind='bar',
        color='steelblue',
        edgecolor='white'
    )
    plt.xlabel('V Gene Segment')
    plt.ylabel('Number of Cells')
    plt.title('Top 20 Most Common V Gene Segments (Beta Chain)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('v_gene_usage.png', dpi=150)
    plt.show()
    print("Saved: v_gene_usage.png")
else:
    print(f"Column '{v_gene_col}' not found. Available IR columns:")
    print([c for c in adata.obs.columns if 'v_call' in c.lower() or 'v_gene' in c.lower()])
```

**Expected output:**

A bar chart showing the frequency of different V gene segments. Some V genes will be much more common than others. This "V gene usage bias" is a well-known phenomenon: certain V genes are physically closer to the J genes in the genome and get used more often during VDJ recombination.

### [DEEP DIVE] Why CDR3 Length Matters

The length of the CDR3 directly affects the shape of the antigen-binding surface. Longer CDR3 loops can reach into deeper pockets on antigens, while shorter ones interact with flatter surfaces. In BCR heavy chains, CDR3 lengths tend to be longer and more variable than in TCR beta chains. Some studies have shown that certain diseases are associated with skewed CDR3 length distributions -- for example, some autoimmune diseases show a narrower CDR3 length distribution, suggesting clonal expansion of specific receptor types.

---

## 14. [DEEP DIVE] Key Databases for TCR/BCR Research

As you progress beyond this bootcamp, you will want to compare your data against curated databases of known receptor sequences. Here are the most important ones in the field.

### VDJdb
**URL:** https://vdjdb.cdr3.net/

A curated database of TCR sequences with **known antigen specificity**. If you have a CDR3 sequence and want to know what antigen it might recognize, VDJdb is the first place to look. Researchers submit experimentally validated TCR-antigen pairs from studies around the world.

**Use case:** "I found a highly expanded clonotype in my cancer patient. Does anyone else's data show the same CDR3, and do we know what antigen it recognizes?"

### IEDB (Immune Epitope Database)
**URL:** https://www.iedb.org/

The largest repository of experimentally characterized immune epitopes (the specific parts of antigens that receptors recognize). It covers both T cell and B cell epitopes across thousands of species and diseases.

**Use case:** "I know my patient's T cells respond to a specific viral peptide. What other epitopes from this virus have been characterized?"

### ImmuneCODE
**URL:** https://clients.adaptivebiotech.com/pub/covid-2020

A massive dataset of TCR sequences from COVID-19 patients, made publicly available by Adaptive Biotechnologies. It was one of the largest coordinated TCR sequencing efforts in history.

**Use case:** "I want to study the T cell response to SARS-CoV-2. ImmuneCODE provides thousands of COVID-associated TCR sequences to analyze."

### TCR3d
**URL:** https://tcr3d.ibbr.umd.edu/

A structural database containing 3D structures of TCR-pMHC complexes (how the TCR physically docks onto the antigen-presenting molecule). Understanding these structures is essential for designing better vaccines and T cell therapies.

**Use case:** "I want to visualize how a specific TCR binds to its target at the atomic level."

### OGRDB (Open Germline Receptor Database)
**URL:** https://ogrdb.airr-community.org/

A reference database for germline (inherited) V, D, and J gene segments. Different people can have slightly different sets of V/D/J genes, and OGRDB curates these population-level differences.

**Use case:** "I want to make sure I am using the correct reference gene segments for my patient's ethnic background when calling V/D/J genes."

These databases will become especially important in Week 4 and the final project, where you will integrate multiple data types and start asking more specific biological questions.

---

## 15. Self-Check

Congratulations on making it through Week 3. Before you move on, work through this self-check to make sure the key concepts have landed.

### Conceptual Questions

Answer these in your own words. If you can explain them clearly, you truly understand the material.

1. **What is a CDR3 region and why is it important?**
   The CDR3 is the most variable part of the immune receptor -- the short amino acid sequence that directly contacts the antigen. It is important because it determines what each individual T or B cell can recognize. In computational analysis, CDR3 sequences are the primary identifier for clonotypes.

2. **What is a clonotype?**
   A clonotype is a group of cells that share the same immune receptor sequence (typically defined by identical CDR3 sequences). Cells in the same clonotype descended from a common ancestor cell.

3. **What does "clonal expansion" mean biologically?**
   When a T cell or B cell recognizes an antigen, it divides rapidly to produce many copies of itself. This process is clonal expansion. In the data, an expanded clonotype is one that appears in many cells -- evidence that the original cell found and responded to a specific target.

4. **What is the difference between repertoire richness and Shannon diversity?**
   Richness simply counts how many unique clonotypes exist. Shannon diversity accounts for both the number of clonotypes AND how evenly they are distributed. A repertoire with 100 equally abundant clonotypes has higher Shannon diversity than one with 100 clonotypes where one dominates 99% of the cells.

5. **In your analysis, which cell clusters showed the most clonal expansion?**
   This answer depends on your specific results. Typically, CD8+ effector and effector memory clusters show the highest expansion, while naive clusters show the least. Look at your `clonal_expansion.png` and cross-tabulation table to answer this for your data.

### Output Files to Verify

Make sure you have generated all of the following files in your working directory:

- [ ] `chain_pairing.png` -- Bar chart of TCR chain pairing quality categories
- [ ] `clonal_expansion.png` -- Stacked bar chart of clonal expansion by cluster
- [ ] `diversity.png` -- Violin plot of Shannon diversity by cluster
- [ ] `umap_clonal.png` -- Side-by-side UMAPs (clusters and clonal expansion)

If you completed the challenge sections, you should also have:

- [ ] `clonotype_overlap.png` -- Heatmap of clonotype sharing between groups
- [ ] `cdr3_length.png` -- Histogram of CDR3 length distribution
- [ ] `v_gene_usage.png` -- Bar chart of V gene segment usage
- [ ] `umap_top_clone.png` -- UMAP highlighting the most expanded clonotype

### Troubleshooting Common Issues

| Problem | Solution |
|---------|----------|
| `ModuleNotFoundError: No module named 'scirpy'` | Run `pip install scirpy` in your terminal with the immunology environment activated |
| `KeyError: 'clone_id'` | You need to run `ir.pp.ir_dist()` and `ir.tl.define_clonotypes()` first |
| `KeyError: 'chain_pairing'` | You need to run `ir.tl.chain_qc()` first |
| `KeyError: 'clonal_expansion'` | You need to run `ir.tl.clonal_expansion()` first -- and clonotypes must be defined before that |
| UMAP looks blank or all one color | Make sure the column you are plotting has non-null values: `print(adata.obs['column_name'].value_counts())` |
| `ir.datasets.wu2020()` gives a download error | Check your internet connection. If the download is blocked, try again later or check scirpy's documentation for alternative download methods |
| Plots are too small or labels overlap | Increase the figure size: change `figsize=(8, 4)` to `figsize=(12, 6)` |

### You Did It

You can now analyze immune receptor repertoires. This is a skill used in cutting-edge cancer immunotherapy research, vaccine development, and infectious disease surveillance. When researchers at companies like Genentech, Adaptive Biotechnologies, or the Broad Institute analyze patient samples to understand T cell responses, they use the exact same tools and techniques you just practiced.

The analyses you performed this week -- clonotype identification, clonal expansion quantification, diversity measurement, and clonotype-to-phenotype mapping -- form the foundation of every immune repertoire study published in journals like Nature, Cell, and Immunity.

---

## 16. Key Papers (Optional Reading)

If you want to go deeper, here are three foundational papers. You do not need to read them to continue the bootcamp, but they provide context for the analyses you just performed.

1. **Dash et al., 2017, Nature** -- "Quantifiable predictive features define epitope-specific T cell receptor repertoires"
   This paper showed that TCR repertoires targeting specific antigens have predictable sequence features. It was one of the first demonstrations that computational methods could predict what a TCR recognizes. This laid the groundwork for modern TCR specificity prediction tools.

2. **Briney et al., 2019, Nature** -- "Commonality despite exceptional diversity in the baseline human antibody repertoire"
   This study sequenced BCR repertoires from multiple healthy individuals at unprecedented depth. They found that while BCR diversity is enormous, there are surprising commonalities between different people -- certain antibody sequences appear again and again. This has implications for vaccine design and antibody therapeutics.

3. **Wu et al., 2020, Nature Medicine** -- "Peripheral T cell expansion predicts tumour infiltration and clinical response"
   This is the source of the dataset you used in this tutorial. The researchers showed that T cell clones found in the tumor are also detectable in peripheral blood, and that tracking these clones in blood samples could predict how well patients respond to immunotherapy. Your clonotype overlap analysis in Section 12 is directly relevant to this finding.

---

## 17. What's Coming Next

In **Week 4: Multi-omics Integration**, you will learn to combine multiple types of data from the same cells:

- **Gene expression** (scRNA-seq) -- which genes are active
- **Receptor sequences** (scTCR-seq/scBCR-seq) -- what each cell recognizes
- **Surface protein levels** (CITE-seq) -- what proteins are on the cell surface
- **Chromatin accessibility** (scATAC-seq) -- which parts of the genome are "open" for use

Integrating these different data modalities gives you a far richer picture of each cell's identity and function than any single measurement alone. You will use tools like MuData and MultiVI to bring it all together.

The receptor analysis skills you built this week are essential for Week 4. The ability to link clonotype identity to gene expression, protein levels, and chromatin state is one of the most powerful techniques in modern immunology research.

See you in [Week 4](../week_4_multi_omics/README.md).

---

## Quick Reference: Key Commands

Here is a one-page summary of the scirpy commands you learned this week, for easy reference:

```python
# --- Setup ---
import scirpy as ir

# --- Load data ---
adata = ir.datasets.wu2020()              # Built-in tutorial dataset

# --- Quality control ---
ir.tl.chain_qc(adata)                     # Assess chain pairing quality

# --- Define clonotypes ---
ir.pp.ir_dist(adata, metric='identity',   # Compute sequence distances
              sequence='nt')
ir.tl.define_clonotypes(adata)            # Group cells into clonotypes

# --- Clonal expansion ---
ir.tl.clonal_expansion(adata)             # Categorize expansion levels
ir.pl.clonal_expansion(adata,             # Plot expansion by group
                        groupby='cluster',
                        clip_at=4,
                        normalize=True)

# --- Diversity ---
ir.tl.alpha_diversity(adata,              # Calculate diversity metrics
                       groupby='cluster',
                       metric='shannon')

# --- Clonotype overlap ---
ir.tl.clonotype_overlap(adata,            # Compare repertoires between groups
                          groupby='patient',
                          metric='jaccard')
ir.pl.clonotype_overlap(adata,            # Visualize overlap as heatmap
                          groupby='patient')
```

---

*This tutorial is part of the [Computational Immunology Bootcamp](../README.md). Go back to [Week 2](../week_2_single_cell_pipeline/README.md) or continue to [Week 4](../week_4_multi_omics/README.md).*
