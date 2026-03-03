# Week 1: Immunology Foundations

**Estimated time: 5-8 hours**

Welcome to Week 1 of the Computational Immunology Bootcamp! Last week you set up your
computing environment -- Python, Miniconda, R, Git, and Jupyter Notebook. This week you
will learn the biology that makes computational immunology possible, and you will write
your first Python code to explore real immune cell data. No biology or programming
background is assumed. Take your time, re-read sections that feel confusing, and remember:
every expert was once a beginner.

---

## Table of Contents

1. [What Is the Immune System?](#1-what-is-the-immune-system)
2. [T Cells and B Cells: The Stars of Adaptive Immunity](#2-t-cells-and-b-cells-the-stars-of-adaptive-immunity)
3. [Antigen Presentation: How the Immune System Shares Information](#3-antigen-presentation-how-the-immune-system-shares-information)
4. [What Is Gene Expression?](#4-what-is-gene-expression)
5. [What Is Single-Cell RNA-seq? (Conceptual)](#5-what-is-single-cell-rna-seq-conceptual)
6. [Your First Python Notebook: Exploring Immune Cell Data](#6-your-first-python-notebook-exploring-immune-cell-data)
7. [Hands-On: Simple Data Exploration](#7-hands-on-simple-data-exploration)
8. [[CHALLENGE] Reading a Scientific Paper](#8-challenge-reading-a-scientific-paper)
9. [Key Vocabulary This Week](#9-key-vocabulary-this-week)
10. [Self-Check](#10-self-check)
11. [What's Coming Next](#11-whats-coming-next)
12. [Key Papers (Optional Reading)](#key-papers-optional-reading)

---

## 1. What Is the Immune System?

### The City Defense Analogy

Imagine your body is a large, walled city. Every day, invaders -- bacteria, viruses,
fungi, parasites -- try to break in. Your immune system is the entire defense
infrastructure of that city: the walls, the guards at the gates, the surveillance cameras,
the police force, and an elite special-operations unit that remembers every criminal it has
ever fought.

The immune system is not a single organ. It is a network of cells, tissues, and molecules
spread throughout your entire body. Its job is to (1) recognize things that do not belong,
(2) destroy them, and (3) remember them so the response is faster next time.

### Innate vs. Adaptive Immunity

The defense system has two major branches:

**Innate immunity** is the part you are born with. Think of it as the city walls, the moat,
and the standing army of guards. It responds to invaders within minutes to hours and
attacks anything that looks foreign -- it does not care exactly *which* bacterium showed
up, only that it is not "self." Innate immunity includes physical barriers (skin, mucus),
chemical barriers (stomach acid, antimicrobial proteins), and a set of fast-responding
immune cells. Innate immunity is broad and immediate, but it does not improve with
experience.

**Adaptive immunity** is the elite special-forces unit. It takes days to mount on a first
encounter, but it is extraordinarily precise: it can tell the difference between two nearly
identical viruses. Most importantly, it has *memory*. After the first fight, adaptive
immune cells remember the invader. If the same pathogen shows up again months or years
later, the adaptive system responds within hours instead of days. This is why vaccines
work -- they train the adaptive system in advance.

### Key Immune Cell Types

The table below lists the major cell types you will encounter repeatedly in this bootcamp.
You do not need to memorize every detail right now, but bookmark this table so you can
refer back to it.

| Cell Type | Branch | Analogy | What It Does |
|-----------|--------|---------|--------------|
| **Neutrophil** | Innate | Front-line foot soldiers | First to arrive at an infection site; engulfs and kills bacteria; short-lived |
| **Macrophage** | Innate | Heavy-duty cleanup crew | Eats (phagocytoses) pathogens and dead cells; presents pieces of what it ate to adaptive cells |
| **Dendritic Cell (DC)** | Innate (bridges to adaptive) | Intelligence officers | Patrols tissues, captures pathogens, travels to lymph nodes to *brief* T cells |
| **NK Cell** | Innate | Undercover assassins | Kills virus-infected cells and tumor cells without needing prior instructions |
| **CD4+ T Cell (Helper T)** | Adaptive | Commanders / Generals | Coordinates the immune response by sending chemical signals (cytokines) that activate other cells |
| **CD8+ T Cell (Killer T)** | Adaptive | Special-ops assassins | Directly kills infected cells by recognizing foreign fragments on their surface |
| **B Cell** | Adaptive | Weapons factory | Produces antibodies -- Y-shaped proteins that tag pathogens for destruction |

> **Quick note on "CD" numbers.** CD stands for "Cluster of Differentiation." These are
> just numbered labels for specific proteins found on the surface of cells. CD4 and CD8 are
> surface proteins that help us distinguish the two major types of T cells. You will see CD
> numbers constantly in immunology -- think of them as name tags.

`[CHECKPOINT 1]` You should now be able to explain, in your own words, the difference
between innate and adaptive immunity, and name at least four types of immune cells.

---

## 2. T Cells and B Cells: The Stars of Adaptive Immunity

T cells and B cells are the two pillars of adaptive immunity. Almost everything in
computational immunology revolves around understanding what these cells are doing, how
diverse they are, and how they respond to threats. Let's look at each in detail.

### T Cells: Born in the Thymus

T cells are born from stem cells in the **bone marrow**, but they grow up and get trained
in a small organ called the **thymus** (located just above your heart -- the "T" in T cell
stands for thymus). In the thymus, young T cells go through a brutal selection process:
any T cell that reacts too strongly to your own body's proteins is killed off. Only the
cells that can tell "self" from "non-self" graduate and enter the bloodstream. This process
is called **thymic selection**, and it is a central topic in immunology.

There are two main types of mature T cells:

- **CD4+ Helper T cells -- the Commanders.** They do not kill anything directly. Instead,
  they recognize threats and release chemical messengers called **cytokines** that tell
  other immune cells what to do. Think of them as generals who issue orders: "B cells,
  start making antibodies!" or "Macrophages, ramp up the inflammation!"

- **CD8+ Killer T cells -- the Assassins.** These cells directly kill your own cells that
  have been infected by a virus or have become cancerous. They physically dock onto the
  sick cell and inject toxic molecules that trigger the cell to self-destruct.

### B Cells: Born in the Bone Marrow

B cells are born *and* trained in the **bone marrow** (the "B" stands for bone marrow).
Their superpower is making **antibodies** -- Y-shaped proteins that are released into the
blood and body fluids. Each antibody is shaped to grab onto one specific molecular pattern
(called an **antigen**). When antibodies stick to a pathogen, they can (1) neutralize it
directly, (2) flag it so macrophages eat it faster, or (3) activate other defense systems.

When a B cell encounters its matching antigen and gets the green light from a helper T
cell, it can multiply rapidly and transform into a **plasma cell** -- an antibody-producing
machine that pumps out thousands of antibody molecules per second. Some B cells instead
become **memory B cells** that persist for years, ready to respond immediately if the same
pathogen returns.

### TCR and BCR: Unique Receptors

Here is one of the most remarkable facts in biology: every single T cell carries a unique
receptor on its surface called the **T Cell Receptor (TCR)**, and every B cell carries a
unique **B Cell Receptor (BCR)**. "Unique" means that your body contains millions of
different T cells, each with a slightly different receptor shape.

Think of it this way: every T cell has a unique **lock**, and it is searching for the one
specific **key** (antigen) that fits. When the key fits the lock, the T cell activates and
mounts an immune response. The same logic applies to B cells and their BCRs.

This raises an obvious question: you only have about 20,000 genes. How can your body
produce millions of unique receptor shapes?

### VDJ Recombination: Shuffling the Deck

The answer is **VDJ recombination**, one of the most elegant tricks in biology. The genes
that encode TCRs and BCRs are not stored as single, complete genes. Instead, they are
stored in segments: **V** (Variable), **D** (Diversity), and **J** (Joining) segments.
There are many versions of each segment.

When a new T cell or B cell is developing, the cell's DNA-editing machinery randomly picks
one V segment, one D segment, and one J segment, and permanently splices them together.
It also adds or removes a few random letters (nucleotides) at the junctions. This
combinatorial process is like shuffling a deck of cards and then cutting the deck at random
positions -- the number of possible combinations is astronomical.

Here is a simplified example:

```
Gene segments in your genome:

  V segments:  V1, V2, V3, V4, V5, ... (dozens of options)
  D segments:  D1, D2, D3, ...          (a handful of options)
  J segments:  J1, J2, J3, J4, ...      (several options)

One cell picks:  V3 -- random letters -- D1 -- random letters -- J4
Another picks:   V5 -- random letters -- D2 -- random letters -- J1

Result: Each cell has a unique receptor sequence.
```

VDJ recombination can theoretically generate more than 10^15 (a quadrillion) different
receptor sequences. This diversity is the reason your adaptive immune system can recognize
virtually any pathogen it might encounter, including ones that have never existed before.

The specific region where the V, D, and J segments join -- including the random insertions
-- is called the **CDR3** (Complementarity-Determining Region 3). The CDR3 is the most
variable part of the receptor and is the main determinant of what antigen the receptor will
bind. When computational immunologists analyze TCR or BCR data, the CDR3 sequence is
usually the star of the show.

`[CHECKPOINT 2]` You should now be able to explain: (a) where T cells and B cells mature,
(b) what CD4 and CD8 T cells do, (c) what antibodies are, and (d) how VDJ recombination
generates receptor diversity.

---

## 3. Antigen Presentation: How the Immune System Shares Information

For adaptive immunity to work, T cells need to "see" what is going on inside cells. But
T cells cannot peek inside other cells directly. Instead, the immune system has evolved an
elegant information-sharing system based on molecules called **MHC** (Major
Histocompatibility Complex).

### MHC Class I: Showing Your ID Badge

Every nucleated cell in your body (that is, nearly every cell except red blood cells)
displays **MHC class I** molecules on its surface. Think of MHC I as a display shelf on
the outside of a cell. The cell constantly chops up some of the proteins it is making
inside, and places small fragments (peptides) onto MHC I molecules for display.

If the cell is healthy, the displayed fragments are all "self" -- normal human proteins.
CD8+ killer T cells patrol the body and glance at these displays. When they see normal
self-peptides, they move on. But if a cell is infected by a virus, some of the displayed
fragments will come from viral proteins. A CD8+ T cell whose receptor matches that viral
fragment will recognize it, dock onto the cell, and kill it.

**Analogy:** Imagine every person in a building is required to wear an ID badge showing
what they have been working on. Security guards (CD8 T cells) walk the halls glancing at
badges. If your badge shows normal work, they pass by. If your badge shows you have been
making counterfeit money (viral proteins), the guard takes you down.

### MHC Class II: Briefing the Commanders

**MHC class II** molecules are only found on a special group of cells called
**antigen-presenting cells (APCs)**. The main APCs are **dendritic cells**, **macrophages**,
and **B cells**. These cells eat (phagocytose) pathogens, chop them up, and display the
fragments on MHC II molecules.

MHC II presentations are read by **CD4+ helper T cells**. When a helper T cell sees a
foreign fragment on MHC II, it activates and starts coordinating the immune response --
sending out cytokines, helping B cells make antibodies, and boosting killer T cell activity.

**Analogy:** Antigen-presenting cells are like intelligence officers who capture an enemy
spy, take their photo and fingerprints, and bring the dossier to the generals (CD4+ T
cells). The generals then decide the battle plan.

### Summary Diagram (Text Version)

```
  INSIDE AN INFECTED CELL                 ANTIGEN-PRESENTING CELL
  +-------------------------+             +-------------------------+
  |  Virus replicating      |             |  Eats pathogen          |
  |  inside the cell        |             |  Chops it up            |
  |                         |             |                         |
  |  Viral protein fragment |             |  Pathogen fragment      |
  |  loaded onto MHC I      |             |  loaded onto MHC II     |
  +--------|----------------+             +--------|----------------+
           |                                       |
           v                                       v
    Displayed on cell surface              Displayed on cell surface
           |                                       |
           v                                       v
    CD8+ Killer T cell sees it          CD4+ Helper T cell sees it
    --> Kills the infected cell          --> Activates immune response
```

`[CHECKPOINT 3]` You should now be able to explain the difference between MHC class I and
MHC class II, and which type of T cell reads each one.

---

## 4. What Is Gene Expression?

Before we can understand single-cell RNA-seq (the technology at the heart of this
bootcamp), we need to understand a foundational concept in biology: **gene expression**.

### DNA -> RNA -> Protein: The Central Dogma

Your DNA is like a massive recipe book stored in the nucleus of every cell. Each **gene**
is one recipe -- instructions for building one specific protein. Proteins are the molecular
machines that do almost everything in your body: they form structures, speed up chemical
reactions, carry signals, and fight infections.

But DNA itself does not leave the nucleus. When a cell needs to build a specific protein,
it makes a temporary copy of the relevant gene. That copy is called **messenger RNA
(mRNA)**. The mRNA travels out of the nucleus and is read by molecular machines called
**ribosomes**, which use the instructions to assemble the protein.

```
DNA  ---[transcription]--->  mRNA  ---[translation]--->  Protein
(recipe book)                (photocopy of one recipe)   (the dish)
```

This flow -- DNA to RNA to Protein -- is called the **Central Dogma of Molecular Biology**.

### Why Gene Expression Matters

Here is the key insight: **every cell in your body has the same DNA** (the same recipe
book), but different cells "read" different recipes. A liver cell turns on liver-specific
genes and turns off brain-specific genes. A T cell turns on immune-related genes and turns
off muscle-specific genes. The set of genes that a cell has turned on (expressed) at any
given moment defines that cell's identity and behavior.

When we say a gene is **"expressed"**, we mean the cell is actively making mRNA copies of
that gene. The more mRNA copies, the more protein the cell is probably making, and the more
"active" that gene is in that cell.

### This Is What Single-Cell RNA-seq Measures

Single-cell RNA-seq captures and counts the mRNA molecules in individual cells. By
measuring which genes are expressed and how strongly, we can figure out:

- What type of cell it is (T cell? B cell? macrophage?)
- What the cell is doing right now (resting? activated? dividing?)
- How the cell differs from its neighbors

`[CHECKPOINT 4]` You should now be able to explain the Central Dogma (DNA -> RNA ->
Protein) and why different cells express different genes even though they share the same
DNA.

---

## 5. What Is Single-Cell RNA-seq? (Conceptual)

### The Old Way: Bulk RNA-seq

Before single-cell technology, scientists would take a tissue sample containing thousands
or millions of cells, grind them all up together, and measure the average gene expression
across all those cells. This is called **bulk RNA-seq**.

**Analogy:** Imagine you have a fruit salad with strawberries, blueberries, mangoes, and
kiwis. You throw the entire salad into a blender and taste the smoothie. You can tell it
is fruity, maybe sweet and tart, but you cannot tell exactly which fruits were in the
salad, their proportions, or what any individual piece of fruit tasted like.

Bulk RNA-seq has the same problem. If your sample contains 60% T cells and 40% B cells,
you get a blended average. You cannot see the individual cell types or rare subpopulations.

### The New Way: Single-Cell RNA-seq (scRNA-seq)

Single-cell RNA-seq measures gene expression in **each cell individually**. Instead of
blending the fruit salad, you taste each piece of fruit one at a time. Now you know
exactly how many strawberries there are, what each one tastes like, and whether any
individual strawberry is different from the rest.

This technology has revolutionized biology. Since about 2015, single-cell methods have
allowed scientists to discover new cell types, track how cells change during development,
and understand how diseases alter immune responses at unprecedented resolution.

### How 10x Genomics Works (Simplified)

The most widely used single-cell RNA-seq platform is made by a company called **10x
Genomics**. Here is a simplified version of how it works:

1. **Cell suspension:** You start with a mixture of individual cells in liquid.
2. **Droplet formation:** A microfluidic chip combines the cells with tiny gel beads in
   oil, forming thousands of nanoliter-sized droplets. Each droplet ideally contains
   **one cell** and **one bead**.
3. **Barcoding:** Each bead is coated with millions of DNA primers that all share the same
   unique **barcode** sequence. When the cell's mRNA is captured by the bead's primers,
   every mRNA molecule from that cell gets tagged with the bead's barcode.
4. **Sequencing:** All the barcoded mRNA is pooled together and sequenced on a standard
   DNA sequencer.
5. **Computational demultiplexing:** A computer reads the barcode on each sequenced
   molecule and assigns it back to the correct cell of origin. The result is a count of
   how many mRNA molecules of each gene were detected in each cell.

### The Data: A Cells-by-Genes Matrix

The end product of a single-cell RNA-seq experiment is a large matrix (spreadsheet):

```
             Gene_A  Gene_B  Gene_C  Gene_D  ...  Gene_20000
Cell_001        0       3       0      12     ...      0
Cell_002        5       0       1       0     ...      2
Cell_003        0       0       0       8     ...      0
   ...
Cell_2700       1       7       0       0     ...      0
```

- **Rows** = individual cells (identified by their barcodes)
- **Columns** = genes (about 20,000 in humans)
- **Values** = the number of mRNA molecules detected for that gene in that cell (called
  "counts" or "UMI counts")

Notice how many zeros there are. Single-cell data is very **sparse** -- most genes in any
given cell have zero detected counts. This is partly because each cell only expresses a
fraction of all genes, and partly because the technology does not capture every single mRNA
molecule (a phenomenon called **dropout**).

This matrix is what you will load and analyze in the hands-on exercises below.

`[CHECKPOINT 5]` You should now be able to explain: (a) the difference between bulk and
single-cell RNA-seq, (b) the basic idea behind how 10x Genomics captures single cells,
and (c) what a cells-by-genes matrix looks like.

---

## 6. Your First Python Notebook: Exploring Immune Cell Data

**[REQUIRED]**

Time to write code! In this exercise you will open a Jupyter Notebook, load a real
single-cell dataset of human immune cells, and explore its basic structure. This dataset
-- called **pbmc3k** -- contains approximately 2,700 **peripheral blood mononuclear cells
(PBMCs)** from a healthy human donor. PBMCs are immune cells isolated from blood, so
this dataset contains a mixture of T cells, B cells, NK cells, monocytes, and dendritic
cells.

### Step 1: Activate Your Environment and Launch Jupyter

Open your terminal and type the following commands, one at a time:

```bash
conda activate immunology
jupyter notebook
```

Your web browser should open with the Jupyter file browser. Navigate to a folder where you
want to keep your work for this bootcamp (for example, create a folder called
`week_1_immunology_foundations` if you have not already).

Click **New -> Python 3** to create a new notebook. Rename it to
**`week1_explore_pbmc3k.ipynb`** by clicking on the title at the top.

### Step 2: Install scanpy (if not already installed)

In the first cell of your notebook, type and run:

```python
# Run this cell only once -- it installs the scanpy package
# After it finishes, you can delete or skip this cell
!pip install scanpy
```

Wait for the installation to complete. You will see output scrolling by as packages are
downloaded. When it finishes, you should see a line near the end that says something like:

```
Successfully installed scanpy-X.X.X ...
```

> **Note:** If you already installed scanpy during Week 0, this cell will simply confirm
> that the requirement is already satisfied. That is fine.

### Step 3: Import Libraries

In a new cell, type and run:

```python
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

print("All imports successful!")
print(f"Scanpy version: {sc.__version__}")
```

**Expected output:**

```
All imports successful!
Scanpy version: 1.x.x
```

(The exact version number may differ -- that is fine.)

### Step 4: Load the PBMC 3k Dataset

In a new cell, type and run:

```python
# Download a small practice dataset of ~2,700 blood immune cells
# This is built into scanpy, so no manual download is needed
adata = sc.datasets.pbmc3k()
print(adata)
```

**Expected output (approximately):**

```
AnnData object with n_obs x n_vars = 2700 x 32738
    var: 'gene_ids'
```

Let's break down what this output means:

- **AnnData object**: This is the data structure that scanpy uses to store single-cell
  data. Think of it as a smart spreadsheet that holds not just the data matrix but also
  metadata about cells and genes.
- **n_obs = 2700**: There are 2,700 observations (cells).
- **n_vars = 32738**: There are 32,738 variables (genes).
- **var: 'gene_ids'**: The gene table has a column called 'gene_ids' that stores Ensembl
  gene identifiers (a standardized naming system).

### Step 5: Explore the Data Shape

In a new cell, type and run:

```python
# How many cells and genes are in this dataset?
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
```

**Expected output:**

```
Number of cells: 2700
Number of genes: 32738
```

### Step 6: Look at Gene Names and Cell Barcodes

In a new cell, type and run:

```python
# Look at some gene names
print("First 10 genes:", list(adata.var_names[:10]))
print()

# Look at some cell barcodes
print("First 5 cells:", list(adata.obs_names[:5]))
```

**Expected output (approximately):**

```
First 10 genes: ['MIR1302-10', 'FAM138A', 'OR4F5', 'RP11-34P13.7', 'RP11-34P13.8', 'AL627309.1', 'RP11-34P13.14', 'RP11-34P13.9', 'AP006222.2', 'RP4-669L17.10']

First 5 cells: ['AAACATACAACCAC-1', 'AAACATTGAGCTAC-1', 'AAACATTGATCAGC-1', 'AAACCGTGCTTCCG-1', 'AAACCGTGTATGCG-1']
```

The gene names are human gene symbols (some look cryptic -- that is normal). The cell
names are barcode sequences assigned by the 10x Genomics instrument, with a `-1` suffix
indicating the sample.

### Step 7: Peek at the Data Matrix

In a new cell, type and run:

```python
# Look at the raw data matrix
print(f"Data matrix type: {type(adata.X)}")
print(f"Data matrix shape: {adata.X.shape}")
```

**Expected output:**

```
Data matrix type: <class 'scipy.sparse.csr_matrix'>
Data matrix shape: (2700, 32738)
```

The data is stored as a **sparse matrix** -- a special format that efficiently stores data
with lots of zeros (which, as we discussed, single-cell data has plenty of). Instead of
storing every single zero, it only stores the non-zero values and their positions.

`[CHECKPOINT 6]` Congratulations! You have successfully loaded a real single-cell dataset
in Python. You know it contains 2,700 cells and over 32,000 genes. You have seen gene
names and cell barcodes. You are officially doing computational biology.

---

## 7. Hands-On: Simple Data Exploration

**[REQUIRED]**

Now let's go further. Continue in the same Jupyter notebook.

### Step 1: How Many Genes Does Each Cell Express?

Not every gene is active in every cell. Let's count how many genes have non-zero expression
in each cell.

In a new cell, type and run:

```python
# Count the number of genes with nonzero expression in each cell
# adata.X is a sparse matrix; indptr marks the start/end of each row's data
genes_per_cell = np.diff(adata.X.indptr)

print(f"Average genes per cell: {genes_per_cell.mean():.0f}")
print(f"Minimum: {genes_per_cell.min()}")
print(f"Maximum: {genes_per_cell.max()}")
print(f"Median: {np.median(genes_per_cell):.0f}")
```

**Expected output (approximately):**

```
Average genes per cell: 817
Minimum: 212
Maximum: 4538
Median: 718
```

(Your exact numbers may vary slightly depending on the scanpy version.)

What does this tell us? On average, each cell expresses about 800 genes out of the 32,000+
in the genome. Some cells express as few as ~200 genes (these might be low-quality cells or
a cell type that is relatively inactive), and some express over 4,000 (these might be large,
highly active cells). This range is typical for 10x Genomics scRNA-seq data.

### Step 2: Visualize with a Histogram

A histogram is one of the most useful plots in data science. It shows the distribution of
a variable -- in this case, how many cells have a given number of expressed genes.

In a new cell, type and run:

```python
plt.figure(figsize=(8, 4))
plt.hist(genes_per_cell, bins=50, color='steelblue', edgecolor='white')
plt.xlabel('Number of Genes Detected')
plt.ylabel('Number of Cells')
plt.title('Distribution of Genes per Cell')
plt.tight_layout()
plt.savefig('genes_per_cell.png', dpi=150)
plt.show()
```

**Expected output:** A histogram should appear below the cell. You should see a
right-skewed distribution -- most cells cluster on the left side (expressing a moderate
number of genes), with a long tail to the right (a smaller number of cells expressing
many genes).

**What does this histogram tell us?**

- The **peak** of the distribution (the tallest bar) tells you the most common number of
  genes per cell -- roughly 500-700 for this dataset.
- The **spread** tells you about cell diversity. Some cells are more active than others.
- Cells in the far left tail (very few genes) might be dead or dying cells -- this is
  something we will learn to filter out in Week 2.
- Cells in the far right tail (very many genes) might be **doublets** -- two cells that
  were accidentally captured in the same droplet, so their mRNA got mixed.

The file `genes_per_cell.png` has been saved to your current directory. You just created
a publication-quality figure from real data.

### Step 3: Check a Specific Gene

Let's look at a specific immune-related gene. **CD3D** is a gene that is expressed
primarily in T cells. Let's see how its expression looks across cells.

In a new cell, type and run:

```python
# Check if CD3D is in the dataset
if 'CD3D' in adata.var_names:
    # Get the expression values for CD3D across all cells
    cd3d_index = list(adata.var_names).index('CD3D')
    cd3d_expression = adata.X[:, cd3d_index].toarray().flatten()

    # How many cells express this gene?
    n_expressing = (cd3d_expression > 0).sum()
    print(f"CD3D is expressed in {n_expressing} out of {adata.n_obs} cells "
          f"({100 * n_expressing / adata.n_obs:.1f}%)")
    print(f"Mean expression (all cells): {cd3d_expression.mean():.2f}")
    print(f"Mean expression (expressing cells only): "
          f"{cd3d_expression[cd3d_expression > 0].mean():.2f}")
else:
    print("CD3D not found in this dataset.")
```

**Expected output (approximately):**

```
CD3D is expressed in 816 out of 2700 cells (30.2%)
Mean expression (all cells): 0.68
Mean expression (expressing cells only): 2.24
```

This makes biological sense. PBMCs are a mixture of many immune cell types, and T cells
make up roughly 40-70% of PBMCs. Not all T cells will have detectable CD3D due to dropout,
so detecting it in ~30% of cells is reasonable. The cells that do express it have a mean
count of about 2 -- meaning on average, each T cell had about 2 mRNA molecules of CD3D
captured.

### Step 4: Compare Two Marker Genes

Let's compare a T cell marker (CD3D) with a B cell marker (CD79A).

In a new cell, type and run:

```python
# Get expression for two marker genes
def get_gene_expression(adata, gene_name):
    """Extract expression values for a single gene."""
    idx = list(adata.var_names).index(gene_name)
    return adata.X[:, idx].toarray().flatten()

cd3d = get_gene_expression(adata, 'CD3D')   # T cell marker
cd79a = get_gene_expression(adata, 'CD79A')  # B cell marker

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

axes[0].hist(cd3d[cd3d > 0], bins=30, color='tomato', edgecolor='white')
axes[0].set_title('CD3D Expression (T cell marker)')
axes[0].set_xlabel('Expression Level (counts)')
axes[0].set_ylabel('Number of Cells')

axes[1].hist(cd79a[cd79a > 0], bins=30, color='cornflowerblue', edgecolor='white')
axes[1].set_title('CD79A Expression (B cell marker)')
axes[1].set_xlabel('Expression Level (counts)')
axes[1].set_ylabel('Number of Cells')

plt.tight_layout()
plt.savefig('marker_genes.png', dpi=150)
plt.show()

print(f"Cells expressing CD3D (T cells):  {(cd3d > 0).sum()}")
print(f"Cells expressing CD79A (B cells): {(cd79a > 0).sum()}")
print(f"Cells expressing both:            {((cd3d > 0) & (cd79a > 0)).sum()}")
print(f"Cells expressing neither:         {((cd3d == 0) & (cd79a == 0)).sum()}")
```

**Expected output:** Two side-by-side histograms and printed counts. You should see that
CD3D is expressed in more cells than CD79A (there are more T cells than B cells in blood),
and very few cells express both markers (because a cell is usually either a T cell or a B
cell, not both). The "neither" group includes NK cells, monocytes, dendritic cells, and
other cell types that express different markers.

`[CHECKPOINT 7]` You just visualized real scientific data from 2,700 human immune cells!
You computed summary statistics, made histograms, and investigated specific marker genes.
These are foundational skills in computational biology.

---

## 8. [CHALLENGE] Reading a Scientific Paper

This section is optional but highly recommended. Learning to read scientific papers is a
critical skill, and the earlier you start, the better.

### How to Read a Paper (The Secret)

Do **not** read a scientific paper from start to finish like a novel. That is a common
mistake. Instead, follow this order:

1. **Title and Abstract** -- Read these carefully. The abstract is a ~200-word summary of
   the entire paper. After reading it, you should know the main question, approach, and
   key finding.

2. **Figures and Figure Legends** -- Skip straight to the figures. In biology, the figures
   tell the story. Read the figure legends (captions) carefully. Try to understand what
   each panel shows.

3. **Introduction (last paragraph)** -- The last paragraph of the introduction almost
   always states the specific goals of the paper.

4. **Results section** -- Now read the results, referring to the figures as you go.

5. **Methods** -- Only read this if you want to understand *how* they did something
   specific.

6. **Discussion** -- Read this to understand the broader implications and limitations.

### Your Challenge: Read the Drop-seq Abstract

Look up this paper:

> Macosko EZ, Basu A, Satija R, et al. "Highly Parallel Genome-wide Expression Profiling
> of Individual Cells Using Nanoliter Droplets." *Cell*, 2015.

You can find it by searching for "Macosko 2015 Cell Drop-seq" on Google Scholar
(scholar.google.com). Read only the **abstract** and look at **Figure 1** (which shows
the Drop-seq technology).

After reading, try to answer these three questions:

1. What problem were the authors trying to solve?
2. What technology did they develop, and what is the basic idea?
3. How many cells did they analyze in their main experiment?

Do not worry if you do not understand every word. The goal is to practice the skill of
extracting key information from a scientific paper. You will get better at this with time.

---

## 9. Key Vocabulary This Week

| Term | Definition |
|------|-----------|
| **Innate immunity** | The branch of the immune system you are born with; fast but non-specific |
| **Adaptive immunity** | The branch that develops specific responses and remembers past infections |
| **Antigen** | Any molecule (usually a protein fragment) that can be recognized by the immune system |
| **Antibody** | A Y-shaped protein made by B cells that binds to a specific antigen |
| **T cell** | An adaptive immune cell that matures in the thymus; includes CD4+ helper and CD8+ killer types |
| **B cell** | An adaptive immune cell that matures in bone marrow and produces antibodies |
| **CD4+ T cell** | Helper T cell; coordinates immune responses by releasing cytokines |
| **CD8+ T cell** | Killer (cytotoxic) T cell; directly destroys infected or cancerous cells |
| **NK cell** | Natural Killer cell; an innate lymphocyte that kills virus-infected and tumor cells |
| **Macrophage** | An innate immune cell that eats pathogens and presents antigens |
| **Dendritic cell** | An innate immune cell that captures antigens and presents them to T cells |
| **TCR** | T Cell Receptor; a unique surface protein on each T cell that recognizes a specific antigen |
| **BCR** | B Cell Receptor; a unique surface protein on each B cell; the membrane-bound form of an antibody |
| **VDJ recombination** | The process of randomly assembling V, D, and J gene segments to create diverse receptors |
| **CDR3** | Complementarity-Determining Region 3; the most variable part of the TCR/BCR, key for antigen binding |
| **MHC class I** | A molecule on all nucleated cells that presents intracellular peptides to CD8+ T cells |
| **MHC class II** | A molecule on antigen-presenting cells that presents extracellular peptides to CD4+ T cells |
| **Gene expression** | The process by which information in a gene is used to make a functional product (usually a protein) |
| **Central Dogma** | The flow of genetic information: DNA -> RNA -> Protein |
| **mRNA** | Messenger RNA; the temporary copy of a gene that is used as a template for protein synthesis |
| **scRNA-seq** | Single-cell RNA sequencing; a technology that measures gene expression in individual cells |
| **10x Genomics** | A company that makes the most widely used scRNA-seq platform |
| **UMI** | Unique Molecular Identifier; a short barcode used to count individual mRNA molecules |
| **AnnData** | The data structure used by scanpy to store single-cell datasets |
| **Sparse matrix** | A data format that efficiently stores matrices with many zero values |
| **PBMC** | Peripheral Blood Mononuclear Cells; immune cells isolated from blood |

---

## 10. Self-Check

### Concepts You Should Be Able to Explain

Go through this checklist. For each item, try to explain it in your own words (out loud or
written down) without looking at the material above. If you cannot, go back and re-read
that section.

- [ ] The difference between innate and adaptive immunity
- [ ] What neutrophils, macrophages, and dendritic cells do (innate)
- [ ] What CD4+ helper T cells and CD8+ killer T cells do
- [ ] What B cells and antibodies do
- [ ] Where T cells and B cells mature
- [ ] What a TCR and BCR are
- [ ] How VDJ recombination generates millions of unique receptors
- [ ] What MHC class I presents and to which T cell type
- [ ] What MHC class II presents and to which T cell type
- [ ] The Central Dogma: DNA -> RNA -> Protein
- [ ] Why different cell types express different genes
- [ ] The difference between bulk RNA-seq and single-cell RNA-seq
- [ ] What a cells-by-genes matrix looks like and what the values represent

### Test Questions

Try to answer these without looking:

1. What is the difference between innate and adaptive immunity?
2. What do CD4+ T cells do? What about CD8+ T cells?
3. What is VDJ recombination and why does it matter?
4. What molecule do B cells produce, and what shape is it?
5. Which MHC class is found on all nucleated cells?
6. What does scRNA-seq measure?
7. How many cells are in the pbmc3k dataset?
8. Approximately how many genes does the average cell in pbmc3k express?
9. Why does single-cell data have so many zeros?
10. What is CD3D a marker for? What about CD79A?

### Coding Skills Checklist

- [ ] Launch Jupyter Notebook from the `immunology` conda environment
- [ ] Import scanpy, matplotlib, and numpy
- [ ] Load the pbmc3k dataset with `sc.datasets.pbmc3k()`
- [ ] Check the number of cells and genes in an AnnData object
- [ ] Calculate basic statistics (mean, min, max) on the data
- [ ] Create a histogram with matplotlib
- [ ] Save a figure to a PNG file
- [ ] Extract expression values for a specific gene

If you checked every box: outstanding work. You now understand the immune system AND can
write Python code to explore real biological data. That is a combination of skills that
very few high school students in the world possess.

---

## 11. What's Coming Next

In **Week 2: Single-Cell Analysis Pipeline**, you will take the pbmc3k dataset and run a
full analysis workflow:

- **Quality control**: filtering out dead cells and doublets
- **Normalization**: making expression values comparable across cells
- **Dimensionality reduction**: compressing 32,000 gene dimensions into 2D for
  visualization (PCA, UMAP)
- **Clustering**: grouping similar cells together automatically
- **Cell type annotation**: figuring out which cluster is T cells, which is B cells, etc.

By the end of Week 2, you will produce a UMAP plot -- the iconic visualization of
single-cell biology -- with colored clusters of different immune cell types. See you there.

---

## Key Papers (Optional Reading)

These are landmark papers in the field. You are not expected to read them in full right now,
but knowing they exist will help you understand the foundations of the technology and data
we will use throughout this bootcamp.

1. **Macosko EZ, Basu A, Satija R, et al.** "Highly Parallel Genome-wide Expression
   Profiling of Individual Cells Using Nanoliter Droplets." *Cell*, 2015;161(5):1202-1214.
   DOI: 10.1016/j.cell.2015.05.002

   *Why it matters:* This paper introduced **Drop-seq**, one of the first droplet-based
   single-cell RNA-seq methods. It demonstrated that thousands of cells could be profiled
   in a single experiment at low cost. The principles described here underlie the 10x
   Genomics platform you will use throughout this bootcamp.

2. **Dominguez Conde C, Xu C, Jarvis LB, et al.** "Cross-tissue immune cell analysis
   reveals tissue-specific features in humans." *Science*, 2022;376(6594):eabl5197.
   DOI: 10.1126/science.abl5197

   *Why it matters:* This paper used single-cell RNA-seq to profile immune cells across
   multiple human tissues, creating a comprehensive **cross-tissue immune cell atlas**. It
   showed how immune cells differ depending on where they reside in the body -- a key
   insight for understanding tissue-specific immunity.
