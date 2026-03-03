# computational immunology bootcamp

# 🧬 Computational Immunology Bootcamp

Welcome to the **4-week introduction to Computational Immunology**.

This bootcamp blends **immunology basics** with **hands-on computational pipelines** for single-cell, TCR/BCR, and multi-omics analysis.

---

## 📅 Weekly Pages

👉 [Week 1 – Setup & Immunology Foundations](about:blank#-week-1--setup--immunology-foundations)

👉 [Week 2 – Single-cell RNA-seq Pipeline](about:blank#-week-2--single-cell-rna-seq-pipeline)

👉 [Week 3 – TCR & BCR Analysis](about:blank#-week-3--tcr--bcr-analysis)

👉 [Week 4 – Multi-omics Integration](about:blank#-week-4--multi-omics-integration)

---

## 📝 Notes & Assignments

Use this space for your **progress tracking** and **custom notes**.

---

# 📄 Week 1 – Setup & Immunology Foundations

**Goal:** Build the right environment and learn basic concepts.

### ✅ Checklist

- [x]  Install Python (Anaconda/Miniconda, Jupyter Notebook)
    - [x]  https://www.python.org/downloads/
    - [x]  https://www.anaconda.com/docs/getting-started/miniconda/main
    - [x]  https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
- [x]  Install R + RStudio
    - [ ]  https://posit.co/download/rstudio-desktop/
- [x]  Connect GitHub & Google Drive
    
    **GitHub** Username: TwngTadpole
    
- [x]  Set up Notion workspace
- [x]  https://phpandmysql.com/extras/command-line-mac/?utm_source=notesre&utm_medium=link
- [x]  Review T cells & B cells basics
- [x]  Understand antigen presentation (MHC I & II)

### 📂 Topics

🔹 Tools setup: conda, pip, renv

🔹 Immunology basics (T vs B cells, antibodies, adaptive vs innate)

### 📖 Resources

- Janeway’s *Immunobiology* (Ch.1–3)
- [Scanpy tutorial](https://scanpy.readthedocs.io/)
- [Seurat tutorial](https://satijalab.org/seurat/)
- Here’s a representative “landmark” paper
    
    ### **Single-cell RNA-seq**
    
    - [x]  **Macosko et al., 2015.**
        
        *Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets.*
        
        *Cell* 161(5):1202–1214.
        
        → One of the pioneering papers introducing **Drop-seq**, which made high-throughput scRNA-seq feasible.
        
    
    ### Immune cell categories
    
    - [x]  **Domínguez Conde et al., 2022.**
    Cross-tissue immune cell analysis reveals tissue-specific features in humans.
    Science 376:eabl5197.
        
        → A comprehensive atlas defining immune cell categories across multiple human tissues using scRNA-seq, establishing a reference for tissue-specific immune variation.
        
    
    ### **TCR (T-cell receptor)**
    
    - [ ]  **Dash et al., 2017.**
        
        *Quantifiable predictive features define epitope-specific T cell receptor repertoires.*
        
        *Nature* 547:89–93.
        
        → A landmark in **TCR repertoire analysis**, introducing predictive models to link TCR sequences with antigen specificity.
        
    
    ### **BCR (B-cell receptor / antibodies)**
    
    - [ ]  **Briney et al., 2019.**
        
        *Commonality despite exceptional diversity in the baseline human antibody repertoire.*
        
        *Nature* 566:393–397.
        
        → A large-scale study of **human BCR repertoires**, revealing both the immense diversity and recurring “public” clonotypes.
        

---

# 📄 Week 2 – Single-cell RNA-seq Pipeline

**Goal:** Learn how to process and analyze scRNA-seq data.

### ✅ Checklist

- [x]  Study 10x Genomics and summarize various strategies along with application examples.
- [x]  QC, normalization, clustering (Seurat/Scanpy)
- [x]  Perform PCA, UMAP, t-SNE
- [x]  Cell Annotation
- [ ]  Make heatmaps and violin plots

### 💎 E**xploring**

- [x]  probabilistic models for single-cell omics (scvi)

### 📂 Topics

🔹 Raw data preprocessing

🔹 Dimensional reduction & clustering

🔹 Visualization

### 📖 Resources

- https://www.sc-best-practices.org/preamble.html 1-21
- https://scanpy-tutorials.readthedocs.io/en/latest/
- https://www.celltypist.org/
- https://scvi-tools.org/
- https://www.kaggle.com/
- https://scikit-learn.org/stable/
- [ ]  [O’Donnell, Timothy J., et al. "Reading the repertoire: Progress in adaptive immune receptor analysis using machine learning." *Cell Systems* 15.12 (2024): 1168-1189.](https://www.cell.com/cell-systems/fulltext/S2405-4712(24)00342-9)
- [ ]  [Weber, Anna, Aurélien Pélissier, and María Rodríguez Martínez. "T-cell receptor binding prediction: A machine learning revolution." *ImmunoInformatics* 15 (2024): 100040.](https://www.sciencedirect.com/science/article/pii/S2667119024000107)

---

# 📄 Week 3 – TCR & BCR Analysis

**Goal:** Explore immune receptor data from single-cell sequencing.

### ✅ Checklist

- [ ]  Review TCR/BCR repertoire concepts (CDR3, clonotypes)
- [x]  Install scRepertoire, Immunarch, scirpy
- [x]  Run clonotype expansion analysis
- [ ]  Calculate diversity metrics
- [x]  Map clonotypes to scRNA-seq clusters

### 📂 Topics

🔹 TCR/BCR repertoire basics

🔹 Clonal expansion

🔹 Diversity & overlap metrics

### 📖 Resources

- [VDJdb](https://vdjdb.cdr3.net/)
- [ImmuneCODE](https://clients.adaptivebiotech.com/pub/covid-2020)
- [TCR3d](https://tcr3d.ibbr.umd.edu/) - UMD
- [TCRdb](https://guolab.wchscu.cn/TCRdb//#/)
- [huARdb](https://huarc.net/v2/database/)
- [OTS](https://opig.stats.ox.ac.uk/webapps/ots)
- [OGRDB](https://ogrdb.airr-community.org/)
- [IEDB](https://www.iedb.org/)

---

# 📄 Week 4 – Multi-omics Integration

**Goal:** Combine modalities for deeper insights.

### ✅ Checklist

- [x]  Review data types (scRNA-seq, TCR/BCR, CITE-seq, ATAC-seq)
- [x]  Install Seurat v5, scvi-tools, ArchR
- [ ]  Perform RNA + protein (CITE-seq) integration
- [ ]  Explore RNA + ATAC integration with ArchR
- [ ]  Case study: integrated immune profiling

### 📂 Topics

🔹 Multi-modal datasets

🔹 Integration frameworks (Seurat, TotalVI, ArchR)

🔹 Case study on immune profiling

### 📖 Resources

- [Seurat multimodal tutorial](https://satijalab.org/seurat/articles/multimodal_reference_mapping.html)
- [scvi-tools docs](https://docs.scvi-tools.org/)
- [ArchR tutorial](https://www.archrproject.com/bookdown/)
- [ ]  [**Steier et al., Nat. Immunol. (2023)](https://urldefense.proofpoint.com/v2/url?u=https-3A__www.nature.com_articles_s41590-2D023-2D01584-2D0&d=DwMGaQ&c=009klHSCxuh5AI1vNQzSO0KGjl4nbi2Q0M1QLJX9BeE&r=dpezMWuR27UYbPumyz_bX_4fXZjeHM0Gll4cQC3py-c&m=jjJJodVx_9-JuQ3YWh4ImioftUyLNK9W3oPjTj5xlLBFJ4_HGq2oU83MKskBnOFQ&s=f8ZzhF0V3CPL68HEhKTCmGGUMOLGwm1VCeeB4WgNlFQ&e=):** RNA-seq with CITE-seq shows thymocytes first initiate CD4+ T cell differentiation, with a second TCR signal determining CD8+ fate via the calcineurin–NFAT–GATA3 pathway.
- [ ]  [**Park et al., Science (2020)](https://urldefense.proofpoint.com/v2/url?u=https-3A__www.science.org_doi_10.1126_science.aay3224&d=DwMGaQ&c=009klHSCxuh5AI1vNQzSO0KGjl4nbi2Q0M1QLJX9BeE&r=dpezMWuR27UYbPumyz_bX_4fXZjeHM0Gll4cQC3py-c&m=jjJJodVx_9-JuQ3YWh4ImioftUyLNK9W3oPjTj5xlLBFJ4_HGq2oU83MKskBnOFQ&s=xLtjF2r4oPi7plw1WrOMQ0u8ofZK5PITQ5pNQkByO-I&e=):** A single-cell RNA-seq atlas of the human thymus across the lifespan, revealing new cell types, TCR recombination biases, and lineage commitment kinetics.
- [ ]  [**Chopp et al., Immunity (2020)](https://urldefense.proofpoint.com/v2/url?u=https-3A__www.sciencedirect.com_science_article_pii_S1074761320304659-3Fvia-253Dihub&d=DwMGaQ&c=009klHSCxuh5AI1vNQzSO0KGjl4nbi2Q0M1QLJX9BeE&r=dpezMWuR27UYbPumyz_bX_4fXZjeHM0Gll4cQC3py-c&m=jjJJodVx_9-JuQ3YWh4ImioftUyLNK9W3oPjTj5xlLBFJ4_HGq2oU83MKskBnOFQ&s=7uJlxRCRQ1Pn8g0LT94x3lkgb7zOHG_PdckG5o6VVmE&e=):** Single-cell transcriptomic and epigenomic profiling uncover asymmetric CD4+/CD8+ lineage emergence and transcriptional networks guiding αβ T cell differentiation.

---

# 🎓 Final Practice – Computational Immunology Bootcamp

A capstone-style, **2-track final practice** that synthesizes Weeks 1–4. You’ll produce a short report, clean notebooks, and shareable figures. Aim for **reproducible code** and **clear reasoning**.

---

## 🧪 Track 1 — Recreate a "Chopp-style" result: RNA + ATAC Integration

**Goal:** Reproduce the core idea of **joint scRNA-seq + scATAC-seq analysis** (as in Chopp et al., *Immunity* 2020) using a public dataset. If you don’t have a thymus dataset, use a well-curated **10x Multiome PBMC** dataset as a stand-in to validate the full pipeline.

### A. Data options

1. **10x Multiome PBMC 10k** (RNA+ATAC in the same cells) – great for full-stack integration
2. **Separate scRNA + scATAC** PBMC datasets – demonstrate label transfer and cross-modality mapping

> If you have access to thymus RNA/ATAC, you can swap in your own data. Otherwise, PBMC is acceptable for method reproducibility.
> 

### B. Analysis milestones (with suggested tools)

1. **QC**
    - RNA: cells × genes, MT% / ribo%, doublet filtering
    - ATAC: TSS enrichment, FRiP, nucleosome signal, blacklist ratio
2. **Normalization & Latent Space**
    - RNA: SCTransform (Seurat) or `sc.pp.normalize_total` + `sc.pp.log1p` (Scanpy)
    - ATAC: TF-IDF + SVD (Signac) or ArchR’s LSI
3. **Integration / Co-embedding**
    - Option A (same cells): **WNN (Seurat v5)** or **totalVI (scvi-tools)**
    - Option B (separate assays): map ATAC→RNA via label transfer using **gene activity** & anchors (Signac/Seurat) or **bridge integration** (scvi-tools + mde)
4. **Cell Type Annotation**
    - Reference mapping (Azimuth/CellTypist/SingleR) + marker-based curation
5. **Regulatory Inference (ATAC)**
    - **Gene activity**: peaks→genes
    - **TF motif/footprinting**: chromVAR / ArchR motif enrichments
    - **Peak2Gene links**: co-accessibility / correlation
6. **Biological Recap**
    - Show lineage structure (e.g., naïve→memory T, B cell maturation, monocyte states)
    - Summarize correspondence between RNA states and accessible chromatin programs

### C. Required outputs

- UMAP(s) colored by modality, cluster, cell type, and batch
- QC panels (per modality)
- Gene activity heatmap for key markers
- TF activity (chromVAR deviations) violin/heatmap
- Peak2Gene link examples for 2–3 canonical markers
- 1–2 **figure panels** that mirror the logic of Chopp et al.: transcriptome-defined states supported by chromatin programs and TF networks

### E. Starter snippets

**Scanpy + scvi-tools (Python)**

```python
import scanpy as sc, scvi
import numpy as np
sc.settings.verbosity = 2

adata = sc.read_10x_mtx("data/raw/pbmc_rna", var_names="gene_symbols", cache=True)
sc.pp.filter_cells(adata, min_genes=500)
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
adata = adata[adata.obs.pct_counts_mt < 15].copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
adata = adata[:, adata.var.highly_variable].copy()
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)

# If multiome (paired RNA+ATAC), consider totalVI
# ad_rna, ad_atac = ... (AnnData objects)
# scvi.model.TOTALVI.setup_anndata(ad_total, batch_key="batch")
# model = scvi.model.TOTALVI(ad_total)
# model.train()
# ad_total.obsm["X_totalvi"] = model.get_latent_representation()

```

---

## 🧬 Track 2 — TCR Specificity Prediction (Kaggle challenge)

**Goal:** Build a baseline→strong model to predict **TCR–epitope specificity**.

### A. Problem framing

- **Input:** TCR sequences (typically CDR3β, sometimes α) + **epitope peptide**
- **Output:** Binary label (binds / does not bind) or probability
- **Metric:** Prioritize **AUROC** and **Average Precision (AP)**; report both

### B. Featurization options

1. **Amino-acid one-hot** (fixed-length with padding/truncation)
2. **K-mer counts** (e.g., k=2–5) → TF-IDF or frequency
3. **Biophysicochemical descriptors** per AA (e.g., Atchley factors, hydropathy, charge)
4. **Substitution matrices** (BLOSUM62) → pairwise/alignment features
5. **Learned embeddings** (if allowed offline): ProtVec, ESM2 small; or your **G2VTCR/graph** approach

> Start simple (1–3), then add learned embeddings if you have them locally.
> 

### C. Baseline→strong model ladder

- **Baseline (fast):** Logistic Regression / Linear SVM on k-mers + AA descriptors
- **Medium:** XGBoost / Random Forest on engineered features
- **Deep:** 1D CNN or small Transformer on sequence tokens (masking + positional encodings)
- **Pair modeling:** Concatenate TCR and epitope embeddings; optionally use **cross-attention**

### D. Cross-validation and splits

- **Avoid leakage**: split by **epitope** or **clone** so identical/near-identical sequences don’t appear in both train/test
- Use **stratified** folds; report mean±sd across folds

### E. Minimal baseline (Python sketch)

```python
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score

# 1) Load your TSV/CSV with columns: cdr3, epitope, label (0/1)
df = pd.read_csv("data/raw/tcr/train.csv")

# 2) Simple joint tokenization: TCR + epitope string
def concat_pair(r):
    return f"T:{r.cdr3} E:{r.epitope}"
text = df.apply(concat_pair, axis=1)
y = df["label"].values

# 3) k-mer vectorizer (k=3)
vec = CountVectorizer(analyzer='char', ngram_range=(3,3))

# 4) Model pipeline
pipe = Pipeline([
    ("vec", vec),
    ("scaler", StandardScaler(with_mean=False)),
    ("clf", LogisticRegression(max_iter=200))
])

# 5) CV
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
aucs, aps = [], []
for tr, te in skf.split(text, y):
    pipe.fit(text.iloc[tr], y[tr])
    p = pipe.predict_proba(text.iloc[te])[:,1]
    aucs.append(roc_auc_score(y[te], p))
    aps.append(average_precision_score(y[te], p))
print(f"AUROC {np.mean(aucs):.3f}±{np.std(aucs):.3f} | AP {np.mean(aps):.3f}±{np.std(aps):.3f}")

```

### F. Stretch goals

- Epitope-conditioned **attention** (TCR tokens attend to epitope tokens)
- **Contrastive pretraining** (positive pairs vs negatives)
- **Calibrated probabilities** (Platt scaling / isotonic)

### G. Reporting

- Curves: ROC, PR; calibration plot (reliability diagram)
- Tables: performance per epitope, per CDR3 length bucket
- Feature inspection: top k-mers, SHAP values for tree models

---

## 🗂 Suggested folder structure

```
project/
  data/
    raw/
    processed/
  notebooks/
    track1_scRNA_ATAC_integration.ipynb
    track2_TCR_specificity_prediction.ipynb
  src/
    io.py
    metrics.py
    models/
  env/
    environment.yml
    requirements.txt
  figs/
  reports/
    final_practice_report.pdf
  README.md

```

---

## 📊 Minimum figure set (both tracks)

- Track 1: RNA/ATAC QC panels; co-embedding UMAP; cell type annotation; TF activity; peak2gene links
- Track 2: ROC & PR curves; calibration; confusion matrix at optimal threshold; top features

---

## ✅ Final report template (headings)

1. **Background & Objectives** (why these questions matter)
2. **Data** (sources, sizes, QC criteria)
3. **Methods** (processing, models, parameters)
4. **Results** (primary figures/tables)
5. **Discussion** (biological insights; limitations)
6. **Reproducibility** (env, seeds, data version)
7. **Next Steps** (what you’d do with more time)

---

## 📝 Hints & Pitfalls

- **Integration:** always verify batch effects; check that biology (not batch) drives the embedding. Use label transfer accuracy and ARI/NMI sanity checks.
- **ATAC:** gene activity is noisy; use multiple evidence lines (motifs, co-accessibility, peak2gene).
- **TCR splits:** prevent epitope/clone leakage; otherwise metrics will be inflated.
- **Metrics:** report AUROC & AUPR + confidence intervals or fold variability.
- **Interpretability:** simple baselines + ablations make your argument stronger.

Good luck—have fun tying together everything from the bootcamp! 🚀