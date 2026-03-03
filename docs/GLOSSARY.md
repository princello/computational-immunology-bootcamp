# Glossary

A plain-language dictionary of every technical term used in this bootcamp. Terms are grouped by topic and listed alphabetically within each group. If a term is not here, try searching the web or asking an AI assistant.

---

## General Biology

| Term | Definition |
|------|-----------|
| **Amino acid** | A building block of proteins. There are 20 standard amino acids, each represented by a single letter (e.g., A = Alanine, C = Cysteine). A protein is a chain of amino acids. |
| **Antigen** | Any molecule that the immune system can recognize and respond to. Often a piece of a virus, bacterium, or abnormal cell. Think of it as the "wanted poster" the immune system is looking for. |
| **Cell** | The basic unit of life. Your body has about 37 trillion cells, each performing specific jobs. |
| **Central Dogma** | The flow of genetic information: DNA is copied into RNA, which is used to make Protein. DNA -> RNA -> Protein. |
| **Chromosome** | A long molecule of DNA packaged tightly. Humans have 23 pairs (46 total). |
| **DNA** | Deoxyribonucleic acid. The instruction manual for building and running an organism. A long molecule made of 4 "letters" (A, T, G, C). |
| **Epitope** | The specific part of an antigen that a receptor (TCR or BCR) binds to. Like the specific notch on a key that a lock grips. |
| **Gene** | A segment of DNA that contains instructions for making a specific protein. Humans have about 20,000 genes. |
| **Gene expression** | The process of turning a gene "on" to produce RNA and ultimately protein. Different cell types express different genes, which is what makes a skin cell different from a blood cell. |
| **Genome** | The complete set of DNA in an organism. The human genome has about 3 billion "letters." |
| **Protein** | A large molecule made of amino acids that performs work in the cell. Enzymes, receptors, antibodies, and structural components are all proteins. |
| **RNA** | Ribonucleic acid. A temporary copy of a gene's instructions. mRNA (messenger RNA) carries the instructions from DNA to the protein-making machinery. |
| **Transcription** | The process of copying DNA into RNA. |
| **Translation** | The process of reading RNA to build a protein. |

---

## Immunology

| Term | Definition |
|------|-----------|
| **Adaptive immunity** | The branch of the immune system that learns and remembers specific threats. Involves T cells and B cells. Takes days to activate but is highly targeted. |
| **Antibody** | A Y-shaped protein made by B cells that binds to a specific antigen. Also called immunoglobulin (Ig). Each antibody binds one specific target. |
| **Antigen presentation** | The process by which cells display pieces of proteins on their surface using MHC molecules, allowing T cells to inspect them. |
| **B cell** | A type of white blood cell that makes antibodies. Born in the Bone marrow. Each B cell produces antibodies with one specific shape. |
| **BCR (B Cell Receptor)** | The receptor on a B cell's surface. It is essentially a membrane-bound antibody. Made of a heavy chain and a light chain. |
| **CD4** | A surface protein found on helper T cells. Used as a marker to identify them. |
| **CD8** | A surface protein found on killer T cells. Used as a marker to identify them. |
| **CDR3** | Complementarity Determining Region 3. The most variable part of a TCR or BCR, responsible for directly contacting the antigen. It is a short amino acid sequence (typically 10-20 characters) and is the "fingerprint" of each receptor. |
| **Clonal expansion** | When an immune cell recognizes a threat, it multiplies rapidly to create an army of identical copies. All copies share the same receptor. |
| **Clonotype** | A group of cells that share the same receptor sequence (same TCR or BCR). They are all descendants of one original cell. |
| **Cytokine** | Small signaling proteins that immune cells use to communicate with each other. Like text messages between cells. |
| **Dendritic cell** | A cell that captures invaders, breaks them into pieces, and presents those pieces to T cells. A key "messenger" between innate and adaptive immunity. |
| **Diversity (repertoire)** | How varied the collection of immune receptors is in a sample. High diversity means many different clonotypes; low diversity means a few clonotypes dominate. |
| **Helper T cell (CD4+)** | A T cell that coordinates the immune response by releasing cytokines to activate other immune cells. The "commander" of the immune army. |
| **Innate immunity** | The branch of the immune system you are born with. Provides immediate, general defense (skin, mucus, inflammation). Does not remember specific threats. |
| **Killer T cell (CD8+)** | A T cell that directly destroys infected or cancerous cells. The "assassin" of the immune system. |
| **Macrophage** | A large cell that engulfs and digests invaders and dead cells. Also presents antigens to T cells. |
| **MHC (Major Histocompatibility Complex)** | Proteins on the cell surface that display antigen fragments. MHC class I is on all cells (for CD8 T cells). MHC class II is on special immune cells (for CD4 T cells). |
| **NK cell (Natural Killer)** | An innate immune cell that kills virus-infected cells and tumor cells without needing prior activation. |
| **Neutrophil** | The most abundant white blood cell. First responder to infection. Part of innate immunity. |
| **Repertoire** | The complete collection of all unique TCRs or BCRs in a person or sample. |
| **T cell** | A type of white blood cell that is central to adaptive immunity. Born in the bone marrow, matures in the Thymus. Main types: helper (CD4+) and killer (CD8+). |
| **TCR (T Cell Receptor)** | The receptor on a T cell's surface that recognizes antigens presented by MHC molecules. Made of an alpha chain and a beta chain. |
| **Thymus** | An organ behind the breastbone where T cells mature and are "trained" to recognize foreign antigens but not attack your own body. |
| **VDJ recombination** | The process by which immune cells shuffle and combine V (Variable), D (Diversity), and J (Joining) gene segments to create unique receptors. This generates the enormous diversity of TCRs and BCRs. |

---

## Single-Cell Biology and Sequencing

| Term | Definition |
|------|-----------|
| **10x Genomics** | A company that makes the most widely used single-cell sequencing technology. Uses tiny gel beads in droplets to capture individual cells. |
| **AnnData** | The data container used by Scanpy (Python). Stores the gene expression matrix plus cell and gene metadata. File extension: `.h5ad`. |
| **Barcode** | A short DNA sequence that uniquely labels each cell in a single-cell experiment. Like a name tag for each cell. |
| **Batch effect** | Technical differences between experiments that are not biological. For example, cells processed on Monday vs. Tuesday might look different due to reagent variability, not biology. |
| **Bulk RNA-seq** | Measuring gene expression from a mixture of many cells at once. You get an average across all cells, losing information about individual cell differences. |
| **Cell annotation** | The process of labeling cells with their identity (e.g., "this cluster is CD4 T cells"). Done using marker genes or automated tools. |
| **CITE-seq** | A technology that measures both RNA and surface proteins from the same cell. Uses antibodies tagged with DNA barcodes. |
| **Clustering** | Grouping similar cells together based on their gene expression patterns. Like sorting a pile of coins by type. |
| **Count matrix** | A table where rows are cells, columns are genes, and values are how many RNA molecules of each gene were detected in each cell. |
| **Dimensionality reduction** | Compressing high-dimensional data (thousands of genes) into fewer dimensions (e.g., 2D for visualization). PCA and UMAP are dimensionality reduction methods. |
| **Doublet** | When two cells are accidentally captured in the same droplet, creating a mixed profile that looks like a single cell. These are usually filtered out. |
| **Dropout** | When a gene is actually active in a cell but was not detected due to technical limitations. Common in single-cell data because each cell has very little RNA. |
| **Feature selection** | Choosing which genes to use for analysis. Typically we select "highly variable genes" -- genes whose expression differs the most between cells. |
| **Highly variable genes (HVGs)** | Genes that show the most variation across cells. These are the most informative for distinguishing cell types. |
| **Leiden algorithm** | A method for finding clusters (communities) in a network of cells. Commonly used in single-cell analysis. |
| **Marker gene** | A gene whose expression is characteristic of a specific cell type. For example, CD3D is a marker for T cells. |
| **Mitochondrial genes** | Genes encoded in the mitochondria (the cell's energy factories). High percentages of mitochondrial RNA usually indicate dying cells. In human data, these gene names start with "MT-". |
| **Normalization** | Adjusting raw counts so that cells with different total RNA amounts can be fairly compared. Like converting test scores to percentages. |
| **PCA (Principal Component Analysis)** | A dimensionality reduction method that finds the main directions of variation in the data. Reduces thousands of genes to ~50 principal components. |
| **QC (Quality Control)** | Filtering out low-quality cells (too few genes, too much mitochondrial RNA) and low-quality genes (detected in too few cells). |
| **scATAC-seq** | Single-cell Assay for Transposase-Accessible Chromatin sequencing. Measures which regions of DNA are "open" (accessible) in each cell. |
| **scRNA-seq** | Single-cell RNA sequencing. Measuring gene expression in individual cells rather than a bulk average. |
| **Seurat** | A popular R package for single-cell analysis (similar to Scanpy for Python). Developed by the Satija Lab at the New York Genome Center. |
| **t-SNE** | t-distributed Stochastic Neighbor Embedding. An older dimensionality reduction method for visualization, largely replaced by UMAP. |
| **UMAP** | Uniform Manifold Approximation and Projection. A method for compressing high-dimensional data into 2D for visualization while preserving the structure. Similar cells appear close together on a UMAP plot. |
| **UMI (Unique Molecular Identifier)** | A random DNA barcode attached to each RNA molecule before amplification. Allows us to count the actual number of original molecules rather than amplified copies. |

---

## Multi-omics and Integration

| Term | Definition |
|------|-----------|
| **ATAC-seq** | Assay for Transposase-Accessible Chromatin. Measures which parts of the DNA are "open" and potentially active. |
| **Chromatin** | The complex of DNA wrapped around histone proteins. "Open" chromatin is accessible for gene expression; "closed" chromatin is silenced. |
| **Chromatin accessibility** | Whether a region of DNA is physically accessible to the cell's machinery. Open = potentially active genes. Closed = silenced genes. |
| **Co-embedding** | Placing cells from different data types into the same coordinate space so they can be compared visually. |
| **Gene activity (ATAC)** | An estimate of gene expression derived from ATAC-seq data by looking at chromatin accessibility near each gene. |
| **Integration** | Combining data from multiple sources or modalities into a unified analysis. |
| **Label transfer** | Using cell type labels from one dataset (e.g., RNA) to annotate cells in another dataset (e.g., ATAC). |
| **LSI (Latent Semantic Indexing)** | A dimensionality reduction method used for ATAC-seq data, similar to PCA but designed for binary/sparse data. |
| **Modality** | A type of measurement. RNA is one modality, protein is another, ATAC is another. Multi-modal means measuring multiple types. |
| **MuData** | A data container (Python) for storing multi-modal data. Like AnnData but with multiple layers for different modalities. |
| **Multiome** | A 10x Genomics product that measures both RNA and ATAC from the same cell in the same experiment. |
| **Peak (ATAC)** | A genomic region with significantly high chromatin accessibility, indicating potential regulatory activity. |
| **Peak-to-gene link** | A statistical association between an accessible chromatin region (peak) and a gene's expression, suggesting the peak regulates the gene. |
| **TF (Transcription Factor)** | A protein that binds to specific DNA sequences and controls gene expression. Like a switch that turns genes on or off. |
| **TF-IDF** | Term Frequency-Inverse Document Frequency. A normalization method used for ATAC-seq data, borrowed from text analysis. |
| **TF motif** | A short DNA sequence pattern that a transcription factor recognizes and binds to. Finding these motifs in open chromatin reveals which TFs may be active. |
| **totalVI** | A deep learning model (from scvi-tools) that jointly models RNA and protein data from CITE-seq experiments. |
| **WNN (Weighted Nearest Neighbors)** | An integration method (from Seurat) that combines multiple modalities by weighting each one differently for each cell. |

---

## Machine Learning and Statistics

| Term | Definition |
|------|-----------|
| **AUROC (Area Under ROC Curve)** | A metric for classification performance. Ranges from 0.5 (random guessing) to 1.0 (perfect). Measures how well the model distinguishes positives from negatives overall. |
| **Average Precision (AP)** | A metric that summarizes the precision-recall curve. More informative than AUROC when classes are imbalanced (many more negatives than positives). |
| **Calibration** | Whether predicted probabilities match actual frequencies. If a model says "80% chance of binding," it should be correct about 80% of the time. |
| **Classification** | Predicting which category something belongs to (e.g., "binds" vs "does not bind"). |
| **Confusion matrix** | A table showing how many predictions were correct (true positives, true negatives) and incorrect (false positives, false negatives). |
| **Cross-validation** | Testing a model by repeatedly splitting data into training and test sets, training on one and testing on the other, to get a reliable performance estimate. |
| **Data leakage** | When information from the test set accidentally leaks into training, making the model appear better than it really is. A major pitfall in TCR prediction. |
| **Embedding** | A numerical representation of data (like a sequence) in a continuous vector space. Similar items have similar embeddings. |
| **Feature** | A measurable property used as input to a model. For TCR prediction, features might include amino acid composition, k-mer counts, or biophysical properties. |
| **Feature engineering** | The process of creating informative features from raw data. Often the most impactful step in building a good model. |
| **Feature importance** | A score indicating how much each feature contributes to the model's predictions. Helps understand what the model learned. |
| **K-mer** | A subsequence of length k. For example, the 3-mers of "CASSF" are "CAS", "ASS", "SSF". Used to turn sequences into numerical features. |
| **Logistic Regression** | A simple classification model that finds a linear boundary between classes. Often a good baseline. |
| **One-hot encoding** | Representing a categorical variable as a binary vector with a 1 in one position and 0s everywhere else. For amino acids: A=[1,0,0,...], C=[0,1,0,...], etc. |
| **Overfitting** | When a model memorizes the training data instead of learning general patterns. Performs well on training data but poorly on new data. |
| **Pipeline** | A sequence of data processing steps chained together. In scikit-learn, pipelines ensure the same transformations are applied consistently. |
| **Precision** | Of all the items the model predicted as positive, what fraction actually were positive? |
| **Random Forest** | An ensemble model that builds many decision trees and averages their predictions. Usually more accurate than a single tree. |
| **Recall** | Of all the actual positives, what fraction did the model correctly identify? |
| **ROC curve** | Receiver Operating Characteristic curve. Plots true positive rate vs false positive rate at different thresholds. |
| **SHAP values** | A method for explaining model predictions by showing how much each feature contributed to each prediction. |
| **Stratified split** | Dividing data into train/test sets while maintaining the same proportion of positive/negative examples in each. |
| **XGBoost** | A powerful gradient-boosted tree model. Often wins machine learning competitions. |

---

## Programming and Tools

| Term | Definition |
|------|-----------|
| **Conda** | A package and environment manager. Installs software and creates isolated environments so different projects don't conflict. |
| **Environment (conda)** | An isolated workspace with its own installed packages. Different projects can use different environments. |
| **Git** | A version control system that tracks changes to your code over time. Like "Track Changes" in Word but much more powerful. |
| **GitHub** | A website for hosting Git repositories. Used to share code, collaborate, and showcase projects. |
| **Jupyter Notebook** | An interactive coding environment that runs in your web browser. You write code in "cells" and see results immediately. Files end in `.ipynb`. |
| **Matplotlib** | A Python library for creating plots and visualizations. |
| **NumPy** | A Python library for numerical computing. Provides fast array operations. |
| **Pandas** | A Python library for data manipulation. Provides DataFrames (like spreadsheets in Python). |
| **pip** | A Python package installer. Like an app store for Python libraries. |
| **R** | A programming language widely used in statistics and bioinformatics. |
| **RStudio** | A user-friendly editor (IDE) for writing R code. |
| **Scanpy** | A Python package for analyzing single-cell data. The Python equivalent of Seurat. |
| **scikit-learn** | A Python library for machine learning. Provides models, metrics, and data processing tools. |
| **Scirpy** | A Python package for analyzing immune receptor (TCR/BCR) data. Works with Scanpy. |
| **scvi-tools** | A collection of deep learning models for single-cell data analysis. Includes scVI, totalVI, MultiVI, and more. |
| **Seaborn** | A Python library for statistical data visualization. Built on top of Matplotlib with nicer defaults. |
| **Terminal** | A text-based interface for interacting with your computer by typing commands. Also called command line, shell, or console. |

---

*Cannot find a term? Search the web for "[term] bioinformatics definition" or ask an AI assistant to explain it.*
