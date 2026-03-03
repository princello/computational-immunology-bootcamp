# Week 0: Setting Up Your Computer

**Time needed:** 2-4 hours
**Goal:** Install all the software you need and verify everything works.

By the end of this week, your computer will be a real bioinformatics workstation -- the same tools that researchers at top universities use every day.

---

## Table of Contents

1. [Meet the Terminal](#1-meet-the-terminal)
2. [Install Miniconda (Python)](#2-install-miniconda-python)
3. [Create Your Bootcamp Environment](#3-create-your-bootcamp-environment)
4. [Install R and RStudio](#4-install-r-and-rstudio)
5. [Set Up Git and GitHub](#5-set-up-git-and-github)
6. [Install Jupyter Notebook](#6-install-jupyter-notebook)
7. [Self-Check](#7-self-check)

---

## 1. Meet the Terminal

The **terminal** (also called "command line" or "shell") is a text-based way to talk to your computer. Instead of clicking icons, you type commands. It looks intimidating at first, but you only need a handful of commands to get started.

### How to Open the Terminal

**Mac:**
1. Press `Cmd + Space` to open Spotlight Search
2. Type `Terminal` and press Enter
3. A window with a blinking cursor appears -- that is your terminal

**Windows:**
1. We will use "Anaconda Prompt" after installing Miniconda (Step 2)
2. Alternatively, install [Windows Terminal](https://apps.microsoft.com/store/detail/windows-terminal/9N0DX20HK701) from the Microsoft Store

**Linux:**
1. Press `Ctrl + Alt + T`

### Essential Terminal Commands

Practice these now. Type each one and press Enter:

```bash
# Print your current location (which folder you are in)
pwd

# List files in the current folder
ls

# Make a new folder called "bootcamp"
mkdir ~/bootcamp

# Move into that folder
cd ~/bootcamp

# Confirm you are in the right place
pwd
# Expected output: something ending in /bootcamp
```

**What the symbols mean:**
- `~` means your home folder (e.g., `/Users/yourname` on Mac)
- `cd` means "change directory" (directory = folder)
- `pwd` means "print working directory" (show where I am)
- `ls` means "list" (show what is here)
- `mkdir` means "make directory" (create a folder)

> **Checkpoint 0.1:** Type `pwd` in your terminal. You should see a path printed. If you see a path, you are ready to continue.

---

## 2. Install Miniconda (Python)

**Miniconda** is a small version of Anaconda. It gives you Python plus a tool called `conda` that manages software packages and environments. Think of it as an app store for scientific computing.

### Why Miniconda and not just Python?

Bioinformatics tools have complex dependencies (other software they need to run). Conda handles all of that automatically so you don't have to.

### Installation Steps

**Mac:**

```bash
# Download Miniconda for Mac
# For Apple Silicon (M1/M2/M3/M4 Macs, most Macs made after 2020):
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh

# For older Intel Macs:
# curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

# Run the installer
bash Miniconda3-latest-MacOSX-arm64.sh
```

When prompted:
- Press Enter to read the license, then type `yes` to accept
- Press Enter to install to the default location
- Type `yes` when asked to initialize Miniconda

Then **close and reopen your terminal**.

**Windows:**

1. Go to https://docs.conda.io/en/latest/miniconda.html
2. Download the Windows 64-bit installer
3. Run the `.exe` file
4. Check "Add Miniconda3 to my PATH environment variable" during installation
5. Click Install, then close and reopen your terminal

**Linux:**

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Follow the same prompts as Mac
```

### Verify Installation

After reopening your terminal:

```bash
conda --version
```

**Expected output:** Something like `conda 24.x.x` (the exact number does not matter).

```bash
python --version
```

**Expected output:** Something like `Python 3.11.x` or `Python 3.12.x`.

> **Checkpoint 0.2:** Both `conda --version` and `python --version` print version numbers without errors.

---

## 3. Create Your Bootcamp Environment

An **environment** is an isolated workspace where you install packages without affecting the rest of your computer. Think of it like having separate toolboxes for different projects.

### Create the Environment

```bash
# Create a new environment named "immunology" with Python 3.11
conda create -n immunology python=3.11 -y
```

The `-y` flag means "yes, go ahead without asking me." This will take 1-2 minutes.

### Activate the Environment

```bash
conda activate immunology
```

**You should see `(immunology)` appear at the beginning of your terminal prompt.** This means you are now "inside" the immunology environment.

> **IMPORTANT:** Every time you open a new terminal to work on this bootcamp, you must run `conda activate immunology` first. If you see `(base)` instead of `(immunology)`, you are in the wrong environment.

### Install Core Python Packages

```bash
# Install the main bioinformatics packages (this may take 5-10 minutes)
pip install scanpy scirpy scvi-tools celltypist leidenalg

# Install data science essentials
pip install pandas numpy matplotlib seaborn scikit-learn jupyter

# Install additional useful tools
pip install anndata mudata session-info
```

Wait for all installations to complete. You will see a lot of text scrolling by -- that is normal.

### Verify the Installations

```bash
python -c "import scanpy; print('Scanpy version:', scanpy.__version__)"
python -c "import scirpy; print('Scirpy version:', scirpy.__version__)"
python -c "import sklearn; print('Scikit-learn version:', sklearn.__version__)"
```

**Expected output:** Three lines showing version numbers, no errors.

> **Checkpoint 0.3:** All three `python -c "import ..."` commands print version numbers without errors.

---

## 4. Install R and RStudio

**R** is another programming language heavily used in bioinformatics. **RStudio** is a user-friendly editor for R (like a word processor for code).

### Install R

**Mac:**
1. Go to https://cran.r-project.org/bin/macosx/
2. Download the `.pkg` file for your Mac type (Apple Silicon or Intel)
3. Double-click to install, follow the prompts

**Windows:**
1. Go to https://cran.r-project.org/bin/windows/base/
2. Click "Download R for Windows"
3. Run the installer, accept all defaults

**Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install r-base
```

### Install RStudio

1. Go to https://posit.co/download/rstudio-desktop/
2. Download the free version for your operating system
3. Install it

### Install R Packages

Open RStudio (find it in your applications). In the **Console** panel at the bottom, type:

```r
# Install Seurat (the main single-cell analysis package for R)
# This can take 10-20 minutes -- be patient
install.packages("Seurat")

# Install other useful packages
install.packages(c("ggplot2", "dplyr", "tidyverse"))

# Install Bioconductor packages (a collection of biology-focused R tools)
install.packages("BiocManager")
BiocManager::install(c("SingleR", "celldex", "scRepertoire"))
```

Each of these will download and install many sub-packages. When asked "Do you want to install from sources the packages which need compilation?", type `no` and press Enter (the pre-built versions are faster to install).

### Verify R Installation

In the RStudio console:

```r
library(Seurat)
packageVersion("Seurat")
```

**Expected output:** A version number like `5.x.x`.

> **Checkpoint 0.4:** You can open RStudio, and `library(Seurat)` runs without errors.

---

## 5. Set Up Git and GitHub

**Git** tracks changes to your code (like "Track Changes" in Word, but much more powerful). **GitHub** is a website where you store and share your code.

### Install Git

**Mac:** Git is usually pre-installed. Check with:
```bash
git --version
```
If not installed, it will prompt you to install Xcode Command Line Tools. Click "Install."

**Windows:** Download from https://git-scm.com/download/win and install with default settings.

### Create a GitHub Account

1. Go to https://github.com
2. Sign up with your email
3. Choose a username (this will be public -- pick something professional)
4. Verify your email

### Configure Git

In your terminal:

```bash
# Replace with YOUR name and email
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

### Create Your Bootcamp Repository

```bash
# Go to your bootcamp folder
cd ~/bootcamp

# Initialize a git repository
git init

# Create a README file
echo "# My Computational Immunology Bootcamp" > README.md

# Stage and commit your first file
git add README.md
git commit -m "Initial commit: start bootcamp"
```

**Expected output:** Something like `[main (root-commit) abc1234] Initial commit: start bootcamp`.

> **Checkpoint 0.5:** Running `git log` shows your first commit.

---

## 6. Install Jupyter Notebook

**Jupyter Notebook** is an interactive coding environment that runs in your web browser. You write code in "cells" and run them one at a time, seeing results immediately. It is the standard tool for data analysis in bioinformatics.

### Verify Installation

You already installed Jupyter in Step 3. Let's test it:

```bash
# Make sure you are in the right environment
conda activate immunology

# Start Jupyter Notebook
jupyter notebook
```

**What should happen:**
1. Your terminal will print some URLs
2. Your web browser will automatically open to a page showing your files
3. This page is the Jupyter dashboard

### Create Your First Notebook

1. In the Jupyter dashboard, click **New** (top right) then **Python 3**
2. A new tab opens with an empty notebook
3. In the first cell, type:

```python
print("Hello, Computational Immunology!")
2 + 2
```

4. Press `Shift + Enter` to run the cell

**Expected output:**
```
Hello, Computational Immunology!
4
```

5. In the next cell, try:

```python
# Let's make sure our bioinformatics tools are working
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

print(f"Scanpy version: {sc.__version__}")
print(f"NumPy version: {np.__version__}")
print("All imports successful!")
```

6. Press `Shift + Enter`

**Expected output:** Version numbers and "All imports successful!"

7. Save your notebook: `Ctrl + S` (or `Cmd + S` on Mac)
8. Name it `test_setup.ipynb`

### How to Stop Jupyter

- Close the browser tab
- Go back to your terminal and press `Ctrl + C`, then type `y` to confirm

> **Checkpoint 0.6:** You can start Jupyter, create a notebook, run a cell with `import scanpy as sc`, and see the output without errors.

---

## 7. Self-Check

Before moving to Week 1, make sure ALL of these are true:

- [ ] You can open a terminal and navigate with `cd`, `ls`, `pwd`
- [ ] `conda --version` prints a version number
- [ ] `conda activate immunology` changes your prompt to show `(immunology)`
- [ ] `python -c "import scanpy"` runs without error
- [ ] `python -c "import scirpy"` runs without error
- [ ] RStudio opens and `library(Seurat)` works
- [ ] `git --version` prints a version number
- [ ] You have a GitHub account
- [ ] `jupyter notebook` opens in your browser
- [ ] You ran your first notebook cell successfully

**If all boxes are checked: Congratulations! You now have a professional bioinformatics workstation -- the same tools used by researchers at places like Columbia, MIT, and the Broad Institute. You are ready for [Week 1](../week_1_immunology_foundations/README.md).**

---

## Your Bootcamp Folder Structure

By now, your `~/bootcamp` folder should look like this:

```
~/bootcamp/
  README.md
  test_setup.ipynb
```

We will add more files each week. Keep everything organized here.

---

## Common Setup Problems

| Problem | Solution |
|---------|----------|
| `conda: command not found` | Close and reopen your terminal. If still broken, re-run the Miniconda installer and type `yes` when asked about initialization. |
| `pip install` says "permission denied" | Make sure you activated your environment first: `conda activate immunology` |
| Jupyter won't open in browser | Copy one of the URLs from the terminal output and paste it into your browser manually |
| R package installation fails | Try restarting RStudio, then run the install command again. Some packages need a C compiler -- on Mac, install Xcode Command Line Tools. |
| `git` not found on Mac | Run `xcode-select --install` in your terminal |
| Everything is very slow | Make sure you have enough free disk space (at least 10 GB). Close other applications. |

For more help, see the [Troubleshooting Guide](../TROUBLESHOOTING.md).
