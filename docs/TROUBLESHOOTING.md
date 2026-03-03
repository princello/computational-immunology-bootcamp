# Troubleshooting Guide

This guide covers the most common errors you will encounter during the bootcamp. Before contacting anyone, work through this checklist:

1. **Read the error message carefully.** The last few lines usually contain the actual problem.
2. **Check this guide** for your specific error.
3. **Google the exact error message** in quotes (e.g., `"ModuleNotFoundError: No module named 'scanpy'"`)
4. **Ask an AI assistant** (ChatGPT or Claude) -- paste your code and the full error message.
5. **Search [Biostars](https://www.biostars.org/) or [Stack Overflow](https://stackoverflow.com/).**

---

## Table of Contents

- [Installation Problems](#installation-problems)
- [Conda and Environment Problems](#conda-and-environment-problems)
- [Jupyter Notebook Problems](#jupyter-notebook-problems)
- [Python Errors](#python-errors)
- [Scanpy and Single-Cell Errors](#scanpy-and-single-cell-errors)
- [R and RStudio Problems](#r-and-rstudio-problems)
- [Git and GitHub Problems](#git-and-github-problems)
- [Data Download Problems](#data-download-problems)
- [Performance and Memory Problems](#performance-and-memory-problems)
- [General Tips](#general-tips)

---

## Installation Problems

### "conda: command not found"

**Cause:** Miniconda was installed but your terminal does not know where to find it.

**Fix:**
1. Close your terminal completely and open a new one
2. If still broken, run:
   ```bash
   # Mac/Linux:
   source ~/miniconda3/bin/activate
   conda init bash   # or: conda init zsh (if you use zsh)
   ```
3. Close and reopen the terminal again
4. If still broken, reinstall Miniconda and make sure to type `yes` when asked "Do you wish to update your shell profile to automatically initialize conda?"

### "pip: command not found"

**Cause:** pip is not installed or your environment is not activated.

**Fix:**
```bash
conda activate immunology
# pip should now work
```

### Package installation hangs or is very slow

**Cause:** Network issues or conda solving dependencies.

**Fix:**
- Be patient -- some installations take 10-20 minutes
- If it hangs for more than 30 minutes, press `Ctrl + C` to cancel and try:
  ```bash
  pip install --no-cache-dir scanpy
  ```
- Try a different network (e.g., switch from Wi-Fi to phone hotspot)

### "ERROR: Could not build wheels for ..."

**Cause:** A package needs to be compiled but your system is missing the required tools.

**Fix (Mac):**
```bash
xcode-select --install
```

**Fix (Windows):**
- Install [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
- Select "C++ build tools" during installation

**Fix (Linux):**
```bash
sudo apt install build-essential python3-dev
```

---

## Conda and Environment Problems

### "(base)" instead of "(immunology)" in my prompt

**Cause:** You forgot to activate the immunology environment.

**Fix:**
```bash
conda activate immunology
```
Do this every time you open a new terminal.

### "EnvironmentNotWritableError"

**Cause:** Permission issue with the conda directory.

**Fix (Mac/Linux):**
```bash
sudo chown -R $USER ~/miniconda3
```

### "PackagesNotFoundError"

**Cause:** The package is not available in the current conda channel.

**Fix:** Use pip instead:
```bash
pip install package_name
```

### I installed a package but Python cannot find it

**Cause:** You installed it in the wrong environment.

**Fix:**
1. Check which environment you are in: look at the prompt (should say `(immunology)`)
2. Activate the right environment: `conda activate immunology`
3. Reinstall: `pip install package_name`
4. Verify: `python -c "import package_name"`

---

## Jupyter Notebook Problems

### Jupyter does not open in my browser

**Cause:** The browser did not launch automatically.

**Fix:** Copy one of the URLs from the terminal output (looks like `http://localhost:8888/?token=...`) and paste it into your browser.

### "Kernel not found" or wrong Python version

**Cause:** Jupyter is using a different Python than your conda environment.

**Fix:**
```bash
conda activate immunology
pip install ipykernel
python -m ipykernel install --user --name immunology --display-name "Python (immunology)"
```
Then in Jupyter, go to Kernel -> Change kernel -> Python (immunology).

### Notebook is not saving

**Cause:** Disk space or permissions issue.

**Fix:**
- Check disk space: `df -h` (Mac/Linux)
- Make sure you are saving to a folder you own (e.g., ~/bootcamp, not a system folder)
- Try File -> Save As and save to a new location

### Cell runs forever (spinning asterisk [*])

**Cause:** The code is taking a long time or is stuck in an infinite loop.

**Fix:**
- Wait a few minutes (some operations are genuinely slow)
- If it seems stuck, go to Kernel -> Interrupt (or press `I I` -- the letter I twice)
- If that does not work, go to Kernel -> Restart (you will lose all variables and need to re-run previous cells)

---

## Python Errors

### "ModuleNotFoundError: No module named 'xyz'"

**Cause:** The package is not installed in your current environment.

**Fix:**
```bash
conda activate immunology
pip install xyz
```
Then restart your Jupyter kernel (Kernel -> Restart) and re-run the import cell.

### "SyntaxError: invalid syntax"

**Cause:** A typo in your code.

**Common causes:**
- Missing colon at the end of `if`, `for`, `def`, `class` statements
- Mismatched parentheses or brackets
- Using a Python 2 syntax in Python 3
- Copying code that has invisible special characters (retype it manually)

### "IndentationError: unexpected indent"

**Cause:** Inconsistent use of spaces and tabs.

**Fix:** Make sure all indentation uses the same number of spaces (4 is standard). In Jupyter, select the problematic lines and press Tab to indent or Shift+Tab to unindent.

### "NameError: name 'xyz' is not defined"

**Cause:** You are trying to use a variable or function that has not been created yet.

**Fix:**
- Did you run all the cells above? (In Jupyter, run cells in order from top to bottom)
- Did you restart the kernel? (If so, re-run all cells from the beginning)
- Check for typos in the variable name

### "TypeError: 'xyz' object is not subscriptable" or similar

**Cause:** You are trying to use an object in a way it does not support.

**Fix:** Check the type of your variable with `type(variable_name)` and make sure you are using the correct syntax for that type.

---

## Scanpy and Single-Cell Errors

### "adata.X is None" or empty results

**Cause:** The data was not loaded correctly, or a previous filtering step removed all cells/genes.

**Fix:**
- Reload the data from scratch
- Check `adata.shape` after each step to make sure you still have cells and genes

### "ValueError: could not convert string to float"

**Cause:** The data matrix contains non-numeric values, often because gene names ended up in the matrix.

**Fix:**
```python
# When loading, make sure gene names go to var_names
adata = sc.read_10x_mtx("path/", var_names="gene_symbols", cache=True)
```

### UMAP/clustering looks wrong (one giant blob or everything scattered)

**Possible causes and fixes:**
1. **Forgot to normalize:** Make sure you ran `sc.pp.normalize_total` and `sc.pp.log1p`
2. **Forgot to select HVGs:** Make sure you selected highly variable genes before PCA
3. **Resolution too low/high:** Try different resolution values in `sc.tl.leiden(adata, resolution=0.5)`
4. **Too few/many PCs:** Try `n_pcs=30` or `n_pcs=40` in `sc.pp.neighbors`
5. **Batch effects:** If mixing datasets, you may need batch correction

### "ImportError: leidenalg is not installed"

**Fix:**
```bash
conda activate immunology
pip install leidenalg
```

### Scanpy plots do not show in Jupyter

**Fix:** Add this at the top of your notebook:
```python
%matplotlib inline
```
Or update matplotlib:
```bash
pip install --upgrade matplotlib
```

### "FutureWarning" or "DeprecationWarning"

**Cause:** You are using a function that will change in a future version.

**Fix:** These are warnings, not errors. Your code still works. You can ignore them or suppress them:
```python
import warnings
warnings.filterwarnings('ignore')
```

---

## R and RStudio Problems

### R package installation fails

**Common fix:**
```r
# Try installing from binary instead of source
install.packages("Seurat", type = "binary")
```

If asked "Do you want to install from sources the package which needs compilation?", type `no`.

### "Error in library(xyz): there is no package called 'xyz'"

**Fix:** Install the package:
```r
install.packages("xyz")
# or for Bioconductor packages:
BiocManager::install("xyz")
```

### RStudio crashes or freezes

**Fix:**
- Close RStudio
- Delete the session file: `rm ~/.local/share/rstudio/sessions/active/*/session-persistent-state` (Mac/Linux)
- Reopen RStudio

---

## Git and GitHub Problems

### "git: command not found"

**Fix (Mac):**
```bash
xcode-select --install
```

**Fix (Windows):** Download and install from https://git-scm.com/download/win

### "fatal: not a git repository"

**Cause:** You are not in a folder that has been initialized with git.

**Fix:**
```bash
cd ~/bootcamp  # or wherever your project is
git init       # only if this is a new project
```

### "error: failed to push some refs"

**Cause:** The remote repository has changes you don't have locally.

**Fix:**
```bash
git pull origin main
# Then try pushing again
git push origin main
```

---

## Data Download Problems

### Download is very slow or keeps failing

**Fix:**
- Try a different network connection
- Use `wget` with resume capability:
  ```bash
  wget -c "URL"  # -c allows resuming interrupted downloads
  ```
- For large files, try downloading during off-peak hours

### "File not found" after download

**Fix:**
- Check where the file was downloaded: `ls ~/Downloads/`
- Move it to the right place: `mv ~/Downloads/filename data/raw/`
- Make sure you are using the correct file path in your code

### Cannot access a dataset URL

**Fix:**
- Some datasets require registration (e.g., 10x Genomics). Create a free account.
- Try a different dataset or use the built-in Scanpy datasets (e.g., `sc.datasets.pbmc3k()`)
- Check if the URL has changed by searching for the dataset name

---

## Performance and Memory Problems

### "MemoryError" or "Killed"

**Cause:** Your computer ran out of RAM.

**Fix:**
- Close other applications (especially web browsers with many tabs)
- Work with a smaller dataset: `adata = adata[:5000, :].copy()` (use only 5000 cells)
- Use sparse matrix operations (Scanpy does this by default)
- If the problem persists, use Google Colab (free, gives you more RAM):
  1. Go to https://colab.research.google.com
  2. Upload your notebook
  3. Install packages with `!pip install scanpy`

### Everything is very slow

**Fix:**
- Check available memory: `python -c "import psutil; print(f'{psutil.virtual_memory().available / 1e9:.1f} GB free')"` (install psutil first: `pip install psutil`)
- Close other applications
- Use fewer cells or genes for testing, then scale up
- Some operations (like scVI training) benefit from a GPU -- consider Google Colab with GPU runtime

---

## General Tips

### How to Read an Error Message

Error messages in Python are called "tracebacks." Read them from **bottom to top**:

```
Traceback (most recent call last):
  File "notebook.py", line 15, in <module>     <-- where the error happened
    adata = sc.read_10x_mtx("wrong_path/")     <-- the line of code
FileNotFoundError: [Errno 2] No such file       <-- THE ACTUAL ERROR
```

The last line is the most important. It tells you what went wrong.

### The Nuclear Option: Start Fresh

If everything is broken and you cannot figure out why:

```bash
# Remove your conda environment
conda deactivate
conda env remove -n immunology

# Recreate it
conda create -n immunology python=3.11 -y
conda activate immunology
pip install scanpy scirpy scvi-tools celltypist leidenalg
pip install pandas numpy matplotlib seaborn scikit-learn jupyter
pip install anndata mudata session-info
```

This takes 10-15 minutes but gives you a clean slate.

### Using Google Colab as a Fallback

If your computer is too slow or you cannot get packages to install, use Google Colab:

1. Go to https://colab.research.google.com
2. Click "New notebook"
3. In the first cell, install your packages:
   ```python
   !pip install scanpy scirpy scvi-tools celltypist
   ```
4. Work through the tutorials just like in Jupyter
5. Your notebooks are saved to Google Drive

Colab is free and gives you a more powerful computer in the cloud. The downside is that your session resets after a few hours of inactivity, so save your work frequently.

---

*Still stuck? That is normal and happens to every programmer. Take a break, come back with fresh eyes, and try again. The solution is almost always simpler than you think.*
