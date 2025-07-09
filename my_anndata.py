# Code to load the converted h5ad file 
# Ethan Chang, Ph.D. Student @ WashU
# ethanc@wustl.edu
# 07/08/2025

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib as mp
from tkinter import filedialog, simpledialog
from scipy.sparse import csr_matrix


# Loads the file
file_path = filedialog.askopenfilename(
    title="Select .h5ad file",
    filetypes=[("h5ad files", "*.h5ad")]
)

if file_path:
    adata = ad.read_h5ad(file_path)
    print(f"Loaded {adata.shape[0]} cells and {adata.shape[1]} genes.")
else:
    print("No file selected.")

print(adata)


# Convert from pandas array to Numpy array ('X_umap" name and Numpy array required by Scanpy)
adata.obsm['X_umap'] = adata.obsm['UMAP'].values


# Prompts users for annotations to plot
colors = []
while True:
    options = "\n".join(adata.obs.columns)

    selected_color = simpledialog.askstring(
        "Choose annotation", 
        f"Available annotations:\n{options}\n\nWhich one would you like to color by?")

    if not selected_color:
        break

    if selected_color not in adata.obs.columns:
        print(f"'{selected_color}' is not a valid annotation.")
    else:
        colors.append(selected_color)

    another = simpledialog.askstring(
        "Add another?",
        "Would you like to add another annotation? (yes/no)")
    
    if another is None or another.lower() != 'yes':
        break

    

# Plots annotations
sc.pl.umap(adata, color=colors)


