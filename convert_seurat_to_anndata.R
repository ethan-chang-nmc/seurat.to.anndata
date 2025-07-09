# Code to convert seurat object to anndata
# Ethan Chang, Ph.D. Student @ WashU
# ethanc@wustl.edu
# 07/08/2025


# Run This block only if you haven't installed zellkonverter before
if (!requireNamespace("BiocManager", quietly=TRUE)) 
{
  install.packages("BiocManager")
}

BiocManager::install("zellkonverter")
###########

library(Seurat)
library(zellkonverter)
library(tcltk)


# Prompt user to select a file and reads file
file_path <- file.choose()
seurat_obj <- readRDS(file_path)


# Check if it's a Seurat object
if (!inherits(seurat_obj, "Seurat")) 
{
  stop("Error: The selected file is not a Seurat object.")
}

print("Seurat object successfully loaded.")


# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj)


# Save as .h5ad
save_path <- tclvalue(tkgetSaveFile(
  filetypes = "{{H5AD Files} {.h5ad}} {{All files} *}"
))

if (save_path != "") 
{
  writeH5AD(sce, save_path)
  message("File successfully saved to: ", save_path)
} 
else 
{
  message("Save operation canceled.")
}


