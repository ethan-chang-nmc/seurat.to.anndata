# seurat.to.anndata
<b>Ethan Chang<b/><br>
<b>Ph.D. Neuroscience @ WashU<b/><br>
<b>Ackerman Lab (rotation), 2025</b><br>

This code was made to convert Seurat (1) files of scRNA data to Anndata (2) using the Zellkonverter library (3). <br>
- convert_seurat_to_anndata.R: Run this script first. This prompts the user to choose a Seurat file and converts it to a h5ad file. This file is downloaded to a user-specified file location. <br>
- environment.yml: Use to get the environment this code was created on. <br>
- my_anndata: Contains functions to load and plot anndata objects <br>

(1) Hao, Y., Stuart, T., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman, A., Srivastava, A., Molla, G., Madad, S., & Satija, R. (2024). Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology, 42(2), 293-304. https://doi.org/10.1038/s41587-023-01767-y <br>
(2) anndata: Annotated data. Isaac Virshup, Sergei Rybakov, Fabian J. Theis, Philipp Angerer, F. Alexander Wolf. JOSS 2024 Sep 16. doi: 10.21105/joss.04371. <br>
(3) Zappia L, Lun A (2025). zellkonverter: Conversion Between scRNA-seq Objects. doi:10.18129/B9.bioc.zellkonverter, R package version 1.18.0, https://bioconductor.org/packages/zellkonverter.
