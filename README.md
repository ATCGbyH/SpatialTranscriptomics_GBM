# Spatial Transcriptomics - Seurat Analysis: Human Glioblastoma (GBM) FFPE Tissue

This project uses 10x Visium CytAssist and Seurat to map gene expression across a human glioblastoma (GBM) tissue section. It compares collagen- and hypoxia-associated regions and asks how similar three *MBP*-rich spot clusters are across their broader expression profiles.

[Open the R notebook](https://github.com/ATCGbyH/SpatialTranscriptomics_GBM/blob/main/R_notebook.pdf) · [View the R Markdown source](GBM_Spatial_Seurat.Rmd)

The notebook contains the complete analysis and figures. For a local preview, download and open `docs/notebook.html` in a browser.

## Dataset

The analysis uses one formalin-fixed, paraffin-embedded (FFPE) GBM section profiled on an 11 mm Visium CytAssist capture area. The filtered Space Ranger 2.0.1 output yielded **18,085 features across 10,878 tissue-covered spots** in Seurat. The dataset is available under CC BY 4.0. [10x Genomics dataset and downloads](https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard)

## Workflow

1. Inspect counts per spot and normalise the spatial assay with `SCTransform()`.
2. Use the first 30 principal components for graph-based clustering and UMAP, then map the **21 spot clusters** onto the H&E image.
3. Rank spatially variable genes using Moran’s I and a mark variogram. Identify cluster markers with `FindAllMarkers()`.
4. Map scores for collagen, hypoxia, MES1- and MES2-associated markers, and a three-gene GBM marker panel (*BCAN*, *PTPRZ1*, *DDR1*).
5. Compare cluster-average SCT residuals using Spearman correlation, both across all spots and within the *MBP*-high subset.

## Findings

**Collagen and hypoxia.** Clusters 5 and 9 showed elevated expression of collagen genes and a cluster-profile correlation of **ρ = 0.571**. Areas with high hypoxia marker scores overlapped cluster 7.

**Myelin-associated regions.** Using the section-wide 80th percentile of *MBP* as the cutoff, **92.0%** of cluster 1, **86.8%** of cluster 14 and **95.4%** of cluster 12 were *MBP*-high. These spots had similar median *PLP1/MOG/MOBP* module scores across the three clusters (2.13–2.14).

| Cluster pair | Spearman ρ, all spots | Spearman ρ, *MBP*-high spots |
| --- | ---: | ---: |
| 1 and 14 | 0.944 | 0.939 |
| 12 and 1 | 0.457 | 0.493 |
| 12 and 14 | 0.434 | 0.468 |

Excluding the four myelin markers yielded similar correlations. Cluster 12 had lower mean SCT residuals for several myeloid and complement-associated genes, including *CD163*, *C1QA* and *VSIG4*. The results describe differences in spot-level expression within one section. Assigning malignant cell identity would require additional evidence, such as histopathology and tumour-associated copy-number profiles. [Myeloid programmes in human glioma](https://doi.org/10.1038/s41586-025-08633-8)

## Reproduce the analysis

1. Use R with Seurat v5 and install the packages loaded by the notebook, including `SeuratData`, `hdf5r`, `Rfast2`, `dplyr`, `ggplot2`, `patchwork`, `pheatmap`, `htmltools` and `vembedr`. Use `rmarkdown` and `knitr` to render the notebook. [Seurat installation guide](https://satijalab.org/seurat/articles/install_v5.html)
2. Download the filtered feature-barcode matrix and spatial imaging archive from 10x Genomics. Place the matrix at `data/filtered_feature_bc_matrix.h5` and extract the images, coordinates and scale factors into `data/spatial/`.
3. Place `data/` beside the R Markdown file or update `data_dir` to point to the input folder. Run the chunks in order.

*Analysis by Hugo Roca*
