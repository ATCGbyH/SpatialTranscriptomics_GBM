# Spatial Transcriptomics of Human Glioblastoma

A Seurat analysis of 10x Visium CytAssist data from a human glioblastoma (GBM) tissue section. I mapped transcriptional clusters onto the H&E image and compared where collagen, hypoxia and other gene programs were expressed.

[Project page](https://atcgbyh.github.io/SpatialTranscriptomics_GBM/) · [Analysis notebook](https://atcgbyh.github.io/SpatialTranscriptomics_GBM/notebook.html) · [R Markdown source](docs/GBM_Spatial_Seurat.Rmd)

![UMAP embedding and H&E tissue image with the same 21 colour-coded spot clusters](docs/assets/spot-clusters.png)

## Key observations

- **21 spot clusters** were identified across 10,878 tissue-covered spots. Cluster 0 was the largest, with 1,262 spots.
- **Clusters 5 and 9** showed elevated collagen-gene expression and had a Spearman correlation of **0.571** across cluster-mean SCTransform residuals. Areas with high hypoxia module scores overlapped **cluster 7**.
- **Clusters 1 and 14** occupied neighbouring areas of the section. They had the most similar cluster-mean profiles among seven clusters selected for high `MBP` expression (**Spearman ρ = 0.944**).

## Dataset

The notebook imports **18,085 features** from one formalin-fixed, paraffin-embedded human glioblastoma section profiled on an 11 mm Visium CytAssist capture area. The source data were processed with Space Ranger 2.0.1 and released under CC BY 4.0. [10x Genomics dataset and downloads](https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard)

## Analysis workflow

1. **Import and normalise.** Inspect total counts per spot and normalise with `SCTransform()`.
2. **Cluster and map spots.** Use the first 30 principal components for graph-based clustering and UMAP, then overlay cluster labels on the tissue image.
3. **Rank spatially variable genes.** Apply Moran’s I to 1,000 highly variable genes and a mark variogram to 10 highly variable genes.
4. **Compare clusters.** Rank positive markers by log₂ fold change and compare cluster-mean SCT residuals using Spearman correlation.
5. **Score marker panels.** Map collagen, hypoxia, MES1-associated, MES2-associated and GBM marker panels using `AddModuleScore()`.

The module scores describe relative expression at spot level in this tissue section.

## View or rerun the notebook

Open `docs/index.html` locally to view the project page, or `docs/notebook.html` for the saved analysis. Select **Code → Download Rmd** in the notebook to get the editable source.

1. Use R with Seurat v5. The source loads `SeuratData`, `hdf5r`, `Rfast2`, `dplyr`, `ggplot2`, `patchwork`, `pheatmap`, `htmltools` and `vembedr`. Use `rmarkdown` and `knitr` to render the notebook. [Seurat installation guide](https://satijalab.org/seurat/articles/install_v5.html)
2. Download the filtered feature-barcode HDF5 matrix and spatial imaging archive from 10x Genomics. Save the matrix as `data/filtered_feature_bc_matrix.h5` and extract the images, coordinates and scale factors into `data/spatial/`.
3. Place `data/` beside the `.Rmd` file or set `data_dir` to its location. Run the notebook chunks in order.

*Analysis by Hugo Roca.*
