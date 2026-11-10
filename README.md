# Spatial Transcriptomics - Seurat Analysis: Human Glioblastoma (GBM) FFPE Tissue

This project uses 10x Visium CytAssist and Seurat to map gene expression across a human glioblastoma (GBM) tissue section. It compares collagen- and hypoxia-associated regions and asks how similar three *MBP*-rich spot clusters are across their broader expression profiles.

[Open the R notebook](https://github.com/ATCGbyH/SpatialTranscriptomics_Seurat_GBM/blob/main/R_notebook.pdf) · [View the R Markdown source](https://github.com/ATCGbyH/SpatialTranscriptomics_Seurat_GBM/blob/main/GBM_Spatial_Seurat.Rmd)

The notebook contains the complete analysis and figures. For a local preview, download and open [HTML R notebook](https://github.com/ATCGbyH/SpatialTranscriptomics_Seurat_GBM/blob/main/notebook.html) in a browser.

<img src="https://i.postimg.cc/D08ZCC9d/GBM-Findings-Summary.png" height="80%" width="80%" alt="Workflow Chart"/>

## Dataset

The analysis uses one formalin-fixed, paraffin-embedded (FFPE) GBM section profiled on an 11 mm Visium CytAssist capture area. The filtered Space Ranger 2.0.1 output yielded **18,085 features across 10,878 tissue-covered spots** in Seurat. 

*H&E and Visium data from [10x Genomics](https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard), licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).*

<img src="https://i.postimg.cc/Pfgp3HqS/GBM-QC-n-Count-Spatial.png" height="80%" width="80%" alt="Workflow Chart"/>

## Workflow

1. Inspect counts per spot and normalise the spatial assay with `SCTransform()`.
2. Use the first 30 principal components for graph-based clustering and UMAP, then map the **21 spot clusters** onto the H&E image.
3. Rank spatially variable genes using Moran’s I and a mark variogram. Identify cluster markers with `FindAllMarkers()`.
4. Map scores for collagen, hypoxia, MES1- and MES2-associated markers, and a three-gene GBM marker panel (*BCAN*, *PTPRZ1*, *DDR1*).
5. Compare cluster-average SCT residuals using Spearman correlation, both across all spots and within the *MBP*-high subset.

## Findings and interpretation

<img src="https://i.postimg.cc/50ndRLLX/GBM-Tissue-to-Clusters.png" height="80%" width="80%" alt="Workflow Chart"/>

Seurat identified **21 clusters across 10,878 tissue-covered Visium spots**. The spatial map shows where the clusters occur on the H&E section, while the UMAP arranges the same spots by gene expression. 

Clusters 1 and 14 occupy neighbouring areas in the upper-left of the tissue and sit close together on the UMAP. Cluster 12 lies mainly on the right side of the section and separates from both on the UMAP.

### Collagen enrichment and hypoxia

Clusters 5 and 9 showed elevated collagen-gene expression and moderate agreement in gene rankings based on cluster-average SCT residuals (Spearman ρ = 0.571). High scores for a hypoxia-associated panel that included *CA9*, *VEGFA* and *HILPDA* overlapped cluster 7.

The collagen signal is consistent with extracellular matrix remodelling. Comparing the collagen and hypoxia maps with the matched H&E image would help assess whether either pattern is associated with necrosis in this GBM section.

<img src="https://i.postimg.cc/jdY26h8J/spatial-gene-programs.png" height="80%" width="80%" alt="Workflow Chart"/>

### Shared myelin signal (Cluster 1 and 14) and divergence in cluster 12

Within the *MBP*-high subset, clusters 1, 14 and 12 had similar median *PLP1/MOG/MOBP* module scores (2.13–2.14). Across the broader set of genes, clusters 1 and 14 showed strong agreement in their gene rankings, while cluster 12 was moderately correlated with both. This pattern persisted after excluding *MBP*, *PLP1*, *MOG* and *MOBP* from the correlation.

<img src="https://i.postimg.cc/hGXKRJyJ/mbp-cluster-correlation.png" height="80%" width="80%" alt="Workflow Chart"/>

| Cluster pair | Spearman ρ, all spots | Spearman ρ, *MBP*-high spots |
| --- | ---: | ---: |
| 1 and 14 | 0.944 | 0.939 |
| 12 and 1 | 0.457 | 0.493 |
| 12 and 14 | 0.434 | 0.468 |

Relative to clusters 1 and 14, cluster 12 had lower mean SCT residuals for monocyte- and macrophage-associated genes (*CD14*, *CD163*) and complement-associated genes (*C1QA*, *VSIG4*). For example, the mean *CD163* residual was −1.67 in cluster 12, compared with +2.64 in cluster 1 and +3.12 in cluster 14. The agreement across these genes is consistent with a stronger myeloid-associated RNA signal in clusters 1 and 14 (Miller et al., 2025).

*KIF5A* showed the opposite pattern. Its mean SCT residual was +1.80 in cluster 12 and negative in both clusters 1 and 14. Given the role of *KIF5A* in axonal transport, this pattern is compatible with a greater neuronal RNA contribution in cluster 12 (Guerra San Juan et al., 2024). Joyce and Lindsay (2024) also reported an inverse association between *KIF5A* expression and glioma grade in a separate cohort analysis.

This could be a good explanation on why these clusters present such high similarity but also, somewhat dissimilar, with a high chance of a GBM tumour slowly forming within cluster 1 and 14, while cluster 12 cells demonstrate to be closely within healthy states.

Another major factor that may indicate this is the highly expression of CD14 and CD163 cells within cluster 1 and 14, indicating direct activity of immune system in this region and possibly part of a forming GBM. 

The three clusters therefore had similar myelin module scores within the *MBP*-high subset, alongside differences in myeloid- and neuronal-associated signals across their broader spot profiles. These differences could reflect variation in the cell types captured by each Visium spot, changes in their transcriptional states, or both.

<img src="https://i.postimg.cc/br5QLVvb/GBM-Cluster-12-Comparison.png" height="80%" width="80%" alt="Workflow Chart"/>

These differences could also merely reflect variation in cell composition or transcriptional state across Visium spots. Histology and tumour-associated copy-number profiles would help assess whether the myelin-rich regions also differ in local tumour infiltration.

## Reproduce the analysis

1. Use R with Seurat v5 and install the packages loaded by the notebook, including `SeuratData`, `hdf5r`, `Rfast2`, `dplyr`, `ggplot2`, `patchwork`, `pheatmap`, `htmltools` and `vembedr`. Use `rmarkdown` and `knitr` to render the notebook. [Seurat installation guide](https://satijalab.org/seurat/articles/install_v5.html)
2. Download the filtered feature-barcode matrix and spatial imaging archive from 10x Genomics. Place the matrix at `data/filtered_feature_bc_matrix.h5` and extract the images, coordinates and scale factors into `data/spatial/`.
3. Place `data/` beside the R Markdown file or update `data_dir` to point to the input folder. Run the chunks in order.

## References

· Joyce, L.J. and Lindsay, A.J. (2024) ‘A systematic computational analysis of the endosomal recycling pathway in glioblastoma’, Biochemistry and Biophysics Reports, 38, p. 101700. doi:10.1016/j.bbrep.2024.101700. 

· Miller, T.E. et al. (2025) ‘Programs, origins and immunomodulatory functions of myeloid cells in glioma’, Nature, 640(8060), pp. 1072–1082. doi:10.1038/s41586-025-08633-8. 

· Guerra San Juan, I. et al. (2024) Kif5a regulates axonal repair and time-dependent axonal transport of SFPQ granules and mitochondria in human motor neurons [Preprint]. doi:10.1101/2024.09.06.611684. 

*Analysis by Hugo Roca*
