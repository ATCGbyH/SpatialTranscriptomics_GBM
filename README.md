# Spatial Transcriptomics - Seurat Analysis: Human Glioblastoma (GBM) FFPE Tissue

This project uses 10x Visium CytAssist and Seurat to map gene expression across a human glioblastoma (GBM) tissue section. It compares collagen- and hypoxia-associated regions and asks how similar three *MBP*-rich spot clusters are across their broader expression profiles.

[Open the R notebook in PDF for visualisation](https://github.com/ATCGbyH/SpatialTranscriptomics_Seurat_GBM/blob/main/GBM_Spatial_Seurat.pdf) · [View the R Markdown source](https://github.com/ATCGbyH/SpatialTranscriptomics_Seurat_GBM/blob/main/GBM_Spatial_Seurat.Rmd)

The notebook contains the complete analysis and figures. For a local preview, download and open [HTML R notebook](https://github.com/ATCGbyH/SpatialTranscriptomics_Seurat_GBM/blob/main/notebook.html) in a browser.

<img src="https://i.postimg.cc/PJHHk6PN/GBM-Findings-Summary.png" height="80%" width="80%" alt="Workflow Chart"/>

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

Clusters 5 and 9 shared elevated expression of *COL1A1*, *COL1A2* and *COL3A1*. Their broader gene rankings showed moderate agreement based on cluster-average SCT residuals (Spearman ρ = 0.571). The collagen profiles still differed within that shared programme. For example, *COL6A3* was detected in **88.7% of cluster 5 spots** and **45.7% of cluster 9 spots** in the marker analysis.

<img src="https://i.postimg.cc/zvtZ93Dw/GBM-Top5-Marker-Heatmap.png" width="80%" alt="Top marker genes across glioblastoma spot clusters"/>

High scores for the hypoxia-associated panel overlapped cluster 7. The individual markers supported this pattern. *CA9* was detected in **80.3% of cluster 7 spots** compared with **7.9% of spots outside cluster 7**. *HILPDA* was detected in **98.0%** and **28.3%**, respectively. The panel also included *CA12*, *ANGPTL4* and *VEGFA*.

<img src="https://i.postimg.cc/jdY26h8J/spatial-gene-programs.png" width="80%" alt="Spatial maps of collagen, hypoxia and other GBM-associated gene programmes"/>

The collagen-associated expression is consistent with extracellular matrix remodelling. Spatial studies have identified collagen-rich regions around vessels in recurrent GBM and linked hypoxia-associated expression to tissue organisation around necrotic cores (Onubogu et al., 2024; Greenwald et al., 2024). Reviewing the matched H&E alongside the score maps would help assess whether the collagen-enriched spots lie near vessels or necrotic borders, and whether the cluster 7 signal overlaps tissue surrounding necrosis.

### Shared myelin signal (Cluster 1 and 14) and divergence in cluster 12

Within the *MBP*-high subset, clusters 1, 14 and 12 had similar median *PLP1/MOG/MOBP* module scores (2.13–2.14). Across the broader set of genes, clusters 1 and 14 showed strong agreement in their gene rankings, while cluster 12 was moderately correlated with both. This pattern persisted after excluding *MBP*, *PLP1*, *MOG* and *MOBP* from the correlation.

<img src="https://i.postimg.cc/hGXKRJyJ/mbp-cluster-correlation.png" height="80%" width="80%" alt="Workflow Chart"/>

| Cluster pair | Spearman ρ, all spots | Spearman ρ, *MBP*-high spots |
| --- | ---: | ---: |
| 1 and 14 | 0.944 | 0.939 |
| 12 and 1 | 0.457 | 0.493 |
| 12 and 14 | 0.434 | 0.468 |

Relative to clusters 1 and 14, cluster 12 had lower mean SCT residuals for monocyte- and macrophage-associated genes (*CD14*, *CD163*) and complement-associated genes (*C1QA*, *VSIG4*). For example, the mean *CD163* residual was −1.67 in cluster 12, compared with +2.64 in cluster 1 and +3.12 in cluster 14. The agreement across these genes is consistent with a stronger myeloid-associated RNA signal in clusters 1 and 14 (Miller et al., 2025).

*KIF5A* showed the opposite pattern. Its mean SCT residual was +1.80 in cluster 12 and negative in both clusters 1 and 14. Given the role of *KIF5A* in axonal transport, this pattern is compatible with a greater neuronal RNA contribution in cluster 12. Joyce and Lindsay (2024) also reported an inverse association between *KIF5A* expression and glioma grade in a separate cohort analysis.

<img src="https://i.postimg.cc/GtW2ZMm5/GBM-Cluster-12-Comparison-Numbered.png" height="80%" width="80%" alt="Workflow Chart"/>

This could be a good explanation on why these clusters may be from the specific cell type, although, with a high chance of a GBM tumour being slowly forming within cluster 1 and 14, while cluster 12 cells demonstrate to be closely within healthy states.

Another major factor that may indicate this is the highly expression of CD14 and CD163 cells within cluster 1 and 14, indicating direct activity of immune system in this region and possibly part of a forming GBM. 

The three clusters therefore had similar myelin module scores within the *MBP*-high subset, alongside differences in myeloid- and neuronal-associated signals across their broader spot profiles. These differences could reflect variation in the cell types captured by each Visium spot, changes in their transcriptional states, or both.

These differences could also merely reflect variation in cell composition or transcriptional state across Visium spots. Histology and tumour-associated copy-number profiles would help assess whether the myelin-rich regions also differ in local tumour infiltration.

## Reproduce the analysis

1. Use R with Seurat v5 and install the packages loaded by the notebook, including `SeuratData`, `hdf5r`, `Rfast2`, `dplyr`, `ggplot2`, `patchwork`, `pheatmap`, `htmltools` and `vembedr`. Use `rmarkdown` and `knitr` to render the notebook. [Seurat installation guide](https://satijalab.org/seurat/articles/install_v5.html)
2. Download the filtered feature-barcode matrix and spatial imaging archive from 10x Genomics. Place the matrix at `data/filtered_feature_bc_matrix.h5` and extract the images, coordinates and scale factors into `data/spatial/`.
3. Place `data/` beside the R Markdown file or update `data_dir` to point to the input folder. Run the chunks in order.

## References

· Joyce, L.J. and Lindsay, A.J. (2024) ‘A systematic computational analysis of the endosomal recycling pathway in glioblastoma’, Biochemistry and Biophysics Reports, 38, p. 101700. https://pubmed.ncbi.nlm.nih.gov/38638676/

· Miller, T.E. et al. (2025) ‘Programs, origins and immunomodulatory functions of myeloid cells in glioma’, Nature, 640(8060), pp. 1072–1082. https://www.nature.com/articles/s41586-025-08633-8

· Greenwald, A.C. et al. (2024). ‘Integrative spatial analysis reveals a multi-layered organization of glioblastoma’. *Cell*, 187(10), pp. 2485–2501.e26. https://doi.org/10.1016/j.cell.2024.03.029
  
· Onubogu, U. et al. (2024). ‘Spatial analysis of recurrent glioblastoma reveals perivascular niche organization’. *JCI Insight*, 9(12), https://insight.jci.org/articles/view/179853

*Analysis by Hugo Roca*
