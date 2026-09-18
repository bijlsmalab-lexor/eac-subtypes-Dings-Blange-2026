# Tumor-intrinsic signatures in esophageal adenocarcinoma

Analysis code accompanying:

> Dings M.P.G.\*, Blangé D.\*, van Schie D.M., *et al.*, van Laarhoven H.W.M.\*, Bijlsma M.F.\*
> **Tumor-intrinsic signatures in esophageal adenocarcinoma associate with cellular phenotypes and responses to therapy.**
> *Nature Communications* (in press). <!-- add DOI on acceptance -->

\* Equal contribution.

This repository contains the R scripts and R Markdown notebooks used to generate the figures in the manuscript. It defines two tumor-intrinsic molecular subtypes of esophageal adenocarcinoma (EAC) — hereafter **IL** and **ML** — and relates them to single-cell phenotypes, regulatory networks, and responses to EGFR/HER2, cMET and IGF1R inhibition.

## Repository structure

Each folder contains the code for the corresponding main or supplementary figure.

| Folder | Contents |
|--------|----------|
| `Figure 1/`  | Subtype discovery / consensus classification of the EAC cohort (`eac.bijlsma.186.s4b.Rmd`, rendered `.nb.html`) |
| `Figure 2/`  | Single-cell subtype signal (`2E.FLOTplot.R`) |
| `Figure 3/`  | Deconvolution of bulk (CPCT) tumors (`deconvolve_CPCT.Rmd`, rendered `.html`) |
| `Figure 5/`  | Network analysis and gastric-cancer classifier signatures (`eac.bijlsma.186.s4b.network-analysis.Rmd`, `gc.classifier.signatures.Rmd`) |
| `Figure 6/`  | Transcription-factor regulon analysis (`TFs.Rmd`) |
| `Figure 7/`  | Subtype-associated drug sensitivity in GDSC/ESCA (`dings.km.gdsc.esca.eac_egfr.R`) |
| `Figure 8/`  | (Phospho)proteome and kinase-activity (INKA) analyses (`Global_proteome.R`, `Kinome.R`, `8C.INKA_viability_plot_spearman.R`) |
| `Figure S2/` | Single-cell processing setup and reference mapping (`setup.R`, `S2B.UMAP_HPCA.R`) |
| `Figure S3/` | Single-cell differential expression and marker plots (`setup.R`, `S3A.volcanoPlot.R`, `S3C.vlnPlots.R`) |

`setup.R` in the supplementary folders loads the single-cell environment shared by several notebooks.

## Requirements

- **R ≥ 4.2** with **Bioconductor** (installed via `BiocManager`).
- Analyses were run in a mix of base-R scripts (`.R`) and R Markdown notebooks (`.Rmd`); rendered `.html`/`.nb.html` outputs are included for reference.

Key packages by domain (see the `library()` calls at the top of each script / `setup.R` for the complete list):

- **Single-cell:** Seurat, SeuratObject, SingleCellExperiment, scater, SingleR, celldex, harmony, clustree, AUCell, MAST, progeny, Nebulosa, SCpubr, scCustomize
- **Clustering / networks:** ConsensusClusterPlus, diceR, WeightedCluster, RTN, igraph, ggraph, MCL
- **Enrichment / genomics:** limma, fgsea, clusterProfiler, org.Hs.eg.db, msigdbr, GEOquery
- **Data handling & plotting:** tidyverse, data.table, ComplexHeatmap, pheatmap, ggplot2 (+ ggpubr, ggrepel, ggbeeswarm, cowplot, patchwork), viridis, RColorBrewer, circlize

A convenience installer for the CRAN/Bioconductor packages:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "Seurat","SingleCellExperiment","scater","SingleR","celldex","harmony",
  "clustree","AUCell","MAST","progeny","Nebulosa",
  "ConsensusClusterPlus","diceR","RTN","igraph","ggraph",
  "limma","fgsea","clusterProfiler","org.Hs.eg.db","msigdbr","GEOquery",
  "tidyverse","data.table","ComplexHeatmap","pheatmap","ggpubr","cowplot","patchwork"
))
```

## Data availability

The scripts read processed data objects that are **not** included in this repository because of size and data-sharing policies. Source data underlying the figures are provided with the published article, and primary datasets are deposited as described in the manuscript's Data Availability statement:

- Single-cell and bulk RNA-seq: `<GEO accession — to add>`
- (Phospho)proteomics: `<PRIDE/ProteomeXchange accession — to add>`
- Public reference datasets: GDSC, TCGA-ESCA, CPCT (see manuscript for accessions)

Update the placeholders above with the final accession numbers before making the repository public.

## Reproducing the figures

1. Obtain the processed data objects (see Data availability) and place them where each script expects them (paths are defined near the top of each script).
2. Install the packages listed above.
3. Run the script or knit the notebook in the corresponding figure folder.

## Citation

Please cite the manuscript above. A `CITATION.cff` can be added once the DOI is assigned.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).

## Contact

Bijlsma Lab, Laboratory of Experimental Oncology and Radiobiology, Amsterdam UMC.
Correspondence: Maarten F. Bijlsma (m.f.bijlsma@amsterdamumc.nl).
