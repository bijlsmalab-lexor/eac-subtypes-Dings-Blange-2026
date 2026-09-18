library(cowplot)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)

setwd("Z:/BijlsmaTeam/DenivanSchie/20241001_MarkEACSubtypingStudy/Deni Rebuttal")

source(file = "setup.R", local = TRUE)

# HPCA
# Load in reference dataset
hpca.se <- celldex::HumanPrimaryCellAtlasData()

# Running SingleR on hpca
pred <- SingleR(test = seurat_integrated@assays$RNA@counts, ref = hpca.se, assay.type.test=1,
                labels = hpca.se$label.main)

# Adding labels to dataset
seurat_integrated$singleR_hpca <- pred$pruned.labels

# # Inspect quality of the predictions
# plotScoreHeatmap(pred)
# plotDeltaDistribution(pred, ncol = 4, dots.on.top = FALSE)

# Generate UMAP colored by predicted cell type
# Define a custom color palette for similar cell types
custom_palette <- c(
  "Endothelial_cells" = "#E31A1C",      # Red for endothelial cells
  "Fibroblasts" = "#CAB2D6",           # Light purple for fibroblasts
  "T_cells" = "#FF7F00",               # Orange for T cells
  "Macrophage" = "#FDBF6F",            # Yellowish-orange for macrophages
  "Epithelial_cells" = "#A6CEE3",      # Light blue for epithelial cells
  "Stem_cells" = "#B2DF8A",            # Light green for stem cells
  "Other" = "#D9D9D9"                  # Grey for miscellaneous/NA
)

p <- DimPlot(object = seurat_integrated, group.by = "singleR_hpca", label = FALSE) + 
  ggtitle("All cells, colored by predicted cell type") +
  theme(aspect.ratio = 1, axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2") + scale_color_manual(values = custom_palette) + labs(colour="Cell type")
p

# Save figure
ggsave(p, file=file.path(figures_folder, "allCells_dimPlot_byPredictedcellType_HPCA.png"), width=8, height=8)