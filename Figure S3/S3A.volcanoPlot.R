library(cowplot)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)

setwd("Z:/BijlsmaTeam/DenivanSchie/20241001_MarkEACSubtypingStudy/Deni Rebuttal")

source(file = "setup.R", local = TRUE)

# Create metadata column indicating Epithelial vs Rest
seurat_integrated$Epithelial <- if_else(grepl("tumor", seurat_integrated$cell_type, ignore.case = TRUE), 
                                   "Epithelial", 
                                   "Rest"
)
table(seurat_integrated$cell_type, seurat_integrated$Epithelial)

# Tumor intrinsic signature
top250genes <- read_excel("data/Table S3 - Top 250 Signature genes.xlsx", skip = 1)[-1]
tumorintrinsic_genes <- c(top250genes$`Sig 1 - Mesenchymal-like`, top250genes$`Sig 3 - Intestinal-like`)
ml_genes <- c(top250genes$`Sig 1 - Mesenchymal-like`)
il_genes <- c(top250genes$`Sig 3 - Intestinal-like`)

# DEG Epithelial vs Rest
Idents(seurat_integrated) <- seurat_integrated$Epithelial

markers_df <- FindAllMarkers(object = seurat_integrated, only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0, test.use="MAST")

# Write the markers to a CSV file
markers_file <- file.path(output_folder, paste0("allCellsmarkers_MAST_epithelialvsrestnocutoff.csv"))
write.csv(markers_df, file = markers_file, row.names = FALSE)

# Check overlap with tumor intrinsic genes and Epithelial cluster
epithelial_up <- markers_df[markers_df$cluster == "Epithelial", ]$gene
overlap_tumorintrinsic_genes <- intersect(epithelial_up, tumorintrinsic_genes)
overlap_il_genes <- intersect(epithelial_up, il_genes)
overlap_ml_genes <- intersect(epithelial_up, ml_genes)

# Get the overlapping markers for epithelial and write to CSV
subset_markers_df <- markers_df[markers_df$gene %in% c(overlap_il_genes, overlap_ml_genes), ]
subset_markers_df$signature <- ifelse(
  subset_markers_df$gene %in% overlap_il_genes, "IL", "ML"
)

# Define the output file path
output_file <- file.path(output_folder, "allCellsmarkers_MAST_epithelialvsrestOverlapnocutoff.csv")

# Save the resulting data frame to a CSV file
write.csv(subset_markers_df, file = output_file, row.names = FALSE)

DEGS <- read.csv("epithelium vs rest/allCellsmarkers_MAST_epithelialvsrestnocutoff.csv")
ML <- read_excel("Z:/BijlsmaTeam/DenivanSchie/20241001_MarkEACSubtypingStudy/Mark Analyses/Underlying datasets/Table S3 - Top 250 Signature genes.xlsx", skip = 1)$`Sig 1 - Mesenchymal-like`
IL <- read_excel("Z:/BijlsmaTeam/DenivanSchie/20241001_MarkEACSubtypingStudy/Mark Analyses/Underlying datasets/Table S3 - Top 250 Signature genes.xlsx", skip = 1)$`Sig 3 - Intestinal-like`

# Add a column to indicate significance (adjust threshold as needed) and signature
DEGS$significance <- ifelse(DEGS$p_val_adj < 0.05 & abs(DEGS$avg_log2FC) > 0.25, "Significant", "Not Significant")

p_cutoff   <- 0.05
fc_cutoff  <- 1.5   

DEGS <- DEGS %>%
  mutate(signature = case_when(
    gene %in% IL ~ "IL",
    gene %in% ML ~ "ML",
    TRUE         ~ NA_character_
  ))

DEGS_epithelial <- DEGS %>%
  filter(cluster == "Epithelial", abs(avg_log2FC) > 0.1, p_val_adj <= p_cutoff) %>%
  mutate(signature = factor(
    ifelse(is.na(signature), "None", signature),
    levels = c("None", "IL", "ML")
  )) %>%
  mutate(
    neglog10_padj = -log10(p_val_adj),
  ) %>% 
  arrange(signature == "None")

max_finite <- max(DEGS_epithelial$neglog10_padj[is.finite(DEGS_epithelial$neglog10_padj)], 
                  na.rm = TRUE)

# replace Inf with max finite
DEGS_epithelial <- DEGS_epithelial %>%
  mutate(
    neglog10_padj = ifelse(is.infinite(neglog10_padj), max_finite, neglog10_padj)
  )

# New df with labeled genes
# genes_to_label <- c("ITGA3", "LGALS4")
genes_to_label <- NULL

DEGS_labels <- DEGS_epithelial %>%
  filter(gene %in% genes_to_label)

# Separate into 2 groups: one with NA ("None"), one with labeled ("IL", "ML")
DEGS_none <- DEGS_epithelial %>% filter(signature == "None")
DEGS_labeled <- DEGS_epithelial %>% filter(signature != "None") %>%
  mutate(signature = factor(signature, levels = c("IL", "ML")))

# Plot
p <- ggplot() +
  geom_point(
    data = DEGS_none,
    aes(x = avg_log2FC, y = neglog10_padj),
    color = "grey70", alpha = 0.2, size = 1
  ) +
  geom_point(
    data = DEGS_labeled,
    aes(x = avg_log2FC, y = neglog10_padj, color = signature),
    alpha = 1, size = 1
  ) +
  # geom_hline(yintercept = -log10(p_cutoff),        linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
  # geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
  scale_color_manual(values = c("IL" = "#4E79A7", "ML" = "#E15759")) +
  theme_classic() +
  labs(title = "", x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme(plot.title = element_text(hjust = 0.5)) + labs(color = "Signature") +
  coord_cartesian(ylim = c(0, max_finite)) +
  geom_label_repel(
    data = DEGS_labels,
    aes(
      x = avg_log2FC,
      y = neglog10_padj,
      label = gene
    ),
    size = 3,
    label.size = 0,              # no label border
    fill = "white",
    segment.color = "black",
    segment.size = 0.4,
    min.segment.length = 0,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = Inf
  )

# Add side labels below the plot area
p_annotated <- ggdraw(p) +
  draw_label("Upregulated in rest", x = 0.05, y = 0.04, hjust = 0, size = 11) +
  draw_label("Upregulated in epithelium", x = 0.95, y = 0.04, hjust = 1, size = 11)

ggsave(p_annotated, height=4,width=5.5,dpi=300, filename="epithelium vs rest/volcanoPlotDEGStumorvsepithelial_Cutoff_0border.pdf", useDingbats=FALSE)
