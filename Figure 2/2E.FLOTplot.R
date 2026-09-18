## -----------------------------
## 1. Load libraries
## -----------------------------
library(Seurat)
library(readxl)
library(ggplot2)
# install.packages("ggExtra") if needed
library(ggExtra)
library(ggside)
library(dplyr)

# Setting directories and loading main dataset
wd <- getwd()
data_folder <- file.path("data/")
output_folder <- file.path("output/")
figures_folder <- file.path("figures")

## -----------------------------
## 2. Read in signature genes and subset Seurat data
## -----------------------------

seurat_integrated <- readRDS(file.path(data_folder, "soupAll_final.rds"))

plan(sequential)  # Disable parallelization

# Redoing variable genes and UMAP dimension reduction
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- FindVariableFeatures(seurat_integrated, selection.method = "vst", nfeatures = 3000, verbose = TRUE)
seurat_integrated <- ScaleData(seurat_integrated, verbose = TRUE, vars.to.regress = NULL)
seurat_integrated <- RunPCA(seurat_integrated, verbose = TRUE)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20, reduction = "pca")

saveRDS(seurat_integrated, file = file.path(data_folder, "seurat_integrated_preClustering.rds"))

data_folder <- "data"  # EDIT ME

seurat_integrated <- readRDS(file.path(data_folder, "seurat_integrated_preClustering.rds"))
epithelial_subset <- subset(
  seurat_integrated,
  subset = highLevelType == "Epithelial" & tissueType == "T"
)

NMF_sigs <- read_excel(
  file.path(data_folder, "Table S3 - Top 250 Signature genes.xlsx"),
  skip = 1
)

# Define each signature vector
NMF_ml         <- NMF_sigs$`Sig 1 - Mesenchymal-like`
NMF_cryptogenic <- NMF_sigs$`Sig 2 - Cryptogenic`
NMF_il         <- NMF_sigs$`Sig 3 - Intestinal-like`
NMF_stromal    <- NMF_sigs$`Sig 4 - Stroma`
NMF_squamous   <- NMF_sigs$`Sig 5 - Squamous`
NMF_immune     <- NMF_sigs$`Sig 6 - Immune`
NMF_hepatic    <- NMF_sigs$`Sig 7 - Hepatic`

# Assemble them into a list for easy iteration
signatures <- list(
  "NMF_ml"         = NMF_ml,
  "NMF_cryptogenic"= NMF_cryptogenic,
  "NMF_il"         = NMF_il,
  "NMF_stromal"    = NMF_stromal,
  "NMF_squamous"   = NMF_squamous,
  "NMF_immune"     = NMF_immune,
  "NMF_hepatic"    = NMF_hepatic
)

## -----------------------------
## 3. Add module scores + Z-scores on the subset
## -----------------------------
for (sig_name in names(signatures)) {
  # Filter signature genes to those present in the subset
  filtered_genes <- intersect(signatures[[sig_name]], rownames(epithelial_subset))
  
  # Add module score to the subset object
  epithelial_subset <- AddModuleScore(
    object   = epithelial_subset,
    features = list(filtered_genes),
    name     = paste0(sig_name, "_Score")
  )
  # By default, this will add columns like "NMF_ml_Score1" to epithelial_subset@meta.data
  
  # Calculate a Z-score (average scaled expression for these genes)
  expr_mat <- GetAssayData(epithelial_subset, assay = "RNA", slot = "data")[filtered_genes, ]
  
  if (!is.null(dim(expr_mat)) && nrow(expr_mat) > 1) {
    z_mat <- scale(t(as.matrix(expr_mat)))  # Scale on transposed expression matrix (cells as rows)
    # Save the rowMeans of scaled values in a new metadata column
    epithelial_subset[[paste0(sig_name, "_Zscore")]] <- rowMeans(z_mat, na.rm = TRUE)
  } else {
    # If only 1 or 0 genes remain after filtering
    epithelial_subset[[paste0(sig_name, "_Zscore")]] <- 0
  }
}

## -----------------------------
## 3. Add module scores + Z-scores on the full dataset
## -----------------------------
for (sig_name in names(signatures)) {
  # Filter signature genes to those present in the subset
  filtered_genes <- intersect(signatures[[sig_name]], rownames(seurat_integrated))
  
  # Add module score to the subset object
  seurat_integrated <- AddModuleScore(
    object   = seurat_integrated,
    features = list(filtered_genes),
    name     = paste0(sig_name, "_Score")
  )
  # By default, this will add columns like "NMF_ml_Score1" to seurat_integrated@meta.data
  
  # Calculate a Z-score (average scaled expression for these genes)
  expr_mat <- GetAssayData(seurat_integrated, assay = "RNA", slot = "data")[filtered_genes, ]
  
  if (!is.null(dim(expr_mat)) && nrow(expr_mat) > 1) {
    z_mat <- scale(t(as.matrix(expr_mat)))  # Scale on transposed expression matrix (cells as rows)
    # Save the rowMeans of scaled values in a new metadata column
    seurat_integrated[[paste0(sig_name, "_Zscore")]] <- rowMeans(z_mat, na.rm = TRUE)
  } else {
    # If only 1 or 0 genes remain after filtering
    seurat_integrated[[paste0(sig_name, "_Zscore")]] <- 0
  }
}

## -----------------------------
## 4. Prepare data for plotting
## -----------------------------
# Example: Plot "NMF_ml_Zscore" on X-axis vs "NMF_il_Zscore" on Y-axis
plot_df <- data.frame(
  cell       = colnames(epithelial_subset),
  Sig1       = epithelial_subset$NMF_ml_Score1,
  Sig3       = epithelial_subset$NMF_il_Score1,
  donor   = epithelial_subset$donor,  # or other grouping variable
  chemo = epithelial_subset$chemo
  )

# Create a new variable combining patient and condition for unique colors
plot_df <- plot_df %>%
  mutate(patient_condition = paste(donor, chemo, sep = "_"))

## -----------------------------
## 5. Create a scatter plot
## -----------------------------
# Get unique donors for each condition
pre_donors <- unique(plot_df$donor[plot_df$chemo == "pre"])
post_donors <- unique(plot_df$donor[plot_df$chemo == "post"])

# Generate color palettes with a broader range
pre_palette <- colorRampPalette(c("lightblue1", "dodgerblue4"))(length(pre_donors))
post_palette <- colorRampPalette(c("indianred1", "firebrick4"))(length(post_donors))

# Map colors to patient-condition combinations
colors <- c(
  setNames(pre_palette, unique(plot_df$donor[plot_df$chemo == "pre"])),
  setNames(post_palette, unique(plot_df$donor[plot_df$chemo == "post"]))
)

# Format legend labels with "Pre" and "Post" titles
legend_labels <- c(
  paste(pre_donors, "(Pre)"),
  paste(post_donors, "(Post)")
)

# Plot
p <- ggplot(plot_df, aes(x = Sig1, y = Sig3)) +
  geom_point(aes(color = donor), alpha = 0.7) +  # Adjust alpha for transparency
  scale_color_manual(
    values = colors,
    name = NULL,  # No overall title
    breaks = names(colors),  # Use names to match colors
    labels = legend_labels  # Custom labels
  ) +
  theme_classic() +
  theme(aspect.ratio = 1,
    legend.position = "bottom",
    legend.box = "vertical",  # Arrange legend vertically
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  guides(
    color = guide_legend(
      title = NULL,
      nrow = 2,  # Two rows to separate Pre and Post
      byrow = TRUE,
      override.aes = list(size = 4)  # Larger legend points
    )
  )

p

# Add marginal boxplots using ggMarginal
p_box <- ggMarginal(
  p,
  type       = "boxplot",
  groupFill  = TRUE,
  size       = 6, 
)

p_box

ggsave(plot=p_box, height=4,width=5.5,dpi=300, filename="dotplot.f1f3signature2_byPatient.pdf", useDingbats=FALSE)

write.csv(plot_df, "f1f3signature_byCell2.csv")

## -----------------------------
## 6. Create a dot plot with fness
## -----------------------------
p <- DotPlot(seurat_integrated, features = c("NMF_ml_Score1","NMF_cryptogenic_Score1","NMF_il_Score1","NMF_stromal_Score1","NMF_squamous_Score1","NMF_immune_Score1", "NMF_hepatic_Score1"), col.min=0) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  coord_flip()
p

ggsave(plot=p,height=5,width=8,dpi=200, filename="dotplot.signatures.pdf", useDingbats=FALSE)

## -----------------------------
## 6. Create a UMAP for different marker genes and cell types
## -----------------------------
# Define markers of interest
markers_list <- list(c(
  "KRT81", "MYC", "CAV1", "HMGA2", "FOSL1", "LIN28B", "CDX2", "KRT20", 
  "LGALS4", "ACTA2", "MYH11", "KRT5", "KRT6", "TP63", "MMP9", "CD83", 
  "CCBE1", "FOXA3", "CLDN14", "BTNL8", "L1CAM", "MUC17"
))


# Function to split the marker list into batches of a given size
split_into_batches <- function(marker_list, batch_size) {
  split(marker_list, ceiling(seq_along(marker_list) / batch_size))
}

# Define the batch size (e.g., 10 markers per batch)
batch_size <- 10
batches <- split_into_batches(unlist(markers_list), batch_size)

# Iterate over batches
for (batch_idx in seq_along(batches)) {
  # Get the current batch of markers
  markers_batch <- batches[[batch_idx]]
  
  # Check if markers exist in the expression data or metadata
  valid_markers <- markers_batch[markers_batch %in% rownames(seurat_integrated@assays$RNA@data)]
  
  if (length(valid_markers) == 0) {
    message("No valid markers found for batch ", batch_idx)
    next  # Skip this batch if no markers are valid
  }
  
  # Generate FeaturePlots for the batch
  feature_plots <- FeaturePlot(
    object = seurat_integrated, 
    features = valid_markers, 
    cols = c("lightgrey", "black"), 
    ncol = 2, 
    pt.size = 0.05  # Adjust point size as needed
  )
  
  # Apply themes to all plots in the batch
  feature_plots <- lapply(feature_plots, function(plot) {
    plot + 
      theme(
        aspect.ratio = 1,  # Ensure square aspect ratio
        axis.title.x = element_text(size = 12),  # Customize axis title size
        axis.title.y = element_text(size = 12),  # Customize axis title size
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank()  # Remove axis ticks
      ) +
      labs(x = "UMAP1", y = "UMAP2")  # Add axis labels
  })
  
  # Combine all FeaturePlots in the batch
  combined_plot <- wrap_plots(feature_plots, ncol = 5) + 
    plot_annotation(
      title = paste("FeaturePlots - Batch", batch_idx),
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
      )
    )
  
  # Print the combined plot for the batch
  print(combined_plot)
  
  # Save the combined plot
  ggsave(
    combined_plot, 
    filename = paste0("featurePlots_batch_", batch_idx, ".png"), 
    width = 20, height = 12, dpi = 300
  )
}

dim_plot <- DimPlot(seurat_integrated) + 
  ggtitle("All cells, colored by cell type") +
  theme(aspect.ratio = 1, axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = "UMAP1", y = "UMAP2")
ggsave(dim_plot, file=file.path("allCells_dimPlot_byCellType.png"), width=8, height=8, dpi = 300)