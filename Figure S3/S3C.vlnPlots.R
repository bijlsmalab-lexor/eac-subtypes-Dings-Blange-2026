library(cowplot)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)

setwd("Z:/BijlsmaTeam/DenivanSchie/20241001_MarkEACSubtypingStudy/Deni Rebuttal")

source(file = "setup.R", local = TRUE)

# Epithelial vs rest
seurat_integrated$Epithelial <- if_else(grepl("tumor", seurat_integrated$cell_type, ignore.case = TRUE), 
                                   "Epithelial", 
                                   "Rest"
)

# Fibroblast vs rest
seurat_integrated$Fibroblast <- if_else(grepl("Fibroblast", seurat_integrated$cell_type, ignore.case = TRUE), 
                                   "Fibroblast", 
                                   "Rest"
)

# Immune vs rest
seurat_integrated$Immune <- if_else(
  grepl("\\bT\\b|\\bMacrophage\\b", seurat_integrated$cell_type, ignore.case = TRUE),
  "Immune",
  "Rest"
)

# Generate pallette
pal_manu <- c(
  Epithelial  = "#1D79B5",
  Fibroblast  = "#73BB68",
  Immune      = "#96C2E0",
  Endothelial = "#ED5456"
)

# Create functions to reduce bloat
cols_for <- function(obj, group.by, highlight, pal=pal_manu, grey="#D3D3D3"){
  lev <- levels(factor(FetchData(obj, group.by)[,1]))
  cols <- setNames(rep(grey, length(lev)), lev)
  keep <- intersect(highlight, names(pal))
  cols[keep] <- pal[keep]
  cols
}

vln_highlight <- function(obj, feature, group.by, cols,
                          pt_alpha=1, pt_size=3, stroke=.4){
  p <- VlnPlot(obj, features=feature, group.by=group.by, pt.size=0) +
    scale_fill_manual(values=cols) + theme(legend.position="none", aspect.ratio=1)
  
  df <- FetchData(obj, vars=c(feature, group.by))
  names(df) <- c("y","grp"); df$grp <- factor(df$grp, levels=names(cols))
  
  p <- p + geom_jitter(data=df, aes(x=grp, y=y, fill=grp),
                       width=.15, size=pt_size, alpha=pt_alpha,
                       shape=21, color="black", stroke=stroke, inherit.aes=FALSE) +
    scale_fill_manual(values=cols)
  p$layers <- rev(p$layers)  # violin above points
  p
}

vln_highlight_multi <- function(obj, features, group.by, highlight,
                                pal=pal_manu, grey="#D3D3D3", ncol=NULL){
  cols <- cols_for(obj, group.by, highlight, pal, grey)
  plots <- lapply(features, function(f) vln_highlight(obj, f, group.by, cols))
  wrap_plots(plots, ncol = ncol %||% length(features))
}

# ---- balance helper: equal # cells per group ----
balance_groups <- function(obj, group.by, n = NULL, seed = 123){
  grp <- FetchData(obj, vars = group.by)[,1]
  names(grp) <- colnames(obj)
  
  # drop NA groups (if any)
  keep_ok <- !is.na(grp)
  grp <- grp[keep_ok]
  obj <- obj[, names(grp)]
  
  tab <- table(grp)
  stopifnot(length(tab) >= 2)
  
  # target per-group size
  if (is.null(n)) {
    n <- rep(min(tab), length(tab))
  } else if (length(n) == 1) {
    n <- rep(as.integer(n), length(tab))
  } else if (length(n) != length(tab)) {
    stop("n must be length 1 or length equal to number of groups")
  }
  names(n) <- names(tab)
  n <- pmin(n, as.integer(tab))   # cap at available per group
  
  set.seed(seed)
  keep_cells <- unlist(lapply(names(tab), function(g){
    cells_g <- names(grp)[grp == g]
    sample(cells_g, n[[g]])
  }))
  obj[, keep_cells]
}

# Plot
## ── genes ──
epi_genes <- c("ITGA3","LGALS4","KRT5","EPCAM")
genes_all <- c(epi_genes, "ACTA2","MMP9","PECAM1")

## ── Epithelial vs Rest: singles + panel (others grey) ──
for(g in epi_genes){
  p <- vln_highlight(seurat_integrated, g, "Epithelial",
                     cols_for(seurat_integrated, "Epithelial", "Epithelial"))
  ggsave(file.path("figures/comparison_compartmentPlots", paste0("vlnPlot_", g, "_Epithelial_pt.pdf")),
         p, width=6, height=5)
}
p_epi <- vln_highlight_multi(seurat_integrated, epi_genes, "Epithelial", "Epithelial")
ggsave("figures/comparison_compartmentPlots/vlnPanel_Epithelial_pt.pdf", p_epi, width=20, height=5)

## ── Fibroblast vs Rest (others grey) ──
p_fb <- vln_highlight(seurat_integrated, "ACTA2", "Fibroblast",
                      cols_for(seurat_integrated, "Fibroblast", "Fibroblast"))
ggsave("figures/comparison_compartmentPlots/vlnPlot_ACTA2_Fibroblast_pt.pdf", p_fb, width=6, height=5)

## ── Immune vs Rest (others grey) ──
p_im <- vln_highlight(seurat_integrated, "MMP9", "Immune",
                      cols_for(seurat_integrated, "Immune", "Immune"))
ggsave("figures/comparison_compartmentPlots/vlnPlot_MMP9_Immune_pt.pdf", p_im, width=6, height=5)

# ## ── Compartment panel (all compartments in manuscript colors) ──
# p_comp_all <- vln_highlight_multi(seurat_integrated, genes_all, "Compartment",
#                                   highlight = names(pal_manu))  # no grey; use palette for all
# ggsave("figures/comparison_compartmentPlots/vlnPlots_compartment_colored.pdf", p_comp_all, width=20, height=5)
# 
# ## ── Compartment panels highlighting one group at a time (others grey) ──
# for(h in names(pal_manu)){
#   p <- vln_highlight_multi(seurat_integrated, genes_all, "Compartment", highlight = h)
#   ggsave(file.path("figures/comparison_compartmentPlots",
#                    paste0("vlnPlots_compartment_highlight_", h, ".pdf")),
#          p, width=20, height=5)
# }

# Balanced groups
## Epithelial vs Rest
epi_eq <- balance_groups(seurat_integrated, "Epithelial")  # or n = 1000
for(g in epi_genes){
  p <- vln_highlight(epi_eq, g, "Epithelial",
                     cols_for(epi_eq, "Epithelial", "Epithelial"))
  ggsave(file.path("figures/comparison_compartmentPlots",
                   paste0("vlnPlot_", g, "_Epithelial_pt_balanced.pdf")),
         p, width=6, height=5)
}
p_epi <- vln_highlight_multi(epi_eq, epi_genes, "Epithelial", "Epithelial")
ggsave("figures/comparison_compartmentPlots/vlnPanel_Epithelial_pt.pdf",
       p_epi, width=20, height=5)

## Fibroblast vs Rest
fb_eq <- balance_groups(seurat_integrated, "Fibroblast")
p_fb <- vln_highlight(fb_eq, "ACTA2", "Fibroblast",
                      cols_for(fb_eq, "Fibroblast", "Fibroblast"))
ggsave("figures/comparison_compartmentPlots/vlnPlot_ACTA2_Fibroblast_pt_balanced.pdf",
       p_fb, width=6, height=5)

## Immune vs Rest
im_eq <- balance_groups(seurat_integrated, "Immune")
p_im <- vln_highlight(im_eq, "MMP9", "Immune",
                      cols_for(im_eq, "Immune", "Immune"))
ggsave("figures/comparison_compartmentPlots/vlnPlot_MMP9_Immune_pt_balanced.pdf",
       p_im, width=6, height=5)

# ============================================================
# Dot-only plots (balanced groups, zeros kept) + statistics
#   - no violin: same styling as the _balanced vln plots
#     (VlnPlot base + geom_jitter), violin layer dropped
#   - ALL cells kept, including expression == 0 (zeros contribute
#     to the statistics and to the point cloud)
#   - Wilcoxon test per gene; stars overlaid on plot, exact
#     p-values / effect sizes exported as a combined stats table
# ============================================================

# Output directory for stats tables
stats_dir <- file.path(output_folder, "comparison_compartmentStats")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

## ── dot-only plot for one feature: same look as _balanced vln plots ──
##    (built on the VlnPlot base like vln_highlight, then the violin
##     layer is dropped so only the styled dots remain; zeros kept)
dot_highlight <- function(obj, feature, group.by, cols,
                          pt_alpha = 1, pt_size = 3, stroke = .4,
                          test = "wilcox.test", stat_label = "p.signif"){
  df <- FetchData(obj, vars = c(feature, group.by))
  names(df) <- c("y", "grp"); df$grp <- factor(df$grp, levels = names(cols))
  df <- df[!is.na(df$grp), ]                    # keep zeros (non-detected cells included)
  
  # same base as vln_highlight so theme / axes match exactly
  p <- VlnPlot(obj, features = feature, group.by = group.by, pt.size = 0) +
    theme(legend.position = "none", aspect.ratio = 1)
  
  # drop the violin layer, then re-root the plot on df so every
  # downstream layer (jitter + the internal ggsignif layer that
  # stat_compare_means adds) reads y/grp instead of the feature name.
  # NOTE: keep fill OUT of the global mapping — a global fill = grp
  # splits the data by group, so ggsignif can't see both groups in one
  # subset and silently draws no bracket. fill lives on the jitter only.
  p$layers  <- p$layers[!vapply(p$layers, function(l) inherits(l$geom, "GeomViolin"),
                                logical(1))]
  p$data    <- df
  p$mapping <- aes(x = grp, y = y)
  
  p <- p + geom_jitter(aes(fill = grp),
                       width = .15, size = pt_size, alpha = pt_alpha,
                       shape = 21, color = "black", stroke = stroke) +
    scale_fill_manual(values = cols)
  
  grps <- levels(droplevels(df$grp))
  comparisons <- if (length(grps) >= 2) combn(grps, 2, simplify = FALSE) else list()
  if (length(comparisons) > 0) {
    p <- p + stat_compare_means(comparisons = comparisons, method = test,
                                label = stat_label, tip.length = .01) +
      # headroom so the significance bracket isn't clipped at the top
      scale_y_continuous(expand = expansion(mult = c(.05, .15)))
  }
  p
}

## ── tidy stats table for one feature (balanced, zeros kept) ──
dot_stats_table <- function(obj, feature, group.by, comparison = group.by, balanced = TRUE){
  df <- FetchData(obj, vars = c(feature, group.by))
  names(df) <- c("y", "grp")
  df <- df[!is.na(df$grp), ]                   # keep zeros (non-detected cells included)
  df$grp <- droplevels(factor(df$grp))
  if (nlevels(df$grp) < 2) return(NULL)
  
  res  <- df %>% rstatix::wilcox_test(y ~ grp) %>% rstatix::add_significance()
  eff  <- df %>% rstatix::wilcox_effsize(y ~ grp)
  res  <- dplyr::left_join(res, eff[, c("group1", "group2", "effsize", "magnitude")],
                           by = c("group1", "group2"))
  
  n_tab  <- as.data.frame(table(df$grp), responseName = "n")
  med    <- tapply(df$y, df$grp, median)
  mn     <- tapply(df$y, df$grp, mean)
  n_str  <- paste(sprintf("%s=%d", n_tab$Var1, n_tab$n), collapse = "; ")
  med_str<- paste(sprintf("%s=%.3f", names(med), med), collapse = "; ")
  mn_str <- paste(sprintf("%s=%.3f", names(mn),  mn),  collapse = "; ")
  
  # direction of the difference: which group is higher (by mean expression),
  # with the per-comparison log2 fold change (pseudocount to keep zeros finite)
  hi  <- ifelse(mn[res$group1] > mn[res$group2], res$group1,
                ifelse(mn[res$group2] > mn[res$group1], res$group2, "tie"))
  l2fc<- log2((mn[res$group1] + 1e-9) / (mn[res$group2] + 1e-9))
  
  res$feature      <- feature
  res$comparison   <- comparison
  res$balanced     <- balanced
  res$n            <- n_str
  res$median       <- med_str
  res$mean         <- mn_str
  res$higher_in    <- unname(hi)
  res$log2FC_g1_g2 <- unname(round(l2fc, 3))
  res
}

## ── run for each compartment comparison ──
dot_specs <- list(
  list(obj = epi_eq, group = "Epithelial", genes = epi_genes),
  list(obj = fb_eq,  group = "Fibroblast", genes = "ACTA2"),
  list(obj = im_eq,  group = "Immune",     genes = "MMP9")
)

stats_all <- list()
for (spec in dot_specs) {
  cols <- cols_for(spec$obj, spec$group, spec$group)
  for (g in spec$genes) {
    # dot plot with stats stars
    p <- dot_highlight(spec$obj, g, spec$group, cols)
    ggsave(file.path("figures/comparison_compartmentPlots",
                     paste0("dotPlot_", g, "_", spec$group, "_balanced.pdf")),
           p, width = 6, height = 5)
    
    # stats table
    tbl <- dot_stats_table(spec$obj, g, spec$group)
    if (!is.null(tbl)) stats_all[[paste(spec$group, g, sep = "_")]] <- tbl
  }
}

## ── combine + export stats output ──
dot_stats_combined <- dplyr::bind_rows(stats_all)
dot_stats_combined <- dplyr::relocate(dot_stats_combined,
                                      comparison, feature, group1, group2,
                                      higher_in, log2FC_g1_g2)
write.csv(dot_stats_combined,
          file.path(stats_dir, "dotPlot_compartment_wilcox_stats.csv"),
          row.names = FALSE)
print(dot_stats_combined)
