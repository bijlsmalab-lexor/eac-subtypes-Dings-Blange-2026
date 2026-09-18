library(ggplot2)
library(ggrepel)
library(pheatmap)

set.seed(648)

ds = read.delim("/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/PrOEF-230602-OPL2052-DB_MD-pTyrIP-Human.Esophageal.Adenocarcinoma.(EAC).Cell-expression-tp231116_global_normCount.txt", row.names = 1)
ds = t(log2(ds+1))

split_results <- lapply(rownames(ds), function(x) unlist(strsplit(x, "_")))
result_matrix <- do.call(rbind, split_results)
result_df <- as.data.frame(result_matrix)
colnames(result_df) = c("type","method","treatment","mol.subtype", "culture.medium", "cell.line.name","replicate.no")
annotations = unique(result_df[,c(4,5,6)])
rownames(annotations) = annotations$cell.line.name

ds = cbind(result_df, ds)

ds_avg <- ds %>%
  filter(treatment != "Tx") %>%
  group_by(cell.line.name) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
ds_avg = as.data.frame(ds_avg)
rownames(ds_avg) = ds_avg$cell.line.name
ds_avg$cell.line.name = NULL

# Perform PCA
# Set a threshold for near-zero variance
variance_threshold <- 1

# Identify constant and near-zero variance columns
low_var_columns <- apply(ds_avg, 2, function(x) var(x) <= variance_threshold)
remove_columns <- which(colMeans(ds_avg) == 0 | low_var_columns)

# Remove identified columns from the dataset
ds_no_low_var <- ds_avg[, -remove_columns]

# Retry PCA on the modified dataset
pca_result <- prcomp(ds_no_low_var, scale. = TRUE)

#pca_result <- prcomp(ds_avg, scale. = TRUE)
# Summary of PCA
summary(pca_result)
# Scree plot to visualize the proportion of variance explained by each principal component
scree_data <- data.frame(
  PC = 1:ncol(pca_result$rotation),
  Proportion_of_Variance = pca_result$sdev^2 / sum(pca_result$sdev^2)
)

ggplot(scree_data, aes(x = PC, y = Proportion_of_Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion of Variance")

# Biplot to visualize the loadings and scores of each observation
biplot_data <- as.data.frame(pca_result$x)
biplot_data = merge(x = biplot_data, y = annotations, by.x = 0, by.y = 3)

pca = ggplot(biplot_data, aes(x = PC1, y = PC2, color = mol.subtype, label = Row.names)) +
  geom_text_repel() +
  geom_point() +
  scale_color_manual(values = c("Intestinal" = "#4E79A7", "Mesenchymal" = "#E15759")) +
  labs(title = "Global proteome",x = paste0("Principal Component 1\n(",round(scree_data$Proportion_of_Variance[1],2)," variance explained)"), y = paste0("Principal Component 2\n(",round(scree_data$Proportion_of_Variance[2],2)," variance explained)")) +
  theme_classic()+ theme(axis.text.x = element_text(color = "black"),
                         axis.text.y = element_text(color = "black"))
pca
pdf("/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/Figure Output/global.proteome.pca.pdf", width = 4.45, height = 3.15)
print(pca)
dev.off()

## Heatmap
# color palette
library(RColorBrewer)
cols <- c("darkblue","blue", "white","red","darkred")
newcol <- colorRampPalette(cols)
col_breaks = seq(from=-3,to=3,by=0.1)
ncols <- length(col_breaks)
cols2 <- newcol(ncols)#apply the function to get 100 colours

ph.colors= list(mol.subtype = c(Intestinal = "#4E79A7", Mesenchymal = "#E15759"))

cluster_distance = "correlation" ## euclidean or correlation

heatmap = pheatmap(t(scale(ds_no_low_var)),
         show_rownames = F,
         show_colnames = T,
         annotation_col = annotations[,1,drop=F],
         annotation_colors = ph.colors,
         breaks=col_breaks,
         cols2,
         cutree_cols = 2,
         cutree_rows = 3,
         clustering_distance_cols = cluster_distance)

pdf(paste0("/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/Figure Output/global.proteome.heatmap.",cluster_distance,".pdf"), height = 4.45, width = 4.45)
print(heatmap)
dev.off()

