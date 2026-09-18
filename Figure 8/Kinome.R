library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)



for (type.of.enrichment in c("IMAC", "pTyr")) {
  for (type.of.analysis in c("Kinome", "INKA")) {
    for (scale in c("linear", "log2")) {
      for (cluster_distance in c("euclidean", "correlation")) {
        
        set.seed(648)
        
        #type.of.enrichment = "pTyr" ## IMAC or pTyr
        #type.of.analysis = "INKA" ## Kinome or INKA
        #scale = "linear" ## linear or log2
        #cluster_distance = "correlation" ## euclidean or correlation
        
        if (type.of.enrichment == "IMAC"){
          if (type.of.analysis == "Kinome"){
            ds = read.delim(
              "/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/PrOEF-230602-OPL2052-DB_MD-pTyrIP-Human.Esophageal.Adenocarcinoma.(EAC).Cell-IMAC-tp231116_Kinase_IMAC.txt",
              row.names = 1
            )
          }
          if (type.of.analysis == "INKA"){
            ds = read.delim(
              "/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/PrOEF-230602-OPL2052-DB_MD-pTyrIP-Human.Esophageal.Adenocarcinoma.(EAC).Cell-IMAC-tp231116_INKA_IMAC.txt",
              row.names = 1
            )
          }
        }
        if (type.of.enrichment == "pTyr"){
          if (type.of.analysis == "Kinome"){  
            ds = read.delim(
              "/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/PrOEF-230602-OPL2052-DB_MD-pTyrIP-Human.Esophageal.Adenocarcinoma.(EAC).Cell-IMAC-tp231115_Kinase_pTyr.txt",
              row.names = 1
            )
          }
          if (type.of.analysis == "INKA"){
            ds = read.delim(
              "/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/PrOEF-230602-OPL2052-DB_MD-pTyrIP-Human.Esophageal.Adenocarcinoma.(EAC).Cell-IMAC-tp231115_INKA_pTyr.txt",
              row.names = 1
            )    
          }
        }
        
        if (scale == "linear") {
          ds = t(ds)  
        }
        if (scale == "log2") {
          ds = t(log2(ds+1))
        }
        
        split_results <- lapply(rownames(ds), function(x) unlist(strsplit(x, "_")))
        result_matrix <- do.call(rbind, split_results)
        result_df <- as.data.frame(result_matrix)
        colnames(result_df) = c("type","sample.no","method","treatment","mol.subtype", "culture.medium", "cell.line.name","replicate.no")
        annotations = unique(result_df[,c(5,6,7)])
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
        
        # Calculate the variance for each protein
        protein_variances <- apply(ds_avg, 2, var)
        
        # Identify the top 100 proteins with the highest variances
        top_proteins <- names(head(sort(protein_variances, decreasing = TRUE), round(0.4*length(protein_variances))))
        
        # Remove non-top proteins from the dataset
        ds_top_proteins <- ds_avg[, top_proteins]
        
        # Retry PCA on the modified dataset
        pca_result <- prcomp(ds_top_proteins, scale. = TRUE)
        
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
          labs(title = paste0(type.of.enrichment, " ", type.of.analysis), x = paste0("Principal Component 1\n(",round(scree_data$Proportion_of_Variance[1],2)," variance explained)"), y = paste0("Principal Component 2\n(",round(scree_data$Proportion_of_Variance[2],2)," variance explained)")) + 
          theme_classic()+ theme(axis.text.x = element_text(color = "black"),
                                 axis.text.y = element_text(color = "black"))
        
        pdf(paste0("/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/Figure Output/",type.of.enrichment,".",type.of.analysis,".",scale,".pca.pdf"), width = 4.45, height = 3.15)
        print(pca)
        dev.off()
        
        ## Heatmap
        # color palette
        library(RColorBrewer)
        cols <- c("darkblue","blue", "white","red","darkred")
        newcol <- colorRampPalette(cols)
        col_breaks = seq(from=-4,to=4,by=0.1)
        ncols <- length(col_breaks)
        cols2 <- newcol(ncols)#apply the function to get 100 colours
        
        ph.colors= list(mol.subtype = c(Intestinal = "#4E79A7", Mesenchymal = "#E15759"))
        
        heatmap = pheatmap(t(scale(ds_top_proteins)),
                           show_rownames = F,
                           show_colnames = T,
                           annotation_col = annotations[,1,drop=F],
                           annotation_colors = ph.colors,
                           breaks=col_breaks,
                           cols2,
                           cutree_cols = 3,
                           cutree_rows = 5,
                           clustering_distance_cols = cluster_distance,
                           main = paste0(type.of.enrichment, " ", type.of.analysis))
        
        pdf(paste0("/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/Figure Output/",type.of.enrichment,".",type.of.analysis,".",scale,".",cluster_distance,".heatmap.pdf"), height = 4.45, width = 4.45)
        print(heatmap)
        dev.off()
        
        if (scale == "linear" && type.of.analysis == "INKA") {
          
          ds_avg$Total <- rowSums(ds_avg)
          # Create a new data frame with reordered rows for plotting
          ds_avg_reordered <- ds_avg[order(ds_avg$Total, decreasing = TRUE), ]
          
          # Create a bar plot
          plot.total <- ggplot(ds_avg_reordered, aes(x = Total, y = reorder(rownames(ds_avg_reordered), Total))) +
            geom_bar(stat = "identity") +
            labs(title = paste0("Sum of ", type.of.enrichment, " ", type.of.analysis), x = "Score", y = "") +
            theme_classic()
          ds_avg$Total = NULL
          
          pdf(paste0("/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/Figure Output/",type.of.enrichment,".",type.of.analysis,".totalScore.pdf"), width = 3, height = 3)
          print(plot.total)
          dev.off()
          
          ##
          # Convert row names to a column
          ds_avg <- ds_avg %>%
            mutate(row_name = row.names(ds_avg))
          ds_avg_reordered = ds_avg[c(3,8,5,6,1,2,4,7),]
          # Reshape data for ggplot
          ds_avg_long <- gather(ds_avg_reordered, key = "variable", value = "value", -row_name)
          
          # Create a list to store individual plots
          plots <- lapply(unique(ds_avg_long$row_name), function(row_name) {
            # Filter data for the specific row
            row_ds_avg <- ds_avg_long[ds_avg_long$row_name == row_name, ]
            
            # Order and select top 3 variables
            top_ds_avg <- row_ds_avg %>%
              arrange(desc(value)) %>%
              head(10)
            
            # Create ggplot for the current row
            ggplot(top_ds_avg, aes(x = value, y = reorder(variable, -value, decreasing = T))) +
              geom_bar(stat = "identity") +
              ggtitle(row_name)+
              xlab("Score") +
              ylab("")+
              theme_classic() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1,color = "black"),
                    axis.text.y = element_text(color = "black"))
            
          })
          
          # Display the individual plots
          print(plots)
          
          grid = cowplot::plot_grid(plotlist = plots, ncol=2, align = "v", axis = "l")
          
          # Arrange individual plots using gridExtra
          pdf(paste0("/Users/markdings/Documents/PhD/Scripts/eac.phosphoprot/batch2/Figure Output/",type.of.enrichment,".",type.of.analysis,".sepScore.pdf"), width = 5, height = 12)
          #gridExtra::grid.arrange(grobs = plots, ncol=2, top = paste0(type.of.enrichment, " ", type.of.analysis))
          print(grid)
          dev.off()
          
        }
      }
    }
  }
}
