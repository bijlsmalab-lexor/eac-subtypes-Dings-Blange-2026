library(data.table)
library(ConsensusClusterPlus)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(reshape2)
require(gridExtra)

cols <- c("darkblue","blue", "white","red","darkred")
newcol <- colorRampPalette(cols)
col_breaks = seq(from=-3,to=3,by=0.1)
ncols <- length(col_breaks)
cols2 <- newcol(ncols)#apply the function to get 100 colours

ph.colors= list(Dings.km_subtype = c(IL = "#4E79A7", ML = "#E15759"),
                Class = c(IL = "#4E79A7", ML = "#E15759"))


ds.expression = as.data.frame(fread("/Users/markdings/Documents/PhD/Scripts/PanGI_Network_analysis/datasets/GDSC/Cell_line_RMA_proc_basalExp.txt"))
ds.expression = ds.expression[(!duplicated(ds.expression$GENE_SYMBOLS)),]
rownames(ds.expression) = ds.expression$GENE_SYMBOLS
ds.expression$GENE_SYMBOLS=NULL
ds.expression$GENE_title=NULL
colnames(ds.expression) = sapply(strsplit(colnames(ds.expression), split='.', fixed=TRUE), function(x) (x[2]))

ds.info = fread("/Users/markdings/Documents/PhD/Scripts/PanGI_Network_analysis/datasets/GDSC/PANCANCER_Info.txt")

## align samples
ds.info$cosmic_id = as.character(ds.info$cosmic_id)

ds.expression = ds.expression[,colnames(ds.expression) %in% ds.info$cosmic_id]
ds.info = ds.info[ds.info$cosmic_id %in% colnames(ds.expression),]
ds.expression = ds.expression[,colnames(ds.expression) %in% ds.info$cosmic_id]
ds.expression = ds.expression[,match(table=colnames(ds.expression),x=ds.info$cosmic_id)]
ds.info = ds.info[match(table=ds.info$cosmic_id,x=colnames(ds.expression)),]

index = ds.info$Tissue_descriptor_3=="oesophagus.eac"

ds.expression = ds.expression[,index]
ds.info = ds.info[index,]
colnames(ds.expression) = ds.info$sample_name

## Classify
## load gene signature
gene.sig1 = read.csv("/Users/markdings/Documents/PhD/Scripts/PanGI_Network_analysis/classification/classification.by.Dings.eacZalm6/signatureGene.feature1.txt")
gene.sig2 = read.csv("/Users/markdings/Documents/PhD/Scripts/PanGI_Network_analysis/classification/classification.by.Dings.eacZalm6/signatureGene.feature3.txt")
gene.sig = rbind(gene.sig1[1:50,2,drop=F], gene.sig2[1:50,2,drop=F])
colnames(gene.sig) = "Genes"
rownames(gene.sig) = gene.sig$Genes
gene.sig$Class = "NULL"
gene.sig$Class[1:50] = "ML"
gene.sig$Class[51:100] = "IL"
## 
expression.set = ds.expression
expression.set = expression.set[toupper(rownames(expression.set)) %in% toupper(gene.sig$Genes),]
gene.sig = gene.sig[toupper(gene.sig$Genes) %in% toupper(rownames(expression.set)), ,drop=F]
expression.set = as.data.frame(t(scale(t(expression.set))))

expression.set = expression.set[match(toupper(gene.sig$Genes), toupper(rownames(expression.set))),]

set.seed(100)
clusterResult = ConsensusClusterPlus(as.matrix(expression.set), maxK = 3, clusterAlg = "km",reps=500, pItem = 0.5)
k=2
clusterK = as.data.frame(clusterResult[[k]]$consensusClass)
colnames(clusterK) = c("Dings.km_subtype")
clusterK$Dings.km_subtype = paste0("subtype", clusterK$Dings.km_subtype)

pheatmap(expression.set[,order(clusterK$Dings.km_subtype)], annotation_col = clusterK,
         breaks=col_breaks,
         annotation_row = gene.sig[,2,drop = F],
         show_rownames = T,
         cluster_rows = F,
         cluster_cols = F,
         cols2)

clusterK$Dings.km_subtype <- gsub('subtype1', 'IL', clusterK$Dings.km_subtype)
clusterK$Dings.km_subtype <- gsub('subtype2', 'ML', clusterK$Dings.km_subtype)

p=pheatmap(expression.set[,order(clusterK$Dings.km_subtype)], annotation_col = clusterK,
           breaks=col_breaks,
           annotation_row = gene.sig[,2,drop = F],
           show_rownames = F,
           cluster_rows = F,
           cluster_cols = F,
           border_color=NA,
           cols2,
           annotation_colors = ph.colors)

#pdf("/Users/markdings/Documents/PhD/Scripts/PanGI_Network_analysis/output.classification/dings/gdsc.esca.eac.kMeans.100genes.heatmap.pdf", width=4.95, height = 3.57)
p
#dev.off()
prop.table(table(clusterK$Dings.km_subtype))
clusterK

ds.info=cbind(ds.info,clusterK)

#write.table(clusterK, file=paste0("/Users/markdings/Documents/PhD/Scripts/PanGI_Network_analysis/output.classification/dings/gdsc.esca.eac.kMeans.100genes.txt"), append = FALSE, sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)

######## Drug data analysis:
ds.gdsc1 = as.data.frame(fread("/Users/markdings/Documents/PhD/Scripts/PanGI_Network_analysis/datasets/GDSC/GDSC1_fitted_dose_response_25Feb20.txt"))
ds.gdsc1 = ds.gdsc1[ds.gdsc1$CELL_LINE_NAME %in% ds.info$sample_name,]
ds.gdsc1 = merge(x=ds.gdsc1, by.x = "CELL_LINE_NAME", y=ds.info, by.y = "sample_name")
ds.gdsc1$DRUG_NAME = as.character(ds.gdsc1$DRUG_NAME)
ds.gdsc1$LN_IC50 = as.numeric(ds.gdsc1$LN_IC50)
ds.gdsc1$AUC = as.numeric(ds.gdsc1$AUC)
ds.gdsc1$Z_SCORE = as.numeric(ds.gdsc1$Z_SCORE)

names(table(as.character(ds.gdsc1$DRUG_NAME))>4)

ds.gdsc1 = ds.gdsc1[ds.gdsc1$DRUG_NAME %in% names(which(table(as.character(ds.gdsc1$DRUG_NAME))>4)),]

stat.test <- ds.gdsc1 %>%
  group_by(DRUG_NAME) %>%
  t_test(Z_SCORE ~ Dings.km_subtype)

drug_subtype_means <- ds.gdsc1 %>%
  group_by(DRUG_NAME, Dings.km_subtype) %>%
  dplyr::summarize(Mean = mean(Z_SCORE, na.rm=TRUE))
drug_subtype_wide =dcast(drug_subtype_means, formula = DRUG_NAME ~Dings.km_subtype, value.var = "Mean")
drug_subtype_wide$delta = drug_subtype_wide$IL - drug_subtype_wide$ML
stat.test = merge(stat.test, drug_subtype_wide, by = "DRUG_NAME")

stat.test$label = NA
                                                                                                                                                    stat.test$DRUG_NAME == "Paclitaxel"|stat.test$DRUG_NAME == "Oxaliplatin"|stat.test$DRUG_NAME == "SN-38"|stat.test$DRUG_NAME == "5-Fluorouracil"]

potential.targets = unique(unlist(strsplit(ds.gdsc1$PUTATIVE_TARGET,", ")))
potential.targets = potential.targets[potential.targets != "EGFR" & potential.targets != "ERBB2" & potential.targets != "ERBB3" & potential.targets != "ERBB4" ]

keep.index = !grepl(paste(potential.targets, collapse='|'), ds.gdsc1$PUTATIVE_TARGET)
keep.drugs = ds.gdsc1[keep.index,]
keep.drugs = keep.drugs %>%
  filter(PATHWAY_NAME == "RTK signaling" | PATHWAY_NAME == "EGFR signaling")
keep.drugs = unique(keep.drugs$DRUG_NAME)

stat.test$label[stat.test$DRUG_NAME %in% keep.drugs] = stat.test$DRUG_NAME[stat.test$DRUG_NAME %in% keep.drugs]
stat.test$highlight = !is.na(stat.test$label)

drug.family = data.frame(DRUG_NAME = ds.gdsc1$DRUG_NAME, PATHWAY_NAME = ds.gdsc1$PATHWAY_NAME)
drug.family = drug.family %>% distinct()
stat.test = merge(x=stat.test, y= drug.family, by.x = "DRUG_NAME", by.y = "DRUG_NAME", all.x = T)
stat.test$PATHWAY_NAME[stat.test$highlight] = "EGFR/ERBB family"
stat.test$sign = stat.test$p < 0.05
#write.table(x= stat.test, "/Users/markdings/Documents/PhD/Scripts/PanGI_Network_analysis/gdsc analysis/stat.test.EAC.gdsc1.pathway.txt")
p1 = ggplot(stat.test, aes(x=delta, y=-log10(p),label = label))+ geom_text_repel() + ggtitle("GDSC1") + 
  geom_point(color = ifelse(stat.test$highlight , "red", "grey50"), alpha=0.95)+theme_classic()

ggplot(stat.test, aes(x=highlight, y=delta)) + geom_boxplot()
beeswarm1 = ggplot(stat.test, aes(y=PATHWAY_NAME, x=delta)) + geom_beeswarm(aes(color=sign)) +stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                                                                                                           geom = "crossbar", color = "red", width = 1) + scale_color_manual(
                                                                                                             breaks = c("FALSE", "TRUE"),
                                                                                                             values=c("black", "black")
                                                                                                           ) + theme_classic() + geom_vline(xintercept = 0)
#pdf("/Users/markdings/Documents/PhD/Scripts/PanGI_Network_analysis/gdsc analysis/stat.test.EAC.gdsc1.pathway.pdf", width = 5, height = 3)
print(beeswarm1)
#dev.off()


