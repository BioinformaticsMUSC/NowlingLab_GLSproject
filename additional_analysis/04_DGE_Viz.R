# Load libraries
suppressPackageStartupMessages({
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(tidyverse)
library(preprocessCore)
library(future.apply)
library(DESeq2)
library(pheatmap)
library(sva)
library(viridis)
library(limma)
library(janitor)
library(UpSetR)
library(ComplexHeatmap)
})

dir.create("dge_viz")

####################
# IL PLOTS         #
####################
tmp1<- readxl::read_xlsx("input/IL24_CTRL24_Sign_DGE_Filtered_padj05_L2FC_.xlsx")
tmp2<- readxl::read_xlsx("input/IL48_CTRL48_Sign_DGE_Filtered_padj05_L2FC_.xlsx")
tmp3<- readxl::read_xlsx("input/IL72_CTRL72_Sign_DGE_Filtered_padj05_L2FC_.xlsx")


tmp1$Class <- "IL24_CTRL24"
tmp2$Class <- "IL48_CTRL48"
tmp3$Class <- "IL72_CTRL72"

tmp_1 <- tmp1 %>%
      select(Gene,Class)

tmp_2 <- tmp2 %>%
      select(Gene,Class)

tmp_3 <- tmp3 %>%
      select(Gene,Class)


tmp_upset <- rbind(tmp_1,tmp_2,tmp_3)

l <- split(as.character(tmp_upset$Gene),tmp_upset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))

pdf("dge_viz/Upset_Plot_Intersection_ILvsCTL.pdf", width = 6, height = 4)
upset(fromList(l),,nsets = 3, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(IL24_CTRL24 = "#89E651", IL48_CTRL48 = "#DB8BE4",IL72_CTRL72 = "#DDA187"), 
    alpha = 0.5))))
dev.off()




# ComplexHeatmap Fold Changes

dge <- rbind(tmp1,tmp2,tmp3)

tab<-table(dge$Class, dge$Direction)

mat<-matrix(nrow=length(unique(dge$Gene)),ncol=3)
rownames(mat)<-unique(dge$Gene)
colnames(mat)<-unique(dge$Class)

for (i in 1:nrow(mat)){
   for (j in 1:ncol(mat)){
     gene_tmp<-dge[which(dge$Gene==rownames(mat)[i]),]
     gene_tmp<-gene_tmp[which(gene_tmp$Class==colnames(mat)[j]),]
     mat[i,j]<-ifelse(nrow(gene_tmp)>0, gene_tmp$logFC,0)
   }
 }

noDup<-tmp1[!duplicated(tmp1$Gene),]
mat<-mat[order(noDup$Class, noDup$logFC),]
mat[ , colSums(is.na(mat))==0]

tab <- tab[rownames(tab) %in% colnames(mat),]

 # Heatmap
 ha2<-HeatmapAnnotation(Class=colnames(mat), 
      col= list(Class=c(
                           "IL24_CTRL24"="#89E651",
                           "IL48_CTRL48"="#DB8BE4",
                           "IL72_CTRL72"="#DDD3DD"
                           )), 
      show_legend=F)
 
 tab <- tab[rownames(tab) %in% colnames(mat),]
 tab <- tab[match(colnames(mat), rownames(tab)),]

# Top annotaiton
 bar1<-HeatmapAnnotation(up=anno_barplot(tab[,2], gp=gpar(fill="red"), axis_param=list(labels=c("","","","")), ylim=c(0,1000)),
                         down=anno_barplot((tab[,1] * -1), gp=gpar(fill="blue"), axis_param=list(labels=c("","","")), ylim=c(-1000,0),
                         show_annotation_name=F,
                         show_legend=F))
 ha<-c( bar1, ha2)


pdf("dge_viz/DGE_heatmap_foldCs_ILvsCTL.pdf", width=2.5, height=6)
Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=circlize::colorRamp2(c(-2,0,2),c("navy", "white", "firebrick3")), 
        show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha)
dev.off()



direction_colors <- c("UpReg" = "orange", "DownReg" = "purple")

 
df_wide <- dge %>%
  select(Gene, logFC, Class) %>%
  pivot_wider(names_from = Class, values_from = logFC)


a <- ggplot(df_wide, aes(x = IL24_CTRL24, y = IL48_CTRL48)) +
  geom_point(alpha = 0.7, color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  stat_cor(method = "pearson", label.x = min(df_wide$IL24_CTRL24, na.rm = TRUE), 
           label.y = max(df_wide$IL48_CTRL48, na.rm = TRUE), size = 5) +
  labs(
    x = "IL-1β vs CTRL (24h)",
    y = "IL-1β vs CTRL (48h)"
  )


b <- ggplot(df_wide, aes(x = IL24_CTRL24, y = IL72_CTRL72)) +
  geom_point(alpha = 0.7, color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  stat_cor(method = "pearson", label.x = min(df_wide$IL24_CTRL24, na.rm = TRUE), 
           label.y = max(df_wide$IL72_CTRL72, na.rm = TRUE), size = 5) +
  labs(
    x = "IL-1β vs CTRL (24h)",
    y = "IL-1β vs CTRL (72h)"
  )


c <- ggplot(df_wide, aes(x = IL48_CTRL48, y = IL72_CTRL72)) +
  geom_point(alpha = 0.7, color = "darkblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  stat_cor(method = "pearson", label.x = min(df_wide$IL48_CTRL48, na.rm = TRUE), 
           label.y = max(df_wide$IL72_CTRL72, na.rm = TRUE), size = 5) +
  labs(
    x = "IL-1β vs CTRL (48h)",
    y = "IL-1β vs CTRL (72h)"
  )


multi_plot <- cowplot::plot_grid(
  a, b, c,
  labels = c("A", "B", "C"),
  ncol = 3,  # Use nrow = 1 if you want a horizontal layout
  align = 'hv'
)

ggsave("dge_viz/logFC_comparison_ILvsCTL.pdf", multi_plot, width = 15, height = 5)



####################
# IL + Elig PLOTS  #
####################
tmp1<- readxl::read_xlsx("input/ILelig24_CTRL24_Sign_DGE_Filtered_padj05_L2FC_.xlsx")
tmp2<- readxl::read_xlsx("input/ILelig48_CTRL48_Sign_DGE_Filtered_padj05_L2FC_.xlsx")
tmp3<- readxl::read_xlsx("input/ILelig72_CTRL72_Sign_DGE_Filtered_padj05_L2FC_.xlsx")


tmp1$Class <- "ILelig24_CTRL24"
tmp2$Class <- "ILelig48_CTRL48"
tmp3$Class <- "ILelig72_CTRL72"

tmp_1 <- tmp1 %>%
      select(Gene,Class)

tmp_2 <- tmp2 %>%
      select(Gene,Class)

tmp_3 <- tmp3 %>%
      select(Gene,Class)


tmp_upset <- rbind(tmp_1,tmp_2,tmp_3)

l <- split(as.character(tmp_upset$Gene),tmp_upset$Class)
Class <- names(l)
ToTGene <- as.numeric(sapply(l, length))
metadata <- as.data.frame(cbind(Class, ToTGene))
names(metadata) <- c("Class", "ToTGene")
metadata$ToTGene <- as.numeric(as.character(metadata$ToTGene))

pdf("dge_viz/Upset_Plot_Intersection_ILeligvsCTL.pdf", width = 6, height = 4)
upset(fromList(l),,nsets = 3, set.metadata = list(data = metadata, plots = list(list(type = "hist", 
    column = "ToTGene", assign = 20), 
    list(type = "matrix_rows", column = "sets", colors = c(ILelig24_CTRL24 = "#89E651", ILelig48_CTRL48 = "#DB8BE4",ILelig72_CTRL72 = "#DDA187"), 
    alpha = 0.5))))
dev.off()




# ComplexHeatmap Fold Changes

dge <- rbind(tmp1,tmp2,tmp3)

tab<-table(dge$Class, dge$Direction)

mat<-matrix(nrow=length(unique(dge$Gene)),ncol=3)
rownames(mat)<-unique(dge$Gene)
colnames(mat)<-unique(dge$Class)

for (i in 1:nrow(mat)){
   for (j in 1:ncol(mat)){
     gene_tmp<-dge[which(dge$Gene==rownames(mat)[i]),]
     gene_tmp<-gene_tmp[which(gene_tmp$Class==colnames(mat)[j]),]
     mat[i,j]<-ifelse(nrow(gene_tmp)>0, gene_tmp$logFC,0)
   }
 }

noDup<-tmp1[!duplicated(tmp1$Gene),]
mat<-mat[order(noDup$Class, noDup$logFC),]
mat[ , colSums(is.na(mat))==0]

tab <- tab[rownames(tab) %in% colnames(mat),]

 # Heatmap
 ha2<-HeatmapAnnotation(Class=colnames(mat), 
      col= list(Class=c(
                           "ILelig24_CTRL24"="#89E651",
                           "ILelig48_CTRL48"="#DB8BE4",
                           "ILelig72_CTRL72"="#DDD3DD"
                           )), 
      show_legend=F)
 
 tab <- tab[rownames(tab) %in% colnames(mat),]
 tab <- tab[match(colnames(mat), rownames(tab)),]

# Top annotaiton
 bar1<-HeatmapAnnotation(up=anno_barplot(tab[,2], gp=gpar(fill="red"), axis_param=list(labels=c("","","","")), ylim=c(0,1000)),
                         down=anno_barplot((tab[,1] * -1), gp=gpar(fill="blue"), axis_param=list(labels=c("","","")), ylim=c(-1000,0),
                         show_annotation_name=F,
                         show_legend=F))
 ha<-c( bar1, ha2)


pdf("dge_viz/DGE_heatmap_foldCs_ILeligvsCTL.pdf", width=2.5, height=6)
Heatmap(mat, cluster_rows=F,show_row_dend=F,cluster_columns=F,col=circlize::colorRamp2(c(1,0,-1),c("navy", "white", "firebrick3")), 
        show_column_names=T, show_row_names=F, column_title=NULL,row_names_gp = gpar(fontsize = 3), top_annotation=ha)
dev.off()



# Compare FC

direction_colors <- c("UpReg" = "orange", "DownReg" = "purple")

 
df_wide <- dge %>%
  select(Gene, logFC, Class) %>%
  pivot_wider(names_from = Class, values_from = logFC)


a <- ggplot(df_wide, aes(x = ILelig24_CTRL24, y = ILelig48_CTRL48)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  stat_cor(method = "pearson", label.x = min(df_wide$ILelig24_CTRL24, na.rm = TRUE), 
           label.y = max(df_wide$ILelig48_CTRL48, na.rm = TRUE), size = 5) +
  labs(
    x = "IL-1β + Eliglustat vs CTRL (24h)",
    y = "IL-1β + Eliglustat vs CTRL (48h)"
  )


b <- ggplot(df_wide, aes(x = ILelig24_CTRL24, y = ILelig72_CTRL72)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  stat_cor(method = "pearson", label.x = min(df_wide$ILelig24_CTRL24, na.rm = TRUE), 
           label.y = max(df_wide$ILelig72_CTRL72, na.rm = TRUE), size = 5) +
  labs(
    x = "IL-1β + Eliglustat vs CTRL (24h)",
    y = "IL-1β + Eliglustat vs CTRL (72h)"
  )


c <- ggplot(df_wide, aes(x = ILelig48_CTRL48, y = ILelig72_CTRL72)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  stat_cor(method = "pearson", label.x = min(df_wide$ILelig24_CTRL24, na.rm = TRUE), 
           label.y = max(df_wide$ILelig72_CTRL72, na.rm = TRUE), size = 5) +
  labs(
    x = "IL-1β + Eliglustat vs CTRL (48h)",
    y = "IL-1β + Eliglustat vs CTRL (72h)"
  )


multi_plot <- plot_grid(
  a, b, c,
  labels = c("A", "B", "C"),
  ncol = 3,  # Use nrow = 1 if you want a horizontal layout
  align = 'hv'
)

ggsave("dge_viz/logFC_comparison_ILeligvsCTL.pdf", multi_plot, width = 15, height = 5)


