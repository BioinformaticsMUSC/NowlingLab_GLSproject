# tammy nowling

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
library(emmeans)
library(broom)
library(janitor)
library(variancePartition)
library(parallel)
library(doParallel)
library(biomaRt)
library(readxl)
library(dplyr)
library(httr2)
library(scToppR)
library(fgsea)
})

setwd("~/biocm/projects/bulks/nowling_4-25/")

source("./scripts/dge/utils/Utils.R")

plan(multisession, workers = 5)

exp <- read.csv("./outs/counts/gene_count.csv")

exp <- exp %>%
        filter(gene_biotype == "protein_coding") %>%
        dplyr::select(-gene_id,-gene_chr,-gene_start,-gene_end,-gene_strand,-gene_length,-gene_biotype,-gene_description,-tf_family) %>%
        group_by(gene_name) %>%
        summarise_each(list(sum)) %>%
        arrange(gene_name) %>%
        column_to_rownames("gene_name")

meta <- read.delim("./outs/counts/meta_il.txt")
rownames(meta) <- meta$Sample_ID
meta$Sample_ID <- NULL

exp <- exp[, colnames(exp) %in% rownames(meta)]
cpm <- future_apply(exp, 2, function(x) x/sum(as.numeric(x)) * 10^6)


meta$Treatment <- as.factor(meta$Treatment)
meta$Time <- as.factor(meta$Time)
meta$Sex <- as.factor(meta$Sex)
meta$Age <- as.factor(meta$Age)

Untreated <- c(); IL1b <- c(); IL1b_Elig <- c()
for(row in c(1:nrow(meta))) {
    if(meta$Treatment[row] == "Untreated") {
        Untreated <- append(Untreated, row)
    }
    else if(meta$Treatment[row] == "IL1b") {
    	IL1b <- append(IL1b, row)
    }
    else {
        IL1b_Elig<-append(IL1b_Elig, row)
    }
}

filter=apply(cpm, 1, function(x) ( all(x[Untreated] >= 0.5) | all(x[IL1b] >= 0.5) | all(x[IL1b_Elig] >= 0.5) ))
count_filt <- exp[filter,]
cpm_filt <- cpm[filter,]

p <- log2(cpm_filt+1)

save(exp, meta, cpm, cpm_filt, count_filt, p, file="./outs/dge/initial_data_prep.RData")

pdf("./outs/dge/PCA.pdf",width=6,height=6,useDingbats=FALSE)
pca.Sample<-prcomp(t(p))
PCi<-data.frame(pca.Sample$x, Treatment=meta$Treatment, Time=meta$Time, ID = rownames(meta))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
           color = "Time",palette=c("red","black","purple"), 
           shape = "Treatment", size = 3, label = "ID")+
     xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
     ylab(paste("PC2 (",round(variance[2],1),"% )"))+
     theme_classic()
dev.off()


# Variance partitioning 
pdf("./outs/dge/Variance_Explained.pdf", width=5, height=4)
var <- VarExp(p,meta,10,FALSE)
plotVarExp(var,"Variance Explained")
dev.off()


cl <- makeCluster(10)
registerDoParallel(cl)
formula <- ~ (1|Treatment) + (1|Time) + (1|Sex) + (1|Age)
varPart <- fitExtractVarPartModel(p, formula, meta)
varPart_ordered <- varPart[order(varPart$Treatment, decreasing=TRUE),]

pdf("./outs/dge/Variance_Explained_ByGene.pdf", width=4, height=6)
plotPercentBars(varPart_ordered[1:50,]) + 
theme_classic()
dev.off()



meta_sva <- meta %>%
            dplyr::select(Sex, Age) %>% #Removing 
            droplevels()

betas<-future_lapply(1:nrow(p), function(x)
            {
                lm(unlist(p[x,])~., data = meta_sva)
            })

residuals<-future_lapply(betas, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
p_regressed <- residuals+matrix(future_apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(p_regressed)<-rownames(p)

write.table(p_regressed,"./outs/dge/expression_regressed.txt",sep="\t",quote=F)



pdf("./outs/dge/PCA_AdjustedForConfound.pdf",width=8,height=8,useDingbats=FALSE)
pca.Sample<-prcomp(t(p_regressed))
PCi<-data.frame(pca.Sample$x, Treatment=meta$Treatment, Time=meta$Time, ID = rownames(meta))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Time",palette=c("red","black","purple"), 
          shape = "Treatment", size = 3,label = "ID")+
xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
ylab(paste("PC2 (",round(variance[2],1),"% )"))+
theme_classic()
dev.off()



model1 <- 'geneExpr ~ Treatment * Time'

# This is a function to get the P-value out of this post-hoc analysis.
# This is based on a linear model
emmeans_lm_p <- function(vectorizedExpression) {
  tmpMetaData <- cbind(meta, data.frame(geneExpr = unname(vectorizedExpression)))
  residuals1 <- lm(model1,data=tmpMetaData)
  pval <- emmeans::emmeans(residuals1,pairwise ~ Treatment * Time, adjust="tukey")$contrasts %>%
          broom::tidy() %>%
          dplyr::select(contrast,adj.p.value) %>% 
          pivot_wider(names_from = contrast, values_from = adj.p.value) %>%
          as.data.frame()
       
}

emmeans_lm_p_fun <- function(vectorizedExpression) {
  tryCatch(emmeans_lm_p(vectorizedExpression))
}


posthoc_pval <- future_apply(p_regressed, 1, emmeans_lm_p_fun) %>% 
                      bind_rows() %>%
                      clean_names() %>%
                      as.data.frame() 

rownames(posthoc_pval) <- rownames(p_regressed)

# Same function but capturing the effect sizes
emmeans_lm_eff <- function(vectorizedExpression) {
  tmpMetaData <- cbind(meta, data.frame(geneExpr = unname(vectorizedExpression)))
  residuals1 <- lm(model1,data=tmpMetaData)
  effect_size <- emmeans::emmeans(residuals1, adjust="tukey",pairwise ~ Treatment * Time)$contrasts %>%
          broom::tidy() %>%
          dplyr::select(contrast,estimate) %>%  
          pivot_wider(names_from = contrast, values_from = estimate) %>%
          as.data.frame()
       
}

emmeans_lm_eff_fun <- function(vectorizedExpression) {
  tryCatch(emmeans_lm_eff(vectorizedExpression))
}


posthoc_eff <- future_apply(p_regressed, 1, emmeans_lm_eff_fun) %>% 
                      bind_rows() %>%
                      clean_names() %>%
                      as.data.frame() 

rownames(posthoc_eff) <- rownames(p_regressed)


write.table(posthoc_pval,"./outs/dge/nowling_Dge_Interaction_Posthoc_Pvalues.txt",sep="\t",quote=F)
write.table(posthoc_eff,"./outs/dge/nowling_Dge_Interaction_Posthoc_EffSize.txt",sep="\t",quote=F)


sign_posthoc_pval <- posthoc_pval %>%
                   filter(if_any(is.numeric, ~ .x < 0.05)) %>%
                   rownames_to_column("Gene")

sign_posthoc_eff <- posthoc_eff[rownames(posthoc_eff) %in% sign_posthoc_pval$Gene,] %>%
                   rownames_to_column("Gene")


openxlsx::write.xlsx(sign_posthoc_pval, 
                     file = "./outs/dge/nowling_Dge_Interaction_Posthoc_Pvalues_Sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")


openxlsx::write.xlsx(sign_posthoc_eff, 
                     file = "./outs/dge/nowling_Dge_Interaction_Posthoc_EffSize_Sign.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

save(posthoc_eff,posthoc_pval,sign_posthoc_pval,sign_posthoc_eff,  file = "./outs/dge/nowling_Dge_Interaction_Posthoc.RData")


save(p_regressed, file = "./outs/dge/p_regressed.RData")

dir.create("./outs/dge/IL24_CTRL24")
dir.create("./outs/dge/IL48_CTRL48")
dir.create("./outs/dge/IL72_CTRL72")
dir.create("./outs/dge/ILelig24_CTRL24")
dir.create("./outs/dge/ILelig48_CTRL48")
dir.create("./outs/dge/ILelig72_CTRL72")

IL24_CTRL24 <- data.frame(Gene = rownames(posthoc_pval), logFC = posthoc_eff$il1b_time24_untreated_time24, FDR = posthoc_pval$il1b_time24_untreated_time24)
IL48_CTRL48 <- data.frame(Gene = rownames(posthoc_pval), logFC = posthoc_eff$il1b_time48_untreated_time48, FDR = posthoc_pval$il1b_time48_untreated_time48)
IL72_CTRL72 <- data.frame(Gene = rownames(posthoc_pval), logFC = posthoc_eff$il1b_time72_untreated_time72, FDR = posthoc_pval$il1b_time72_untreated_time72)
ILelig24_CTRL24 <- data.frame(Gene = rownames(posthoc_pval), logFC = posthoc_eff$il1b_elig_time24_untreated_time24, FDR = posthoc_pval$il1b_elig_time24_untreated_time24)
ILelig48_CTRL48 <- data.frame(Gene = rownames(posthoc_pval), logFC = posthoc_eff$il1b_elig_time48_untreated_time48, FDR = posthoc_pval$il1b_elig_time48_untreated_time48)
ILelig72_CTRL72 <- data.frame(Gene = rownames(posthoc_pval), logFC = posthoc_eff$il1b_elig_time72_untreated_time72, FDR = posthoc_pval$il1b_elig_time72_untreated_time72)

IL24_CTRL24_Sign <- IL24_CTRL24 %>% filter(FDR < 0.05, abs(logFC) > 0.3)
IL48_CTRL48_Sign <- IL48_CTRL48 %>% filter(FDR < 0.05, abs(logFC) > 0.3)
IL72_CTRL72_Sign <- IL72_CTRL72 %>% filter(FDR < 0.05, abs(logFC) > 0.3)
ILelig24_CTRL24_Sign <- ILelig24_CTRL24 %>% filter(FDR < 0.05, abs(logFC) > 0.3)
ILelig48_CTRL48_Sign <- ILelig48_CTRL48 %>% filter(FDR < 0.05, abs(logFC) > 0.3)
ILelig72_CTRL72_Sign <- ILelig72_CTRL72 %>% filter(FDR < 0.05, abs(logFC) > 0.3)



save(IL24_CTRL24, IL48_CTRL48, IL72_CTRL72,
     IL24_CTRL24_Sign, IL48_CTRL48_Sign, IL72_CTRL72_Sign,
     ILelig24_CTRL24, ILelig48_CTRL48, ILelig72_CTRL72,
     ILelig24_CTRL24_Sign, ILelig48_CTRL48_Sign, ILelig72_CTRL72_Sign, 
     file = "./outs/dge/nowling_Dge_Defined.RData")

### box plots ###

# IL24
IL24_CTRL24_Sign_df <- IL24_CTRL24_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")

openxlsx::write.xlsx(IL24_CTRL24_Sign_df,
                     file = "./outs/dge/IL24_CTRL24/IL24_CTRL24_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL24_CTRL24_Sign_df_top_labelled <- IL24_CTRL24_Sign_df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(FDR) %>%
  top_n(n = 10, wt = LOG)

IL24_CTRL24_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 24, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 24, ]$Treatment)


# IL48
IL48_CTRL48_Sign_df <- IL48_CTRL48_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")

openxlsx::write.xlsx(IL48_CTRL48_Sign_df,
                     file = "./outs/dge/IL48_CTRL48/IL48_CTRL48_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL48_CTRL48_Sign_df_top_labelled <- IL48_CTRL48_Sign_df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(FDR) %>%
  top_n(n = 10, wt = LOG)

IL48_CTRL48_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 48, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 48, ]$Treatment)


# IL72
IL72_CTRL72_Sign_df <- IL72_CTRL72_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")

openxlsx::write.xlsx(IL72_CTRL72_Sign_df,
                     file = "./outs/dge/IL72_CTRL72/IL72_CTRL72_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL72_CTRL72_Sign_df_top_labelled <- IL72_CTRL72_Sign_df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(FDR) %>%
  top_n(n = 10, wt = LOG)

IL72_CTRL72_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 72, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 72, ]$Treatment)


# ILelig24
ILelig24_CTRL24_Sign_df <- ILelig24_CTRL24_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")

openxlsx::write.xlsx(ILelig24_CTRL24_Sign_df,
                     file = "./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

ILelig24_CTRL24_Sign_df_top_labelled <- ILelig24_CTRL24_Sign_df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(FDR) %>%
  top_n(n = 10, wt = LOG)

ILelig24_CTRL24_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 24, ]), Treatment=meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 24, ]$Treatment)



# ILelig48
ILelig48_CTRL48_Sign_df <- ILelig48_CTRL48_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")

openxlsx::write.xlsx(ILelig48_CTRL48_Sign_df,
                     file = "./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

ILelig48_CTRL48_Sign_df_top_labelled <- ILelig48_CTRL48_Sign_df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(FDR) %>%
  top_n(n = 10, wt = LOG)

ILelig48_CTRL48_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 48, ]), Treatment=meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 48, ]$Treatment)


# ILelig72
ILelig72_CTRL72_Sign_df <- ILelig72_CTRL72_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")

openxlsx::write.xlsx(ILelig72_CTRL72_Sign_df,
                     file = "./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

ILelig72_CTRL72_Sign_df_top_labelled <- ILelig72_CTRL72_Sign_df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(FDR) %>%
  top_n(n = 10, wt = LOG)

ILelig72_CTRL72_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 72, ]), Treatment=meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 72, ]$Treatment)


#  boxplots

# IL24
mat <- p_regressed[rownames(p_regressed) %in% IL24_CTRL24_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL24_CTRL24_coldata$sample] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = IL24_CTRL24_coldata$Treatment) %>%
  pivot_longer(!Treatment, names_to = "Gene", values_to="Exp")


pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_Sign_03_Boxplots_TopGenes.pdf",width=8,height=5,useDingbats=FALSE)
ggboxplot(mat, "Treatment", "Exp", fill = "Treatment",
          palette = c("blue", "red")) +
  xlab("") +
  ylab("log2(Expression Adjusted)") +
  theme_classic() +
  facet_wrap(.~Gene,scales="free",ncol=5,nrow=5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


# IL48
mat <- p_regressed[rownames(p_regressed) %in% IL48_CTRL48_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL48_CTRL48_coldata$sample] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = IL48_CTRL48_coldata$Treatment) %>%
  pivot_longer(!Treatment, names_to = "Gene", values_to="Exp")


pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_Sign_03_Boxplots_TopGenes.pdf",width=8,height=5,useDingbats=FALSE)
ggboxplot(mat, "Treatment", "Exp", fill = "Treatment",
          palette = c("blue", "red")) +
  xlab("") +
  ylab("log2(Expression Adjusted)") +
  theme_classic() +
  facet_wrap(.~Gene,scales="free",ncol=5,nrow=5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


# IL72
mat <- p_regressed[rownames(p_regressed) %in% IL72_CTRL72_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL72_CTRL72_coldata$sample] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = IL72_CTRL72_coldata$Treatment) %>%
  pivot_longer(!Treatment, names_to = "Gene", values_to="Exp")


pdf("./outs/dge/IL72_CTRL72/IL72_CTRL72_Sign_03_Boxplots_TopGenes.pdf",width=8,height=5,useDingbats=FALSE)
ggboxplot(mat, "Treatment", "Exp", fill = "Treatment",
          palette = c("blue", "red")) +
  xlab("") +
  ylab("log2(Expression Adjusted)") +
  theme_classic() +
  facet_wrap(.~Gene,scales="free",ncol=5,nrow=5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


# ILelig24
mat <- p_regressed[rownames(p_regressed) %in% ILelig24_CTRL24_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% ILelig24_CTRL24_coldata$sample] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = ILelig24_CTRL24_coldata$Treatment) %>%
  pivot_longer(!Treatment, names_to = "Gene", values_to="Exp")


pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_Sign_03_Boxplots_TopGenes.pdf",width=8,height=5,useDingbats=FALSE)
ggboxplot(mat, "Treatment", "Exp", fill = "Treatment",
          palette = c("blue", "red")) +
  xlab("") +
  ylab("log2(Expression Adjusted)") +
  theme_classic() +
  facet_wrap(.~Gene,scales="free",ncol=5,nrow=5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


# ILelig48
mat <- p_regressed[rownames(p_regressed) %in% ILelig48_CTRL48_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% ILelig48_CTRL48_coldata$sample] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = ILelig48_CTRL48_coldata$Treatment) %>%
  pivot_longer(!Treatment, names_to = "Gene", values_to="Exp")


pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_Sign_03_Boxplots_TopGenes.pdf",width=8,height=5,useDingbats=FALSE)
ggboxplot(mat, "Treatment", "Exp", fill = "Treatment",
          palette = c("blue", "red")) +
  xlab("") +
  ylab("log2(Expression Adjusted)") +
  theme_classic() +
  facet_wrap(.~Gene,scales="free",ncol=5,nrow=5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


# ILelig72
mat <- p_regressed[rownames(p_regressed) %in% ILelig72_CTRL72_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% ILelig72_CTRL72_coldata$sample] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = ILelig72_CTRL72_coldata$Treatment) %>%
  pivot_longer(!Treatment, names_to = "Gene", values_to="Exp")


pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_Sign_03_Boxplots_TopGenes.pdf",width=8,height=5,useDingbats=FALSE)
ggboxplot(mat, "Treatment", "Exp", fill = "Treatment",
          palette = c("blue", "red")) +
  xlab("") +
  ylab("log2(Expression Adjusted)") +
  theme_classic() +
  facet_wrap(.~Gene,scales="free",ncol=5,nrow=5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()



### volc plots ###

# IL24
vol <- IL24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_CTRL24_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_no_text_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  #geom_text_repel(data = IL24_CTRL24_Sign_df_top_labelled,
  #                mapping = aes(label = Gene),
  #                size = 5,
  #                box.padding = unit(0.4, "lines"),
  #                point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

vol_sig <- IL24_CTRL24_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_Sign_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol_sig,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_CTRL24_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# IL48
vol <- IL48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_CTRL48_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_no_text_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  #geom_text_repel(data = IL48_CTRL48_Sign_df_top_labelled,
  #                mapping = aes(label = Gene),
  #                size = 5,
  #                box.padding = unit(0.4, "lines"),
  #                point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

vol_sig <- IL48_CTRL48_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_Sign_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol_sig,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_CTRL48_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# ILelig24
vol <- ILelig24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig24_CTRL24_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_no_text_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  #geom_text_repel(data = ILelig24_CTRL24_Sign_df_top_labelled,
  #                mapping = aes(label = Gene),
  #                size = 5,
  #                box.padding = unit(0.4, "lines"),
  #                point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

vol_sig <- ILelig24_CTRL24_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_Sign_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol_sig,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig24_CTRL24_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# ILelig48
vol <- ILelig48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig48_CTRL48_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_no_text_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  #geom_text_repel(data = ILelig48_CTRL48_Sign_df_top_labelled,
  #                mapping = aes(label = Gene),
  #                size = 5,
  #                box.padding = unit(0.4, "lines"),
  #                point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

vol_sig <- ILelig48_CTRL48_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_Sign_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol_sig,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig48_CTRL48_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# ILelig72
vol <- ILelig72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig72_CTRL72_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_no_text_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  #geom_text_repel(data = ILelig72_CTRL72_Sign_df_top_labelled,
  #                mapping = aes(label = Gene),
  #                size = 5,
  #                box.padding = unit(0.4, "lines"),
  #                point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

vol_sig <- ILelig72_CTRL72_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_Sign_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol_sig,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig72_CTRL72_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


###  heat maps  ###

# IL24
allmat <- p_regressed[rownames(p_regressed) %in% IL24_CTRL24_Sign$Gene, colnames(p_regressed) %in% IL24_CTRL24_coldata$sample]
hmat <- p_regressed[rownames(p_regressed) %in% IL24_CTRL24_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL24_CTRL24_coldata$sample]

anno <- IL24_CTRL24_coldata
Treatment <- c("#FFA500", "#BF40BF")
names(Treatment) <- c("IL1b", "Untreated")
anno_colors <- list(Treatment = Treatment)
rownames(anno) <- IL24_CTRL24_coldata$sample

pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_HeatPlot_Top20.pdf",width=8,height=10)
pheatmap(hmat,scale="row",show_rownames = T,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()

pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_HeatPlot_All.pdf",width=8,height=10)
pheatmap(allmat,scale="row",show_rownames = F,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()


#IL48
allmat <- p_regressed[rownames(p_regressed) %in% IL48_CTRL48_Sign$Gene, colnames(p_regressed) %in% IL48_CTRL48_coldata$sample]
hmat <- p_regressed[rownames(p_regressed) %in% IL48_CTRL48_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL48_CTRL48_coldata$sample]

anno <- IL48_CTRL48_coldata
Treatment <- c("#FFA500", "#BF40BF")
names(Treatment) <- c("IL1b", "Untreated")
anno_colors <- list(Treatment = Treatment)
rownames(anno) <- IL48_CTRL48_coldata$sample

pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_HeatPlot_Top20.pdf",width=8,height=10)
pheatmap(hmat,scale="row",show_rownames = T,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()

pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_HeatPlot_All.pdf",width=8,height=10)
pheatmap(allmat,scale="row",show_rownames = F,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()


#IL72
allmat <- p_regressed[rownames(p_regressed) %in% IL72_CTRL72_Sign$Gene, colnames(p_regressed) %in% IL72_CTRL72_coldata$sample]
hmat <- p_regressed[rownames(p_regressed) %in% IL72_CTRL72_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL72_CTRL72_coldata$sample]

anno <- IL72_CTRL72_coldata
Treatment <- c("#FFA500", "#BF40BF")
names(Treatment) <- c("IL1b", "Untreated")
anno_colors <- list(Treatment = Treatment)
rownames(anno) <- IL72_CTRL72_coldata$sample

pdf("./outs/dge/IL72_CTRL72/IL72_CTRL72_HeatPlot_Top20.pdf",width=8,height=10)
pheatmap(hmat,scale="row",show_rownames = T,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()

pdf("./outs/dge/IL72_CTRL72/IL72_CTRL72_HeatPlot_All.pdf",width=8,height=10)
pheatmap(allmat,scale="row",show_rownames = F,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()


# ILelig24
allmat <- p_regressed[rownames(p_regressed) %in% ILelig24_CTRL24_Sign$Gene, colnames(p_regressed) %in% ILelig24_CTRL24_coldata$sample]
hmat <- p_regressed[rownames(p_regressed) %in% ILelig24_CTRL24_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% ILelig24_CTRL24_coldata$sample]

anno <- ILelig24_CTRL24_coldata
Treatment <- c("#FFA500", "#BF40BF")
names(Treatment) <- c("IL1b_Elig", "Untreated")
anno_colors <- list(Treatment = Treatment)
rownames(anno) <- ILelig24_CTRL24_coldata$sample

pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_HeatPlot_Top20.pdf",width=8,height=10)
pheatmap(hmat,scale="row",show_rownames = T,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_HeatPlot_All.pdf",width=8,height=10)
pheatmap(allmat,scale="row",show_rownames = F,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()


# ILelig48
allmat <- p_regressed[rownames(p_regressed) %in% ILelig48_CTRL48_Sign$Gene, colnames(p_regressed) %in% ILelig48_CTRL48_coldata$sample]
hmat <- p_regressed[rownames(p_regressed) %in% ILelig48_CTRL48_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% ILelig48_CTRL48_coldata$sample]

anno <- ILelig48_CTRL48_coldata
Treatment <- c("#FFA500", "#BF40BF")
names(Treatment) <- c("IL1b_Elig", "Untreated")
anno_colors <- list(Treatment = Treatment)
rownames(anno) <- ILelig48_CTRL48_coldata$sample

pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_HeatPlot_Top20.pdf",width=8,height=10)
pheatmap(hmat,scale="row",show_rownames = T,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_HeatPlot_All.pdf",width=8,height=10)
pheatmap(allmat,scale="row",show_rownames = F,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()


# ILelig72
allmat <- p_regressed[rownames(p_regressed) %in% ILelig72_CTRL72_Sign$Gene, colnames(p_regressed) %in% ILelig72_CTRL72_coldata$sample]
hmat <- p_regressed[rownames(p_regressed) %in% ILelig72_CTRL72_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% ILelig72_CTRL72_coldata$sample]

anno <- ILelig72_CTRL72_coldata
Treatment <- c("#FFA500", "#BF40BF")
names(Treatment) <- c("IL1b_Elig", "Untreated")
anno_colors <- list(Treatment = Treatment)
rownames(anno) <- ILelig72_CTRL72_coldata$sample

pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_HeatPlot_Top20.pdf",width=8,height=10)
pheatmap(hmat,scale="row",show_rownames = T,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_HeatPlot_All.pdf",width=8,height=10)
pheatmap(allmat,scale="row",show_rownames = F,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()



###  toppgene  ###

# IL24
IL24_CTRL24_Sign_df$p_val<-sign_posthoc_pval[sign_posthoc_pval$Gene %in% IL24_CTRL24_Sign_df$Gene, ]$il1b_time24_untreated_time24
IL24_CTRL24_Sign_df$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL24_CTRL24_toppgene<-toppFun(IL24_CTRL24_Sign_df, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
IL24_CTRL24_toppgene$Cluster<-'bulk'

save(IL24_CTRL24_toppgene, file="./outs/dge/IL24_CTRL24/IL24_CTRL24_toppgene_queries.RData")

pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(IL24_CTRL24_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(IL24_CTRL24_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(IL24_CTRL24_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(IL24_CTRL24_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# IL48
IL48_CTRL48_Sign_df$p_val<-sign_posthoc_pval[sign_posthoc_pval$Gene %in% IL48_CTRL48_Sign_df$Gene, ]$il1b_time48_untreated_time48
IL48_CTRL48_Sign_df$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL48_CTRL48_toppgene<-toppFun(IL48_CTRL48_Sign_df, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
IL48_CTRL48_toppgene$Cluster<-'bulk'

save(IL48_CTRL48_toppgene, file="./outs/dge/IL48_CTRL48/IL48_CTRL48_toppgene_queries.RData")

pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(IL48_CTRL48_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(IL48_CTRL48_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(IL48_CTRL48_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(IL48_CTRL48_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# IL72
IL72_CTRL72_Sign_df$p_val<-sign_posthoc_pval[sign_posthoc_pval$Gene %in% IL72_CTRL72_Sign_df$Gene, ]$il1b_time72_untreated_time72
IL72_CTRL72_Sign_df$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL72_CTRL72_toppgene<-toppFun(IL72_CTRL72_Sign_df, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
IL72_CTRL72_toppgene$Cluster<-'bulk'

save(IL72_CTRL72_toppgene, file="./outs/dge/IL72_CTRL72/IL72_CTRL72_toppgene_queries.RData")

pdf("./outs/dge/IL72_CTRL72/IL72_CTRL72_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(IL72_CTRL72_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL72_CTRL72/IL72_CTRL72_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(IL72_CTRL72_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL72_CTRL72/IL72_CTRL72_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(IL72_CTRL72_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL72_CTRL72/IL72_CTRL72_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(IL72_CTRL72_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# ILelig24
ILelig24_CTRL24_Sign_df$p_val<-sign_posthoc_pval[sign_posthoc_pval$Gene %in% ILelig24_CTRL24_Sign_df$Gene, ]$il1b_elig_time72_untreated_time72
ILelig24_CTRL24_Sign_df$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
ILelig24_CTRL24_toppgene<-toppFun(ILelig24_CTRL24_Sign_df, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
ILelig24_CTRL24_toppgene$Cluster<-'bulk'

save(ILelig24_CTRL24_toppgene, file="./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_toppgene_queries.RData")

pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig24_CTRL24_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig24_CTRL24_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig24_CTRL24_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(ILelig24_CTRL24_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# ILelig48
ILelig48_CTRL48_Sign_df$p_val<-sign_posthoc_pval[sign_posthoc_pval$Gene %in% ILelig48_CTRL48_Sign_df$Gene, ]$il1b_elig_time72_untreated_time72
ILelig48_CTRL48_Sign_df$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
ILelig48_CTRL48_toppgene<-toppFun(ILelig48_CTRL48_Sign_df, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
ILelig48_CTRL48_toppgene$Cluster<-'bulk'

save(ILelig48_CTRL48_toppgene, file="./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_toppgene_queries.RData")

pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig48_CTRL48_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig48_CTRL48_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig48_CTRL48_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(ILelig48_CTRL48_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# ILelig72
ILelig72_CTRL72_Sign_df$p_val<-sign_posthoc_pval[sign_posthoc_pval$Gene %in% ILelig72_CTRL72_Sign_df$Gene, ]$il1b_elig_time72_untreated_time72
ILelig72_CTRL72_Sign_df$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
ILelig72_CTRL72_toppgene<-toppFun(ILelig72_CTRL72_Sign_df, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
ILelig72_CTRL72_toppgene$Cluster<-'bulk'

save(ILelig72_CTRL72_toppgene, file="./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_toppgene_queries.RData")

pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig72_CTRL72_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig72_CTRL72_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig72_CTRL72_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(ILelig72_CTRL72_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()



# Boxplot all genes
p_regressed_df <- as.data.frame(p_regressed)
mat <-  p_regressed_df[rownames(p_regressed_df)%in% sign_posthoc_pval$Gene,] %>%
        rownames_to_column("Gene") %>%
        reshape2::melt() %>%
        as.data.frame()
tmp <- merge(mat,meta,by.x="variable",by.y="row.names") %>%
        mutate(Gene == as.factor(Gene)) %>%
        mutate(Treatment = fct_relevel(Treatment, c("Untreated","IL1b", "IL1b_Elig")), Time = fct_relevel(Time, c("24","48","72")))
dir.create("./outs/dge/gene_boxplots/")
cl <- colors(distinct = TRUE)
set.seed(15887) # to set random generator seed
cols <- c("#454B87","#BD3106","#FFA500")
doPlot = function(sel_name)
{
    df = subset(tmp, Gene == sel_name)
    PLOT= ggboxplot(df, "Treatment", "value", color = "Time",
              palette = cols,
              outlier.shape = NA) +
      xlab("")+
      ylab("Gene Exp")+
        theme_classic() +
        rotate_x_text(angle = 45)
    print(PLOT)
    ggsave(sprintf("./outs/dge/gene_boxplots/%s.pdf", sel_name),width=4,height=4)
 }
lapply(unique(tmp$Gene), doPlot)



pathways.hallmark <- gmtPathways("./scripts/dge/utils/h.all.v7.5.symbols.gmt")


# IL24

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL24_CTRL24$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL24_CTRL24$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL24_CTRL24, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_all_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=6)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL24_CTRL24_Sign_df$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL24_CTRL24_Sign_df$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL24_CTRL24_Sign_df, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL24_CTRL24/IL24_CTRL24_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()


# IL48

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL48_CTRL48$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL48_CTRL48$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL48_CTRL48, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_all_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=6)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL48_CTRL48_Sign_df$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL48_CTRL48_Sign_df$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL48_CTRL48_Sign_df, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL48_CTRL48/IL48_CTRL48_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()


# IL72

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL72_CTRL72$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL72_CTRL72$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL72_CTRL72, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL72_CTRL72/IL72_CTRL72_all_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=6)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL72_CTRL72_Sign_df$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL72_CTRL72_Sign_df$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL72_CTRL72_Sign_df, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL72_CTRL72/IL72_CTRL72_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()


# ILelig24

# Create a named vector with fold changes as rank and run GSEA
ranks <- ILelig24_CTRL24$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- ILelig24_CTRL24$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(ILelig24_CTRL24, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_all_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=6)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()

# Create a named vector with fold changes as rank and run GSEA
ranks <- ILelig24_CTRL24_Sign_df$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- ILelig24_CTRL24_Sign_df$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(ILelig24_CTRL24_Sign_df, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/ILelig24_CTRL24/ILelig24_CTRL24_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()


# ILelig48

# Create a named vector with fold changes as rank and run GSEA
ranks <- ILelig48_CTRL48$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- ILelig48_CTRL48$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(ILelig48_CTRL48, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_all_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=6)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()


# Create a named vector with fold changes as rank and run GSEA
ranks <- ILelig48_CTRL48_Sign_df$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- ILelig48_CTRL48_Sign_df$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(ILelig48_CTRL48_Sign_df, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/ILelig48_CTRL48/ILelig48_CTRL48_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()


# ILelig72

# Create a named vector with fold changes as rank and run GSEA
ranks <- ILelig72_CTRL72$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- ILelig72_CTRL72$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(ILelig72_CTRL72, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_all_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=6)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()


# Create a named vector with fold changes as rank and run GSEA
ranks <- ILelig72_CTRL72_Sign_df$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- ILelig72_CTRL72_Sign_df$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(ILelig72_CTRL72_Sign_df, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/ILelig72_CTRL72/ILelig72_CTRL72_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()






#################


dir.create("./outs/dge/IL24_ILelig24")
dir.create("./outs/dge/IL48_ILelig48")
dir.create("./outs/dge/IL72_ILelig72")



IL24_elig24 <- data.frame(Gene = rownames(posthoc_pval), logFC = posthoc_eff$il1b_time24_il1b_elig_time24, FDR = posthoc_pval$il1b_time24_il1b_elig_time24)
IL48_elig48 <- data.frame(Gene = rownames(posthoc_pval), logFC = posthoc_eff$il1b_time48_il1b_elig_time48, FDR = posthoc_pval$il1b_time48_il1b_elig_time48)
IL72_elig72 <- data.frame(Gene = rownames(posthoc_pval), logFC = posthoc_eff$il1b_time72_il1b_elig_time72, FDR = posthoc_pval$il1b_time72_il1b_elig_time72)

IL24_elig24_Sign <- IL24_elig24 %>% filter(FDR < 0.05, abs(logFC) > 0.3)
IL48_elig48_Sign <- IL48_elig48 %>% filter(FDR < 0.05, abs(logFC) > 0.3)
IL72_elig72_Sign <- IL72_elig72 %>% filter(FDR < 0.05, abs(logFC) > 0.3)

save(IL24_CTRL24, IL48_CTRL48, IL72_CTRL72,
     IL24_CTRL24_Sign, IL48_CTRL48_Sign, IL72_CTRL72_Sign,
     ILelig24_CTRL24, ILelig48_CTRL48, ILelig72_CTRL72,
     ILelig24_CTRL24_Sign, ILelig48_CTRL48_Sign, ILelig72_CTRL72_Sign,
     IL24_elig24, IL48_elig48, IL72_elig72,
     IL24_elig24_Sign, IL48_elig48_Sign, IL72_elig72_Sign,
     file = "./outs/dge/nowling_Dge_Defined.RData")

### box plots ###

#
IL24_elig24_Sign_df <- IL24_elig24_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")

openxlsx::write.xlsx(IL24_elig24_Sign_df,
                     file = "./outs/dge/IL24_ILelig24/IL24_ILelig24_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL24_elig24_Sign_df_top_labelled <- IL24_elig24_Sign_df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(FDR) %>%
  top_n(n = 10, wt = LOG)

IL24_elig24_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 24, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 24, ]$Treatment)

#
IL48_elig48_Sign_df <- IL48_elig48_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")

openxlsx::write.xlsx(IL48_elig48_Sign_df,
                     file = "./outs/dge/IL48_ILelig48/IL48_ILelig48_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL48_elig48_Sign_df_top_labelled <- IL48_elig48_Sign_df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(FDR) %>%
  top_n(n = 10, wt = LOG)

IL48_elig48_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 48, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 48, ]$Treatment)

#
IL72_elig72_Sign_df <- IL72_elig72_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")

openxlsx::write.xlsx(IL72_elig72_Sign_df,
                     file = "./outs/dge/IL72_ILelig72/IL72_ILelig72_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL72_elig72_Sign_df_top_labelled <- IL72_elig72_Sign_df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(FDR) %>%
  top_n(n = 10, wt = LOG)

IL72_elig72_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 72, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 72, ]$Treatment)

#  boxplots

# 
mat <- p_regressed[rownames(p_regressed) %in% IL24_elig24_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL24_elig24_coldata$sample] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = IL24_elig24_coldata$Treatment) %>%
  pivot_longer(!Treatment, names_to = "Gene", values_to="Exp")


pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_Sign_03_Boxplots_TopGenes.pdf",width=8,height=5,useDingbats=FALSE)
ggboxplot(mat, "Treatment", "Exp", fill = "Treatment",
          palette = c("blue", "red")) +
  xlab("") +
  ylab("log2(Expression Adjusted)") +
  theme_classic() +
  facet_wrap(.~Gene,scales="free",ncol=5,nrow=5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

# 
mat <- p_regressed[rownames(p_regressed) %in% IL48_elig48_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL48_elig48_coldata$sample] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = IL24_elig24_coldata$Treatment) %>%
  pivot_longer(!Treatment, names_to = "Gene", values_to="Exp")


pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_Sign_03_Boxplots_TopGenes.pdf",width=8,height=5,useDingbats=FALSE)
ggboxplot(mat, "Treatment", "Exp", fill = "Treatment",
          palette = c("blue", "red")) +
  xlab("") +
  ylab("log2(Expression Adjusted)") +
  theme_classic() +
  facet_wrap(.~Gene,scales="free",ncol=5,nrow=5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

# 
mat <- p_regressed[rownames(p_regressed) %in% IL72_elig72_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL72_elig72_coldata$sample] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Treatment = IL72_elig72_coldata$Treatment) %>%
  pivot_longer(!Treatment, names_to = "Gene", values_to="Exp")


pdf("./outs/dge/IL72_ILelig72/IL72_Ilelig72_Sign_03_Boxplots_TopGenes.pdf",width=8,height=5,useDingbats=FALSE)
ggboxplot(mat, "Treatment", "Exp", fill = "Treatment",
          palette = c("blue", "red")) +
  xlab("") +
  ylab("log2(Expression Adjusted)") +
  theme_classic() +
  facet_wrap(.~Gene,scales="free",ncol=5,nrow=5) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()


### volc plots ###

# 
vol <- IL24_elig24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_elig24_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_no_text_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  #geom_text_repel(data = IL24_elig24_Sign_df_top_labelled,
  #                mapping = aes(label = Gene),
  #                size = 5,
  #                box.padding = unit(0.4, "lines"),
  #                point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

vol_sig <- IL24_elig24_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_Sign_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol_sig,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_elig24_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# 
vol <- IL48_elig48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_elig48_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_no_text_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  #geom_text_repel(data = IL48_elig48_Sign_df_top_labelled,
  #                mapping = aes(label = Gene),
  #                size = 5,
  #                box.padding = unit(0.4, "lines"),
  #                point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

vol_sig <- IL48_elig48_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_Sign_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol_sig,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_elig48_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# 
vol <- IL72_elig72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL72_elig72_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_no_text_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  #geom_text_repel(data = IL72_elig72_Sign_df_top_labelled,
  #                mapping = aes(label = Gene),
  #                size = 5,
  #                box.padding = unit(0.4, "lines"),
  #                point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

vol_sig <- IL72_elig72_Sign %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg"))

pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_Sign_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol_sig,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL72_elig72_Sign_df_top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


###  heat maps  ###

# IL24
allmat <- p_regressed[rownames(p_regressed) %in% IL24_elig24_Sign$Gene, colnames(p_regressed) %in% IL24_elig24_coldata$sample]
hmat <- p_regressed[rownames(p_regressed) %in% IL24_elig24_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL24_elig24_coldata$sample]

anno <- IL24_elig24_coldata
Treatment <- c("#FFA500", "#BF40BF")
names(Treatment) <- c("IL1b", "IL1b_Elig")
anno_colors <- list(Treatment = Treatment)
rownames(anno) <- IL24_elig24_coldata$sample

pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_HeatPlot_Top20.pdf",width=8,height=10)
pheatmap(hmat,scale="row",show_rownames = T,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()

pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_HeatPlot_All.pdf",width=8,height=10)
pheatmap(allmat,scale="row",show_rownames = F,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()


#
allmat <- p_regressed[rownames(p_regressed) %in% IL48_elig48_Sign$Gene, colnames(p_regressed) %in% IL48_elig48_coldata$sample]
hmat <- p_regressed[rownames(p_regressed) %in% IL48_elig48_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL48_elig48_coldata$sample]

anno <- IL48_elig48_coldata
Treatment <- c("#FFA500", "#BF40BF")
names(Treatment) <- c("IL1b", "IL1b_Elig")
anno_colors <- list(Treatment = Treatment)
rownames(anno) <- IL48_elig48_coldata$sample

pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_HeatPlot_Top20.pdf",width=8,height=10)
pheatmap(hmat,scale="row",show_rownames = T,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()

pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_HeatPlot_All.pdf",width=8,height=10)
pheatmap(allmat,scale="row",show_rownames = F,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()


#
allmat <- p_regressed[rownames(p_regressed) %in% IL72_elig72_Sign$Gene, colnames(p_regressed) %in% IL72_elig72_coldata$sample]
hmat <- p_regressed[rownames(p_regressed) %in% IL72_elig72_Sign_df_top_labelled$Gene, colnames(p_regressed) %in% IL72_elig72_coldata$sample]

anno <- IL72_elig72_coldata
Treatment <- c("#FFA500", "#BF40BF")
names(Treatment) <- c("IL1b", "IL1b_Elig")
anno_colors <- list(Treatment = Treatment)
rownames(anno) <- IL72_elig72_coldata$sample

pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_HeatPlot_Top20.pdf",width=8,height=10)
pheatmap(hmat,scale="row",show_rownames = T,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()

pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_HeatPlot_All.pdf",width=8,height=10)
pheatmap(allmat,scale="row",show_rownames = F,
         annotation_col=anno, annotation_colors = anno_colors)
dev.off()



###  toppgene  ###

# IL24
IL24_elig24_Sign_df$p_val<-sign_posthoc_pval[sign_posthoc_pval$Gene %in% IL24_elig24_Sign_df$Gene, ]$il1b_time24_il1b_elig_time24
IL24_elig24_Sign_df$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL24_elig24_toppgene<-toppFun(IL24_elig24_Sign_df, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
IL24_elig24_toppgene$Cluster<-'bulk'

save(IL24_elig24_toppgene, file="./outs/dge/IL24_ILelig24/IL24_ILelig24_toppgene_queries.RData")

pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(IL24_elig24_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(IL24_elig24_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(IL24_elig24_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(IL24_elig24_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()



# 
IL48_elig48_Sign_df$p_val<-sign_posthoc_pval[sign_posthoc_pval$Gene %in% IL48_elig48_Sign_df$Gene, ]$il1b_time48_il1b_elig_time48
IL48_elig48_Sign_df$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL48_elig48_toppgene<-toppFun(IL48_elig48_Sign_df, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
IL48_elig48_toppgene$Cluster<-'bulk'

save(IL48_elig48_toppgene, file="./outs/dge/IL48_ILelig48/IL48_ILelig48_toppgene_queries.RData")

pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(IL48_elig48_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(IL48_elig48_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(IL48_elig48_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(IL48_elig48_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()



# 
IL72_elig72_Sign_df$p_val<-sign_posthoc_pval[sign_posthoc_pval$Gene %in% IL72_elig72_Sign_df$Gene, ]$il1b_time72_il1b_elig_time72
IL72_elig72_Sign_df$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL72_elig72_toppgene<-toppFun(IL72_elig72_Sign_df, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
IL72_elig72_toppgene$Cluster<-'bulk'

save(IL72_elig72_toppgene, file="./outs/dge/IL72_ILelig72/IL72_ILelig72_toppgene_queries.RData")

pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(IL72_elig72_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(IL72_elig72_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(IL72_elig72_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(IL72_elig72_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()



#

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL24_elig24$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL24_elig24$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL24_elig24, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_all_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=6)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL24_elig24_Sign_df$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL24_elig24_Sign_df$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL24_elig24_Sign_df, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL24_ILelig24/IL24_ILelig24_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()


#

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL48_elig48$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL48_elig48$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL48_elig48, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_all_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=6)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL48_elig48_Sign_df$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL48_elig48_Sign_df$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL48_elig48_Sign_df, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL48_ILelig48/IL48_ILelig48_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()


#

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL72_elig72$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL72_elig72$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL72_elig72, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_all_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=6)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()

# Create a named vector with fold changes as rank and run GSEA
ranks <- IL72_elig72_Sign_df$logFC # take the FOLD CHANGE of all the genes
names(ranks) <- IL72_elig72_Sign_df$Gene
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, minSize = 10)
# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  dplyr::rename(Gene = "SYMBOL") %>%
  inner_join(IL72_elig72_Sign_df, by="Gene")
filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)
tmp <- fgseaResTidy
tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)
pdf("./outs/dge/IL72_ILelig72/IL72_ILelig72_GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
ggbarplot(tmp %>% filter(ES > 0),
          x = "pathway", y = "NES",
          fill = "adjPvalue",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = c("#000004FF","#CD4071FF"),            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Enrichment Score",
          legend.title = "FDR < 0.05",
          rotate = TRUE,
          ggtheme = theme_minimal()
          )
dev.off()



##### GENES OF INTEREST ###############


dir.create("./outs/dge/IL24_CTRL24/GOI")
dir.create("./outs/dge/IL48_CTRL48/GOI")
dir.create("./outs/dge/IL72_CTRL72/GOI")
dir.create("./outs/dge/ILelig24_CTRL24/GOI")
dir.create("./outs/dge/ILelig48_CTRL48/GOI")
dir.create("./outs/dge/ILelig72_CTRL72/GOI")
dir.create("./outs/dge/IL24_ILelig24/GOI")
dir.create("./outs/dge/IL48_ILelig48/GOI")
dir.create("./outs/dge/IL72_ILelig72/GOI")




#GOI <- read_xlsx("./outs/counts/Genes_of_Interest.xlsx")
#GOI <- unique(append(GOI$Gene, c("ERN1", "XBP1", "HSPA13")))

GOI<-c("ATF5","C3","CANX","CCL2","CEBPD","CXCL8","DNAJB11","DNAJC10","DNAJC3","ERLIN1","ERN1","GSN",
"HSPA13","IL6","INF2","MAP2K3","MGST1","NFE2L1","NFKB2","OSBPL10","P4HA2","PDIA4","PTPN12",
"RCN1","RPN1","SLC39A14","SOD2","TMX1","TNIP1","VAPA","XBP1")

### box plots ###

# IL24
IL24_CTRL24_Sign_df_GOI <- IL24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI)

openxlsx::write.xlsx(IL24_CTRL24_Sign_df_GOI,
                     file = "./outs/dge/IL24_CTRL24/GOI/GOI_IL24_CTRL24_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL24_CTRL24_Sign_df_GOI_top_labelled <- IL24_CTRL24_Sign_df_GOI %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)

#IL24_CTRL24_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 24, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 24, ]$Treatment)


# IL48
IL48_CTRL48_Sign_df_GOI <- IL48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI)

openxlsx::write.xlsx(IL48_CTRL48_Sign_df_GOI,
                     file = "./outs/dge/IL48_CTRL48/GOI/GOI_IL48_CTRL48_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL48_CTRL48_Sign_df_GOI_top_labelled <- IL48_CTRL48_Sign_df_GOI %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)

#IL48_CTRL48_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 48, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 48, ]$Treatment)


# IL72
IL72_CTRL72_Sign_df_GOI <- IL72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI)

openxlsx::write.xlsx(IL72_CTRL72_Sign_df_GOI,
                     file = "./outs/dge/IL72_CTRL72/GOI/GOI_IL72_CTRL72_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL72_CTRL72_Sign_df_GOI_top_labelled <- IL72_CTRL72_Sign_df_GOI %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)

#IL72_CTRL72_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 72, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'Untreated') & meta$Time == 72, ]$Treatment)


# ILelig24
ILelig24_CTRL24_Sign_df_GOI <- ILelig24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI)

openxlsx::write.xlsx(ILelig24_CTRL24_Sign_df_GOI,
                     file = "./outs/dge/ILelig24_CTRL24/GOI/GOI_ILelig24_CTRL24_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

ILelig24_CTRL24_Sign_df_GOI_top_labelled <- ILelig24_CTRL24_Sign_df_GOI %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)

#ILelig24_CTRL24_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 24, ]), Treatment=meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 24, ]$Treatment)



# ILelig48
ILelig48_CTRL48_Sign_df_GOI <- ILelig48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI)

openxlsx::write.xlsx(ILelig48_CTRL48_Sign_df_GOI,
                     file = "./outs/dge/ILelig48_CTRL48/GOI/GOI_ILelig48_CTRL48_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

ILelig48_CTRL48_Sign_df_GOI_top_labelled <- ILelig48_CTRL48_Sign_df_GOI %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)

#ILelig48_CTRL48_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 48, ]), Treatment=meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 48, ]$Treatment)


# ILelig72
ILelig72_CTRL72_Sign_df_GOI <- ILelig72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI)

openxlsx::write.xlsx(ILelig72_CTRL72_Sign_df_GOI,
                     file = "./outs/dge/ILelig72_CTRL72/GOI/GOI_ILelig72_CTRL72_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

ILelig72_CTRL72_Sign_df_GOI_top_labelled <- ILelig72_CTRL72_Sign_df_GOI %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)

#ILelig72_CTRL72_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 72, ]), Treatment=meta[meta$Treatment %in% c('IL1b_Elig', 'Untreated') & meta$Time == 72, ]$Treatment)


#
IL24_elig24_Sign_df_GOI <- IL24_elig24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI)

openxlsx::write.xlsx(IL24_elig24_Sign_df_GOI,
                     file = "./outs/dge/IL24_ILelig24/GOI/GOI_IL24_ILelig24_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL24_elig24_Sign_df_GOI_top_labelled <- IL24_elig24_Sign_df_GOI %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)

#IL24_elig24_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 24, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 24, ]$Treatment)

#
IL48_elig48_Sign_df_GOI <- IL48_elig48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI)

openxlsx::write.xlsx(IL48_elig48_Sign_df_GOI,
                     file = "./outs/dge/IL48_ILelig48/GOI/GOI_IL48_ILelig48_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL48_elig48_Sign_df_GOI_top_labelled <- IL48_elig48_Sign_df_GOI %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)

#IL48_elig48_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 48, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 48, ]$Treatment)

#
IL72_elig72_Sign_df_GOI <- IL72_elig72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI)

openxlsx::write.xlsx(IL72_elig72_Sign_df_GOI,
                     file = "./outs/dge/IL72_ILelig72/GOI/GOI_IL72_ILelig72_Sign_DGE_Filtered_padj05_L2FC_.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

IL72_elig72_Sign_df_GOI_top_labelled <- IL72_elig72_Sign_df_GOI %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)

#IL72_elig72_coldata<-data.frame(sample=rownames(meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 72, ]), Treatment=meta[meta$Treatment %in% c('IL1b', 'IL1b_Elig') & meta$Time == 72, ]$Treatment)


### volc plots ###

# EnhancedVolcano(
#     vol, 
#     lab=vol$Gene,
#     x="logFC", 
#     y="FDR", 
#     selectLab=c("ATF5","C3","CANX","CCL2","CEBPD","CXCL8","DNAJB11","DNAJC10","DNAJC3","ERLIN1","ERN1","GSN",
#                 "HSPA13","IL6","INF2","MAP2K3","MGST1","NFE2L1","NFKB2","OSBPL10","P4HA2","PDIA4","PTPN12",
#                 "RCN1","RPN1","SLC39A14","SOD2","TMX1","TNIP1","VAPA","XBP1"), 
#     pCutoff = 0.05,
#     FCcutoff = 0.3,
#     pointSize = 1.0,
#     labSize = 4.0,
#     colAlpha = 1,
#     col=c('darkgrey','darkgrey','darkgrey', "#BF40BF"),
#     legendPosition = "None",
#     drawConnectors = TRUE,
#     widthConnectors = 0.75)

# IL24
vol <- IL24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) #%>%
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL24_CTRL24/GOI/GOI_IL24_CTRL24_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_CTRL24_Sign_df_GOI_top_labelled,
                  mapping = aes(label = Gene),
                  size = 3,  force=20,#size = 2.5, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+ #unit(0.2,"lines")
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# IL48
vol <- IL48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) #%>%
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL48_CTRL48/GOI/GOI_IL48_CTRL48_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_CTRL48_Sign_df_GOI_top_labelled,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# IL72
vol <- IL72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) #%>%
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL72_CTRL72/GOI/GOI_IL72_CTRL72_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL72_CTRL72_Sign_df_GOI_top_labelled,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# ILelig24
vol <- ILelig24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) #%>%
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/ILelig24_CTRL24/GOI/GOI_ILelig24_CTRL24_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig24_CTRL24_Sign_df_GOI_top_labelled,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()



# ILelig48
vol <- ILelig48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) #%>%
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/ILelig48_CTRL48/GOI/GOI_ILelig48_CTRL48_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig48_CTRL48_Sign_df_GOI_top_labelled,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


# ILelig72
vol <- ILelig72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) #%>%
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/ILelig72_CTRL72/GOI/GOI_ILelig72_CTRL72_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig72_CTRL72_Sign_df_GOI_top_labelled,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()



# 
vol <- IL24_elig24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) #%>%
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL24_ILelig24/GOI/GOI_IL24_ILelig24_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_elig24_Sign_df_GOI_top_labelled,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()




# 
vol <- IL48_elig48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) #%>%
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL48_ILelig48/GOI/GOI_IL48_ILelig48_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_elig48_Sign_df_GOI_top_labelled,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


#
vol <- IL72_elig72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) #%>%
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL72_ILelig72/GOI/GOI_IL72_ILelig72_04_Volcano_Plot.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL72_elig72_Sign_df_GOI_top_labelled,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()


###  toppgene  ###

# IL24
IL24_CTRL24_Sign_df_GOI$p_val<-posthoc_pval[rownames_to_column(posthoc_pval, "Gene")$Gene %in% IL24_CTRL24_Sign_df_GOI$Gene, ]$il1b_time24_untreated_time24
IL24_CTRL24_Sign_df_GOI$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL24_CTRL24_GOI_toppgene<-toppFun(IL24_CTRL24_Sign_df_GOI, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
IL24_CTRL24_GOI_toppgene$Cluster<-'bulk'

save(IL24_CTRL24_GOI_toppgene, file="./outs/dge/IL24_CTRL24/GOI/GOI_IL24_CTRL24_toppgene_queries.RData")

pdf("./outs/dge/IL24_CTRL24/GOI/GOI_IL24_CTRL24_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(IL24_CTRL24_GOI_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL24_CTRL24/GOI/GOI_IL24_CTRL24_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(IL24_CTRL24_GOI_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL24_CTRL24/GOI/GOI_IL24_CTRL24_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(IL24_CTRL24_GOI_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL24_CTRL24/GOI/GOI_IL24_CTRL24_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(IL24_CTRL24_GOI_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# IL48
IL48_CTRL48_Sign_df_GOI$p_val<-posthoc_pval[rownames_to_column(posthoc_pval, "Gene")$Gene %in% IL48_CTRL48_Sign_df_GOI$Gene, ]$il1b_time48_untreated_time48
IL48_CTRL48_Sign_df_GOI$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL48_CTRL48_GOI_toppgene<-toppFun(IL48_CTRL48_Sign_df_GOI, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
IL48_CTRL48_GOI_toppgene$Cluster<-'bulk'

save(IL48_CTRL48_GOI_toppgene, file="./outs/dge/IL48_CTRL48/GOI/GOI_IL48_CTRL48_toppgene_queries.RData")

pdf("./outs/dge/IL48_CTRL48/GOI/GOI_IL48_CTRL48_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(IL48_CTRL48_GOI_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL48_CTRL48/GOI/GOI_IL48_CTRL48_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(IL48_CTRL48_GOI_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL48_CTRL48/GOI/GOI_IL48_CTRL48_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(IL48_CTRL48_GOI_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL48_CTRL48/GOI/GOI_IL48_CTRL48_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(IL48_CTRL48_GOI_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# IL72
IL72_CTRL72_Sign_df_GOI$p_val<-posthoc_pval[rownames_to_column(posthoc_pval, "Gene")$Gene %in% IL72_CTRL72_Sign_df_GOI$Gene, ]$il1b_time72_untreated_time72
IL72_CTRL72_Sign_df_GOI$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL72_CTRL72_GOI_toppgene<-toppFun(IL72_CTRL72_Sign_df_GOI, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
IL72_CTRL72_GOI_toppgene$Cluster<-'bulk'

save(IL72_CTRL72_GOI_toppgene, file="./outs/dge/IL72_CTRL72/GOI/GOI_IL72_CTRL72_toppgene_queries.RData")

pdf("./outs/dge/IL72_CTRL72/GOI/GOI_IL72_CTRL72_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(IL72_CTRL72_GOI_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL72_CTRL72/GOI/GOI_IL72_CTRL72_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(IL72_CTRL72_GOI_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL72_CTRL72/GOI/GOI_IL72_CTRL72_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(IL72_CTRL72_GOI_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/IL72_CTRL72/GOI/GOI_IL72_CTRL72_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(IL72_CTRL72_GOI_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# ILelig24
ILelig24_CTRL24_Sign_df_GOI$p_val<-posthoc_pval[rownames_to_column(posthoc_pval, "Gene")$Gene %in% ILelig24_CTRL24_Sign_df_GOI$Gene, ]$il1b_elig_time72_untreated_time72
ILelig24_CTRL24_Sign_df_GOI$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
ILelig24_CTRL24_GOI_toppgene<-toppFun(ILelig24_CTRL24_Sign_df_GOI, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
ILelig24_CTRL24_GOI_toppgene$Cluster<-'bulk'

save(ILelig24_CTRL24_GOI_toppgene, file="./outs/dge/ILelig24_CTRL24/GOI/GOI_ILelig24_CTRL24_toppgene_queries.RData")

pdf("./outs/dge/ILelig24_CTRL24/GOI/GOI_ILelig24_CTRL24_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig24_CTRL24_GOI_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/GOI/GOI_ILelig24_CTRL24_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig24_CTRL24_GOI_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/GOI/GOI_ILelig24_CTRL24_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig24_CTRL24_GOI_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/GOI/GOI_ILelig24_CTRL24_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(ILelig24_CTRL24_GOI_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# ILelig48
ILelig48_CTRL48_Sign_df_GOI$p_val<-posthoc_pval[rownames_to_column(posthoc_pval, "Gene")$Gene %in% ILelig48_CTRL48_Sign_df_GOI$Gene, ]$il1b_elig_time72_untreated_time72
ILelig48_CTRL48_Sign_df_GOI$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
ILelig48_CTRL48_GOI_toppgene<-toppFun(ILelig48_CTRL48_Sign_df_GOI, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
ILelig48_CTRL48_GOI_toppgene$Cluster<-'bulk'

save(ILelig48_CTRL48_GOI_toppgene, file="./outs/dge/ILelig48_CTRL48/GOI/GOI_ILelig48_CTRL48_toppgene_queries.RData")

pdf("./outs/dge/ILelig48_CTRL48/GOI/GOI_ILelig48_CTRL48_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig48_CTRL48_GOI_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/GOI/GOI_ILelig48_CTRL48_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig48_CTRL48_GOI_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/GOI/GOI_ILelig48_CTRL48_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig48_CTRL48_GOI_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/GOI/GOI_ILelig48_CTRL48_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(ILelig48_CTRL48_GOI_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# ILelig72
ILelig72_CTRL72_Sign_df_GOI$p_val<-posthoc_pval[rownames_to_column(posthoc_pval, "Gene")$Gene %in% ILelig72_CTRL72_Sign_df_GOI$Gene, ]$il1b_elig_time72_untreated_time72
ILelig72_CTRL72_Sign_df_GOI$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
ILelig72_CTRL72_GOI_toppgene<-toppFun(ILelig72_CTRL72_Sign_df_GOI, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")
ILelig72_CTRL72_GOI_toppgene$Cluster<-'bulk'

save(ILelig72_CTRL72_GOI_toppgene, file="./outs/dge/ILelig72_CTRL72/GOI/GOI_ILelig72_CTRL72_toppgene_queries.RData")

pdf("./outs/dge/ILelig72_CTRL72/GOI/GOI_ILelig72_CTRL72_molecular_function_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig72_CTRL72_GOI_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/GOI/GOI_ILelig72_CTRL72_biological_process_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig72_CTRL72_GOI_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/GOI/GOI_ILelig72_CTRL72_cellular_component_gene_ont.pdf",width=8,height=8)
toppPlot(ILelig72_CTRL72_GOI_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/GOI/GOI_ILelig72_CTRL72_pathway_gene_ont.pdf",width=16,height=8)
toppPlot(ILelig72_CTRL72_GOI_toppgene, category="Pathway", clusters=c("bulk"))
dev.off()


# IL24
IL24_elig24_Sign_df_GOI$p_val<-posthoc_pval[rownames_to_column(posthoc_pval, "Gene")$Gene %in% IL24_elig24_Sign_df_GOI$Gene, ]$il1b_time24_il1b_elig_time24
IL24_elig24_Sign_df_GOI$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL24_elig24_GOI_toppgene<-toppFun(IL24_elig24_Sign_df_GOI, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")

# NO RESULTS --^

# IL24_elig24_GOI_toppgene$Cluster<-'bulk'

# save(IL24_elig24_GOI_toppgene, file="./outs/dge/IL24_ILelig24/GOI/GOI_IL24_ILelig24_toppgene_queries.RData")

# pdf("./outs/dge/IL24_ILelig24/GOI/GOI_IL24_ILelig24_molecular_function_gene_ont.pdf",width=8,height=8)
# toppPlot(IL24_elig24_GOI_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
# dev.off()

# pdf("./outs/dge/IL24_ILelig24/GOI/GOI_IL24_ILelig24_biological_process_gene_ont.pdf",width=8,height=8)
# toppPlot(IL24_elig24_GOI_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
# dev.off()

# pdf("./outs/dge/IL24_ILelig24/GOI/GOI_IL24_ILelig24_cellular_component_gene_ont.pdf",width=8,height=8)
# toppPlot(IL24_elig24_GOI_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
# dev.off()

# pdf("./outs/dge/IL24_ILelig24/GOI/GOI_IL24_ILelig24_pathway_gene_ont.pdf",width=16,height=8)
# toppPlot(IL24_elig24_GOI_toppgene, category="Pathway", clusters=c("bulk"))
# dev.off()

# 
IL48_elig48_Sign_df_GOI$p_val<-posthoc_pval[rownames_to_column(posthoc_pval, "Gene")$Gene %in% IL48_elig48_Sign_df_GOI$Gene, ]$il1b_time48_il1b_elig_time48
IL48_elig48_Sign_df_GOI$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL48_elig48_GOI_toppgene<-toppFun(IL48_elig48_Sign_df_GOI, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")

# NO RESULTS --^

# IL48_elig48_GOI_toppgene$Cluster<-'bulk'

# save(IL48_elig48_GOI_toppgene, file="./outs/dge/IL48_ILelig48/GOI/GOI_IL48_ILelig48_toppgene_queries.RData")

# pdf("./outs/dge/IL48_ILelig48/GOI/GOI_IL48_ILelig48_molecular_function_gene_ont.pdf",width=8,height=8)
# toppPlot(IL48_elig48_GOI_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
# dev.off()

# pdf("./outs/dge/IL48_ILelig48/GOI/GOI_IL48_ILelig48_biological_process_gene_ont.pdf",width=8,height=8)
# toppPlot(IL48_elig48_GOI_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
# dev.off()

# pdf("./outs/dge/IL48_ILelig48/GOI/GOI_IL48_ILelig48_cellular_component_gene_ont.pdf",width=8,height=8)
# toppPlot(IL48_elig48_GOI_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
# dev.off()

# pdf("./outs/dge/IL48_ILelig48/GOI/GOI_IL48_ILelig48_pathway_gene_ont.pdf",width=16,height=8)
# toppPlot(IL48_elig48_GOI_toppgene, category="Pathway", clusters=c("bulk"))
# dev.off()

# 
IL72_elig72_Sign_df_GOI$p_val<-posthoc_pval[rownames_to_column(posthoc_pval, "Gene")$Gene %in% IL72_elig72_Sign_df_GOI$Gene, ]$il1b_time72_il1b_elig_time72
IL72_elig72_Sign_df_GOI$cluster<-"bulk"

#httr::set_config(config(ssl_verifypeer=0)) # for whatever reason you get an curl SSL error about self-signed certs if you don't set this
#Sys.setenv(CURL_CA_BUNDLE="/private/etc/ssl/cacert.pem")
IL72_elig72_GOI_toppgene<-toppFun(IL72_elig72_Sign_df_GOI, gene_col="Gene", p_val_col="p_val", logFC_col="logFC", cluster_col="cluster")

# NO RESULTS --^

# IL72_elig72_GOI_toppgene$Cluster<-'bulk'

# save(IL72_elig72_GOI_toppgene, file="./outs/dge/IL72_ILelig72/GOI/GOI_IL72_ILelig72_toppgene_queries.RData")

# pdf("./outs/dge/IL72_ILelig72/GOI/GOI_IL72_ILelig72_molecular_function_gene_ont.pdf",width=8,height=8)
# toppPlot(IL72_elig72_GOI_toppgene, category="GeneOntologyMolecularFunction", clusters=c("bulk"))
# dev.off()

# pdf("./outs/dge/IL72_ILelig72/GOI/GOI_IL72_ILelig72_biological_process_gene_ont.pdf",width=8,height=8)
# toppPlot(IL72_elig72_GOI_toppgene, category="GeneOntologyBiologicalProcess", clusters=c("bulk"))
# dev.off()

# pdf("./outs/dge/IL72_ILelig72/GOI/GOI_IL72_ILelig72_cellular_component_gene_ont.pdf",width=8,height=8)
# toppPlot(IL72_elig72_GOI_toppgene, category="GeneOntologyCellularComponent", clusters=c("bulk"))
# dev.off()

# pdf("./outs/dge/IL72_ILelig72/GOI/GOI_IL72_ILelig72_pathway_gene_ont.pdf",width=16,height=8)
# toppPlot(IL72_elig72_GOI_toppgene, category="Pathway", clusters=c("bulk"))
# dev.off()



##### GENES OF INTEREST 2 ###############
  
dir.create("./outs/dge/IL24_CTRL24/GOI_2")
dir.create("./outs/dge/IL48_CTRL48/GOI_2")
dir.create("./outs/dge/IL72_CTRL72/GOI_2")
dir.create("./outs/dge/ILelig24_CTRL24/GOI_2")
dir.create("./outs/dge/ILelig48_CTRL48/GOI_2")
dir.create("./outs/dge/ILelig72_CTRL72/GOI_2")
dir.create("./outs/dge/IL24_ILelig24/GOI_2")
dir.create("./outs/dge/IL48_ILelig48/GOI_2")
dir.create("./outs/dge/IL72_ILelig72/GOI_2")

#GOI <- read_xlsx("./outs/counts/Genes_of_Interest.xlsx")
#GOI <- unique(append(GOI$Gene, c("ERN1", "XBP1", "HSPA13")))

GOI_2<-c("CXCL8", "SOD2", "ERN1", "HSPA13", "CANX", "DNAJC3", "NFE2L1")

### box plots ###

# IL24
IL24_CTRL24_Sign_df_GOI_2 <- IL24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_2)


IL24_CTRL24_Sign_df_GOI_top_labelled_2 <- IL24_CTRL24_Sign_df_GOI_2 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# IL48
IL48_CTRL48_Sign_df_GOI_2 <- IL48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_2)


IL48_CTRL48_Sign_df_GOI_top_labelled_2 <- IL48_CTRL48_Sign_df_GOI_2 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# IL72
IL72_CTRL72_Sign_df_GOI_2 <- IL72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_2)


IL72_CTRL72_Sign_df_GOI_top_labelled_2 <- IL72_CTRL72_Sign_df_GOI_2 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# ILelig24
ILelig24_CTRL24_Sign_df_GOI_2 <- ILelig24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_2)


ILelig24_CTRL24_Sign_df_GOI_top_labelled_2 <- ILelig24_CTRL24_Sign_df_GOI_2 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# ILelig48
ILelig48_CTRL48_Sign_df_GOI_2 <- ILelig48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_2)


ILelig48_CTRL48_Sign_df_GOI_top_labelled_2 <- ILelig48_CTRL48_Sign_df_GOI_2 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# ILelig72
ILelig72_CTRL72_Sign_df_GOI_2 <- ILelig72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_2)


ILelig72_CTRL72_Sign_df_GOI_top_labelled_2 <- ILelig72_CTRL72_Sign_df_GOI_2 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


#
IL24_elig24_Sign_df_GOI_2 <- IL24_elig24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_2)


IL24_elig24_Sign_df_GOI_top_labelled_2 <- IL24_elig24_Sign_df_GOI_2 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


#
IL48_elig48_Sign_df_GOI_2 <- IL48_elig48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_2)


IL48_elig48_Sign_df_GOI_top_labelled_2 <- IL48_elig48_Sign_df_GOI_2 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


#
IL72_elig72_Sign_df_GOI_2 <- IL72_elig72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_2)


IL72_elig72_Sign_df_GOI_top_labelled_2 <- IL72_elig72_Sign_df_GOI_2 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)



### volc plots ###

# EnhancedVolcano(
#     vol, 
#     lab=vol$Gene,
#     x="logFC", 
#     y="FDR", 
#     selectLab=c("ATF5","C3","CANX","CCL2","CEBPD","CXCL8","DNAJB11","DNAJC10","DNAJC3","ERLIN1","ERN1","GSN",
#                 "HSPA13","IL6","INF2","MAP2K3","MGST1","NFE2L1","NFKB2","OSBPL10","P4HA2","PDIA4","PTPN12",
#                 "RCN1","RPN1","SLC39A14","SOD2","TMX1","TNIP1","VAPA","XBP1"), 
#     pCutoff = 0.05,
#     FCcutoff = 0.3,
#     pointSize = 1.0,
#     labSize = 4.0,
#     colAlpha = 1,
#     col=c('darkgrey','darkgrey','darkgrey', "#BF40BF"),
#     legendPosition = "None",
#     drawConnectors = TRUE,
#     widthConnectors = 0.75)

# IL24
vol <- IL24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL24_CTRL24_Sign_df_GOI_top_labelled_2$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL24_CTRL24/GOI_2/GOI_IL24_CTRL24_04_Volcano_Plot_2.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_CTRL24_Sign_df_GOI_top_labelled_2,
                  mapping = aes(label = Gene),
                  size = 3,  force=20,#size = 2.5, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+ #unit(0.2,"lines")
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL24_CTRL24/GOI_2/GOI_IL24_CTRL24_04_Volcano_Plot_2_2.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#FF7777", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#7DA8E6", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "darkred", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL24_CTRL24_Sign_df_GOI_top_labelled_2, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$")  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


# IL48
vol <- IL48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL48_CTRL48_Sign_df_GOI_top_labelled_2$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL48_CTRL48/GOI_2/GOI_IL48_CTRL48_04_Volcano_Plot_2.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_CTRL48_Sign_df_GOI_top_labelled_2,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL48_CTRL48/GOI_2/GOI_IL48_CTRL48_04_Volcano_Plot_2_2.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#FF7777", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#7DA8E6", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "darkred", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL48_CTRL48_Sign_df_GOI_top_labelled_2, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$")  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


# IL72
vol <- IL72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL72_CTRL72_Sign_df_GOI_top_labelled_2$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL72_CTRL72/GOI_2/GOI_IL72_CTRL72_04_Volcano_Plot_2.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL72_CTRL72_Sign_df_GOI_top_labelled_2,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL72_CTRL72/GOI_2/GOI_IL72_CTRL72_04_Volcano_Plot_2_2.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#FF7777", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#7DA8E6", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "darkred", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL72_CTRL72_Sign_df_GOI_top_labelled_2, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$")  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()

# ILelig24
vol <- ILelig24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% ILelig24_CTRL24_Sign_df_GOI_top_labelled_2$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/ILelig24_CTRL24/GOI_2/GOI_ILelig24_CTRL24_04_Volcano_Plot_2.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig24_CTRL24_Sign_df_GOI_top_labelled_2,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/GOI_2/GOI_ILelig24_CTRL24_04_Volcano_Plot_2_2.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#FF7777", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#7DA8E6", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "darkred", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = ILelig24_CTRL24_Sign_df_GOI_top_labelled_2, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$")  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()




# ILelig48
vol <- ILelig48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% ILelig48_CTRL48_Sign_df_GOI_top_labelled_2$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/ILelig48_CTRL48/GOI_2/GOI_ILelig48_CTRL48_04_Volcano_Plot_2.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig48_CTRL48_Sign_df_GOI_top_labelled_2,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/GOI_2/GOI_ILelig48_CTRL48_04_Volcano_Plot_2_2.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#FF7777", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#7DA8E6", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "darkred", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = ILelig48_CTRL48_Sign_df_GOI_top_labelled_2, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$")  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


# ILelig72
vol <- ILelig72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% ILelig72_CTRL72_Sign_df_GOI_top_labelled_2$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/ILelig72_CTRL72/GOI_2/GOI_ILelig72_CTRL72_04_Volcano_Plot_2.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig72_CTRL72_Sign_df_GOI_top_labelled_2,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/GOI_2/GOI_ILelig72_CTRL72_04_Volcano_Plot_2_2.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#FF7777", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#7DA8E6", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "darkred", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = ILelig72_CTRL72_Sign_df_GOI_top_labelled_2, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$")  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()

# 
vol <- IL24_elig24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL24_elig24_Sign_df_GOI_top_labelled_2$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL24_ILelig24/GOI_2/GOI_IL24_ILelig24_04_Volcano_Plot_2.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_elig24_Sign_df_GOI_top_labelled_2,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL24_ILelig24/GOI_2/GOI_IL24_ILelig24_04_Volcano_Plot_2_2.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#FF7777", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#7DA8E6", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "darkred", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL24_elig24_Sign_df_GOI_top_labelled_2, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$")  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


# 
vol <- IL48_elig48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL48_elig48_Sign_df_GOI_top_labelled_2$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL48_ILelig48/GOI_2/GOI_IL48_ILelig48_04_Volcano_Plot_2.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_elig48_Sign_df_GOI_top_labelled_2,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL48_ILelig48/GOI_2/GOI_IL48_ILelig48_04_Volcano_Plot_2_2.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#FF7777", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#7DA8E6", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "darkred", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL48_elig48_Sign_df_GOI_top_labelled_2, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$")  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


#
vol <- IL72_elig72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL72_elig72_Sign_df_GOI_top_labelled_2$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL72_ILelig72/GOI_2/GOI_IL72_ILelig72_04_Volcano_Plot_2.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL72_elig72_Sign_df_GOI_top_labelled_2,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL72_ILelig72/GOI_2/GOI_IL72_ILelig72_04_Volcano_Plot_2_2.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#FF7777", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#7DA8E6", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "darkred", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL72_elig72_Sign_df_GOI_top_labelled_2, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$")  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$") +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()




##### GENES OF INTEREST 3 ###############
  
dir.create("./outs/dge/IL24_CTRL24/GOI_3")
dir.create("./outs/dge/IL48_CTRL48/GOI_3")
dir.create("./outs/dge/IL72_CTRL72/GOI_3")
dir.create("./outs/dge/ILelig24_CTRL24/GOI_3")
dir.create("./outs/dge/ILelig48_CTRL48/GOI_3")
dir.create("./outs/dge/ILelig72_CTRL72/GOI_3")
dir.create("./outs/dge/IL24_ILelig24/GOI_3")
dir.create("./outs/dge/IL48_ILelig48/GOI_3")
dir.create("./outs/dge/IL72_ILelig72/GOI_3")


GOI_3<-c("CXCL8", "SOD2", "ERN1", "HSPA13", "CANX", "DNAJC3", "NFE2L1", "CXCL5", "CXCL8",
         "CXCL5", "CXCL6", "CXCL1", "CCL2", "C3", "IL6", "SOD2", "CANX", "DNAJC3", "ERN1",
         "HSPA13")

GOI_3 <- unique(GOI_3)

### box plots ###

# IL24
IL24_CTRL24_Sign_df_GOI_3 <- IL24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_3)


IL24_CTRL24_Sign_df_GOI_top_labelled_3 <- IL24_CTRL24_Sign_df_GOI_3 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# IL48
IL48_CTRL48_Sign_df_GOI_3 <- IL48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_3)


IL48_CTRL48_Sign_df_GOI_top_labelled_3 <- IL48_CTRL48_Sign_df_GOI_3 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# IL72
IL72_CTRL72_Sign_df_GOI_3 <- IL72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_3)


IL72_CTRL72_Sign_df_GOI_top_labelled_3 <- IL72_CTRL72_Sign_df_GOI_3 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# ILelig24
ILelig24_CTRL24_Sign_df_GOI_3 <- ILelig24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_3)


ILelig24_CTRL24_Sign_df_GOI_top_labelled_3 <- ILelig24_CTRL24_Sign_df_GOI_3 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# ILelig48
ILelig48_CTRL48_Sign_df_GOI_3 <- ILelig48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_3)


ILelig48_CTRL48_Sign_df_GOI_top_labelled_3 <- ILelig48_CTRL48_Sign_df_GOI_3 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


# ILelig72
ILelig72_CTRL72_Sign_df_GOI_3 <- ILelig72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_3)


ILelig72_CTRL72_Sign_df_GOI_top_labelled_3 <- ILelig72_CTRL72_Sign_df_GOI_3 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


#
IL24_elig24_Sign_df_GOI_3 <- IL24_elig24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_3)


IL24_elig24_Sign_df_GOI_top_labelled_3 <- IL24_elig24_Sign_df_GOI_3 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


#
IL48_elig48_Sign_df_GOI_3 <- IL48_elig48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_3)


IL48_elig48_Sign_df_GOI_top_labelled_3 <- IL48_elig48_Sign_df_GOI_3 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)


#
IL72_elig72_Sign_df_GOI_3 <- IL72_elig72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Gene %in% GOI_3)


IL72_elig72_Sign_df_GOI_top_labelled_3 <- IL72_elig72_Sign_df_GOI_3 %>%
  group_by(Direction) %>%
  #na.omit() %>%
  arrange(FDR) #%>%
  #top_n(n = 10, wt = LOG)



### volc plots ###

# EnhancedVolcano(
#     vol, 
#     lab=vol$Gene,
#     x="logFC", 
#     y="FDR", 
#     selectLab=c("ATF5","C3","CANX","CCL2","CEBPD","CXCL8","DNAJB11","DNAJC10","DNAJC3","ERLIN1","ERN1","GSN",
#                 "HSPA13","IL6","INF2","MAP2K3","MGST1","NFE2L1","NFKB2","OSBPL10","P4HA2","PDIA4","PTPN12",
#                 "RCN1","RPN1","SLC39A14","SOD2","TMX1","TNIP1","VAPA","XBP1"), 
#     pCutoff = 0.05,
#     FCcutoff = 0.3,
#     pointSize = 1.0,
#     labSize = 4.0,
#     colAlpha = 1,
#     col=c('darkgrey','darkgrey','darkgrey', "#BF40BF"),
#     legendPosition = "None",
#     drawConnectors = TRUE,
#     widthConnectors = 0.75)

# IL24
vol <- IL24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL24_CTRL24_Sign_df_GOI_top_labelled_3$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL24_CTRL24/GOI_3/GOI_IL24_CTRL24_04_Volcano_Plot_3.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_CTRL24_Sign_df_GOI_top_labelled_3,
                  mapping = aes(label = Gene),
                  size = 3,  force=20,#size = 2.5, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+ #unit(0.2,"lines")
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL24_CTRL24/GOI_3/GOI_IL24_CTRL24_04_Volcano_Plot_3_3.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#BF40BF", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#FFA500", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "purple4", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL24_CTRL24_Sign_df_GOI_top_labelled_3, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$",fontsize=10)  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$",fontsize=10) +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


# IL48
vol <- IL48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL48_CTRL48_Sign_df_GOI_top_labelled_3$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL48_CTRL48/GOI_3/GOI_IL48_CTRL48_04_Volcano_Plot_3.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_CTRL48_Sign_df_GOI_top_labelled_3,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL48_CTRL48/GOI_3/GOI_IL48_CTRL48_04_Volcano_Plot_3_3.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#BF40BF", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#FFA500", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "purple4", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL48_CTRL48_Sign_df_GOI_top_labelled_3, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$", fontsize=10)  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$", fontsize=10) +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


# IL72
vol <- IL72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL72_CTRL72_Sign_df_GOI_top_labelled_3$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL72_CTRL72/GOI_3/GOI_IL72_CTRL72_04_Volcano_Plot_3.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL72_CTRL72_Sign_df_GOI_top_labelled_3,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL72_CTRL72/GOI_3/GOI_IL72_CTRL72_04_Volcano_Plot_3_3.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#BF40BF", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#FFA500", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "purple4", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL72_CTRL72_Sign_df_GOI_top_labelled_3, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$",fontsize=10)  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$",fontsize=10) +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()

# ILelig24
vol <- ILelig24_CTRL24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% ILelig24_CTRL24_Sign_df_GOI_top_labelled_3$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/ILelig24_CTRL24/GOI_3/GOI_ILelig24_CTRL24_04_Volcano_Plot_3.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig24_CTRL24_Sign_df_GOI_top_labelled_3,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/ILelig24_CTRL24/GOI_3/GOI_ILelig24_CTRL24_04_Volcano_Plot_3_3.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#BF40BF", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#FFA500", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "purple4", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = ILelig24_CTRL24_Sign_df_GOI_top_labelled_3, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$",fontsize=10)  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$",fontsize=10) +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()




# ILelig48
vol <- ILelig48_CTRL48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% ILelig48_CTRL48_Sign_df_GOI_top_labelled_3$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/ILelig48_CTRL48/GOI_3/GOI_ILelig48_CTRL48_04_Volcano_Plot_3.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig48_CTRL48_Sign_df_GOI_top_labelled_3,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/ILelig48_CTRL48/GOI_3/GOI_ILelig48_CTRL48_04_Volcano_Plot_3_3.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#BF40BF", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#FFA500", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "purple4", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = ILelig48_CTRL48_Sign_df_GOI_top_labelled_3, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$",fontsize=10)  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$",fontsize=10) +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


# ILelig72
vol <- ILelig72_CTRL72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% ILelig72_CTRL72_Sign_df_GOI_top_labelled_3$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/ILelig72_CTRL72/GOI_3/GOI_ILelig72_CTRL72_04_Volcano_Plot_3.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = ILelig72_CTRL72_Sign_df_GOI_top_labelled_3,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/ILelig72_CTRL72/GOI_3/GOI_ILelig72_CTRL72_04_Volcano_Plot_3_3.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#BF40BF", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#FFA500", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "purple4", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = ILelig72_CTRL72_Sign_df_GOI_top_labelled_3, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$",fontsize=10)  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$",fontsize=10) +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()

# 
vol <- IL24_elig24 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL24_elig24_Sign_df_GOI_top_labelled_3$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL24_ILelig24/GOI_3/GOI_IL24_ILelig24_04_Volcano_Plot_3.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL24_elig24_Sign_df_GOI_top_labelled_3,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL24_ILelig24/GOI_3/GOI_IL24_ILelig24_04_Volcano_Plot_3_3.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#BF40BF", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#FFA500", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "purple4", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL24_elig24_Sign_df_GOI_top_labelled_3, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$",fontsize=10)  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$",fontsize=10) +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


# 
vol <- IL48_elig48 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL48_elig48_Sign_df_GOI_top_labelled_3$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL48_ILelig48/GOI_3/GOI_IL48_ILelig48_04_Volcano_Plot_3.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL48_elig48_Sign_df_GOI_top_labelled_3,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL48_ILelig48/GOI_3/GOI_IL48_ILelig48_04_Volcano_Plot_3_3.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#BF40BF", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#FFA500", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "purple4", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL48_elig48_Sign_df_GOI_top_labelled_3, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$",fontsize=10)  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$",fontsize=10) +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()


#
vol <- IL72_elig72 %>%
  mutate(LOG = -log10(FDR), ABS = abs(logFC)) %>%
  mutate(Threshold = if_else(FDR < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(logFC > 0.3 & FDR < 0.05 ~ "UpReg",
                               logFC < -0.3 & FDR < 0.05 ~ "DownReg")) %>%
  mutate(candidate = abs(logFC) > 0.3 & FDR < 0.05) %>%
  mutate(coolgenes = Gene %in% IL72_elig72_Sign_df_GOI_top_labelled_3$Gene)
  #dplyr::filter(Gene %in% GOI)

pdf("./outs/dge/IL72_ILelig72/GOI_3/GOI_IL72_ILelig72_04_Volcano_Plot_3.pdf",width=6,height=6,useDingbats=FALSE)

ggscatter(vol,
          x = "logFC",
          y = "LOG",
          color = "Direction",#color = "Threshold",
          palette=c("#FFA500", "#BF40BF"),
          size = 2,max.overlaps = Inf,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = IL72_elig72_Sign_df_GOI_top_labelled_3,
                  mapping = aes(label = Gene),
                  size = 3, force=20,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,16.0) + xlim(-10.0,+10.0)
dev.off()

pdf("./outs/dge/IL72_ILelig72/GOI_3/GOI_IL72_ILelig72_04_Volcano_Plot_3_3.pdf",width=6,height=6,useDingbats=FALSE)
vol %>%
     as.data.frame() %>%
     tidyplot(x = logFC, y = LOG)  %>%
  add_data_points(data = filter_rows(!candidate),color = "lightgrey", rasterize = TRUE, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "UpReg"),color = "#BF40BF", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(candidate, Direction == "DownReg"),color = "#FFA500", alpha = 0.5, size=0.2)  %>%
  add_data_points(data = filter_rows(coolgenes),color = "purple4", alpha = 1, size=1)  %>%
  add_reference_lines(x = c(-0.3, 0.3), y = -log10(0.05))  %>%
  add_data_labels_repel(data = IL72_elig72_Sign_df_GOI_top_labelled_3, label = Gene,
                        color = "#000000", min.segment.length = 0, background = TRUE)  %>%
  adjust_x_axis_title("$Log[2]~fold~change$",fontsize=10)  %>%
  adjust_y_axis_title("$-Log[10]~italic(P)~adjusted$",fontsize=10) +
      theme(legend.position="none")+
      ylim(0,16) + xlim(-10,+10)
dev.off()
