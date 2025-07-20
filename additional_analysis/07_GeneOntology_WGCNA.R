
suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(BiocParallel)
library(ggpubr)
library(magrittr)
library(broom)
library(data.table)
library(cowplot)
library(BiocSingular)
library(clusterProfiler)
library(future.apply)
library(enrichR)
library(scToppR)
library(fgsea)
library(tidyplots)
})

options(future.globals.maxSize = +Inf)
plan(multisession, workers=10)

dir.create("wgcna_output/functional_enrichment/")
mod <- read.table("wgcna_output/ModuleOutput.txt",header=T)


p <- read.table("wgcna_output/modTraitP.txt")
rownames(p) <- gsub("ME","",rownames(p))
cor <- read.table("wgcna_output/modTraitCor.txt")
rownames(cor) <- gsub("ME","",rownames(cor))


p_filt <- p %>% 
      rownames_to_column("Module") %>%
      filter(TreatmentIL1b < 0.05 | TreatmentIL1b_Elig < 0.05 | TreatmentUntreated < 0.05)

cor_filt <- cor[rownames(cor) %in% p_filt$Module,colnames(cor)%in%colnames(p_filt)] %>% 
      rownames_to_column("Module")



l <- split(mod, mod$ModuleColor)

GOI <- list()
GeneOnto <- list()

for(i in 1:length(l)){
GOI[[i]] <- bitr(as.character(l[[i]]$Gene),  
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)

GeneOnto[[i]] <- enrichGO(gene = unique(GOI[[i]]$ENTREZID), 
                     keyType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff  = 1, 
                     qvalueCutoff = 1, 
                     readable = TRUE)

openxlsx::write.xlsx(as.data.frame(GeneOnto[[i]]), 
                     file = sprintf("wgcna_output/functional_enrichment/%s_GO.xlsx", names(l)[[i]]), 
                     colNames = TRUE,
                     rowNames = TRUE, 
                     borders = "columns",
                     sheetName="GeneOnto", 
                     overwrite = TRUE)

}


cur_result <- list()
for(i in 1:length(l)){
    cur_result[[i]] <- GeneOnto[[i]] %>% as.data.frame() %>% remove_rownames()
    cur_result[[i]]$Module <- as.character(unique(mod$ModuleColor))[[i]]
    collapsed_output <- bind_rows(cur_result) %>% as.data.frame()
}


collapsed_output %>%
  write.csv(file='wgcna_output/functional_enrichment/ENRICHR_Modules_GO_terms.csv')


input_bub <- collapsed_output %>% 
    filter(Module %in% cor_filt$Module) %>% 
    mutate(log = -log10(pvalue)) %>%
    group_by(Module) %>%
    top_n(2,log) %>%
    mutate(Term2 = gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", Description, perl=TRUE)) %>%    
    mutate(Term2 = as.factor(Term2), Module = as.factor(Module))


colors <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF","#E16462FF", "#FCA636FF", "#F0F921FF")

pdf("wgcna_output/functional_enrichment/ENRICHR_HDAC5c_SpecificMod_GO_bubblechart.pdf", width = 15, height = 6)
ggballoonplot(input_bub, x = "Module", y = "Term2",
              size = "Count", fill = "log") +
   scale_fill_gradientn(colors = colors) +
     scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  guides(size = FALSE) + 
  coord_flip()
dev.off()



#######################
#######################




dir.create("wgcna_output/functional_enrichment_scToppR/")

mod <- read.table("wgcna_output/ModuleOutput.txt",header=T)
l <- split(mod$Gene, mod$ModuleColor)


toppSalm <-  toppFun(
       l$salmon,
       type = "marker_list",
       topp_categories = NULL,
       cluster_col = "ModuleColor",
       gene_col = "Gene",
       p_val_col = NULL,
       logFC_col = NULL,
       min_genes = 5,
       max_genes = 5000,
       max_results = 50
     )

toppSalm %>%
  write.csv(file='wgcna_output/functional_enrichment_scToppR/SALMON_TOPPGENE.csv')


toppSteel <-  toppFun(
       l$steelblue,
       type = "marker_list",
       topp_categories = NULL,
       cluster_col = "ModuleColor",
       gene_col = "Gene",
       p_val_col = NULL,
       logFC_col = NULL,
       min_genes = 5,
       max_genes = 5000,
       max_results = 50
     )

toppFlo <-  toppFun(
       l$floralwhite,
       type = "marker_list",
       topp_categories = NULL,
       cluster_col = "ModuleColor",
       gene_col = "Gene",
       p_val_col = NULL,
       logFC_col = NULL,
       min_genes = 5,
       max_genes = 5000,
       max_results = 50
     )

toppFlo %>%
  write.csv(file='wgcna_output/functional_enrichment_scToppR/FLORALWHITE_TOPPGENE.csv')

toppSalm$Cluster <- "Salmon"
toppSteel$Cluster <- "Steelblue"
toppFlo$Cluster <- "Floralwhite"
toppData <- rbind(toppSalm,toppSteel,toppFlo)


toppPlot(
       toppData,
       category = "GeneOntologyMolecularFunction",
       num_terms = 5,
       p_val_adj = "BH",
       p_val_display = "log",
       save = TRUE,
       save_dir = "wgcna_output/functional_enrichment_scToppR",
       width = 5,
       height = 3
     )

toppPlot(
       toppData,
       category = "GeneOntologyBiologicalProcess",
       num_terms = 5,
       p_val_adj = "BH",
       p_val_display = "log",
       save = TRUE,
       save_dir = "wgcna_output/functional_enrichment_scToppR",
       width = 5,
       height = 3
     )


toppPlot(
       toppData,
       category = "GeneOntologyCellularComponent",
       num_terms = 5,
       p_val_adj = "BH",
       p_val_display = "log",
       save = TRUE,
       save_dir = "wgcna_output/functional_enrichment_scToppR",
       width = 5,
       height = 3
     )

toppPlot(
       toppData,
       category = "Drug",
       num_terms = 5,
       p_val_adj = "BH",
       p_val_display = "log",
       save = TRUE,
       save_dir = "wgcna_output/functional_enrichment_scToppR",
       width = 5,
       height = 3
     )

toppPlot(
       toppData,
       category = "Disease",
       num_terms = 5,
       p_val_adj = "BH",
       p_val_display = "log",
       save = TRUE,
       save_dir = "wgcna_output/functional_enrichment_scToppR",
       width = 5,
       height = 3
     )


toppPlot(
       toppData,
       category = "Pathway",
       num_terms = 5,
       p_val_adj = "BH",
       p_val_display = "log",
       save = TRUE,
       save_dir = "wgcna_output/functional_enrichment_scToppR",
       width = 5,
       height = 3
     )


# GSEA Salmon and Floralwhite

dir.create("wgcna_output/functional_enrichment_GSEA/")

mod <- read.table("wgcna_output/ModuleOutput.txt",header=T)
l <- split(mod, mod$ModuleColor)
df <- l$salmon %>%
       as.data.frame() %>%
       arrange(desc(kWithin))


pathways.hallmark <- gmtPathways("utils/h.all.v7.5.symbols.gmt")

# Create a named vector with fold changes as rank and run GSEA
ranks <- df$kWithin
names(ranks) <- df$Gene

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks,minSize = 10)

# Tidy result
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)


gene.in.pathway <- pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest(cols = c(SYMBOL)) %>% 
  rename(Gene = "SYMBOL") %>% 
  inner_join(df, by="Gene")


filt <-  fgseaResTidy %>%
            filter(padj < 0.05, ES > 0)

tmp <- fgseaResTidy

tmp$adjPvalue <- ifelse(tmp$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")

tmp$pathway <- gsub("HALLMARK_","",tmp$pathway)

pdf("wgcna_output/functional_enrichment_GSEA/GSEA_enrichment_Hallmark_Barplot.pdf", width=6, height=5)
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



# viz scatter for sctoppR

tmp1   <- toppData %>%
              filter(Cluster == "Salmon" & str_detect(Category, "GeneOntology")) %>%
              mutate(LOG = -log10(PValue)) %>%
              group_by(Category) %>%
              slice_max(order_by = LOG, n = 2) %>%
              ungroup() %>%
              as.data.frame()


a <- tmp1 |>
  tidyplot(x = LOG, y = GenesInTermInQuery, color = Category) %>%  
    add_data_points(alpha = 0.4) %>%
    add_data_labels_repel(data = tmp1, Name)  %>%
    adjust_y_axis_title("$Numeber~of~Genes$") %>%
    adjust_x_axis_title("$-Log[10]~italic(P)~adjusted$") +
    theme(legend.position = "none") +
    ylim(0, 60) + xlim(8, 15)


tmp2   <- toppData %>%
              filter(Cluster == "Floralwhite" & str_detect(Category, "GeneOntology")) %>%
              mutate(LOG = -log10(PValue)) %>%
              group_by(Category) %>%
              slice_max(order_by = LOG, n = 2) %>%
              ungroup() %>%
              as.data.frame()

b <- tmp2 |>
  tidyplot(x = LOG, y = GenesInTermInQuery, color = Category) %>%  
    add_data_points(alpha = 0.4) %>%
    add_data_labels_repel(data = tmp2, Name)  %>%
    adjust_y_axis_title("$Numeber~of~Genes$") %>%
    adjust_x_axis_title("$-Log[10]~italic(P)~adjusted$") +
    theme(legend.position = "bottom") +
    ylim(0, 60) + xlim(5, 15)



multi_plot <- plot_grid(
  a, b,
  labels = c("A", "B"),
  ncol = 2,  # Use nrow = 1 if you want a horizontal layout
  align = 'hv'
)

ggsave("wgcna_output/functional_enrichment_scToppR/multi_category_plot.pdf", multi_plot, width = 15, height = 5)


