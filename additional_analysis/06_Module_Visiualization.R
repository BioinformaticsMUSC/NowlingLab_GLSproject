# WGCNA
suppressPackageStartupMessages({
library(tidyverse)
library(WGCNA)
library(cluster)
library(ggplot2)
library(ggpubr)
library(igraph)
library(reshape2)
library(clusterProfiler)
library(enrichR)
library(cowplot)
library(tidygraph)
library(ggraph)
})
enableWGCNAThreads()


# Module Viz
dir.create("wgcna_output/modules_visualizations/")

load("input/p_regressed.RData")

load("input/initial_data_prep.RData")


# Create a meta file
pd <- meta %>%
     select(-Sex,-Age) %>% 
     mutate(Time = as.numeric(Time),Treatment = as.factor(Treatment)) 

# Start WGCNA analysis

datTraits <- model.matrix(~ .-1, data = pd) %>%
                        as.data.frame()

datExpr <- t(log2(cpm_filt+1))
logCPM <- log2(cpm_filt+1)

exp <- logCPM

mod <- read.table("wgcna_output/ModuleOutput.txt",header=T)



# Hubs 
hubs <- mod %>% 
            group_by(ModuleColor) %>% 
            top_n(n = 10, wt = kWithin) %>%
            as.data.frame()

hub_list <- split(hubs, hubs$ModuleColor)


tmp <- list()
adjMat <- list()
modColors <- list()
valueList <- list()
colorList <- list()
g1 <- list()
layoutFR <- list()

for(i in 1:length(hub_list)){
tmp[[i]] <- exp[rownames(exp) %in% hub_list[[i]]$Gene,]
adjMat[[i]] <- bicor(t(tmp[[i]]))
adjMat[[i]][abs(adjMat[[i]])<0.5] <- 0
adjMat[[i]] <- abs(adjMat[[i]])
adjMat[[i]] <- adjMat[[i]][match(hub_list[[i]]$Gene,rownames(adjMat[[i]])), match(hub_list[[i]]$Gene,colnames(adjMat[[i]]))]
modColors[[i]] <- as.matrix(t(hub_list[[i]]$ModuleColor))
valueList[[i]] <- lapply(1:ncol(modColors[[i]]), function(x) as.numeric(!is.na(modColors[[i]][,x])))
colorList[[i]] <- lapply(1:ncol(modColors[[i]]), function(x) modColors[[i]][,x])
g1[[i]] <- graph.adjacency(as.matrix(adjMat[[i]]),mode="undirected",weighted=TRUE,diag=FALSE)
g1[[i]] <- delete.vertices(g1[[i]],which(degree(g1[[i]])<1))
layoutFR[[i]] <- layout_with_lgl(g1[[i]],maxiter = 500)

pdf(paste("wgcna_output/modules_visualizations/", names(hub_list)[[i]], "_igraph.pdf",sep = ""),width=3.5,height=3.5,useDingbats=FALSE)
plot.igraph(g1[[i]],
            rescale=TRUE,
            vertex.label.dist=1,
            vertex.size=15,
            vertex.label.color="black",
            vertex.color=as.character(hub_list[[i]]$ModuleColor),
            vertex.label.cex=0.8,
            vertex.frame.color="black",
            layout=layoutFR[[i]],
            edge.color=adjustcolor("grey", alpha.f = .3))
dev.off()
}


# GGnet
dir.create("wgcna_output/MODULES_CYTOSCAPE")

files = list.files(path="wgcna_output/Cyto",pattern = 'CytoEdge',full.names =TRUE)
samples <- gsub("wgcna_output/Cyto/CytoEdge", "", files )
samples <- gsub(".txt", "", samples )
c = lapply(files, read.table,sep="\t",fill=TRUE,header=T)
names(myfiles) <- samples

list_mod <- within(myfiles, rm(grey)) 

list_mod_top <- list()
for(i in 1:length(list_mod)){
pdf(paste("wgcna_output/MODULES_CYTOSCAPE/", names(list_mod)[[i]], "_cytoscape.pdf",sep = ""),width=8,height=8,useDingbats=FALSE)
list_mod_top[[i]] <- list_mod[[i]] %>% 
                        top_n(n = 200, wt = weight) %>%
                        as.data.frame() %>%
                        ggnet::ggnet2(label=TRUE, size = "degree", size.cut = 20,alpha = 0.5, color=names(list_mod)[[i]]) +
                        guides(color = FALSE, size = FALSE)
print(list_mod_top[[i]]) 
dev.off()
}

# Only for the two module of interest

salmon <- myfiles$salmon %>%
            select(fromNode, toNode, weight)
dge <- readxl::read_excel("input/ILelig24_CTRL24_Sign_DGE_Filtered_padj05_L2FC_.xlsx")



salmon_filt <- salmon %>% 
        filter(fromNode %in% dge$Gene) %>%
        top_n(n = 100, wt = weight) %>%
        as.data.frame()

graph <- graph_from_edgelist(as.matrix(salmon_filt[,1:2]), directed = FALSE) %>%
            as_tbl_graph()


fc <- dge %>%
        filter(Gene %in% data.frame(graph)$name) %>%
        group_by(Gene) %>%
        summarise(logFC = mean(logFC, na.rm = TRUE))


pdf("wgcna_output/MODULES_CYTOSCAPE/SALMON_IGRAPH_Network.pdf",width=10,height=10,useDingbats=FALSE)
graph %>% 
  activate("nodes") %>% 
  mutate(degree = centrality_degree(),
         betweenness = centrality_betweenness()) %>% 
      full_join(fc, by = c("name" = "Gene")) %>%
ggraph(layout = 'kk') + 
  geom_edge_link0(edge_color = "gray")+
  geom_node_point(shape = 21, aes(size=degree,fill=logFC))+
  geom_node_text(aes(label = name, size = degree), repel = TRUE)+
  #scale_fill_gradient(low="blue",high="red")+
  #scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, breaks = breaks) +
  scale_fill_gradientn(colours = c("white", "salmon", "salmon4"),values = scales::rescale(c(0, 5))) +
  scale_size(range = c(4, 10), guide = "none") +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white', base_family = 'Helvetica') +
  theme(legend.position = "bottom")
dev.off()



steelblue <- myfiles$steelblue %>%
            select(fromNode, toNode, weight)
dge <- readxl::read_excel("input/ILelig24_CTRL24_Sign_DGE_Filtered_padj05_L2FC_.xlsx")



steelblue_filt <- steelblue %>% 
        filter(fromNode %in% dge$Gene) %>%
        top_n(n = 100, wt = weight) %>%
        as.data.frame()

graph <- graph_from_edgelist(as.matrix(steelblue_filt[,1:2]), directed = FALSE) %>%
            as_tbl_graph()


fc <- dge %>%
        filter(Gene %in% data.frame(graph)$name) %>%
        group_by(Gene) %>%
        summarise(logFC = mean(logFC, na.rm = TRUE))


pdf("wgcna_output/MODULES_CYTOSCAPE/STEELBLUE_IGRAPH_Network.pdf",width=10,height=10,useDingbats=FALSE)
graph %>% 
  activate("nodes") %>% 
  mutate(degree = centrality_degree(),
         betweenness = centrality_betweenness()) %>% 
      full_join(fc, by = c("name" = "Gene")) %>%
ggraph(layout = 'kk') + 
  geom_edge_link0(edge_color = "gray")+
  geom_node_point(shape = 21, aes(size=degree,fill=logFC))+
  geom_node_text(aes(label = name, size = degree), repel = TRUE)+
  #scale_fill_gradient(low="blue",high="red")+
  #scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, breaks = breaks) +
  scale_fill_gradientn(colours = c("white", "orange", "darkorange"),values = scales::rescale(c(0, 5))) +
  scale_size(range = c(4, 10), guide = "none") +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white', base_family = 'Helvetica') +
  theme(legend.position = "bottom")
dev.off()



floralwhite <- myfiles$floralwhite %>%
            select(fromNode, toNode, weight)
dge <- readxl::read_excel("input/ILelig24_CTRL24_Sign_DGE_Filtered_padj05_L2FC_.xlsx")



floralwhite_filt <- floralwhite %>% 
        filter(fromNode %in% dge$Gene) %>%
        top_n(n = 100, wt = weight) %>%
        as.data.frame()

graph <- graph_from_edgelist(as.matrix(floralwhite_filt[,1:2]), directed = FALSE) %>%
            as_tbl_graph()


fc <- dge %>%
        filter(Gene %in% data.frame(graph)$name) %>%
        group_by(Gene) %>%
        summarise(logFC = mean(logFC, na.rm = TRUE))


pdf("wgcna_output/MODULES_CYTOSCAPE/FLORALWHITE_IGRAPH_Network.pdf",width=10,height=10,useDingbats=FALSE)
graph %>% 
  activate("nodes") %>% 
  mutate(degree = centrality_degree(),
         betweenness = centrality_betweenness()) %>% 
      full_join(fc, by = c("name" = "Gene")) %>%
ggraph(layout = 'kk') + 
  geom_edge_link0(edge_color = "gray")+
  geom_node_point(shape = 21, aes(size=degree,fill=logFC))+
  geom_node_text(aes(label = name, size = degree), repel = TRUE)+
  #scale_fill_gradient(low="blue",high="red")+
  #scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, breaks = breaks) +
  scale_fill_gradientn(colours = c("white", "green", "darkgreen"),values = scales::rescale(c(0, 5))) +
  scale_size(range = c(4, 10), guide = "none") +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white', base_family = 'Helvetica') +
  theme(legend.position = "bottom")
dev.off()









# Traits
p <- read.table("wgcna_output/modTraitP.txt")
rownames(p) <- gsub("ME","",rownames(p))
cor <- read.table("wgcna_output/modTraitCor.txt")
rownames(cor) <- gsub("ME","",rownames(cor))

p_filt <- p %>% 
      rownames_to_column("Module") %>%
      filter(TreatmentIL1b < 0.05 | TreatmentIL1b_Elig < 0.05 | TreatmentUntreated < 0.05)

cor_filt <- cor[rownames(cor) %in% p_filt$Module,colnames(cor)%in%colnames(p_filt)] %>% 
      rownames_to_column("Module")

p <- melt(p_filt)
cor <- melt(cor_filt)
p$cor <- -1*cor$value
p$log <- -log10(p$value)
p$Direction <- ifelse(p$cor > 0,"Pos","Neg") %>% as.factor()
#p$Module <- factor(p$Module,levels=paste("WM",1:26,sep=""))
#p$variable <- gsub("_LR","",p$variable)

p$log[p$log < 1.3] <- NA



p$abs<- ifelse(is.na(p$log), p$log, p$cor) %>% abs()


pdf("wgcna_output/Module_Traits_Bubble.pdf",width=8,height=2)
ggscatter(
        p, 
                  x = "variable",
                  y = "Module",
                  size="abs",
                  color="Direction",
                  palette=c("blue","red"),
                  alpha = 0.8,
                  xlab = "",
            ylab = "",) +
                  theme_minimal() +
                  rotate_x_text(angle = 45)+
        coord_flip()
dev.off()


# Module Eigengene
tab=read.table("wgcna_output/Matrix_module_correlation.txt")
df=melt(tab)
df$variable=gsub("ME","",df$variable)

coolmod=p_filt$Module
df=df[df$variable %in% coolmod,]
df$Time <- pd$Time

pdf("wgcna_output/Module_Eigengene_CoolMod.pdf",width=8,height=6,useDingbats=FALSE)
ggplot(df, aes(x=Treatment, y=value, fill=Treatment,group=Treatment)) +
geom_pointrange(mapping = aes(x = Treatment, y = value,group=Treatment,colour=Treatment), position=position_dodge(width=0.5),
stat = "summary",
fun.ymin = function(z) {quantile(z,0.25)},
fun.ymax = function(z) {quantile(z,0.75)},
fun.y = median)+
theme_classic()+
geom_hline(yintercept = 0,linetype="dashed")+
scale_colour_manual(values=c("blue", "orange","purple","green"))+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
labs(title="Module EigenGene",x="", y = "Eigengene")+
ylim(-0.3,+0.3) +
facet_wrap(.~variable)
dev.off()

pdf("wgcna_output/Module_Eigengene_CoolMod_Time.pdf",width=8,height=8,useDingbats=FALSE)
ggscatter(df, x = "Time", y = "value",
          add = "loess",                         # Add regression line
          conf.int = FALSE,                          # Add confidence interval
          color = "Treatment", palette = "jco",           # Color by groups "cyl"
          shape = "Treatment"                             # Change point shape by groups "cyl"
          )+
theme_classic() +
facet_wrap(.~variable)
dev.off()


