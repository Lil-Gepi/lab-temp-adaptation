# ── minimal R recipe: read → graph → plot ────────────────────────────────
library(data.table)
library(tidygraph)
library(ggraph)
library(ggrepel)

# 1. read FlyBase TSV (skips “##” comment lines)
dt <- fread("~/Downloads/gene_genetic_interactions_fb_2025_02.tsv",
  sep = "\t",  col.names = c("geneA_symbol","geneA_FBgn","geneB_symbol","geneB_FBgn","type","pub"))[,-6]
gene.universe <- unlist(read.table("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Script/Gene.universe.txt", header = F))
dt <- dt[geneA_FBgn%in%gene.universe&geneB_FBgn%in% gene.universe,]
lab.gene <- unlist(read.table("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Script/Lab.gene.txt", header = F))
temp.gene <- unlist(read.table("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Script/Lab.gene.txt", header = F))
inter.gene <- unlist(read.table("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Script/Lab.gene.txt", header = F))


dt <- dt[geneA_FBgn%in%lab.gene | geneB_FBgn%in% temp.gene| geneB_FBgn%in% inter.gene,]

# 2. build an undirected graph, add degree
g <- as_tbl_graph(dt[, .(geneA_symbol, geneB_symbol, type)], directed = T) 
g
ggraph(g, layout = "auto") +
  geom_edge_link(aes(color = type), alpha = 0.4) +
  geom_node_point(size = 1.2, colour = "steelblue") +
  theme_void()


ggraph(g, layout = 'linear', circular = TRUE) + 
  geom_edge_arc(aes(colour = factor(type))) + 
  coord_fixed()
