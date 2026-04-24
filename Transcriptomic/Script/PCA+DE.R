rm(list=ls())
gc()
####step0: import required libraries and customized function####
library(data.table)
setDTthreads(percent = 80) 
library(edgeR)
library(ggplot2);library(ggrepel);library(cowplot)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(ExpressionNormalizationWorkflow);library(scales)
setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Script/")

## read in data, but keep only the relevant samples of the manuscript ----
counts <- read.csv("../Data/RNAseq_count_table_28182010_AncColdHot.csv",stringsAsFactors = F,row.names = 1)
colnames(counts)
setcolorder(counts, neworder = sort(colnames(counts)))
colnames(counts)
counts <- counts[apply(cpm(counts),1,function(x) !sum(x<1)>=2),]#filtered for lowly expressed genes

cge <- sub("_.*", "", colnames(counts))
evo <- sub("^[^_]+_([A-Z]).*", "\\1", colnames(counts))
group <- paste0(cge,evo)
replicate <- sub("^[^_]+_[A-Z]([0-9]).*", "\\1", colnames(counts))
replicate[which(evo == "A")] <- "0"
bio.sample <- sub("^(([^_]+_[^_]+))_.*", "\\1", colnames(counts))
y <- DGEList(counts=counts, group = group)
y <- calcNormFactors(y, method = "TMM")
#break here
### PCA plot ----
pca <- prcomp(t(cpm(y, log = T)), center = T, scale. = F)
coord <- as.data.frame(pca$x)
coord$cge <- cge
coord$evo <- evo
# Define a custom color scale for your population
pop_colors <- c("A" = "forestgreen", "C" = "steelblue", "H" = "maroon")
pop_labels <- c("A" = "Ancestral", "C" = "Cold-evolved", "H" = "Hot-evolved")

p1 <- ggplot(coord, aes(x = PC1, y = PC2, color = evo)) +
  # Points as normal dots
  geom_point(shape = 16, size = 2) +
  # Text labels not in legend
  geom_text_repel(aes(label = rownames(coord)), show.legend = FALSE, alpha = 0.8) +
  scale_color_manual(name = "Population",values = pop_colors,labels = pop_labels) +
  labs(x = paste0("PC1 (", round(100*pca$sdev[1]^2 / sum(pca$sdev^2), 2), "%)"),
       y = paste0("PC2 (", round(100*pca$sdev[2]^2 / sum(pca$sdev^2), 2), "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  )

p2 <- ggplot(coord, aes(x = PC3, y = PC2, color = evo)) +
  # Points as normal dots
  geom_point(shape = 16, size = 2) +
  # Text labels not in legend
  geom_text_repel(aes(label = rownames(coord)),  show.legend = FALSE, alpha = 0.8) +
  scale_color_manual(name = "Population", values = pop_colors,labels = pop_labels) +
  labs(x = paste0("PC3 (", round(100*pca$sdev[3]^2 / sum(pca$sdev^2), 2), "%)"),
    y = paste0("PC2 (", round(100*pca$sdev[2]^2 / sum(pca$sdev^2), 2), "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  )
# Remove individual legends
p1_no_legend <- p1 + theme(legend.position = "none")
p2_no_legend <- p2 + theme(legend.position = "none")
# Extract legend from one plot (they share the same scale, so either is fine)
legend <- get_legend(p1)
# Combine the two plots side by side
plots_combined <- plot_grid(p1_no_legend, p2_no_legend, ncol = 2, align = "hv")
# Add the legend on the right
combined_with_legend <- plot_grid(plots_combined, legend, rel_widths = c(1, 0.15))
# Add an overall title
final_plot <- plot_grid(
  ggdraw() + draw_label("PCA on gene expression data", fontface = 'bold', x = 0.5),
  combined_with_legend,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Save or display the final figure
png("../Plot/PCA.png", width = 10, height = 5, unit = "in", res = 600)
print(final_plot)
dev.off()




### PVCA ----
#prior to PVCA, we want to see how variance is cumulatively increasing given more PCs considered
# Calculate variance explained by each PC
# var_explained and cum_var from your prcomp() object
var_explained <- pca$sdev^2 / sum(pca$sdev^2)
cum_var <- cumsum(var_explained)

df <- data.frame(PC = 1:length(cum_var),CumulativeVariance = cum_var)

ggplot(df, aes(x = PC, y = CumulativeVariance)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "red") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
  labs(
    title = "Cumulative Variance Explained by Principal Components",
    x = "Principal Component",
    y = "Cumulative Variance Explained"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

## real PVCA starts now
meta.table <- data.frame(cge = cge, evo = evo, row.names = rownames(pca$x))
meta.table$cge_evo <- interaction(meta.table$cge, meta.table$evo, drop = TRUE)
colnames(meta.table)[3] <- "cge:evo"
annot <- data.frame(
  labelDescription = c("Factor levels", "Factor levels", "Factor levels"),
  row.names = c("cge", "evo", "cge:evo")
)
annot_factors <- AnnotatedDataFrame(data = meta.table, varMetadata = annot)

expr.set <- ExpressionSet(
  assayData = as.matrix(cpm(y, log = TRUE)),
  phenoData = annot_factors
)

pvca_res <- pvcAnaly(expr.set, pct_threshold = 0.5, batch_factors = c("cge", "evo", "cge:evo"))
str(pvca_res)

pvca_df <- data.frame(Effect = as.vector(pvca_res$label),var = as.numeric(as.vector(pvca_res$dat)))
pvca_df$Effect <- c("Env:Evo","Evo","Env","Residual")
pvca_df$Effect <- factor(pvca_df$Effect, levels = c("Env", "Evo","Env:Evo","Residual"))

png("../Plot/PVCA.png", width = 5, height = 7, unit = "in", res = 600)

ggplot(pvca_df, aes(x = Effect, y = var)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = scales::percent(var, accuracy = 1)),vjust = -0.5,size = 4) +
  scale_y_continuous(labels = percent_format(accuracy = 1),limits = c(0, 1)) +
  labs(title = "Principal Variant Component Analysis",
    x = NULL,y = "Proportion of Variance Explained") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold")
  )
dev.off()


### DE analysis ----
ModelDesign <- model.matrix(~0+group, list(levels(group)))
colnames(ModelDesign) <- sub("^group", "", colnames(ModelDesign))

y <- estimateDisp(y, ModelDesign, robust=TRUE)
plotBCV(y, xlim = c(2,15))
fit <- glmFit(y, ModelDesign)
mycontrast <- makeContrasts(evo.hot = (hotcgeH + hotcgeC)/2 - hotcgeA,
                            evo.cold = (coldcgeH + coldcgeC)/2 - coldcgeA,
                            hc.hot = hotcgeH - hotcgeC,
                            hc.cold = coldcgeH - coldcgeC,
                            h.hot = hotcgeH - hotcgeA,
                            c.hot = hotcgeC - hotcgeA,
                            h.cold = coldcgeH - coldcgeA,
                            c.cold = coldcgeC - coldcgeA,
                            levels = ModelDesign)
evo.hot <- glmLRT(fit, contrast = mycontrast[,"evo.hot"])$table
evo.hot$p.adjust <- p.adjust(evo.hot$PValue,method = "fdr")
setDT(evo.hot, keep.rownames = T)

result <- copy(evo.hot[,c(1,3)])
setnames(result, 1,"gene")

for (contr in colnames(mycontrast)) {
  print(contr)
  temp <- glmLRT(fit, contrast = mycontrast[,contr])$table
  temp$p.adjust <- p.adjust(temp$PValue,method = "fdr")
  temp_df <- data.frame(logFC = temp$logFC,fdr = temp$p.adjust)
  colnames(temp_df) <- c(paste0("logFC.", contr), paste0("fdr.", contr))
  result <- cbind(result, temp_df)}
rm(temp, temp_df)

###### VennDiagram ----
library(VennDiagram)
venn.diagram(list(hot = result[fdr.evo.hot<0.05,gene],
                  cold = result[fdr.evo.cold<0.05,gene],
                  hc.hot = result[fdr.hc.hot<0.05, gene], 
                  hc.cold = result[fdr.hc.cold<0.05,gene]),
             disable.logging = T, filename = "../Plot/VennDiagram.png")

### lab vs. temp proportion ----
short.list.hot <- result[fdr.evo.hot < 0.05 | fdr.hc.hot < 0.05, .(gene, logFC.evo.hot, logFC.hc.hot, fdr.evo.hot, fdr.hc.hot)]
short.list.hot[fdr.evo.hot < 0.05 & fdr.hc.hot >= 0.05, sig := "Lab.sig"]
short.list.hot[fdr.evo.hot >= 0.05 & fdr.hc.hot < 0.05, sig := "Temp.sig"]
short.list.hot[fdr.evo.hot < 0.05 & fdr.hc.hot < 0.05, sig := "Inter.sig"]
short.list.cold <- result[fdr.evo.cold < 0.05 | fdr.hc.cold < 0.05, .(gene, logFC.evo.cold, logFC.hc.cold, fdr.evo.cold, fdr.hc.cold)]
short.list.cold[fdr.evo.cold < 0.05 & fdr.hc.cold >= 0.05, sig := "Lab.sig"]
short.list.cold[fdr.evo.cold >= 0.05 & fdr.hc.cold < 0.05, sig := "Temp.sig"]
short.list.cold[fdr.evo.cold < 0.05 & fdr.hc.cold < 0.05, sig := "Inter.sig"]
setnames(short.list.hot, old = c("logFC.evo.hot","logFC.hc.hot"), new = c("logFC.evo", "logFC.hc"))
setnames(short.list.cold, old = c("logFC.evo.cold","logFC.hc.cold"), new = c("logFC.evo", "logFC.hc"))
short.list.hot[,cge := "hotcge"]
short.list.cold[,cge := "coldcge"]


## draw the logFC_cold vs logFC_hot plot, highlighting the trhee categories of genes of our interest. ----
logfc.dt.hotcge <- result[,.(gene, logFC.h.hot, logFC.c.hot)]
logfc.dt.coldcge <- result[,.(gene, logFC.h.cold, logFC.c.cold)]
logfc.dt.hotcge <- merge(logfc.dt.hotcge, short.list[cge == "hotcge",.(gene, sig)], by = "gene", all.x = T)
logfc.dt.coldcge <- merge(logfc.dt.coldcge, short.list[cge == "coldcge",.(gene, sig)], by = "gene", all.x = T)
logfc.dt.hotcge$cge <- "hotcge";logfc.dt.coldcge$cge <- "coldcge"
setnames(logfc.dt.hotcge, 2:3, c("logFC.Hot", "logFC.Cold"))
setnames(logfc.dt.coldcge, 2:3, c("logFC.Hot", "logFC.Cold"))
logfc.dt <- rbind(logfc.dt.coldcge, logfc.dt.hotcge)
# Define factor levels to control legend order
logfc.dt[, sig := factor(sig, levels = c("Lab.sig", "Temp.sig", "Inter.sig"))]

bg <- data.frame(
  cge  = c("coldcge", "hotcge"),      
  xmin = -Inf, xmax =  Inf,           
  ymin = -Inf, ymax =  Inf)

png("../Plot/DE.LogFC.png", width = 8, height = 6, units = "in", res = 450)
ggplot(logfc.dt, aes(x = logFC.Hot, y = logFC.Cold)) +
  facet_wrap(~cge)+
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  scale_fill_manual(values = c("coldcge" = "lightblue",
                               "hotcge" = "lightpink")) +
  geom_point(data = logfc.dt[is.na(sig)], color = "grey", size = 1) +
  geom_point(data = logfc.dt[!is.na(sig)], aes(color = sig), size = 1.8, alpha = 0.7) +
  scale_color_manual(
    values = c("Lab.sig" = "purple","Temp.sig" = "yellow","Inter.sig" = "orange"),
    name = "Gene Category",labels = c("Lab-selected", "Temperature-selected", "Interaction")) +
  theme_minimal(base_size = 13) +
  guides(fill = "none")+
  labs(x = "logFC (Hot.evo - Anc)",y = "logFC (Cold.evo - Anc)",
    title = "Differential Expression: Cold vs Hot Evolution") +
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5))
dev.off()
### continue with proportion analysis----
short.list <- rbind(short.list.hot[,c(1,2,3,6,7)], short.list.cold[,c(1,2,3,6,7)])
rm(short.list.cold, short.list.hot)
## to polarize the FC to be positibe on lab response, while the temp response *-1
short.list[logFC.evo < 0, `:=`(logFC.evo = -1*logFC.evo, logFC.hc = -1*logFC.hc)]
short.list[logFC.hc > 0, temp.prop := logFC.hc/(logFC.evo+logFC.hc)]
short.list[logFC.hc < 0, temp.prop := (-1*logFC.hc)/(logFC.evo-logFC.hc)]
short.list[, lab.prop := 1 - temp.prop]

plot.dt <- short.list[,.(gene, cge,sig, lab.prop, temp.prop)]
setorder(plot.dt, lab.prop)
plot.dt <- melt(plot.dt, id.vars = c("gene","cge","sig"), variable.name = "type", value.name = "proportion")
plot.dt[,gene.cge := paste(gene, cge, sep = ".")]
plot.dt[, type:= factor( type, levels = c("temp.prop", "lab.prop"))]
order_vec <- unique(plot.dt[type == "lab.prop"][order(-proportion), gene.cge])
plot.dt[ , gene.cge := factor(gene.cge, levels = order_vec) ]
plot.dt[,sig := factor(sig, levels = c("Lab.sig","Inter.sig","Temp.sig"))]

png("../Plot/GeneExpression.proportion.lab.temp.FC.png", width = 8, height = 6, units = "in", res = 600)
ggplot(data = plot.dt, aes(x = gene.cge, y = proportion, fill = type))+
  geom_col(width = 1.01) +
  facet_grid(~sig, scale = "free_x", space = "free_x", drop = T)+
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL,
       title = "Temp- vs. Lab-driven gene expression changes",
       fill = "Type") +
  scale_fill_manual(values = c("lab.prop"  = "purple",
                               "temp.prop" = "gold"),
                    name   = "Type",
                    labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x        = element_blank(),strip.clip    = "off", 
        panel.spacing.x = unit(0.2, "lines"),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0))
dev.off()


### indistinguishing between cges ----
cpm.mat <- cpm(y, normalized.lib.sizes = T, log = T, prior.count = 1)

hotcge.lab <- result[fdr.evo.hot<0.05 & fdr.hc.hot>0.05,gene]
hotcge.temp <- result[fdr.hc.hot<0.05 & fdr.evo.hot>0.05, gene]
hotcge.inter <- result[fdr.evo.hot<0.05 & fdr.hc.hot < 0.05, gene]

coldcge.lab <- result[fdr.evo.cold<0.05 & fdr.hc.cold>0.05,gene]
coldcge.temp <- result[fdr.hc.cold<0.05 & fdr.evo.cold>0.05, gene]
coldcge.inter <- result[fdr.evo.cold<0.05 & fdr.hc.cold < 0.05, gene]

hot.list <- list(lab = hotcge.lab,
            temp = hotcge.temp,
            inter = hotcge.inter)

cold.list <- list(lab = coldcge.lab,
             temp = coldcge.temp,
             inter = coldcge.inter)

### first hot----
mat <- cpm.mat[rownames(cpm.mat)%in%unlist(hot.list),15:26,drop=FALSE]
mat <- t(scale(t(mat),center=TRUE,scale=TRUE));mat[is.na(mat)]<-0
cat <- factor(sapply(rownames(mat), function(g) {
  paste(names(hot.list)[vapply(hot.list, function(v) g %in% v, logical(1))], collapse=";")
}), levels = c("lab", "temp", "inter"))
hc.r <- hclust(dist(mat,"euclidean"),"complete");hc.c <- hclust(dist(t(mat),"euclidean"),"complete")
mat <- mat[hc.r$order,hc.c$order];cat <- cat[hc.r$order]
dt <- as.data.table(mat);dt[,gene:=factor(rownames(mat),levels=rownames(mat))];dt[,category:=cat]
dt <- melt(dt,id.vars=c("gene","category"),variable.name="sample",value.name="z")
dt[,sample:=factor(sample,levels=colnames(mat))]
png("../Plot/DEG_category_general_plot/Heatmap.hotcgeDE.png", width = 8, height = 6, units =  "in", res = 400)
ggplot(dt,aes(sample,gene,fill=z))+
  geom_tile()+
  scale_fill_gradient2(low="steelblue",mid="white",high="tomato",midpoint=0)+
  facet_grid(category~.,scales="free_y",space="free_y",switch="y")+
  theme_minimal(base_size=12)+
  ggtitle("DE genes detected in the hotcge")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.clip    = "off", 
        panel.spacing.y=unit(0.2,"lines"),
        strip.background = element_rect(fill="white",colour="grey",linewidth=.6),  # ← add border
        strip.text.y     = element_text(size=10,face="bold")     )
dev.off()
### then cold----
mat <- cpm.mat[rownames(cpm.mat) %in% unlist(cold.list), 1:14, drop = FALSE]
mat <- t(scale(t(mat),center=TRUE,scale=TRUE));mat[is.na(mat)]<-0
cat <- factor(sapply(rownames(mat), function(g) {
  paste(names(cold.list)[vapply(cold.list, function(v) g %in% v, logical(1))], collapse=";")
}), levels = c("lab", "temp", "inter"))
hc.r <- hclust(dist(mat,"euclidean"),"complete");hc.c <- hclust(dist(t(mat),"euclidean"),"complete")
mat <- mat[hc.r$order,hc.c$order];cat <- cat[hc.r$order]
dt <- as.data.table(mat);dt[,gene:=factor(rownames(mat),levels=rownames(mat))];dt[,category:=cat]
dt <- melt(dt,id.vars=c("gene","category"),variable.name="sample",value.name="z")
dt[,sample:=factor(sample,levels=colnames(mat))]

png("../Plot/DEG_category_general_plot/Heatmap.coldcgeDE.png", width = 8, height = 6, units = "in", res = 400)
ggplot(dt,aes(sample,gene,fill=z))+
  geom_tile()+
  scale_fill_gradient2(low="steelblue",mid="white",high="tomato",midpoint=0)+
  facet_grid(category~.,scales="free_y",space="free_y",switch="y")+
  theme_minimal(base_size=12)+
  ggtitle("DE genes detected in the coldcge")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        strip.clip    = "off", 
        panel.spacing.y=unit(0.2,"lines"),
        strip.background = element_rect(fill="white",colour="grey",linewidth=.6),  # ← add border
        strip.text.y     = element_text(size=10,face="bold")     )
dev.off()
## GO enrichment on teh three categories of DE genes ----
library(topGO);library(org.Dm.eg.db);library(AnnotationDbi);library(patchwork)

gene2GO <- AnnotationDbi::select(org.Dm.eg.db,
                                 keys   = unique(result$gene),
                                 keytype= "FLYBASE",
                                 columns= "GO")
gene2GO <- gene2GO[!is.na(gene2GO$GO), ]
gene2GO <- split(gene2GO$GO, gene2GO$FLYBASE)

universe <- names(gene2GO)

run_topgo <- function(target, title){
  geneList <- factor(as.integer(universe %in% target)); names(geneList) <- universe
  obj <- new("topGOdata", ontology="BP", allGenes=geneList, nodeSize=10,
             annot=annFUN.gene2GO, gene2GO=gene2GO)
  tbl <- GenTable(obj, p = runTest(obj, algorithm="weight01", statistic="fisher"),
                  topNodes=500)
  tbl$GeneRatio <- as.numeric(tbl$Significant)/as.numeric(tbl$Annotated)
  tbl <- tbl[order(as.numeric(tbl$p)),][1:7,]
  tbl$Term <- factor(str_wrap(tbl$Term,35), levels = rev(str_wrap(tbl$Term,35)))
  ggplot(tbl,aes(GeneRatio,Term,size=Significant,colour=as.numeric(p)))+
    geom_point()+scale_colour_continuous(type="viridis",name="p")+
    labs(title=title,x="Gene ratio",y=NULL,size="Count")+
    theme_bw(base_size=14)+theme(panel.grid.major.y=element_blank())
}

plots <- list(
  run_topgo(hotcge.lab,  "hot lab genes"),
  run_topgo(hotcge.temp, "hot temp genes"),
  run_topgo(hotcge.inter,"hot inter genes"),
  run_topgo(coldcge.lab,  "cold lab genes"),
  run_topgo(coldcge.temp, "cold temp genes"),
  run_topgo(coldcge.inter,"cold inter genes")
)

png("../Plot/DEG_category_general_plot/GO.enrich.temp.lab.inter.png", width = 12, height = 7, units = "in", res = 400)
wrap_plots(plots,ncol=3)+
  plot_annotation(title="GO over-representation (topGO)")
dev.off()
png("../Plot/DEG_category_general_plot/GO.enrich.temp.lab.inter.png", width = 14, height = 7, units = "in", res = 600)
leg   <- cowplot::get_legend(plots[[4]] + theme(legend.position="right"))
plots <- lapply(plots, \(p) p + theme(legend.position="none"))
grid  <- plot_grid(plotlist=plots, ncol=3, align="hv")
final <- plot_grid(grid, leg, ncol=2, rel_widths=c(1,0.12))
print(final)
dev.off()

### old scripts below ----
hotcge.lab <- result[fdr.evo.hot<0.05 & fdr.hc.hot>0.05,gene]
hotcge.temp <- result[fdr.hc.hot<0.05 & fdr.evo.hot>0.05, gene]
hotcge.inter <- result[fdr.evo.hot<0.05 & fdr.hc.hot < 0.05, gene]

coldcge.lab <- result[fdr.evo.cold<0.05 & fdr.hc.cold>0.05,gene]
coldcge.temp <- result[fdr.hc.cold<0.05 & fdr.evo.cold>0.05, gene]
coldcge.inter <- result[fdr.evo.cold<0.05 & fdr.hc.cold < 0.05, gene]

type.list <- list(hotcge.lab = hotcge.lab,
                  hotcge.temp=hotcge.temp,
                  hotcge.inter=hotcge.inter,
                  coldcge.lab=coldcge.lab,
                  coldcge.temp=coldcge.temp,
                  coldcge.inter=coldcge.inter)

cpm.mat <- cpm(y, normalized.lib.sizes = T, log = T, prior.count = 1)
for (type in names(type.list)) {
  dir.create(path = paste0("../Plot/DEG_category_plot/", type,"/"))
  temp.genelist <- type.list[[type]]
  temp.dat <- cpm.mat[rownames(cpm.mat) %in% temp.genelist,]
  temp.dat <- setDT(as.data.frame(temp.dat), keep.rownames = T)
  long_data <- melt(temp.dat, id.vars = "rn", variable.name = "sample", value.name = "expression")
  long_data[, c("environment", "replicate") := tstrsplit(sample, "_")[1:2]]
  long_data[, population := substr(replicate, 1, 1)]
  long_data[, replicate := substr(replicate, 2, 2)]
  for (gene in unique(long_data$rn)) {
    p<- ggplot(long_data[rn == gene,], aes(x = population, y = expression, color = population)) +
      facet_wrap(~environment)+
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(size  =3,width = 0.2, shape = factor(replicate))+
      labs(title = paste0("Gene ", gene, ", categorized as ", type), x = "Population", y = "logCPM (Expression Intensity)",shape = "Replicate") +
      scale_color_manual(values = c("forestgreen","steelblue","maroon"))+
      theme_minimal()
    png(filename = paste0("../Plot/DEG_category_plot/", type,"/",gene,".png"), width = 7, height = 5, units = "in", res = 350)
    print(p)
    dev.off()
  }
}

for (type in names(type.list)) {
  temp.genelist <- type.list[[type]]
  temp.dat <- cpm.mat[rownames(cpm.mat) %in% temp.genelist, ]
  annotation_col <- data.frame(Group = group)
  rownames(annotation_col) <- colnames(temp.dat)
  annotation_colors <- list(
    Group = c(
      "coldcgeA" = "forestgreen",
      "coldcgeC" = "steelblue",
      "coldcgeH" = "maroon",
      "hotcgeA"  = "lightgreen",
      "hotcgeC"  = "lightblue",
      "hotcgeH"  = "lightcoral"  
    )
  )
  # Save the plot as a PNG file
  png(filename = paste0("../plot/DEG_category_general_plot/", type, "_pheatmap.png"),
      width = 8, height = 6, units = "in", res = 600)
  
  pheatmap(temp.dat,scale = "row",cluster_cols = T,cutree_cols = 4,
    cluster_rows = TRUE,show_rownames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    border_color = NA,main = paste0("Heatmap for category: ", type))
  dev.off()
}











### PC loading GSEA ----
library(clusterProfiler)
library(org.Dm.eg.db)
loadings <- pca$rotation  # genes x PCs
pc_list <- list(PC1 = loadings[, "PC1"],
                PC2 = loadings[, "PC2"],
                PC3 = loadings[, "PC3"])
gsea_results <- list()

for (pc_name in names(pc_list)) {
  message("Running GSEA for ", pc_name, "...")
  # Sort loadings in decreasing order
  pc_loadings <- pc_list[[pc_name]]
  pc_loadings <- sort(pc_loadings, decreasing = TRUE)
  flybase_ids <- names(pc_loadings)
  # This returns a data frame with columns FLYBASE and ENTREZID.
  flybase_to_entrez <- bitr(
    flybase_ids, 
    fromType = "FLYBASE",     # the input ID type
    toType   = "ENTREZID",    # the ID type we want
    OrgDb    = org.Dm.eg.db)
  # Some FlyBase IDs might not map to Entrez IDs, so subset to only mapped IDs
  matched_idx <- flybase_to_entrez$FLYBASE %in% flybase_ids
  flybase_to_entrez <- flybase_to_entrez[matched_idx, ]
  pc_loadings <- pc_loadings[flybase_to_entrez$FLYBASE]
  names(pc_loadings) <- flybase_to_entrez$ENTREZID
  gsea <- gseKEGG(
    geneList     = pc_loadings,
    organism     = "dme",  
    keyType      = "ncbi-geneid",
    minGSSize    = 20,  
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    exponent = 1.1,
    verbose      = FALSE)
  gsea_results[[pc_name]] <- gsea
}

p1 <- dotplot(gsea_results[["PC1"]], showCategory = 15) + ggtitle("GSEA for PC1")+ xlim(0,1)+theme(legend.position = "right")
p2 <- dotplot(gsea_results[["PC2"]], showCategory = 15) + ggtitle("GSEA for PC2")+ xlim(0,1)+ theme(legend.position = "right")
p3 <- dotplot(gsea_results[["PC3"]], showCategory = 15) + ggtitle("GSEA for PC3")+ xlim(0,1)+ theme(legend.position = "right")

plots_combined <- plot_grid(  p1,p2,p3, ncol = 3,rel_widths = c(1,1,1),  align = "v")
png(filename = "../Plot/PC.GSEA.png", width = 18, height = 12, units = "in",res = 600)
print(plots_combined)
dev.off()


### find interesting Pathways and highlight the DE genes ----
ModelDesign <- model.matrix(~0+group, list(levels(group)))
colnames(ModelDesign) <- sub("^group", "", colnames(ModelDesign))

y <- estimateDisp(y, ModelDesign, robust=TRUE)
plotBCV(y, xlim = c(2,15))
fit <- glmFit(y, ModelDesign)
mycontrast <- makeContrasts(hotcge.hot = hotcgeH - hotcgeA,
                            hotcge.cold = hotcgeC - hotcgeA,
                            coldcge.hot = coldcgeH - coldcgeA,
                            coldcge.cold = coldcgeC - coldcgeA,
                            levels = ModelDesign)
hotcge.hot <- glmLRT(fit, contrast = mycontrast[,"hotcge.hot"])$table
hotcge.hot$p.adjust <- p.adjust(hotcge.hot$PValue,method = "fdr")
setDT(hotcge.hot, keep.rownames = T)

result <- copy(hotcge.hot[,c(1,3)])
setnames(result, 1,"gene")

for (contr in colnames(mycontrast)) {
  print(contr)
  temp <- glmLRT(fit, contrast = mycontrast[,contr])$table
  temp$p.adjust <- p.adjust(temp$PValue,method = "fdr")
  temp_df <- data.frame(logFC = temp$logFC,fdr = temp$p.adjust)
  colnames(temp_df) <- c(paste0("logFC.", contr), paste0("fdr.", contr))
  result <- cbind(result, temp_df)
}

hotcge.hot <- result[fdr.hotcge.hot<0.05 ,gene]
hotcge.cold <- result[fdr.hotcge.cold<0.05 , gene]
coldcge.hot <- result[fdr.coldcge.hot<0.05,gene]
coldcge.cold <- result[fdr.coldcge.cold<0.05, gene]

type.list <- list(hotcge.hot = hotcge.hot,
                  hotcge.cold=hotcge.cold,
                  coldcge.hot=coldcge.hot,
                  coldcge.cold=coldcge.cold)

library(pathview)
flybase_to_entrez <- bitr(
  result$gene,
  fromType = "FLYBASE",
  toType   = "ENTREZID",
  OrgDb    = org.Dm.eg.db
)
result_merged <- merge(
  x = result,
  y = flybase_to_entrez,
  by.x = "gene",
  by.y = "FLYBASE",
  all.x = TRUE
)
# Remove rows that failed to map (ENTREZID=NA)
result_merged <- subset(result_merged, !is.na(ENTREZID))

logfc.data <- as.matrix(result_merged[,c(3,5,7,9)]); rownames(logfc.data) <- result_merged$ENTREZID
# A helper function to run pathview
run_pathview <- function(gene_vector, pathway = pathway_id, suffix) {
  pathview(
    gene.data   = gene_vector,
    pathway.id  = pathway,
    species     = "dme",        # for D. melanogaster
    kegg.native = F,
    out.suffix  = suffix,
    limit       = list(gene = c(-1, 1)) # color scale from -3 to +3
  )
}

for (pathway_id in gsea_results[["PC1"]]@result$ID) {
  run_pathview(logfc.data, pathway = pathway_id, "PC1")
} 
for (pathway_id in gsea_results[["PC2"]]@result$ID) {
  run_pathview(logfc.data, pathway = pathway_id, "PC2")
} 
for (pathway_id in gsea_results[["PC3"]]@result$ID) {
  run_pathview(logfc.data, pathway = pathway_id, "PC3")
} 
run_pathview(logfc.data, pathway = "dme04081", "Hormone")
pv.hotCold  <- run_pathview(geneVector_hotCold,  suffix="hotcge.cold")
pv.coldHot  <- run_pathview(geneVector_coldHot,  suffix="coldcge.hot")
pv.coldCold <- run_pathview(geneVector_coldCold, suffix="coldcge.cold")


gseaplot(gsea_results[["PC2"]], geneSetID = "dme04141")




# 
# 
# ### Double-checking, are these enrichment legit? do we get the same from DE genes? ----
# ModelDesign <- model.matrix(~0+group, list(levels(group)))
# colnames(ModelDesign) <- sub("^group", "", colnames(ModelDesign))
# 
# y <- estimateDisp(y, ModelDesign, robust=TRUE)
# plotBCV(y, xlim = c(2,15))
# fit <- glmFit(y, ModelDesign)
# mycontrast <- makeContrasts(lab.effect = (hotcgeH + hotcgeC + coldcgeH + coldcgeC)/4 - (hotcgeA + coldcgeA)/2,
#                             temp.effect = (hotcgeH + coldcgeH)/2 - (hotcgeC + coldcgeC)/2, 
#                             cge.effect = (hotcgeH + hotcgeC + hotcgeA)/3 - (coldcgeH + coldcgeC + coldcgeA)/3,
#                             levels = ModelDesign)
# # Step 2: Loop through contrasts and run GSEA
# contrast_names <- colnames(mycontrast)
# gsea_results_de <- list()
# for (contrast_name in contrast_names) {
#   message("Running DE-GSEA for ", contrast_name, "...")
#   de_table <- glmLRT(fit, contrast = mycontrast[, contrast_name])$table
#   de_table$p.adjust <- p.adjust(de_table$PValue, method = "fdr")
#   de_table$effect <- de_table$logFC 
#   # Named vector for GSEA
#   gene_vector <- de_table$effect
#   names(gene_vector) <- rownames(de_table)
#   gene_vector <- sort(gene_vector, decreasing = TRUE)
#   plot(gene_vector)
#   # Map FlyBase to Entrez
#   flybase_to_entrez <- bitr(names(gene_vector), fromType = "FLYBASE", toType = "ENTREZID", OrgDb = org.Dm.eg.db)
#   gene_vector <- gene_vector[flybase_to_entrez$FLYBASE]
#   names(gene_vector) <- flybase_to_entrez$ENTREZID
#   # Run GSEA
#   gsea_de <- gseKEGG(
#     geneList     = gene_vector,
#     organism     = "dme",
#     keyType      = "ncbi-geneid",
#     minGSSize    = 20,
#     maxGSSize    = 500,
#     pvalueCutoff = 0.05,
#     pAdjustMethod= "BH",
#     exponent     = 1,
#     verbose      = FALSE)
#   # Store result
#   gsea_results_de[[contrast_name]] <- gsea_de
# }
# 
# p1 <- dotplot(gsea_results_de[["cge.effect"]], showCategory = 15) + ggtitle("GSEA for CGE effect")+ xlim(0,1)+theme(legend.position = "right")
# p2 <- dotplot(gsea_results_de[["lab.effect"]], showCategory = 15) + ggtitle("GSEA for Lab effect")+ xlim(0,1)+ theme(legend.position = "right")
# p3 <- dotplot(gsea_results_de[["temp.effect"]], showCategory = 15) + ggtitle("GSEA for Temp effect")+ xlim(0,1)+ theme(legend.position = "right")
# 
# plots_combined <- plot_grid(p1,p2,p3, ncol = 3,rel_widths = c(1,1,1),  align = "v")
# png(filename = "../Plot/DE.GSEA.png", width = 18, height = 12, units = "in",res = 600)
# print(plots_combined)
# dev.off()
# 
# 
# gseaplot(gsea_results_de[["cge.effect"]], geneSetID = "dme03010")
# 
# 
# ### Double-checking, with originally interesting DE genes ----
# ModelDesign <- model.matrix(~0+group, list(levels(group)))
# colnames(ModelDesign) <- sub("^group", "", colnames(ModelDesign))
# 
# y <- estimateDisp(y, ModelDesign, robust=TRUE)
# plotBCV(y, xlim = c(2,15))
# fit <- glmFit(y, ModelDesign)
# mycontrast <- makeContrasts(hotcge.common = (hotcgeH + hotcgeC )/2 - hotcgeA,
#                             hotcge.diff = hotcgeH  - hotcgeC,
#                             coldcge.common = (coldcgeH + coldcgeC )/2 - coldcgeA,
#                             coldcge.diff = coldcgeH  - coldcgeC,
#                             levels = ModelDesign)
# 
# # Step 2: Loop through contrasts and run GSEA
# contrast_names <- colnames(mycontrast)
# gsea_results_de <- list()
# for (contrast_name in contrast_names) {
#   message("Running DE-GSEA for ", contrast_name, "...")
#   de_table <- glmLRT(fit, contrast = mycontrast[, contrast_name])$table
#   de_table$p.adjust <- p.adjust(de_table$PValue, method = "fdr")
#   de_table$effect <- de_table$logFC 
#   # Named vector for GSEA
#   gene_vector <- de_table$effect
#   names(gene_vector) <- rownames(de_table)
#   gene_vector <- sort(gene_vector, decreasing = TRUE)
#   plot(gene_vector)
#   # Map FlyBase to Entrez
#   flybase_to_entrez <- bitr(names(gene_vector), fromType = "FLYBASE", toType = "ENTREZID", OrgDb = org.Dm.eg.db)
#   gene_vector <- gene_vector[flybase_to_entrez$FLYBASE]
#   names(gene_vector) <- flybase_to_entrez$ENTREZID
#   # Run GSEA
#   gsea_de <- gseKEGG(
#     geneList     = gene_vector,
#     organism     = "dme",
#     keyType      = "ncbi-geneid",
#     minGSSize    = 20,
#     maxGSSize    = 500,
#     pvalueCutoff = 0.05,
#     pAdjustMethod= "BH",
#     exponent     = 1,
#     verbose      = FALSE)
#   # Store result
#   gsea_results_de[[contrast_name]] <- gsea_de
# }
# 
# p1 <- dotplot(gsea_results_de[["hotcge.common"]], showCategory = 15) + ggtitle("GSEA for hotcge lab effect")+ xlim(0,1)+theme(legend.position = "right")
# p2 <- dotplot(gsea_results_de[["hotcge.diff"]], showCategory = 15) + ggtitle("GSEA for hotcge temp effect")+ xlim(0,1)+ theme(legend.position = "right")
# p3 <- dotplot(gsea_results_de[["coldcge.common"]], showCategory = 15) + ggtitle("GSEA for coldcge lab effect")+ xlim(0,1)+ theme(legend.position = "right")
# p4 <- dotplot(gsea_results_de[["coldcge.diff"]], showCategory = 15) + ggtitle("GSEA for coldcge temp effect")+ xlim(0,1)+ theme(legend.position = "right")
# 
# plots_combined <- plot_grid(p1,p2,p3,p4, ncol = 2,  align = "v")
# png(filename = "../Plot/DE.categories.GSEA.png", width = 18, height = 18, units = "in",res = 600)
# print(plots_combined)
# dev.off()
# 
# 
# gseaplot(gsea_results_de[["cge.effect"]], geneSetID = "dme03010")
