rm(list=ls())
gc()

library(data.table)
setDTthreads(percent = 80) 
library(edgeR)
library(ggplot2);library(ggrepel);library(cowplot)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(ExpressionNormalizationWorkflow);library(scales)
library(biomaRt)
setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Script/")

## 1. Read in data ----
counts <- read.csv("../Data/RNAseq_count_table_28182010_AncColdHot.csv",stringsAsFactors = F,row.names = 1)[,-c(6,13)]
colnames(counts)
setcolorder(counts, neworder = sort(colnames(counts)))
colnames(counts)
setDT(counts, keep.rownames = T)
colnames(counts)[1] <- "gene"

### Update the gene IDs to the latest
latest_FB_ID <- read.delim("./fbgn_annotation_ID_current.tsv",header = T,stringsAsFactors = F,sep="\t")
for(i in 1:dim(latest_FB_ID)[1]){
  counts[gene%in%strsplit2(latest_FB_ID[i,4],","), gene:= latest_FB_ID[i,3]]
}
counts <- counts[gene%in%latest_FB_ID$primary_FBgn.,]# keep only annotated genes, for down stream analysis already

counts <- counts[, lapply(.SD, sum), by = gene] # some genes are later merged in flybase, so we also merge their counts here
setDF(counts)
rownames(counts) <- counts$gene; counts <- counts[,2:ncol(counts)]
counts <- counts[apply(cpm(counts),1,function(x) sum(x<1)<=5),]#filtered for lowly expressed genes

### define some variables
cge <- sub("_.*", "", colnames(counts))
evo <- sub("^[^_]+_([A-Z]).*", "\\1", colnames(counts))
group <- paste0(cge,evo)
replicate <- sub("^[^_]+_[A-Z]([0-9]).*", "\\1", colnames(counts))
replicate[which(evo == "A")] <- "0"
bio.sample <- sub("^(([^_]+_[^_]+))_.*", "\\1", colnames(counts))
y <- DGEList(counts=counts, group = group)
y <- calcNormFactors(y, method = "TMM")

## 2. PCA ----
pca <- prcomp(t(cpm(y, log = T)), center = T, scale. = F)
coord <- as.data.frame(pca$x)
coord$cge <- cge
coord$evo <- evo

pop_colors <- c("A" = "forestgreen", "C" = "steelblue", "H" = "maroon")
pop_labels <- c("A" = "Ancestral", "C" = "Cold-evolved", "H" = "Hot-evolved")

p1 <- ggplot(coord, aes(x = PC1, y = PC2, color = evo)) +
  geom_point(shape = 16, size = 2) +
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
  geom_point(shape = 16, size = 2) +
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
p1_no_legend <- p1 + theme(legend.position = "none")
p2_no_legend <- p2 + theme(legend.position = "none") # remove indicidual legend, cuz they ugly
legend <- get_legend(p1) # take one legend out so I can print it just once

plots_combined <- plot_grid(p1_no_legend, p2_no_legend, ncol = 2, align = "hv")
combined_with_legend <- plot_grid(plots_combined, legend, rel_widths = c(1, 0.15))
final_plot <- plot_grid(
  ggdraw() + draw_label("PCA on gene expression data", x = 0.5),
  combined_with_legend, ncol = 1 , rel_heights = c(0.1, 1))

png("../Plot/PCA.png", width = 10, height = 5, unit = "in", res = 600)
print(final_plot)
dev.off()


## 3. DE analysis ----
ModelDesign <- model.matrix( ~ cge*evo)
colnames(ModelDesign)
y <- DGEList(counts=counts)
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, ModelDesign)
plotBCV(y)
fit <- glmFit(y, ModelDesign)

### 3.1 detect & remove genes with significant cge:evo effect----
interaction_terms <- grep("cgehotcge:evo*", colnames(ModelDesign))
lrt_interaction <- glmLRT(fit, coef=interaction_terms)
interaction_fdr <- p.adjust(lrt_interaction$table$PValue, method="fdr")
sum(interaction_fdr < 0.2)

significant_interaction <- interaction_fdr < 0.2

no_interaction_genes <- rownames(y)[!significant_interaction]
y_filtered <- y[no_interaction_genes, , keep.lib.sizes=FALSE]
fit_filtered <- glmFit(y_filtered, ModelDesign)

rm(lrt_interaction,significant_interaction, y, fit,  interaction_fdr, interaction_terms)

### 3.2 hot evolution ----
lrt_hot <- glmLRT(fit_filtered, coef=4)
lrt_hot <-lrt_hot$table
lrt_hot$fdr <- p.adjust(lrt_hot$PValue, method="fdr")
setDT(lrt_hot, keep.rownames = T)
setkey(lrt_hot, rn)
lrt_hot[fdr<0.05, .N]
lrt_hot <- lrt_hot[,.(rn,logFC, fdr)]
setnames(lrt_hot, 1:3, c("gene","logFC.hot","fdr.hot"))

result <- copy(lrt_hot)
### 3.3 cold evolution ----
lrt_cold <- glmLRT(fit_filtered, coef=3) 
lrt_cold <-lrt_cold$table
lrt_cold$fdr <- p.adjust(lrt_cold$PValue, method="fdr")
setDT(lrt_cold, keep.rownames = T)
setkey(lrt_cold, rn)
lrt_cold[fdr<0.05, .N]
lrt_cold <- lrt_cold[,.(rn,logFC, fdr)]
setnames(lrt_cold, 1:3, c("gene","logFC.cold","fdr.cold"))

result <- merge(result, lrt_cold )

### 3.4 lab adaptive ----
lrt_lab <- glmLRT(fit_filtered, contrast = c(0,0,0.5,0.5,0,0)) 
lrt_lab <-lrt_lab$table
lrt_lab$fdr <- p.adjust(lrt_lab$PValue, method="fdr")
setDT(lrt_lab, keep.rownames = T)
setkey(lrt_lab, rn)
lrt_lab[fdr<0.05, .N]
lrt_lab <- lrt_lab[,.(rn,logFC, fdr)]
setnames(lrt_lab, 1:3, c("gene","logFC.lab","fdr.lab"))
result <- merge(result, lrt_lab )

### 3.5 temp adaptive ----
lrt_temp <- glmLRT(fit_filtered, contrast = c(0,0,-1,1,0,0)) 
lrt_temp <-lrt_temp$table
lrt_temp$fdr <- p.adjust(lrt_temp$PValue, method="fdr")
setDT(lrt_temp, keep.rownames = T)
setkey(lrt_temp, rn)
lrt_temp[fdr<0.05, .N]
lrt_temp <- lrt_temp[,.(rn,logFC, fdr)]
setnames(lrt_temp, 1:3, c("gene","logFC.temp","fdr.temp"))
result <- merge(result, lrt_temp )

### 3.6 Env effect ----
lrt_env <- glmLRT(fit_filtered, contrast = c(0,1,0,0,0,0)) 
lrt_env <-lrt_env$table
lrt_env$fdr <- p.adjust(lrt_env$PValue, method="fdr")
setDT(lrt_env, keep.rownames = T)
setkey(lrt_env, rn)
lrt_env[fdr<0.05, .N]
lrt_env <- lrt_env[,.(rn,logFC, fdr)]
setnames(lrt_env, 1:3, c("gene","logFC.env","fdr.env"))
result <- merge(result, lrt_env )

## 4. VennDiagram ----
venn.list <- list(
  hot  = result[fdr.hot < 0.05, gene],
  cold = result[fdr.cold < 0.05, gene],
  lab  = result[fdr.lab < 0.05, gene],
  temp = result[fdr.temp < 0.05, gene]
)
library(ggVennDiagram)
png(filename = "../Plot/VennDiagram.hot.cold.lab.temp.png", width = 8, height = 6, units = "in", res = 400)
ggVennDiagram(venn.list, label_alpha = 0.2, set_size =  6) +
  scale_fill_gradient(low = "white", high = "chartreuse4") +
  theme_void() +
  theme(legend.position = "none")
dev.off()

rm(venn.list)
## 5. Define Lab Temp DEGs ----
short.list <- copy(result)
short.list[fdr.lab < 0.05 & fdr.temp >= 0.05, sig := "Lab DEGs"]
short.list[(fdr.lab >= 0.05 & fdr.temp < 0.05), sig := "Temp DEGs"]
short.list[fdr.lab < 0.05 & fdr.temp < 0.05 , sig := "Lab X Temp DEGs"]
saveRDS(short.list, file = "./short.list.RDS")
## 6. Detailed expression plot ----
cpm.mat <- cpm(y_filtered, normalized.lib.sizes = T, log = T, prior.count = 1)
saveRDS(cpm.mat, file = "./cpm.mat.RDS")
DE.list <- list(lab = short.list[sig == "Lab DEGs",gene],
                temp = short.list[sig == "Temp DEGs",gene],
                inter = short.list[sig == "Lab X Temp DEGs",gene])
fwrite(x = list(rownames(cpm(y_filtered))), file = "./Gene.universe.txt", quote = F, sep = "\t")
fwrite(x = list(DE.list$lab), file = "./Lab.gene.txt", quote = F, sep = "\t")
fwrite(x = list(DE.list$temp), file = "./Temp.gene.txt", quote = F, sep = "\t")
fwrite(x = list(DE.list$inter), file = "./Inter.gene.txt", quote = F, sep = "\t")

# for (type in names(DE.list)) {
#   dir.create(path = paste0("../Plot/DEG_category_plot/", type,"/"))
#   temp.genelist <- DE.list[[type]]
#   temp.dat <- cpm.mat[rownames(cpm.mat) %in% temp.genelist,]
#   temp.dat <- setDT(as.data.frame(temp.dat), keep.rownames = T)
#   long_data <- melt(temp.dat, id.vars = "rn", variable.name = "sample", value.name = "expression")
#   long_data[, c("environment", "replicate") := tstrsplit(sample, "_")[1:2]]
#   long_data[, population := substr(replicate, 1, 1)]
#   long_data[, replicate := substr(replicate, 2, 2)]
#   for (gene in unique(long_data$rn)) {
#     p<- ggplot(long_data[rn == gene,], aes(x = population, y = expression, color = population)) +
#       facet_wrap(~environment)+
#       geom_boxplot(outlier.shape = NA)+
#       labs(title = paste0("Gene ", gene, ", categorized as ", type), x = "Population", y = "logCPM (Expression Intensity)") +
#       scale_color_manual(values = c("forestgreen","steelblue","maroon"))+
#       theme_minimal()
#     png(filename = paste0("../Plot/DEG_category_plot/", type,"/",gene,".png"), width = 7, height = 5, units = "in", res = 350)
#     print(p)
#     dev.off()
#   }
# }
# dev.off()

rm(lrt_cold, lrt_env, lrt_hot, lrt_lab, lrt_temp)

## 7. LogFC plot ----
logfc.dt <- copy(short.list)
logfc.dt[, sig := factor(sig, levels = c("Lab DEGs", "Temp DEGs", "Lab X Temp DEGs"))]

png("../Plot/DE.LogFC.png", width = 10, height = 5, units = "in", res = 450)

ggplot(logfc.dt, aes(x = logFC.hot, y = logFC.cold)) +
  geom_point(data = logfc.dt[is.na(sig)], color = "grey", size = 1) +
  geom_point(data = logfc.dt[!is.na(sig)], aes(color = sig), size = 1.8, alpha = 0.7) +
  scale_color_manual(
    values = c("Lab DEGs" = "purple","Temp DEGs" = "yellow","Lab X Temp DEGs" = "orange"),
    name = "Gene Category",labels = c(paste0("Lab-selected (",logfc.dt[sig == "Lab DEGs", .N],")"), 
                                      paste0("Temp-selected (",logfc.dt[sig == "Temp DEGs", .N],")"), 
                                      paste0("Lab X Temp-selected (",logfc.dt[sig == "Lab X Temp DEGs", .N],")"))) +
  theme_minimal(base_size = 13) +
  guides(fill = "none")+
  ylim(-2,3)+
  xlim(-2,3)+
  labs(x = "logFC (Hot.evo - Anc)",y = "logFC (Cold.evo - Anc)",
       title = "Differential Expression: Cold vs Hot Evolution") +
  theme(legend.position = "right",plot.title = element_text(hjust = 0.5))
dev.off()


##  8. Fig. 5 ----
library(ggplot2)
library(cowplot)
theme_manuscript <- theme_minimal(base_size = 18, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    legend.title = element_text(face = "bold", size = 22)
  )

## Panels A and B: remove legends from both plots
pA <- p1 + 
  theme_manuscript + 
  theme(legend.position = "none")   # PC1 vs PC2

pB <- p2 + 
  theme_manuscript + 
  theme(legend.position = "none")   # PC3 vs PC2

## Extract a shared legend for A and B, positioned at the bottom
legend_AB <- get_legend(
  p1 +
    theme_manuscript +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    ) +
    guides(color = guide_legend(nrow = 1))
)

## Build panel C
pC <- ggplot(logfc.dt, aes(x = logFC.hot, y = logFC.cold)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey", linetype = "dashed") +
  geom_abline(slope = -1, intercept = 0, colour = "grey", linetype = "dashed") +
  geom_point(data = logfc.dt[is.na(sig)], color = "grey", size = 1) +
  geom_point(data = logfc.dt[!is.na(sig)], aes(color = sig), size = 1.8, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Lab DEGs" = "purple",
      "Temp DEGs" = "yellow",
      "Lab X Temp DEGs" = "orange"
    ),
    name = "Gene Category",
    labels = c(
      paste0("Lab adaptive (", logfc.dt[sig == "Lab DEGs", .N], ")"),
      paste0("Temp adaptive (", logfc.dt[sig == "Temp DEGs", .N], ")"),
      paste0("Lab X Temp adaptive (", logfc.dt[sig == "Lab X Temp DEGs", .N], ")")
    )
  ) +
  xlim(-2, 4) +
  ylim(-2, 4) +
  labs(x = "logFC (Hot.evo - Anc)", y = "logFC (Cold.evo - Anc)") +
  theme_manuscript +
  theme(
    legend.position = c(1, 0.2),
    legend.justification = c("right", "bottom")
  )

## Combine A and B first
panel_AB <- plot_grid(
  pA, pB,
  nrow = 1,
  labels = c("A", "B"),
  label_size = 20,
  label_fontface = "bold",
  align = "h",
  axis = "tb"
)

## Put shared legend underneath A and B
panel_AB_with_legend <- plot_grid(
  panel_AB, legend_AB,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

## Add label to panel C separately
panel_C <- ggdraw(pC) +
  draw_plot_label(
    label = "C",
    x = 0.02, y = 0.98,
    hjust = 0, vjust = 1,
    size = 20,
    fontface = "bold"
  )

## Final assembly: (A+B+shared legend) on the left, C on the right
final_fig <- plot_grid(
  panel_AB_with_legend, panel_C,
  nrow = 1,
  rel_widths = c(1.7, 1),
  align = "h",
  axis = "tb"
)

ggsave(
  "../Plot/Figure5_PCA_DE_ABC.png",
  final_fig,
  width = 18,
  height = 7,
  units = "in",
  dpi = 600
)

## 9. Purine synthesis pathway ----

library(org.Dm.eg.db)
library(AnnotationDbi)
library(KEGGREST)

enzymes_list <- c(
  "2.4.2.14", "6.3.4.13", "2.1.2.2", "6.3.5.3", "6.3.3.1",
  "4.1.1.21", "6.3.2.6", "4.3.2.2", "3.5.4.10"
)

ec_mapping <- as.data.table(
  AnnotationDbi::select(
    org.Dm.eg.db,
    keys = enzymes_list,
    keytype = "ENZYME",
    columns = c("FLYBASE", "SYMBOL")
  )
)
# ec_mapping <- ec_mapping[SYMBOL != "Prat",]

setnames(ec_mapping, c("ENZYME", "FLYBASE", "SYMBOL"), c("ec", "gene", "symbol"))
ec_mapping <- unique(ec_mapping[!is.na(gene)])

gene_metadata <- ec_mapping[, .(
  ec_label = paste(sort(unique(ec)), collapse = ", "),
  symbol = symbol[1]
), by = gene]

gene_metadata[, display_label := paste0(symbol, " (", ec_label, ")")]

setDT(short.list)

purine.meta <- merge(short.list, gene_metadata, by = "gene")

## here we filter out Prat, the retro copy of Prat2 because it mostly is expressed in 
## preadult stages, our RNAseq data is from adult male flies so it is not so relevant.
## It is also lowly expressed and not significantly different in any comparison. 
purine.meta <- purine.meta[symbol != "Prat",]

label_levels <- sort(unique(purine.meta$display_label))
purine.meta[, x := match(display_label, label_levels)]

purine.ann <- purine.meta[, .(
  x = x[1],
  y = max(logFC.hot, logFC.cold, na.rm = TRUE) + 0.035,
  sig_label = fcase(
    fdr.temp < 0.001, "***",
    fdr.temp < 0.01,  "**",
    fdr.temp < 0.05,  "*",
    default = "n.s."
  )
), by = display_label]

purine.ann[, `:=`(
  y_tick = y - 0.015,
  y_text = y + 0.015
)]

ymax <- max(0.5, purine.ann$y_text, na.rm = TRUE) + 0.03

purine.meta <- melt(
  purine.meta,
  id.vars = c("gene", "display_label", "x"),
  measure.vars = c("logFC.hot", "logFC.cold"),
  variable.name = "condition",
  value.name = "logFC"
)

purine.meta[, condition := fifelse(condition == "logFC.hot", "Hot-evolved", "Cold-evolved")]

png("../Plot/Purine_gene_expression.png",
    width = 6, height = 12, units = "in", res = 450)

ggplot(purine.meta, aes(x = x, y = logFC, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_segment(data = purine.ann,
               aes(x = x - 0.35, xend = x + 0.35, y = y, yend = y),
               inherit.aes = FALSE) +
  geom_segment(data = purine.ann,
               aes(x = x - 0.35, xend = x - 0.35, y = y, yend = y_tick),
               inherit.aes = FALSE) +
  geom_segment(data = purine.ann,
               aes(x = x + 0.35, xend = x + 0.35, y = y, yend = y_tick),
               inherit.aes = FALSE) +
  geom_text(data = purine.ann,
            aes(x = x, y = y_text, label = sig_label),
            size = 6, inherit.aes = FALSE) +
  scale_x_continuous(
    breaks = seq_along(label_levels),
    labels = label_levels,
    expand = expansion(add = 0.6)
  ) +
  coord_cartesian(ylim = c(-0.2, ymax)) +
  labs(
    y = "logFC against ancestral state",
    x = "Gene (EC Number)",
    fill = "Population"
  ) +
  scale_fill_manual(values = c(
    "Hot-evolved" = "maroon",
    "Cold-evolved" = "steelblue"
  )) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position = "bottom",
    text = element_text(size = 20)
  )

dev.off()

## 9.1 plot out the raw expression value of Prat and Prat2 to found the removal of Prat2----

library(data.table)
library(ggplot2)

genes <- c(Prat = "FBgn0004901", Prat2 = "FBgn0041194",
           Pfas = "FBgn0000052", Gart = "FBgn0000053",
           Paics = "FBgn0020513", Adsl = "FBgn0038467",
           CG11089 = "FBgn0039241")

dt <- as.data.table(as.table(cpm.mat[genes, ]))
setnames(dt, c("gene_id", "sample", "expr"))
dt[, gene := factor(names(genes)[match(gene_id, genes)])]
dt[, cge := factor(fifelse(grepl("^coldcge", sample), "Cold CGE", "Hot CGE"), levels = c("Cold CGE", "Hot CGE"))]
dt[, pop := substr(sub(".*_", "", sample), 1, 1)]
dt[, Evo := factor(fifelse(pop == "A", "Ancestral", fifelse(pop == "C", "Cold-evolved", "Hot-evolved")), levels = c("Ancestral", "Cold-evolved", "Hot-evolved"))]

evo_cols <- c("Ancestral" = "forestgreen", "Cold-evolved" = "steelblue", "Hot-evolved" = "maroon")

pd <- position_dodge(width = 0.55)

p <- ggplot(dt, aes(gene, expr, color = Evo)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55), size = 2.7, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", position = pd, size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = pd, width = 0.18, linewidth = 0.6) +
  facet_wrap(~ cge, nrow = 1) +
  scale_color_manual(values = evo_cols) +
  labs(x = NULL, y = "Expression (Counts per million)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom", legend.title = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"))

p

ggsave("../Plot/SuppFig_purine_syn_expression.png", p, width = 12, height = 6, dpi = 600)
## 10. Purine salvage pathway ----

library(org.Dm.eg.db)
library(AnnotationDbi)

salvage_genes <- c(
  "FBgn0035348",
  "FBgn0000109",
  "FBgn0037661",
  "FBgn0036337",
  "FBgn0034898",
  "FBgn0052626"
)

gene_metadata <- as.data.table(
  AnnotationDbi::select(
    org.Dm.eg.db,
    keys = salvage_genes,
    keytype = "FLYBASE",
    columns = "SYMBOL"
  )
)

setnames(gene_metadata, c("FLYBASE", "SYMBOL"), c("gene", "symbol"))
gene_metadata <- unique(gene_metadata)
gene_metadata[, order_id := match(gene, salvage_genes)]
gene_metadata[, display_label := paste0(symbol, " (", gene, ")")]

setDT(short.list)

purine.salvage <- merge(short.list, gene_metadata, by = "gene")

label_levels <- unique(purine.salvage[order(order_id), display_label])
purine.salvage[, x := match(display_label, label_levels)]

purine.salvage.ann <- purine.salvage[, .(
  x = x[1],
  y = max(0, logFC.hot, logFC.cold, na.rm = TRUE) + 0.035,
  sig_label = fcase(
    fdr.temp < 0.001, "***",
    fdr.temp < 0.01,  "**",
    fdr.temp < 0.05,  "*",
    default = "n.s."
  )
), by = display_label]

purine.salvage.ann[, `:=`(
  y_tick = y - 0.015,
  y_text = y + 0.015
)]

ymax <- max(0.5, purine.salvage.ann$y_text, na.rm = TRUE) + 0.03

purine.salvage <- melt(
  purine.salvage,
  id.vars = c("gene", "display_label", "x"),
  measure.vars = c("logFC.hot", "logFC.cold"),
  variable.name = "condition",
  value.name = "logFC"
)

purine.salvage[, condition := fifelse(condition == "logFC.hot", "Hot-evolved", "Cold-evolved")]

png("../Plot/Suppl.Purine_salvage_gene_expression.png",
    width = 6, height = 12, units = "in", res = 450)

ggplot(purine.salvage, aes(x = x, y = logFC, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_segment(data = purine.salvage.ann,
               aes(x = x - 0.35, xend = x + 0.35, y = y, yend = y),
               inherit.aes = FALSE) +
  geom_segment(data = purine.salvage.ann,
               aes(x = x - 0.35, xend = x - 0.35, y = y, yend = y_tick),
               inherit.aes = FALSE) +
  geom_segment(data = purine.salvage.ann,
               aes(x = x + 0.35, xend = x + 0.35, y = y, yend = y_tick),
               inherit.aes = FALSE) +
  geom_text(data = purine.salvage.ann,
            aes(x = x, y = y_text, label = sig_label),
            size = 6, inherit.aes = FALSE) +
  scale_x_continuous(
    breaks = seq_along(label_levels),
    labels = label_levels,
    expand = expansion(add = 0.6)
  ) +
  coord_cartesian(ylim = c(-0.2, ymax)) +
  labs(
    y = "logFC against ancestral state",
    x = "Gene",
    fill = "Population"
  ) +
  scale_fill_manual(values = c(
    "Hot-evolved" = "maroon",
    "Cold-evolved" = "steelblue"
  )) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, angle = 75, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.position = "bottom",
    text = element_text(size = 20)
  )

dev.off()

