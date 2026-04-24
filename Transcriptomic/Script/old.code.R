temp.dt.4.plot <- cpm(y)[significant_interaction,]
temp.dt.4.plot <- setDT(data.frame(temp.dt.4.plot), keep.rownames = T)
temp.dt.4.plot <- melt(temp.dt.4.plot, id.vars = "rn", variable.name = "sample", value.name = "expression")
setnames(temp.dt.4.plot, "rn", "gene")

temp.dt.4.plot[, env := fifelse(grepl("^coldcge", sample), "coldcge", "hotcge")]
temp.dt.4.plot[, group := toupper(sub("^.*_([ACH])[0-9_]*$", "\\1", sample))]
table(temp.dt.4.plot$group)
table(temp.dt.4.plot$env)
temp.dt.4.plot[, group_label := factor(group, levels = c("A", "C", "H"),
                                       labels = c("Ancestral", "Cold evolved", "Hot evolved"))]
temp.dt.4.plot[, env := factor(env, levels = c("coldcge", "hotcge"))]
pdf("../Plot/sig_evo_cge_genes.pdf", width = 8, height = 5)
genes <- unique(temp.dt.4.plot$gene)
for (g in genes) {
  p <- ggplot(temp.dt.4.plot[gene == g], aes(x = group_label, y = log(expression))) +
    geom_boxplot(outlier.shape = NA, width = 0.5,fill = "grey90") +
    geom_jitter(width = 0.2, size = 3, aes(color = group_label)) +
    facet_wrap(~env, nrow = 1) +
    labs(title = g, x = "Group", y = "CPM (log scale)") +
    scale_y_continuous(trans = 'log1p') +
    scale_color_manual(values = c("forestgreen","steelblue","maroon"))+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  print(p)
}
dev.off()



### 8.5 heatmap ----
for (type in names(DE.list)) {
  temp.genelist <- DE.list[[type]]
  temp.dat <- cpm.mat[rownames(cpm.mat) %in% temp.genelist, ]
  scaled_dat <- t(scale(t(temp.dat)))
  annotation_col <- data.frame(Group = group)
  rownames(annotation_col) <- colnames(temp.dat)
  annotation_colors <- list(
    Group = c(
      "coldcgeA" = "forestgreen",
      "coldcgeC" = "steelblue",
      "coldcgeH" = "maroon",
      "hotcgeA"  = "lightgreen",
      "hotcgeC"  = "lightblue",
      "hotcgeH"  = "lightcoral")
  )
  png(filename = paste0("../plot/Heatmap_", type, ".png"),
      width = 8, height = 6, units = "in", res = 600)
  pheatmap(
    mat = scaled_dat,
    cluster_rows = T,  # Use manual clustering
    cluster_cols = FALSE,
    show_rownames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    border_color = NA,
    main = paste0("Heatmap for category: ", type)
  )
  dev.off()
}


### 10. Proportion analysis----
short.list[, temp.prop := abs(.SD$logFC.hot-.SD$logFC.cold)/max(abs(.SD$logFC.hot), abs(.SD$logFC.cold)), by = gene ]
short.list[temp.prop>1, temp.prop := 1]
short.list[, lab.prop := 1 - temp.prop]

plot.dt <- short.list[fdr.lab< 0.05 | fdr.temp<0.05,.(gene, sig, lab.prop, temp.prop)]
setorder(plot.dt, lab.prop)
plot.dt <- melt(plot.dt, id.vars = c("gene","sig"), variable.name = "type", value.name = "proportion")
plot.dt[, type:= factor( type, levels = c("temp.prop", "lab.prop"))]
order_vec <- unique(plot.dt[type == "lab.prop"][order(-proportion), gene])
plot.dt[ , gene := factor(gene, levels = order_vec) ]
plot.dt[,sig := factor(sig, levels = c("Lab DEGs", "Lab X Temp DEGs", "Temp DEGs"))]

png("../Plot/GeneExpression.proportion.png", width = 8, height = 6, units = "in", res = 600)
ggplot(data = plot.dt, aes(x = gene, y = proportion, fill = type)) +
  geom_col(width = 1.00) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL,y= "Proportion",
    title = "Temp- vs. Lab-driven gene expression changes",
    fill = "Type") +
  scale_fill_manual(
    values = c("lab.prop"  = "purple", "temp.prop" = "gold"),
    name   = "",
    labels = c("lab.prop" = "Lab-driven", "temp.prop" = "Temp-driven")
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_blank(),
    legend.position = "bottom",
    plot.caption = element_text(size = 10, vjust = 5, hjust = 0)
  )
dev.off()

rm(plot.dt, order_vec)

## 11. GO enrichment on the three categories of DE genes ----
library(topGO);library(org.Dm.eg.db);library(AnnotationDbi);library(patchwork);library(stringr); library(tidytext)

gene2GO <- AnnotationDbi::select(org.Dm.eg.db, keys   = unique(short.list$gene),
                                 keytype= "FLYBASE",columns= "GO")
gene2GO <- gene2GO[!is.na(gene2GO$GO), ]
gene2GO <- split(gene2GO$GO, gene2GO$FLYBASE)
universe <- names(gene2GO)

run_topgo <- function(target, ont){
  target   <- intersect(target, universe)
  geneList <- factor(as.integer(universe %in% target)); names(geneList) <- universe
  obj <- new("topGOdata", ontology=ont, allGenes=geneList, nodeSize=10,
             annot=annFUN.gene2GO, gene2GO=gene2GO)
  tbl <- GenTable(obj, p = runTest(obj, algorithm="weight01", statistic="fisher"),
                  topNodes=15, orderBy = "GeneRatio", numChar=100)
  setDT(tbl)
  tbl[,GeneRatio := as.numeric(Significant)/as.numeric(Annotated)]
  tbl <- tbl[order(as.numeric(p), decreasing = F),]
  tbl[p=="<1e-30" | p == "< 1e-30", p:= 1e-30][,p := as.numeric(p)]
  tbl[, Term := str_trim(gsub("\\s*\\([^)]*\\)", "", Term))]
  tbl$Term <- factor(str_wrap(tbl$Term,40), levels = rev(str_wrap(tbl$Term,40)))
  tbl <- tbl[p<0.05, ]
  # tbl <- na.omit(tbl)
  return(tbl)
}

lab.bp <- run_topgo(target = short.list[sig == "Lab DEGs",gene],ont = "BP")[,`:=`(DEG = "Lab", ont = "BP")]
lab.mf <- run_topgo(target = short.list[sig == "Lab DEGs",gene],ont = "MF")[,`:=`(DEG = "Lab", ont = "MF")]
lab.cc <- run_topgo(target = short.list[sig == "Lab DEGs",gene],ont = "CC")[,`:=`(DEG = "Lab", ont = "CC")]

temp.bp <- run_topgo(target = short.list[sig == "Temp DEGs",gene],ont = "BP")[,`:=`(DEG = "Temp", ont = "BP")]
temp.mf <- run_topgo(target = short.list[sig == "Temp DEGs",gene],ont = "MF")[,`:=`(DEG = "Temp", ont = "MF")]
temp.cc <- run_topgo(target = short.list[sig == "Temp DEGs",gene],ont = "CC")[,`:=`(DEG = "Temp", ont = "CC")]

inter.bp <- run_topgo(target = short.list[sig == "Lab X Temp DEGs",gene],ont = "BP")[,`:=`(DEG = "Lab X Temp", ont = "BP")]
inter.mf <- run_topgo(target = short.list[sig == "Lab X Temp DEGs",gene],ont = "MF")[,`:=`(DEG = "Lab X Temp", ont = "MF")]
inter.cc <- run_topgo(target = short.list[sig == "Lab X Temp DEGs",gene],ont = "CC")[,`:=`(DEG = "Lab X Temp", ont = "CC")]

CGE.bp <- run_topgo(target = short.list[fdr.env<0.05,gene],ont = "BP")[,`:=`(DEG = "CGE", ont = "BP")]
CGE.mf <- run_topgo(target = short.list[fdr.env<0.05,gene],ont = "MF")[,`:=`(DEG = "CGE", ont = "MF")]
CGE.cc <- run_topgo(target = short.list[fdr.env<0.05,gene],ont = "CC")[,`:=`(DEG = "CGE", ont = "CC")]

# plotting the results
library(ggplot2); library(tidytext); library(viridis); library(scales)
## 11.1 BP ----
GO.enrich.bp <- rbind(lab.bp, temp.bp, inter.bp, CGE.bp)
GO.enrich.bp[, DEG := factor(DEG, levels = c("Lab", "Temp", "Lab X Temp", "CGE"))]
sz_breaks  <- pretty(GO.enrich.bp$Significant, 4)
col_breaks <- pretty(-log10(GO.enrich.bp$p),   5)

png("../Plot/GO.enrich.BP.png", 16, 12, units = "in", res = 600)
ggplot(GO.enrich.bp,
       aes(GeneRatio,
           reorder_within(Term, GeneRatio, DEG),
           size   = Significant,
           colour = -log10(p))) +
  facet_wrap(~DEG, scales = "free", nrow = 2) +
  geom_point(shape = 15) +
  scale_y_reordered() +
  scale_colour_viridis_c(option = "plasma",
                         breaks  = col_breaks,
                         name    = expression(-log[10](p))) +
  scale_size_continuous(breaks = sz_breaks,
                        range  = c(1.5, 8),
                        name   = "Gene count") +
  labs(x = "Gene ratio",
       y = NULL,
       title = "Gene Ontology Enrichment: Biological Process") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.y = element_blank())
dev.off()

## 11.2 MF ----
GO.enrich.mf <- rbind(lab.mf, temp.mf, inter.mf, CGE.mf)
GO.enrich.mf[, DEG := factor(DEG, levels = c("Lab", "Temp", "Lab X Temp", "CGE"))]
sz_breaks  <- pretty(GO.enrich.mf$Significant, 4)
col_breaks <- pretty(-log10(GO.enrich.mf$p),   5)

png("../Plot/GO.enrich.MF.png", 16, 12, units = "in", res = 600)
ggplot(GO.enrich.mf,
       aes(GeneRatio,
           reorder_within(Term, GeneRatio, DEG),
           size   = Significant,
           colour = -log10(p))) +
  facet_wrap(~DEG, scales = "free", nrow = 2) +
  geom_point(shape = 15) +
  scale_y_reordered() +
  scale_colour_viridis_c(option = "plasma",
                         breaks  = col_breaks,
                         name    = expression(-log[10](p))) +
  scale_size_continuous(breaks = sz_breaks,
                        range  = c(1.5, 8),
                        name   = "Gene count") +
  labs(x = "Gene ratio",
       y = NULL,
       title = "Gene Ontology Enrichment: Molecular Function") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.y = element_blank())
dev.off()

## 11.3 CC ----
GO.enrich.cc <- rbind(lab.cc, temp.cc, inter.cc, CGE.cc)
GO.enrich.cc[, DEG := factor(DEG, levels = c("Lab", "Temp", "Lab X Temp", "CGE"))]
sz_breaks  <- pretty(GO.enrich.cc$Significant, 4)
col_breaks <- pretty(-log10(GO.enrich.cc$p),   5)

png("../Plot/GO.enrich.CC.png", 16, 12, units = "in", res = 600)
ggplot(GO.enrich.cc,
       aes(GeneRatio,
           reorder_within(Term, GeneRatio, DEG),
           size   = Significant,
           colour = -log10(p))) +
  facet_wrap(~DEG, scales = "free", nrow = 2) +
  geom_point(shape = 15) +
  scale_y_reordered() +
  scale_colour_viridis_c(option = "plasma",
                         breaks  = col_breaks,
                         name    = expression(-log[10](p))) +
  scale_size_continuous(breaks = sz_breaks,
                        range  = c(1.5, 8),
                        name   = "Gene count") +
  labs(x = "Gene ratio",
       y = NULL,
       title = "Gene Ontology Enrichment: Cellular Component") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.y = element_blank())
dev.off()

rm(run_topgo, GO.enrich.bp,sz_breaks, col_breaks, GO.enrich.cc, GO.enrich.mf, universe,lab.bp, lab.mf, lab.cc, temp.bp, temp.mf, temp.cc, inter.bp, inter.mf, inter.cc, CGE.bp, CGE.mf, CGE.cc)

### 12. tissue enrichment ----
cont_table <- function(query, background, classifier) {
  ## Keep only genes that belong to the background universe
  query      <- intersect(query,      background)
  classifier <- intersect(classifier, background)

  ## Logical flags over the universe (= background)
  in_q <- background %in% query
  in_c <- background %in% classifier

  matrix(c(sum(in_q & in_c),        # both
           sum(in_q & !in_c),       # query only
           sum(!in_q & in_c),       # classifier only
           sum(!in_q & !in_c)),     # neither
         nrow = 2, byrow = TRUE,
         dimnames = list(Classified = c("Yes", "No"),
                         Query       = c("Yes", "No")))
}
library(readxl)
flyatlas2 <- as.data.table(read_excel(path = "./FlyAtlas2_gene_data_2025.xlsx", sheet = "FPKMs 2025"))
cols_interest <- c("FBgn",grep(" M$", names(flyatlas2), value = TRUE),"Testis", "Accessory Gland")
cols_interest <- cols_interest[!grepl("SD", cols_interest)]

flyatlas2 <- flyatlas2[, ..cols_interest]
flyatlas2 <- flyatlas2[`Whole M` != 0]

expr_cols <- setdiff(names(flyatlas2), c("FBgn", "Whole M"))
expr_cols <- expr_cols[!grepl("SD", expr_cols)]

flyatlas2[ ,(expr_cols) := lapply(.SD,\(x) log2( x / `Whole M` )), .SDcols = expr_cols]
flyatlas2
background <- rownames(y_filtered$counts)
tissue_sets <- lapply(expr_cols,\(col) {flyatlas2[!is.na(get(col)) & get(col) > 1 ,FBgn]}) # log₂FC > 1  (≥ 2-fold)
names(tissue_sets) <- expr_cols
# Initialize lists to store p-values, odds ratios, and enriched gene IDs
p.val <- odds <- tissue_enriched_ID <- vector("list", length(tissue_sets))
names(p.val) <- names(odds) <- names(tissue_enriched_ID) <- names(tissue_sets)

for (tissue in names(tissue_sets)) {
  classifier <- tissue_sets[[tissue]]

  p.val[[tissue]]  <- numeric(length(DE.list))
  odds[[tissue]]   <- numeric(length(DE.list))
  tmp_hits         <- vector("list", length(DE.list))

  for (k in seq_along(DE.list)) {
    query <- DE.list[[k]]
    ft    <- fisher.test(
      cont_table(query, background, classifier),
      alternative = "greater"
    )

    p.val[[tissue]][k] <- ft$p.value
    odds [[tissue]][k] <- unname(ft$estimate)
    tmp_hits[[k]]      <- intersect(query, classifier)
  }

  names(tmp_hits)                <- names(DE.list)
  tissue_enriched_ID[[tissue]]   <- tmp_hits
}

for (k in seq_along(DE.list)) {
  pv      <- sapply(p.val, `[`, k)        # grab kth column
  adj.pv  <- p.adjust(pv, method = "BH")  # Benjamini-Hochberg

  for (tissue in names(p.val))            # put back
    p.val[[tissue]][k] <- adj.pv[tissue]
}

# Reshape the data to a format suitable for ggplot2
plot_data <- data.frame()

for (i in 1:length(DE.list)) {
  for (j in 1:length(p.val)) {
    temp_data <- data.frame(
      Tissue = names(p.val)[j],
      DEG_List = names(DE.list)[i],
      P_Value = p.val[[j]][i],
      Odds = odds[[j]][i]
    )
    plot_data <- rbind(plot_data, temp_data)
  }
}

plot_data$Alpha_Category <- cut(plot_data$P_Value,
                                breaks = c(0, 0.001, 0.005, 0.01, 0.05, 1),
                                labels = c("< 0.001", "0.001 - 0.005", "0.005 - 0.01", "0.01 - 0.05", "> 0.05"),
                                include.lowest = TRUE)

# Assign transparency levels based on the p-value categories
plot_data$Alpha <- factor(plot_data$Alpha_Category,
                          levels = c("< 0.001", "0.001 - 0.005", "0.005 - 0.01", "0.01 - 0.05", "> 0.05"))

# Shorten the alpha legend label and modify x-axis labels
tissue_labels <- c("Hd" = "Head", "Ey" = "Eye", "Br" = "Brain",
                   "Tg" = "Thoracoabdominal Ganglion", "Cr" = "Crop",
                   "Mg" = "Midgut", "Hg" = "Hindgut", "Tu" = "Malpighian Tubules",
                   "Sg" = "Salivary Gland", "Cs" = "Carcass", "Ag" = "Accessory Glands",
                   "Ts" = "Testis",  "Ap" = "Fat Body")

plot_data$Tissue <- sub("_M", "", plot_data$Tissue)
plot_data$DEG_List <- factor(plot_data$DEG_List, levels = c("lab", "temp", "inter", "env"))
# Create the updated barplot with shorter label and updated tissue names
png("../Plot/Tissue.enrich.png", width = 16, height = 8, units = "in", res = 400)
ggplot(plot_data, aes(x = Tissue, y = log2(Odds+0.1), fill = -log10(P_Value))) +
  geom_bar(stat = "identity", width = 0.75) +
  facet_wrap(~DEG_List, scales = "free_x", nrow = 1) +
  scale_fill_gradient(low = "white", high = "darkblue") +
  labs(title = "Tissue Enrichment for DEGs",
       x = "Tissue",
       y = "Log Odds Ratio",
       fill = "-log10(p-value)") +
  coord_flip()+
  geom_hline(yintercept = 0, linewidth = 0.8, linetype = "dashed",color = "red", alpha = 0.4)+
  scale_x_discrete(labels = tissue_labels) +
  theme_minimal() +
  theme(text = element_text(face = "bold", size = 16),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),legend.title = element_text(face = "bold"),
  )
dev.off()







### 14. Individual gene focus ----

short.list[gene == "FBgn0193894",]
## 14.1 LysP (Lysozyme P) ----
# this enzyme is found to be natural antibiotic in mammals,
#  but it was only excreted in the Dmel to digest microbs/food.
short.list[gene == "FBgn0004429",] ## <-- LAB gene

## 14.2 Anp (Andropin) also known as Cec ----
####Andropin (Anp) encodes an antibacterial peptide expressed in the male genital tract.
short.list[gene == "FBgn0000094",] ## <-- LAB gene

## 14.3 DptA (Diptericin A) ----
#### It's an antibacterial peptide with activity against Gram-negative bacteria
#### It is expressed in the fat body during the systemic immune response and in various epithelia.
short.list[gene == "FBgn0004240",] ## <-- LAB gene

## 14.4 DptB (Diptericin B) ----
#### It's an antibacterial peptide with activity against Gram-negative bacteria
#### It is expressed in the fat body during the systemic immune response and in various epithelia.
short.list[gene == "FBgn0034407",] ## <-- LAB gene

## 14.5 Def (Defensin) ----
#### Defensin (Def) encodes an antibacterial peptide with activity against Gram-positive bacteria.
#### It is induced in the fat body during the systemic immune response and is expressed in various epithelia.
short.list[gene == "FBgn0010385",] ## <-- non-sig

## 14.6 Mtk (Metchnikowin) ----
#### Metchnikowin (Mtk) encodes an antifungal peptide that is secreted from the fat body
####  during the systemic immune response, and is produced by various epithelia.
short.list[gene == "FBgn0014865",] ## <-- non-sig


short.list[gene == "FBgn0029913",]
### 14.7 all FB defense against Gram-negative bacteria genes----
Gram_negative_defense_gene <- unlist(fread(file = "./Defense_gram_negative_genes.txt", header = F))
Gram_negative_defense_gene <- short.list[gene %in% Gram_negative_defense_gene,][!is.na(sig),]
Gram_negative_defense_gene$immunity <- "Gram-Negative"

### 14.7 all FB defense against Gram-positive bacteria genes----
Gram_positive_defense_gene <- unlist(fread(file = "./Defense_gram_positive_genes.txt", header = F))
Gram_positive_defense_gene <- short.list[gene %in% Gram_positive_defense_gene,][!is.na(sig),]
Gram_positive_defense_gene$immunity <- "Gram-Positive"
Defense_genes <- rbind(Gram_negative_defense_gene, Gram_positive_defense_gene)

Defense_genes[gene =="FBgn0014865", immunity := "Gram-Both"]
Defense_genes[gene =="FBgn0039102", immunity := "Gram-Both"]
Defense_genes[,.N, by = .(immunity, sig)]
Defense_genes <- Defense_genes[!duplicated(gene),]
saveRDS(Defense_genes, file = "./Defense_genes.RDS")
Defense_genes <- readRDS("./Defense_genes.RDS")



Defense_genes[, significant := fdr.hot < 0.05 | fdr.cold < 0.05]
Defense_genes[, gene_label := ifelse(significant, paste(gene, "*"), gene)]
gene_order <- Defense_genes$gene_label
Defense_genes[, gene_label := factor(gene_label, levels = gene_order)]
long_dt <- melt(Defense_genes, id.vars = c("gene_label","sig","immunity"), measure.vars = c("logFC.hot", "logFC.cold"),
                variable.name = "condition", value.name = "logFC")
long_dt[, condition := ifelse(condition == "logFC.hot", "Hot", "Cold")]
ggplot(long_dt, aes(x = gene_label, y = logFC, fill = condition)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_minimal() +
  facet_wrap(~immunity, scales = "free_x")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(title = "Log Fold Change for Hot and Cold Conditions by Gene",
       x = "Gene",
       y = "logFC",
       fill = "Condition") +
  scale_fill_manual(values = c("Hot" = "red", "Cold" = "blue"))



##  Check whether all genes in Tissue sets are changing in the same direction in the hot-evolution----
carcass.genes <- unlist(tissue_sets["Carcass M"])
cpm.mat <- cpm(y_filtered, normalized.lib.sizes = T, log = T, prior.count = 1)
cpm.mat <- cpm.mat[rownames(cpm.mat) %in% carcass.genes,]
temp.scale <- t(apply(cpm.mat,1,scale, center = T, scale = T));colnames(temp.scale) <- colnames(cpm.mat)
temp.dat <- setDT(as.data.frame(temp.scale), keep.rownames = T)
long_data <- melt(temp.dat, id.vars = "rn", variable.name = "sample", value.name = "expression")
long_data[, c("environment", "replicate") := tstrsplit(sample, "e_")]
long_data[,environment := paste0(environment, "e")]
long_data[, population := substr(replicate, 1, 1)]
long_data[, replicate := substr(replicate, 2, 2)]
long_avg <- long_data[,.(mean(expression)), by = .(rn, environment,population)]
setnames(long_avg, 4, "expression")
png(filename = "../Plot/carcass.expression.unpola.png", width = 7, height = 5, units = "in", res = 350)
ggplot(long_avg, aes(x = population, y = expression, group = rn, color = population)) +
  geom_line(color = "grey40", alpha = 0.1, linewidth = 1.2)+
  facet_wrap(~environment)+
  geom_point(inherit.aes = F, mapping = aes(x = population, y = expression, color = population, alpha = 0.4))+
  # geom_jitter(size  =2, width = 0.2)+
  labs(title = "General expression pattern of Carcass marker genes",
       x = "Population", y = "logCPM (Expression Intensity)") +
  scale_color_manual(values = c("forestgreen","steelblue","maroon"))+
  theme_minimal()
dev.off()

### okay, we may need to polarize the expression to get a clearer picture
pop_means <- long_data[ ,.(pop_mean = mean(expression, na.rm = TRUE)),by = .(rn, environment,population)]
genes_to_flip <- pop_means[ ,.(flip = pop_mean[environment == "coldcge" & population == "A"] > mean(pop_mean[environment == "coldcge" & population != "A"])),by = rn][flip == TRUE, rn]
long_data[rn %in% genes_to_flip, expression := -expression]
long_avg <- long_data[,.(mean(expression)), by = .(rn, environment,population)]
setnames(long_avg, 4, "expression")
png(filename = "../Plot/carcass.expression.pola.png", width = 7, height = 5, units = "in", res = 350)
ggplot(long_avg, aes(x = population, y = expression, group = rn, color = population)) +
  geom_line(color = "grey40", alpha = 0.1, linewidth = 1.2)+
  facet_wrap(~environment)+
  geom_point(inherit.aes = F, mapping = aes(x = population, y = expression, color = population, alpha = 0.4))+
  # geom_jitter(size  =2, width = 0.2)+
  labs(title = "General expression pattern of Carcass marker genes",
       x = "Population", y = "logCPM (Expression Intensity)") +
  scale_color_manual(values = c("forestgreen","steelblue","maroon"))+
  theme_minimal()
dev.off()

## okay, maybe mroe specific Carcass genes have clearer pattern
tissue_sets <- lapply(expr_cols,\(col) {flyatlas2[!is.na(get(col)) & get(col) > 2 ,FBgn]}) # log₂FC > 2  (≥ 4-fold !!!)
names(tissue_sets) <- expr_cols
carcass.genes <- unlist(tissue_sets["Carcass M"])
cpm.mat <- cpm(y_filtered, normalized.lib.sizes = T, log = T, prior.count = 1)
cpm.mat <- cpm.mat[rownames(cpm.mat) %in% carcass.genes,]
temp.scale <- t(apply(cpm.mat,1,scale, center = T, scale = T));colnames(temp.scale) <- colnames(cpm.mat)
temp.dat <- setDT(as.data.frame(temp.scale), keep.rownames = T)
long_data <- melt(temp.dat, id.vars = "rn", variable.name = "sample", value.name = "expression")
long_data[, c("environment", "replicate") := tstrsplit(sample, "e_")]
long_data[,environment := paste0(environment, "e")]
long_data[, population := substr(replicate, 1, 1)]
long_data[, replicate := substr(replicate, 2, 2)]
long_avg <- long_data[,.(mean(expression)), by = .(rn, environment,population)]
setnames(long_avg, 4, "expression")
png(filename = "../Plot/carcass.expression.unpola.logFC.greater2.png", width = 7, height = 5, units = "in", res = 350)
ggplot(long_avg, aes(x = population, y = expression, group = rn, color = population)) +
  geom_line(color = "grey40", alpha = 0.3, linewidth = 1.2)+
  facet_wrap(~environment)+
  geom_point(inherit.aes = F, mapping = aes(x = population, y = expression, color = population, alpha = 0.7))+
  # geom_jitter(size  =2, width = 0.2)+
  labs(title = "General expression pattern of top 20 Carcass marker genes",
       x = "Population", y = "logCPM (Expression Intensity)") +
  scale_color_manual(values = c("forestgreen","steelblue","maroon"))+
  theme_minimal()
dev.off()


## 10. downstream gene list look up ----
mart <- useMart(
  biomart = "metazoa_mart",
  dataset = "dmelanogaster_eg_gene",
  host = "https://metazoa.ensembl.org"
)

## 10.1 glycolysis ----
glyco <- getBM(attributes = c("ensembl_gene_id", 
                              "external_gene_name", 
                              "flybase_gene_id",
                              "go_id",
                              "name_1006"),
               filters = "go",
               values = "GO:0006096", mart = mart)
short.list[gene %in% unique(glyco$ensembl_gene_id),]
glyco.expr <- cpm.mat[unique(glyco$ensembl_gene_id),]
glyco.expr <- setDT(as.data.frame(glyco.expr),keep.rownames = T)

glyco.expr.long <- melt(data = glyco.expr,id.vars = "rn",
                        variable.name = "sample",
                        value.name = "expr")
glyco.expr.long[, c("cge", "tmp") := tstrsplit(sample, "_", fixed = TRUE)]
glyco.expr.long[, evo := sub("^([A-Z]).*$", "\\1", tmp)]
glyco.expr.long[, rep := as.integer(sub("^[A-Z]", "", tmp))]
glyco.expr.long[, tmp := NULL]
glyco.expr.av <- glyco.expr.long[,.(avg.expr = mean(expr)), by = .(rn,cge, evo)]

ggplot(data = glyco.expr.av, aes(x = cge, y = avg.expr, colour = evo))+
  geom_boxplot()
## 10.2 pentose phosphate ----
pentose.p <- getBM(
  attributes = c("ensembl_gene_id", 
                 "external_gene_name", 
                 "flybase_gene_id",
                 "go_id",
                 "name_1006"),
  filters = "go",
  values = "GO:0006098",
  mart = mart
)

short.list[gene %in% unique(pentose.p$ensembl_gene_id),]
pentose.p.expr <- cpm.mat[unique(pentose.p$ensembl_gene_id),]
pentose.p.expr <- setDT(as.data.frame(pentose.p.expr),keep.rownames = T)

pentose.p.expr.long <- melt(data = pentose.p.expr,id.vars = "rn",
                            variable.name = "sample",
                            value.name = "expr")
pentose.p.expr.long[, c("cge", "tmp") := tstrsplit(sample, "_", fixed = TRUE)]
pentose.p.expr.long[, evo := sub("^([A-Z]).*$", "\\1", tmp)]
pentose.p.expr.long[, rep := as.integer(sub("^[A-Z]", "", tmp))]
pentose.p.expr.long[, tmp := NULL]
pentose.p.expr.av <- pentose.p.expr.long[,.(avg.expr = mean(expr)), by = .(rn,cge, evo)]

ggplot(data = pentose.p.expr.av, aes(x = cge, y = avg.expr, colour = evo))+
  geom_boxplot()

## 10.3 purine biosynthesis ----
purine.syn <- getBM(
  attributes = c("ensembl_gene_id", 
                 "external_gene_name", 
                 "flybase_gene_id",
                 "go_id",
                 "name_1006"),
  filters = "go",
  values = "GO:0006164",
  mart = mart
)

short.list[gene %in% unique(purine.syn$ensembl_gene_id),]
purine.syn.expr <- cpm.mat[unique(purine.syn$ensembl_gene_id),]
purine.syn.expr <- setDT(as.data.frame(purine.syn.expr),keep.rownames = T)

purine.syn.expr.long <- melt(data = purine.syn.expr,id.vars = "rn",
                             variable.name = "sample",
                             value.name = "expr")
purine.syn.expr.long[, c("cge", "tmp") := tstrsplit(sample, "_", fixed = TRUE)]
purine.syn.expr.long[, evo := sub("^([A-Z]).*$", "\\1", tmp)]
purine.syn.expr.long[, rep := as.integer(sub("^[A-Z]", "", tmp))]
purine.syn.expr.long[, tmp := NULL]
purine.syn.expr.av <- purine.syn.expr.long[,.(avg.expr = mean(expr)), by = .(rn,cge, evo)]

ggplot(data = purine.syn.expr.av, aes(x = cge, y = avg.expr, colour = evo))+
  geom_jitter()

## 10.4 Ribosome biogenesis ----
ribosome.syn <- getBM(
  attributes = c("ensembl_gene_id", 
                 "external_gene_name", 
                 "flybase_gene_id",
                 "go_id",
                 "name_1006"),
  filters = "go",
  values = "GO:0042254",
  mart = mart
)

short.list[gene %in% unique(ribosome.syn$ensembl_gene_id),]

## 11 use KEGG instead of GO ----
### 11.1 purine biosynthesis ----
library(KEGGREST)

ppp <- keggGet("dme00030")[[1]]$GENE
ppp_fbgn <- ppp[seq(1, length(ppp), 2)]

ppp_fbgn
