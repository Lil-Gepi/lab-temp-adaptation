rm(list=ls())
gc()

library(data.table)
setDTthreads(percent = 80) 
library(edgeR)
library(ggplot2)
library(KEGGREST)
library(org.Dm.eg.db)
library(AnnotationDbi)

cpm.mat <- readRDS(file = "./cpm.mat.RDS")
## 1. Setup: KEGG → FBgn mapping helpers----
# KEGG pathway -> CG IDs (e.g., "CG10202") + descriptions
kegg_pathway_genes <- function(path_id){
  x <- keggGet(path_id)[[1]]
  g <- x$GENE
  dt <- data.table(
    kegg_gene = g[seq(1, length(g), 2)],
    desc      = g[seq(2, length(g), 2)]
  )
  dt[, cg := sub("^dme:", "", kegg_gene)]
  dt[, cg := sub("^Dmel_", "", cg)]
  dt
}

# Map CG -> FBgn (plus optional symbol)
cg_to_fbgn <- function(cg_vec){
  fbgn <- AnnotationDbi::mapIds(
    org.Dm.eg.db,
    keys = unique(cg_vec),
    keytype = "FLYBASECG",
    column = "FLYBASE",
    multiVals = "first"
  )
  data.table(cg = names(fbgn), fbgn = unname(fbgn))
}

fbgn_to_symbol <- function(fbgn_vec){
  sym <- AnnotationDbi::mapIds(
    org.Dm.eg.db,
    keys = unique(fbgn_vec),
    keytype = "FLYBASE",
    column = "SYMBOL",
    multiVals = "first"
  )
  data.table(fbgn = names(sym), symbol = unname(sym))
}

## 2 Pull KEGG pathways I care about----
pw_ids <- c(
  Purine = "dme00230",
  Glycolysis = "dme00010",
  TCA = "dme00020"
)

pw_dt <- rbindlist(lapply(names(pw_ids), function(nm){
  dt <- kegg_pathway_genes(pw_ids[[nm]])
  dt[, pathway := nm]
  dt
}))


map_dt <- cg_to_fbgn(pw_dt$cg)
pw_dt <- merge(pw_dt, map_dt, by = "cg", all.x = TRUE)

# only keep genes that map + exist in your matrix
pw_dt <- pw_dt[!is.na(fbgn)]
pw_dt <- pw_dt[fbgn %in% rownames(cpm.mat)]

# add symbols for plotting
sym_dt <- fbgn_to_symbol(pw_dt$fbgn)
pw_dt <- merge(pw_dt, sym_dt, by = "fbgn", all.x = TRUE)
pw_dt[, gene_label := fifelse(is.na(symbol) | symbol == "", fbgn, paste0(symbol, " (", fbgn, ")"))]

pur <- pw_dt[pathway == "Purine"]

# crude but effective keyword sets (edit/expand as needed)
de_novo_keys <- c(
  "amidophosphoribosyltransferase", "GAR", "FGAM", "AIR", "SAICAR", "AICAR",
  "IMP", "adenylosuccinate", "IMP dehydrogenase", "GMP synthetase"
)

salvage_keys <- c(
  "phosphoribosyltransferase", # HGPRT/APRT-like
  "adenosine kinase",
  "purine nucleoside phosphorylase",
  "xanthine", "hypoxanthine", "guanine", "adenine",
  "nucleotidase", "nucleosidase"
)

pur[, purine_class := "Other"]
pur[grepl(paste(de_novo_keys, collapse="|"), desc, ignore.case = TRUE), purine_class := "De novo"]
pur[grepl(paste(salvage_keys, collapse="|"), collapse="|", ignore.case = TRUE)]  # (ignore; just to show)
pur[grepl(paste(salvage_keys, collapse="|"), desc, ignore.case = TRUE), purine_class := "Salvage"]

table(pur$purine_class)

expr_long_for_fbgn <- function(fbgn_vec, label_dt){
  m <- cpm.mat[fbgn_vec, , drop = FALSE]
  dt <- as.data.table(m, keep.rownames = "fbgn")
  dt <- melt(dt, id.vars = "fbgn", variable.name = "sample", value.name = "expr")
  dt[, c("cge", "tmp") := tstrsplit(sample, "_", fixed = TRUE)]
  dt[, evo := sub("^([A-Z]).*$", "\\1", tmp)]
  dt[, rep := as.integer(sub("^[A-Z]", "", tmp))]
  dt[, tmp := NULL]
  dt <- merge(dt, label_dt, by = "fbgn", all.x = TRUE)
  dt
}

library(ggplot2)

make_avg_plot <- function(dt_long, title){
  avg <- dt_long[, .(avg.expr = scale(expr)), by = .(fbgn, gene_label, cge, evo)]
  ggplot(avg, aes(x = cge, y = avg.expr, color = evo)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
    theme_minimal(base_size = 12) +
    labs(title = title, x = "Environment", y = "Mean logCPM (per gene)") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# Glycolysis
gly_fbgn <- unique(pw_dt[pathway == "Glycolysis", fbgn])
gly_long <- expr_long_for_fbgn(gly_fbgn, unique(pw_dt[pathway=="Glycolysis", .(fbgn, gene_label)]))
p_gly <- make_avg_plot(gly_long, "Glycolysis / Gluconeogenesis (KEGG dme00010)")

# TCA
tca_fbgn <- unique(pw_dt[pathway == "TCA", fbgn])
tca_long <- expr_long_for_fbgn(tca_fbgn, unique(pw_dt[pathway=="TCA", .(fbgn, gene_label)]))
p_tca <- make_avg_plot(tca_long, "Citrate cycle / TCA (KEGG dme00020)")

# Purine de novo
pur_dn_fbgn <- unique(pur[purine_class=="De novo", fbgn])
pur_dn_long <- expr_long_for_fbgn(pur_dn_fbgn, unique(pur[purine_class=="De novo", .(fbgn, gene_label)]))
p_pur_dn <- make_avg_plot(pur_dn_long, "Purine metabolism — De novo subset (from KEGG dme00230)")

# Purine salvage
pur_sal_fbgn <- unique(pur[purine_class=="Salvage", fbgn])
pur_sal_long <- expr_long_for_fbgn(pur_sal_fbgn, unique(pur[purine_class=="Salvage", .(fbgn, gene_label)]))
p_pur_sal <- make_avg_plot(pur_sal_long, "Purine metabolism — Salvage subset (from KEGG dme00230)")

