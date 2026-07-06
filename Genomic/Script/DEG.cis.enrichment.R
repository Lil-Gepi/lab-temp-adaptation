## Locus plots for all DEG categories ----
## Output:
## 1 PDF per DEG category with candidate SNPs
## 1 PDF per DEG category without candidate SNPs
## 1 combined TSV summary
##
## Gene symbols are added with org.Dm.eg.db when available.

rm(list = ls())

library(data.table)
library(ggplot2)

## 1. Paths ----
genomic_dir <- "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Genomic/Script"
rna_dir     <- "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Script"
gtf_file    <- "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Data/dsimM252v1.2.gtf"
outdir      <- "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Genomic/Plot/DEG_locus_position_AF"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


## 2. Load data ----
af     <- readRDS(file.path(genomic_dir, "af.covfilter.RDS"))
p.cold <- readRDS(file.path(genomic_dir, "p.cmh.cold.RDS"))
p.hot  <- readRDS(file.path(genomic_dir, "p.cmh.hot.RDS"))
short.list <- readRDS(file.path(rna_dir, "short.list.RDS"))

setDT(af)
setDT(p.cold)
setDT(p.hot)
setDT(short.list)

setnames(p.cold, names(p.cold)[3], "q_cold")
setnames(p.hot,  names(p.hot)[3],  "q_hot")

af[, CHROM := as.character(CHROM)]
p.cold[, CHROM := as.character(CHROM)]
p.hot[, CHROM := as.character(CHROM)]

setkey(af, CHROM, POS)


## 3. Candidate SNP table ----
cand <- merge(
  p.cold[, .(CHROM, POS, q_cold)],
  p.hot[,  .(CHROM, POS, q_hot)],
  by = c("CHROM", "POS"),
  all = TRUE
)

cand[, candidate := fifelse(q_cold < 0.05 & q_hot < 0.05, "Shared",
                            fifelse(q_cold < 0.05, "Cold",
                                    fifelse(q_hot < 0.05, "Hot", "None")))]

cand[is.na(candidate), candidate := "None"]
setkey(cand, CHROM, POS)


## 4. DEG sets ----
DE.list <- list(
  Lab_DEGs      = short.list[sig == "Lab DEGs", gene],
  Temp_DEGs     = short.list[sig == "Temp DEGs", gene],
  LabXTemp_DEGs = short.list[sig == "Lab X Temp DEGs", gene],
  Hot_DEGs      = short.list[fdr.hot < 0.05, gene],
  Cold_DEGs     = short.list[fdr.cold < 0.05, gene]
)


## 5. Read GTF ----
gtf_lines <- readLines(gtf_file, warn = FALSE)
gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]
gtf_lines <- gtf_lines[nzchar(gtf_lines)]

gtf <- fread(
  text = paste(gtf_lines, collapse = "\n"),
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE
)

setnames(
  gtf,
  c("CHROM", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes")
)

gtf[, CHROM := sub("^chr", "", as.character(CHROM))]
gtf[, gene := sub('.*gene_id "([^"]+)".*', "\\1", attributes)]

gtf <- gtf[CHROM %in% c("X", "2L", "2R", "3L", "3R")]
gtf <- gtf[grepl("^FBgn", gene)]


## 6. Map FBgn to clean gene symbols ----
gene_symbol <- data.table(gene = unique(gtf$gene), symbol = unique(gtf$gene))

if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
    requireNamespace("org.Dm.eg.db", quietly = TRUE)) {
  
  gene_symbol <- as.data.table(
    AnnotationDbi::select(
      org.Dm.eg.db::org.Dm.eg.db,
      keys = unique(gtf$gene),
      keytype = "FLYBASE",
      columns = "SYMBOL"
    )
  )
  
  setnames(gene_symbol, c("FLYBASE", "SYMBOL"), c("gene", "symbol"))
  gene_symbol <- gene_symbol[!is.na(symbol)]
  gene_symbol <- unique(gene_symbol)
  gene_symbol <- gene_symbol[, .(symbol = symbol[1]), by = gene]
  
} else {
  message("AnnotationDbi/org.Dm.eg.db not installed. Gene symbols will stay as FBgn IDs.")
}


## 7. Gene ranges ----
gene_ranges <- gtf[, .(
  CHROM = unique(CHROM)[1],
  strand = unique(strand)[1],
  start = min(as.integer(start), na.rm = TRUE),
  end   = max(as.integer(end), na.rm = TRUE)
), by = gene]

gene_ranges <- merge(gene_ranges, gene_symbol, by = "gene", all.x = TRUE)
gene_ranges[is.na(symbol) | symbol == "", symbol := gene]

gene_ranges <- gene_ranges[!is.na(start) & !is.na(end)]
setkey(gene_ranges, gene)


## 8. Cis settings ----
cis_upstream <- 5000
cis_downstream <- 1000


## 9. Build summary table first ----
all_summary <- list()
k <- 1

for (set_name in names(DE.list)) {
  
  message("Summarizing ", set_name)
  
  genes <- intersect(DE.list[[set_name]], gene_ranges$gene)
  
  for (gene_id in genes) {
    
    g <- gene_ranges[gene == gene_id]
    if (nrow(g) == 0) next
    
    if (g$strand == "+") {
      cis_start <- max(1, g$start - cis_upstream)
      cis_end   <- g$end + cis_downstream
      promoter_start <- max(1, g$start - cis_upstream)
      promoter_end   <- g$start - 1
    } else {
      cis_start <- max(1, g$start - cis_downstream)
      cis_end   <- g$end + cis_upstream
      promoter_start <- g$end + 1
      promoter_end   <- g$end + cis_upstream
    }
    
    region <- af[CHROM == g$CHROM & POS >= cis_start & POS <= cis_end]
    
    if (nrow(region) == 0) {
      all_summary[[k]] <- data.table(
        DEG_set = set_name,
        gene = gene_id,
        symbol = g$symbol,
        CHROM = g$CHROM,
        strand = g$strand,
        gene_start = g$start,
        gene_end = g$end,
        cis_start = cis_start,
        cis_end = cis_end,
        n_snps = 0,
        n_cold = 0,
        n_hot = 0,
        n_shared = 0,
        n_snp_cis = 0,
        n_cold_cis = 0,
        n_hot_cis = 0,
        n_shared_cis = 0,
        n_snp_body = 0,
        n_cold_body = 0,
        n_hot_body = 0,
        n_shared_body = 0,
        has_candidate = FALSE
      )
      k <- k + 1
      next
    }
    
    region <- merge(region, cand, by = c("CHROM", "POS"), all.x = TRUE)
    region[is.na(candidate), candidate := "None"]
    
    region[, in_body := POS >= g$start & POS <= g$end]
    region[, in_cis_only := !in_body]
    
    all_summary[[k]] <- data.table(
      DEG_set = set_name,
      gene = gene_id,
      symbol = g$symbol,
      CHROM = g$CHROM,
      strand = g$strand,
      gene_start = g$start,
      gene_end = g$end,
      cis_start = cis_start,
      cis_end = cis_end,
      
      n_snps = nrow(region),
      n_cold = region[candidate %in% c("Cold", "Shared"), .N],
      n_hot = region[candidate %in% c("Hot", "Shared"), .N],
      n_shared = region[candidate == "Shared", .N],
      
      n_snp_cis = region[in_cis_only == TRUE, .N],
      n_cold_cis = region[in_cis_only == TRUE & candidate %in% c("Cold", "Shared"), .N],
      n_hot_cis = region[in_cis_only == TRUE & candidate %in% c("Hot", "Shared"), .N],
      n_shared_cis = region[in_cis_only == TRUE & candidate == "Shared", .N],
      
      n_snp_body = region[in_body == TRUE, .N],
      n_cold_body = region[in_body == TRUE & candidate %in% c("Cold", "Shared"), .N],
      n_hot_body = region[in_body == TRUE & candidate %in% c("Hot", "Shared"), .N],
      n_shared_body = region[in_body == TRUE & candidate == "Shared", .N],
      
      has_candidate = any(region$candidate != "None")
    )
    
    k <- k + 1
  }
}

all_summary <- rbindlist(all_summary, fill = TRUE)

fwrite(
  all_summary,
  file.path(outdir, "All_DEG_sets_locus_summary.tsv"),
  sep = "\t"
)


## 10. Plot PDFs split by candidate presence ----
for (set_name in names(DE.list)) {
  
  for (candidate_group in c("with_candidate", "without_candidate")) {
    
    message("Plotting ", set_name, " / ", candidate_group)
    
    if (candidate_group == "with_candidate") {
      genes_to_plot <- all_summary[DEG_set == set_name & has_candidate == TRUE, gene]
    } else {
      genes_to_plot <- all_summary[DEG_set == set_name & has_candidate == FALSE, gene]
    }
    
    if (length(genes_to_plot) == 0) {
      message("No genes for ", set_name, " / ", candidate_group, ". Skipping PDF.")
      next
    }
    
    pdf(
      file = file.path(outdir, paste0(set_name, "_", candidate_group, "_locus_AF.pdf")),
      width = 10,
      height = 4.5,
      onefile = TRUE
    )
    
    for (gene_id in genes_to_plot) {
      
      g <- gene_ranges[gene == gene_id]
      if (nrow(g) == 0) next
      
      if (g$strand == "+") {
        cis_start <- max(1, g$start - cis_upstream)
        cis_end   <- g$end + cis_downstream
        promoter_start <- max(1, g$start - cis_upstream)
        promoter_end   <- g$start - 1
      } else {
        cis_start <- max(1, g$start - cis_downstream)
        cis_end   <- g$end + cis_upstream
        promoter_start <- g$end + 1
        promoter_end   <- g$end + cis_upstream
      }
      
      region <- af[CHROM == g$CHROM & POS >= cis_start & POS <= cis_end]
      if (nrow(region) == 0) next
      
      region <- merge(region, cand, by = c("CHROM", "POS"), all.x = TRUE)
      region[is.na(candidate), candidate := "None"]
      
      region[, F0 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^F0_r", names(region), value = TRUE)]
      region[, Cold_F90 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^Cold_F90_r", names(region), value = TRUE)]
      region[, Hot_F90 := rowMeans(.SD, na.rm = TRUE), .SDcols = grep("^Hot_F90_r", names(region), value = TRUE)]
      
      ## Polarize to rising allele on average
      region[, evolved_mean := (Cold_F90 + Hot_F90) / 2]
      region[evolved_mean < F0, `:=`(
        F0 = 1 - F0,
        Cold_F90 = 1 - Cold_F90,
        Hot_F90 = 1 - Hot_F90
      )]

      # ## Polarize to the minor allele in the ancestral state
      # ## If ancestral AF > 0.5, flip all populations to 1 - AF
      # region[F0 > 0.5, `:=`(
      #   F0 = 1 - F0,
      #   Cold_F90 = 1 - Cold_F90,
      #   Hot_F90 = 1 - Hot_F90
      # )]
      # 
      xpos <- data.table(
        population = c("F0", "Cold_F90", "Hot_F90"),
        xoff = c(-6, 0, 6)
      )
      
      plot_dt <- melt(
        region[, .(CHROM, POS, candidate, F0, Cold_F90, Hot_F90)],
        id.vars = c("CHROM", "POS", "candidate"),
        variable.name = "population",
        value.name = "AF"
      )
      
      plot_dt <- merge(plot_dt, xpos, by = "population")
      plot_dt[, POS_plot := POS + xoff]
      
      plot_dt[, population := factor(population, levels = c("F0", "Cold_F90", "Hot_F90"))]
      plot_dt[, candidate := factor(candidate, levels = c("None", "Cold", "Hot", "Shared"))]
      
      plot_dt[candidate == "None", color_group := "None"]
      plot_dt[candidate != "None" & population == "F0", color_group := "F0"]
      plot_dt[candidate != "None" & population == "Cold_F90", color_group := "Cold_F90"]
      plot_dt[candidate != "None" & population == "Hot_F90", color_group := "Hot_F90"]
      plot_dt[, color_group := factor(color_group, levels = c("None", "F0", "Cold_F90", "Hot_F90"))]
      
      exons <- unique(
        gtf[type == "exon" & gene == gene_id, .(start = as.integer(start), end = as.integer(end))]
      )
      
      if (nrow(exons) > 0) {
        setorder(exons, start, end)
        exons[, run_end := cummax(end)]
        exons[, prev_end := c(-Inf, head(run_end, -1))]
        exons[, grp := cumsum(start > prev_end + 1)]
        exons <- exons[, .(start = min(start), end = max(end)), by = grp]
        exons[, c("grp", "run_end", "prev_end") := NULL]
      }
      
      gene_label <- if (is.na(g$symbol) || g$symbol == gene_id) {
        gene_id
      } else {
        paste0(gene_id, " (", g$symbol, ")")
      }
      
      track_ymin <- -0.12
      track_ymax <- -0.08
      gene_ymin  <- -0.108
      gene_ymax  <- -0.092
      
      p <- ggplot() +
        
        annotate("rect",
                 xmin = cis_start, xmax = cis_end,
                 ymin = track_ymin, ymax = track_ymax,
                 fill = "grey88", color = NA) +
        
        annotate("rect",
                 xmin = promoter_start, xmax = promoter_end,
                 ymin = track_ymin, ymax = track_ymax,
                 fill = "khaki", color = NA) +
        
        annotate("rect",
                 xmin = g$start, xmax = g$end,
                 ymin = gene_ymin, ymax = gene_ymax,
                 fill = "grey30", color = NA) +
        
        geom_rect(
          data = exons,
          aes(xmin = start, xmax = end, ymin = track_ymin, ymax = track_ymax),
          inherit.aes = FALSE,
          fill = "black",
          color = NA
        ) +
        
        geom_point(
          data = plot_dt[candidate == "None"],
          aes(x = POS_plot, y = AF, shape = candidate),
          color = "grey70",
          size = 2.0,
          alpha = 0.9
        ) +
        
        geom_point(
          data = plot_dt[candidate != "None"],
          aes(x = POS_plot, y = AF, color = color_group, shape = candidate),
          size = 2.6,
          alpha = 0.95
        ) +
        
        scale_color_manual(
          values = c(
            "None" = "grey70",
            "F0" = "forestgreen",
            "Cold_F90" = "steelblue",
            "Hot_F90" = "maroon"
          )
        ) +
        
        scale_shape_manual(
          values = c(
            "None" = 16,
            "Cold" = 17,
            "Hot" = 15,
            "Shared" = 18
          )
        ) +
        
        scale_y_continuous(
          limits = c(-0.14, 1.02),
          breaks = c(0, 0.25, 0.5, 0.75, 1)
        ) +
        
        theme_bw(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          legend.position = "right"
        ) +
        
        labs(
          title = paste0(gene_label, " | ", g$CHROM, ":", cis_start, "-", cis_end),
          subtitle = paste0(
            "Gene body: ", g$start, "-", g$end,
            " | strand: ", g$strand,
            " | exons = black, intron span = dark grey, promoter-side cis = khaki"
          ),
          x = paste0("Genomic position on ", g$CHROM),
          y = "Polarized allele frequency",
          color = "Population",
          shape = "Candidate status"
        )
      
      print(p)
    }
    
    dev.off()
  }
}

cat("Done. Outputs saved to:\n", outdir, "\n")