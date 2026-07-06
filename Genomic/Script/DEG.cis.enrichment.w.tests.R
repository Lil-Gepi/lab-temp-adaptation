## DEG cis enrichment test ----
rm(list = ls())
gc()

library(data.table)

## 1. Paths ----
genomic_dir <- "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Genomic/Script"
rna_dir     <- "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Script"
gtf_file <- normalizePath(path.expand("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Transcriptomic/Data/dsimM252v1.2.gtf"))
outdir      <- "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Genomic/Result/DEG_cis_enrichment"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


## 2. Load data ----
af     <- readRDS(file.path(genomic_dir, "af.covfilter.RDS"))
p.cold <- readRDS(file.path(genomic_dir, "p.cmh.cold.RDS"))
p.hot  <- readRDS(file.path(genomic_dir, "p.cmh.hot.RDS"))
short.list <- readRDS(file.path(rna_dir, "short.list.RDS"))

setDT(af); setDT(p.cold); setDT(p.hot); setDT(short.list)

setnames(p.cold, names(p.cold)[3], "q_cold")
setnames(p.hot,  names(p.hot)[3],  "q_hot")

af[, CHROM := as.character(CHROM)]
p.cold[, CHROM := as.character(CHROM)]
p.hot[, CHROM := as.character(CHROM)]


## 3. Candidate SNP table ----
snp <- merge(
  unique(af[, .(CHROM, POS)]),
  p.cold[, .(CHROM, POS, q_cold)],
  by = c("CHROM", "POS"),
  all.x = TRUE
)

snp <- merge(
  snp,
  p.hot[, .(CHROM, POS, q_hot)],
  by = c("CHROM", "POS"),
  all.x = TRUE
)

snp[, cold_cand := !is.na(q_cold) & q_cold < 0.05]
snp[, hot_cand  := !is.na(q_hot)  & q_hot  < 0.05]
snp[, `:=`(snp_start = POS, snp_end = POS)]


## 4. Read GTF and get gene ranges ----
gtf <- fread(
  cmd = paste("grep -v '^#'", shQuote(gtf_file)),
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE
)

setnames(gtf, c("CHROM", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

gtf[, CHROM := sub("^chr", "", as.character(CHROM))]
gtf[, gene := sub('.*gene_id "([^"]+)".*', "\\1", attributes)]
gtf <- gtf[CHROM %in% c("X", "2L", "2R", "3L", "3R") & grepl("^FBgn", gene)]

gene_ranges <- gtf[, .(
  CHROM = unique(CHROM)[1],
  strand = unique(strand)[1],
  start = min(as.integer(start)),
  end = max(as.integer(end))
), by = gene]


## 5. Expressed gene universe and DEG labels ----
genes <- short.list[, .(
  fdr.cold = min(fdr.cold, na.rm = TRUE),
  fdr.hot  = min(fdr.hot,  na.rm = TRUE)
), by = gene]

genes[is.infinite(fdr.cold), fdr.cold := NA_real_]
genes[is.infinite(fdr.hot),  fdr.hot  := NA_real_]

genes <- merge(genes, gene_ranges, by = "gene")
genes[, cold_DEG := !is.na(fdr.cold) & fdr.cold < 0.05]
genes[, hot_DEG  := !is.na(fdr.hot)  & fdr.hot  < 0.05]


## 6. Make cis-only regions, excluding gene body ----
cis_upstream <- 5000
cis_downstream <- 1000

genes[, `:=`(
  left_start  = fifelse(strand == "+", pmax(1, start - cis_upstream), pmax(1, start - cis_downstream)),
  left_end    = start - 1,
  right_start = end + 1,
  right_end   = fifelse(strand == "+", end + cis_downstream, end + cis_upstream)
)]

cis <- rbindlist(list(
  genes[left_start <= left_end, .(gene, CHROM, cold_DEG, hot_DEG, cis_start = left_start, cis_end = left_end)],
  genes[right_start <= right_end, .(gene, CHROM, cold_DEG, hot_DEG, cis_start = right_start, cis_end = right_end)]
))


## 7. Count SNPs per gene cis region ----
setkey(cis, CHROM, cis_start, cis_end)
setkey(snp, CHROM, snp_start, snp_end)

ov <- foverlaps(
  snp, cis,
  by.x = c("CHROM", "snp_start", "snp_end"),
  by.y = c("CHROM", "cis_start", "cis_end"),
  nomatch = 0L
)

cis_count <- ov[, .(
  n_snp_cis  = .N,
  n_cold_cis = sum(cold_cand),
  n_hot_cis  = sum(hot_cand)
), by = .(gene, CHROM, cold_DEG, hot_DEG)]

cis_count <- merge(
  unique(genes[, .(gene, CHROM, cold_DEG, hot_DEG)]),
  cis_count,
  by = c("gene", "CHROM", "cold_DEG", "hot_DEG"),
  all.x = TRUE
)

cis_count[is.na(n_snp_cis), `:=`(n_snp_cis = 0L, n_cold_cis = 0L, n_hot_cis = 0L)]
cis_count[, cold_rate := fifelse(n_snp_cis > 0, n_cold_cis / n_snp_cis, NA_real_)]
cis_count[, hot_rate  := fifelse(n_snp_cis > 0, n_hot_cis  / n_snp_cis, NA_real_)]

fwrite(cis_count, file.path(outdir, "DEG_cis_candidate_counts_per_gene.tsv"), sep = "\t")



## 8.1 Binary gene-level permutation: >=1 cis candidate SNP ----
set.seed(42)
n_perm <- 10000

x <- copy(cis_count[n_snp_cis > 0])
x[, cold_has_cis := n_cold_cis >= 1]
x[, hot_has_cis  := n_hot_cis  >= 1]

res <- list()
null <- list()

for (tt in c("cold", "hot")) {
  
  deg_col <- paste0(tt, "_DEG")
  has_col <- paste0(tt, "_has_cis")
  
  deg <- x[get(deg_col) == TRUE]
  
  obs_n <- sum(deg[[has_col]])
  obs_rate <- obs_n / nrow(deg)
  rand_rate <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    rand <- x[sample.int(.N, nrow(deg), replace = T)]
    rand_rate[i] <- sum(rand[[has_col]]) / nrow(rand)
  }
  
  res[[tt]] <- data.table(
    test = paste0(tt, " DEGs with >=1 ", tt, " candidate SNP in cis"),
    n_DEG = nrow(deg),
    observed_genes_with_candidate = obs_n,
    observed_gene_rate = obs_rate,
    mean_random_gene_rate = mean(rand_rate),
    sd_random_gene_rate = sd(rand_rate),
    fold_enrichment = obs_rate / mean(rand_rate),
    p_empirical = (sum(rand_rate >= obs_rate) + 1) / (n_perm + 1)
  )
  
  null[[tt]] <- data.table(
    test = paste0(tt, " DEGs with >=1 ", tt, " candidate SNP in cis"),
    permutation = seq_len(n_perm),
    random_gene_rate = rand_rate,
    observed_gene_rate = obs_rate
  )
}

res <- rbindlist(res)
res[, p_BH := p.adjust(p_empirical, method = "BH")]
null <- rbindlist(null)

fwrite(res, file.path(outdir, "DEG_cis_candidate_binary_gene_enrichment.tsv"), sep = "\t")
fwrite(null, file.path(outdir, "DEG_cis_candidate_binary_gene_enrichment_null.tsv"), sep = "\t")
print(res)

plot_lab <- res[, .(
  test,
  observed_gene_rate,
  label = paste0(
    "obs = ", round(observed_gene_rate, 3),
    "\nnull = ", round(mean_random_gene_rate, 3),
    "\nBH p = ", signif(p_BH, 3)
  )
)]

(p1 <- ggplot(null, aes(x = random_gene_rate)) +
  geom_histogram(bins = 53, color = "white") +
  geom_vline(data = plot_lab, aes(xintercept = observed_gene_rate), linewidth = 0.8) +
  geom_text(data = plot_lab, aes(x = observed_gene_rate, y = Inf, label = label),
            vjust = 1.2, hjust = -0.05, size = 3.2) +
  facet_wrap(~ test, scales = "free_x") +
  theme_classic(base_size = 12) +
  labs(
    x = "Random gene-set rate",
    y = "Permutation count",
    title = "DEG cis enrichment permutation test: >=1 candidate SNP"
  )
)
ggsave( "../Plot/DEG_cis_candidate_binary_gene_enrichment_null.png",
       p1, width = 8, height = 4.5)

## 8.2 Binary gene-level permutation: >=2 cis candidate SNPs ----
set.seed(44)
n_perm <- 10000

x <- copy(cis_count[n_snp_cis > 0])
x[, cold_has_cis := n_cold_cis >= 2]
x[, hot_has_cis  := n_hot_cis  >= 2]

res <- list()
null <- list()

for (tt in c("cold", "hot")) {
  
  deg_col <- paste0(tt, "_DEG")
  has_col <- paste0(tt, "_has_cis")
  
  deg <- x[get(deg_col) == TRUE]
  
  obs_n <- sum(deg[[has_col]])
  obs_rate <- obs_n / nrow(deg)
  rand_rate <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    rand <- x[sample.int(.N, nrow(deg), replace = T)]
    rand_rate[i] <- sum(rand[[has_col]]) / nrow(rand)
  }
  
  res[[tt]] <- data.table(
    test = paste0(tt, " DEGs with >=2 ", tt, " candidate SNPs in cis"),
    n_DEG = nrow(deg),
    observed_genes_with_candidate = obs_n,
    observed_gene_rate = obs_rate,
    mean_random_gene_rate = mean(rand_rate),
    sd_random_gene_rate = sd(rand_rate),
    fold_enrichment = obs_rate / mean(rand_rate),
    p_empirical = (sum(rand_rate >= obs_rate) + 1) / (n_perm + 1)
  )
  
  null[[tt]] <- data.table(
    test = paste0(tt, " DEGs with >=2 ", tt, " candidate SNPs in cis"),
    permutation = seq_len(n_perm),
    random_gene_rate = rand_rate,
    observed_gene_rate = obs_rate
  )
}

res <- rbindlist(res)
res[, p_BH := p.adjust(p_empirical, method = "BH")]
null <- rbindlist(null)

fwrite(res, file.path(outdir, "DEG_cis_candidate_binary_gene_enrichment_2SNPs.tsv"), sep = "\t")
fwrite(null, file.path(outdir, "DEG_cis_candidate_binary_gene_enrichment_2SNPs_null.tsv"), sep = "\t")
print(res)

plot_lab <- res[, .(
  test,
  observed_gene_rate,
  label = paste0(
    "obs = ", round(observed_gene_rate, 3),
    "\nnull = ", round(mean_random_gene_rate, 3),
    "\nBH p = ", signif(p_BH, 3)
  )
)]

(p2 <- ggplot(null, aes(x = random_gene_rate)) +
  geom_histogram(bins = 45, color = "white") +
  geom_vline(data = plot_lab, aes(xintercept = observed_gene_rate), linewidth = 0.8) +
  geom_text(data = plot_lab, aes(x = observed_gene_rate, y = Inf, label = label),
            vjust = 1.2, hjust = -0.05, size = 3.2) +
  facet_wrap(~ test, scales = "free_x") +
  theme_classic(base_size = 12) +
  labs(
    x = "Random gene-set rate",
    y = "Permutation count",
    title = "DEG cis enrichment permutation test: >=2 candidate SNPs"
  )
)
ggsave("../Plot/DEG_cis_candidate_binary_gene_enrichment_2SNPs_null.png",
       p2, width = 8, height = 4.5)
