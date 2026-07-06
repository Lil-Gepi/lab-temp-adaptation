## load in packages----
setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Genomic/Script/")
rm(list = ls())
library(data.table)
setDTthreads(percent = 80)
getDTthreads()
library(matrixStats)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggExtra)
source("readBCF.R")
source("RemoveMulti.R")

## 0. Data Prep----
af <- read_bcf("../Data/Dsim.FL.ACH.piled.AF1090filtered.vcf.gz", what = "AF")
xf <- read_bcf("../Data/Dsim.FL.ACH.piled.AF1090filtered.vcf.gz", what = "XF")
dp <- read_bcf("../Data/Dsim.FL.ACH.piled.AF1090filtered.vcf.gz", what = "DP")
af[, 6:175 := lapply(.SD, as.numeric), .SDcols = 6:175]
af <- prune_multi(af,dp,how="most_freq")
xf <- prune_multi(xf,dp,how="most_freq")

colnames(af)

colnames(af)[6:175] <- c(
  paste0("F0_r",1:10),
  paste0("Cold_F", rep(c(10, 20, 30, 40, 50, 60, 70, 90), each = 10),"_r", rep(1:10, times = 8) ),
  paste0("Hot_F", rep(c(10, 20, 30, 40, 50, 60, 70, 90), each = 10),"_r", rep(1:10, times = 8) )
)

colnames(xf)[6:175] <- c(
  paste0("F0_r",1:10),
  paste0("Cold_F", rep(c(10, 20, 30, 40, 50, 60, 70, 90), each = 10),"_r", rep(1:10, times = 8) ),
  paste0("Hot_F", rep(c(10, 20, 30, 40, 50, 60, 70, 90), each = 10),"_r", rep(1:10, times = 8) )
)

colnames(dp)[3:172] <- c(
  paste0("F0_r",1:10),
  paste0("Cold_F", rep(c(10, 20, 30, 40, 50, 60, 70, 90), each = 10),"_r", rep(1:10, times = 8) ),
  paste0("Hot_F", rep(c(10, 20, 30, 40, 50, 60, 70, 90), each = 10),"_r", rep(1:10, times = 8) )
)

setkey(af, CHROM, POS )
setkey(xf, CHROM, POS)
setkey(dp, CHROM, POS)

saveRDS(af, file = "./af.RDS")
saveRDS(xf, file = "./xf.RDS")
saveRDS(dp, file = "./dp.RDS")
## filter out lowly covered SNPs
af <- readRDS("af.RDS")
xf <- readRDS("xf.RDS")
dp <- readRDS("dp.RDS")

(low.threshold <- dp[, sapply(.SD, quantile, probs = 0.01, na.rm = TRUE), .SDcols = 3:172])
keep <- rep(TRUE, nrow(dp))
for (sample in colnames(dp[,3:172])) {
  sample.threshold <- low.threshold[[paste0(sample,".1%")]]
  value   <- dp[,sample, with = F]
  keep <- keep & !is.na(value) & value >= sample.threshold
}

(dp[,.N] - sum(keep) )/dp[,.N]

dp <- dp[as.vector(keep)]
af <- af[as.vector(keep)]
xf <- xf[as.vector(keep)]
rm(value, keep, low.threshold, sample, sample.threshold, prune_multi, read_bcf)
saveRDS(af, file = "./af.covfilter.RDS")
saveRDS(xf, file = "./xf.covfilter.RDS")
saveRDS(dp, file = "./dp.covfilter.RDS")

## 1. PCA ----
af <- readRDS("af.covfilter.RDS")
xf <- readRDS("xf.covfilter.RDS")
dp <- readRDS("dp.covfilter.RDS")

### 1.1 PCA on hot & cold ----
pca.dat <- af[,6:175]
setDT(pca.dat)
pca.dat[, var_non_zero := rowVars(as.matrix(.SD)), .SDcols = names(pca.dat)]

pca.dat <- pca.dat[var_non_zero !=0, 1:170]
pca.dat <- t(as.matrix(2*asin(sqrt(pca.dat))))
pca.dat<-na.omit(pca.dat)#remove missing data

pca.res <- prcomp(pca.dat,retx = T, center = T, scale. =T)

generation <- c(rep(c(0,10,20,30,40,50,60,70,90), each = 10),
                rep(c(  10,20,30,40,50,60,70,90), each = 10))
pop <- c(rep("Anc", 10), rep("Cold", 80),rep("Hot", 80))
repl <- c(rep(1:10, times = 17))

pca_df <- as.data.frame(pca.res$x)
pca_df$generation <- factor(generation)   # make it a column (not just in the parent env)
pca_df$pop        <- factor(pop, levels = c("Anc", "Cold","Hot"))
pca_df$repl       <- factor(repl)

additional_shapes <- c(0:4, 6:10)

ve_pct <- round(100 * (pca.res$sdev^2 / sum(pca.res$sdev^2)), 1)

gen_levels <- levels(pca_df$generation)
alpha_vals <- seq(0.3, 1, length.out = length(gen_levels))
names(alpha_vals) <- gen_levels

plot_pcs <- function(i, j) {
  ggplot(pca_df,
         aes_string(x = paste0("PC", i), y = paste0("PC", j),
                    color = "pop", shape = "repl", alpha = "generation")
  ) +
    geom_point(size = 2.6) +
    scale_shape_manual(values = additional_shapes) +
    scale_color_manual(values = c("Anc" = "forestgreen", "Cold" = "steelblue", "Hot" = "maroon")) +
    scale_alpha_manual(values = alpha_vals) +
    theme_minimal() +
    labs(
      x = sprintf("PC%d (%.1f%%)", i, ve_pct[i]),
      y = sprintf("PC%d (%.1f%%)", j, ve_pct[j]),
      color = "Population",
      alpha = "Generation",
      shape = "Replicate",
      title = sprintf("PCA: PC%d vs PC%d", i, j)
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10),
      axis.title   = element_text(size = 13),
      axis.text    = element_text(size = 11),
      plot.title   = element_text(size = 14, face = "bold")
    )
}
#### PC1 vs the others ...
k <- 12  
pdf("../Plot/pca_PC1_vs_others.pdf", width = 7, height = 7)
for (j in 2:k) {
  print(plot_pcs(1, j))
}
dev.off()

#### Consecutive pairs: (PC1,PC2), (PC3,PC4), ...
k_pairs <- 12  
pair_idx <- split(1:k_pairs, ceiling((1:k_pairs)/2))  # rough helper; we’ll filter to pairs
pair_idx <- lapply(seq(1, k_pairs, by = 2), function(s) c(s, s + 1))

pdf("../Plot/pca_consecutive_pairs.pdf", width = 7, height = 7)
for (pair in pair_idx) {
  if (max(pair) <= length(ve_pct)) {
    print(plot_pcs(pair[1], pair[2]))
  }
}
dev.off()

pca.ach.pc12 <- plot_pcs(1, 2)

### 1.2 PCA on hot alone ----
pca.dat <- af[,c(6:15, 96:175)]
setDT(pca.dat)
pca.dat[, var_non_zero := rowVars(as.matrix(.SD)), .SDcols = names(pca.dat)]

pca.dat <- pca.dat[var_non_zero !=0, 1:90]
pca.dat <- t(as.matrix(2*asin(sqrt(pca.dat))))
pca.dat<-na.omit(pca.dat)#remove missing data

pca.res <- prcomp(pca.dat,retx = T, center = T, scale. =T)

generation <- c(rep(c(0,10,20,30,40,50,60,70,90), each = 10))
pop <- c(rep("Anc", 10), rep("Hot", 80))
repl <- c(rep(1:10, times = 9))

pca_df <- as.data.frame(pca.res$x)
pca_df$generation <- factor(generation)   # make it a column (not just in the parent env)
pca_df$pop        <- factor(pop, levels = c("Anc", "Hot"))
pca_df$repl       <- factor(repl)

additional_shapes <- c(0:4, 6:10)

ve_pct <- round(100 * (pca.res$sdev^2 / sum(pca.res$sdev^2)), 1)

gen_levels <- levels(pca_df$generation)
alpha_vals <- seq(0.3, 1, length.out = length(gen_levels))
names(alpha_vals) <- gen_levels

pca.ah.pc12 <- plot_pcs(1, 2)

### 1.3 PCA on cold alone ----
pca.dat <- af[,c(6:95)]
setDT(pca.dat)
pca.dat[, var_non_zero := rowVars(as.matrix(.SD)), .SDcols = names(pca.dat)]

pca.dat <- pca.dat[var_non_zero !=0, 1:90]
pca.dat <- t(as.matrix(2*asin(sqrt(pca.dat))))
pca.dat<-na.omit(pca.dat)

pca.res <- prcomp(pca.dat,retx = T, center = T, scale. =T)

generation <- c(rep(c(0,10,20,30,40,50,60,70,90), each = 10))
pop <- c(rep("Anc", 10), rep("Cold", 80))
repl <- c(rep(1:10, times = 9))

pca_df <- as.data.frame(pca.res$x)
pca_df$generation <- factor(generation)   # make it a column (not just in the parent env)
pca_df$pop        <- factor(pop, levels = c("Anc", "Cold"))
pca_df$repl       <- factor(repl)

additional_shapes <- c(0:4, 6:10)

ve_pct <- round(100 * (pca.res$sdev^2 / sum(pca.res$sdev^2)), 1)

gen_levels <- levels(pca_df$generation)
alpha_vals <- seq(0.3, 1, length.out = length(gen_levels))
names(alpha_vals) <- gen_levels

pca.ac.pc12 <- plot_pcs(1, 2)


save(list = c("pca.ach.pc12", "pca.ac.pc12", "pca.ah.pc12"), file = "pca.PC12.RDS")

## clean up 
rm(additional_shapes, alpha_vals, gen_levels, j, k, k_pairs, pair, ve, ve_pct, pair_idx, pca_df, pca.dat, pca.res)
load("pca.PC12.RDS")
### 1.4 plot three PCA plots together ----
library(cowplot)

theme_set(theme_cowplot())

shared_legend <- get_legend(
  pca.ach.pc12 + theme(legend.position = "right")
)

pca.ach.pc12_noleg <- pca.ach.pc12 + theme(legend.position = "none") + ggtitle("")
pca.ac.pc12_noleg  <- pca.ac.pc12  + theme(legend.position = "none")+ ggtitle("")
pca.ah.pc12_noleg  <- pca.ah.pc12  + theme(legend.position = "none")+ ggtitle("")

row2 <- plot_grid(
  pca.ac.pc12_noleg,
  pca.ah.pc12_noleg,
  ncol = 2,
  align = "hv",
  axis  = "tblr",
  labels = c("B", "C")
)

main_plot <- plot_grid(
  pca.ach.pc12_noleg,
  row2,
  ncol = 1,
  rel_heights = c(1.75, 1),
  labels = c("A", "")
)

final_figure <- plot_grid(
  main_plot,
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.15)
)

final_figure
ggsave("../Plot/Figure3.PCA.png", final_figure, width = 10, height = 8, dpi = 300)

rm(generation, pop, repl, plot_pcs,pca.ac.pc12, pca.ac.pc12_noleg, pca.ach.pc12, pca.ach.pc12_noleg, pca.ah.pc12, pca.ah.pc12_noleg, row2, shared_legend, main_plot, final_figure)

## 2. Ne estimation ----
library(poolSeq)
xf.start.end <- xf[,c(1:15,86:95,166:175)]
setkey(xf.start.end, CHROM, POS)
dp.start.end <- dp[,c(1:12,83:92,163:172)]
setkey(dp.start.end, CHROM, POS)


ne_estimates <- NULL
nb_SNPs <- 500
nb_rounds <- 1000
times <- c(0,90)

summary(dp.start.end)

for (chrom in c("X","2L","2R","3L","3R")) {
  for(j in seq_len(nb_rounds)){
    set.seed(seed = j)
    cov_trial <- sample_n(dp.start.end[CHROM == chrom, ], size = nb_SNPs)
    (SNP <- cov_trial[,.(CHROM,POS)])
    setkey(SNP, CHROM, POS)
    (af_trial <- xf.start.end[CHROM== chrom& POS %in% SNP$POS,])
    for (evo in c("Cold","Hot")) {
      for (replicate in 1:10) {
        pref.i <- paste0("F0_r", replicate) #time point i
        pref.j <- paste0(evo,"_F",times[2],"_r", replicate) #time point i
        
        pi <- unlist(af_trial[,pref.i, with =F])
        pj <- unlist(af_trial[,pref.j, with =F])
        covi <- unlist(cov_trial[,pref.i, with =F])
        covj <- unlist(cov_trial[,pref.j, with =F])
        ne <- estimateNe(p0 = pi, pt = pj, cov0 = covi, covt = covj, t = times[2]-times[1], 
                         ploidy = 2, truncAF = 0.1, method = "P.alt.2step.planI", poolSize = c(625, 625), Ncensus = 1250)
        ne_estimates <- rbind(ne_estimates, data.frame(chrom = chrom, evo = evo, replicate = replicate, start = times[1], end = times[2], trial = j,  ne = ne))
      }
    }
  }
}

ne_estimates <- na.omit(ne_estimates)
is.data.table(ne_estimates)
setDT(ne_estimates)
ne_median <- ne_estimates[,lapply(.SD, FUN = median), by = .(chrom, evo, replicate, start, end)]
ne_median[,`:=`(trial = NULL, ne = round(ne+0.5))]

## quick plot
png(filename = "../Plot/Ne.Replicate.png", width = 7, height = 4, units = "in", res =350)
ggplot(data = ne_median, aes(x = evo, y = ne, color = evo, label = replicate)) +
  geom_text(size = 3, position = position_jitter(width = 0.4,seed = 42)) +
  facet_wrap(~chrom, nrow = 1) +
  scale_color_manual(values = c("Hot" = "maroon", "Cold" = "steelblue")) +
  xlab("") + 
  ylab("Ne estimate") +
  theme_bw()+
  theme(legend.position = "")
dev.off()

png(filename = "../Plot/Ne.png", width = 7, height = 4, units = "in", res =350)
ggplot(data = ne_median, aes(x = evo, y = ne, color = evo, label = replicate)) +
  geom_boxplot()+
  geom_jitter(position = position_jitter(width = 0.4, seed = 42), size = 2)+
  facet_wrap(~chrom, nrow = 1) +
  scale_color_manual(values = c("Hot" = "maroon", "Cold" = "steelblue")) +
  xlab("") + 
  ylab("Ne estimate") +
  theme_bw()+
  theme(legend.position = "")
dev.off()


saveRDS(ne_estimates, "ne_estimate.RDS")

rm(af_trial, cov_trial, dp.start.end,SNP, xf.start.end, ne_estimates, ne_median, chrom, covj, covi, evo, j, nb_rounds, nb_SNPs, ne, pi, pj, pref.i, pref.j, replicate, times)
## 3. CMH ----
ne_estimates <- readRDS("ne_estimate.RDS")

ne_median <- ne_estimates[,lapply(.SD, FUN = median), by = .(chrom, evo, replicate, start, end)]
ne_median[,`:=`(trial = NULL, ne = round(ne+0.5))]

library(ACER)
(focal_cols <- c(paste0("F0_r", 1:10), paste0("Cold_F", rep(c(1:7*10, 90), each = 10), "_r", rep(1:10, times = 8))))
af.cold <- as.matrix(xf[,focal_cols, with = F])
dp.cold <- as.matrix(dp[,focal_cols, with = F])
Ne <- ne_median[evo == "Cold" , .(med = median(ne)), by = replicate]
p.cmh.cold <- adapted.cmh.test(freq = af.cold,coverage = dp.cold,gen = c(0:7*10L, 90),mincov = 10,repl = 1:10,Ne = Ne$med,poolSize = rep(625, 90),IntGen = T,order = 1)
p.cmh.cold <- p.adjust(p.cmh.cold, method = "fdr")
rm(af.cold,dp.cold, Ne)
p.cmh.cold <- cbind(xf[,.(CHROM,POS)], p.cmh.cold)

(focal_cols <- c(paste0("F0_r", 1:10), paste0("Hot_F", rep(c(1:7*10, 90), each = 10), "_r", rep(1:10, times = 8))))
af.hot <- as.matrix(xf[,focal_cols, with = F])
dp.hot <- as.matrix(dp[,focal_cols, with = F])
Ne <- ne_median[evo == "Hot" , .(med = median(ne)), by = replicate]
p.cmh.hot <- adapted.cmh.test(freq = af.hot,coverage = dp.hot,gen = c(0:7*10L, 90),mincov = 10,repl = 1:10,Ne = Ne$med,poolSize = rep(625, 90),IntGen = T,order = 1)
p.cmh.hot <- p.adjust(p.cmh.hot, method = "fdr")
rm(af.hot,dp.hot, Ne)
p.cmh.hot <- cbind(xf[,.(CHROM,POS)], p.cmh.hot)

saveRDS(p.cmh.cold, "p.cmh.cold.RDS")
saveRDS(p.cmh.hot, "p.cmh.hot.RDS")


# Below is plotting ----
xf <- readRDS("xf.covfilter.RDS")
p.cmh.cold <- readRDS("p.cmh.cold.RDS")
p.cmh.hot  <- readRDS("p.cmh.hot.RDS")

setDT(xf)
setkeyv(xf, c("CHROM", "POS"))

cold.090.target <- p.cmh.cold[p.cmh.cold < 0.05, .(CHROM, POS)]
hot.090.target  <- p.cmh.hot[p.cmh.hot < 0.05, .(CHROM, POS)]

library(grid)
library(cowplot)
library(ggVennDiagram)
library(ggExtra)

theme_bigfig <- theme_minimal(base_size = 18, base_family = "sans") +
  theme(
    plot.title      = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title      = element_text(size = 22),
    axis.text       = element_text(size = 18, colour = "black"),
    strip.text      = element_text(size = 18, face = "plain"),
    legend.title    = element_text(size = 18, face = "bold"),
    legend.text     = element_text(size = 18),
    plot.caption    = element_text(size = 16, hjust = 0.5),
    legend.key.size = unit(0.8, "cm"),
    panel.spacing   = unit(0.7, "lines"),
    panel.grid.minor = element_blank()
  )


cold.snp <- cbind(
  cold.090.target,
  cold.afc = apply(xf[cold.090.target, 86:95, with = FALSE], 1, mean) -
    apply(xf[cold.090.target, 6:15, with = FALSE], 1, mean),
  hot.afc = apply(xf[cold.090.target, 166:175, with = FALSE], 1, mean) -
    apply(xf[cold.090.target, 6:15, with = FALSE], 1, mean)
) |> as.data.table()


hot.snp <- cbind(
  hot.090.target,
  cold.afc = apply(xf[hot.090.target, 86:95, with = FALSE], 1, mean) -
    apply(xf[hot.090.target, 6:15, with = FALSE], 1, mean),
  hot.afc = apply(xf[hot.090.target, 166:175, with = FALSE], 1, mean) -
    apply(xf[hot.090.target, 6:15, with = FALSE], 1, mean)
) |> as.data.table()
## 5. Venn diagram ----
hot_up  <- unique(unlist(hot.snp[hot.afc > 0, .(paste0(CHROM, "_", POS))]))
hot_dn  <- unique(unlist(hot.snp[hot.afc < 0, .(paste0(CHROM, "_", POS))]))
cold_up <- unique(unlist(cold.snp[cold.afc > 0, .(paste0(CHROM, "_", POS))]))
cold_dn <- unique(unlist(cold.snp[cold.afc < 0, .(paste0(CHROM, "_", POS))]))

venn_obj <- list(
  "Hot Increasing"  = hot_up,
  "Hot Decreasing"  = hot_dn,
  "Cold Increasing" = cold_up,
  "Cold Decreasing" = cold_dn
)

p_venn <- ggVennDiagram(
  x = venn_obj,
  label = "both",
  label_color = "grey5",
  set_color = c("maroon", "maroon", "steelblue", "steelblue"),
  set_size = 5,label_size = 5,
  label_alpha = 0.15
) +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  theme_void(base_family = "sans") +
  coord_fixed(clip = "off") +
  theme(
    text            = element_text(size = 20),
    legend.position = "none",
    plot.margin     = margin(30, 30, 30, 30)
  )

png("../Plot/Figure4.VennDiagram.png", width = 8, height = 8, units = "in", res = 600)
print(p_venn)
dev.off()

## 4. AFC rank plot ----
cold.snp[cold.afc < 0, `:=`(
  cold.afc = abs(cold.afc),
  hot.afc  = -1 * hot.afc
)]


hot.snp[hot.afc < 0, `:=`(
  hot.afc  = abs(hot.afc),
  cold.afc = -1 * cold.afc
)]

## 4.1 cold ----
ranks_cold <- cold.snp[order(cold.afc, decreasing = TRUE), .(CHROM, POS)][, rank := .I]

cold.snp.long <- melt(
  cold.snp,
  id.vars = c("CHROM", "POS"),
  measure.vars = c("cold.afc", "hot.afc"),
  variable.name = "Source",
  value.name = "AFC"
)

cold.snp.long <- merge(cold.snp.long, ranks_cold, by = c("CHROM", "POS"))
cold.snp.long[, Source := factor(Source, levels = c("hot.afc", "cold.afc"))]
setorder(cold.snp.long, Source)

c <- ggplot(cold.snp.long, aes(x = rank, y = AFC, color = Source)) +
  geom_point(size = 1.5, alpha = 0.1) +
  scale_color_manual(
    values = c(cold.afc = "steelblue", hot.afc = "maroon"),
    labels = c(cold.afc = "AFC in Cold (focal)", hot.afc = "AFC in Hot")
  ) +
  labs(
    x = "SNP rank by mean cold AFC",
    y = "Mean AFC across replicates",
    color = NULL
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  coord_cartesian(ylim = c(-0.3, 0.6)) +
  theme_bigfig +
  theme(
    legend.position      = c(0.85, 1),
    legend.justification = c("right", "top"),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    axis.line            = element_blank(),
    plot.background      = element_rect(fill = "white", color = NA),
    plot.margin          = margin(12, 12, 25, 12)
  )

p_cold <- ggMarginal(
  c,
  type = "density",
  margins = "y",
  groupColour = TRUE,
  groupFill = TRUE,
  size = 5
)

png("../Plot/Figure4.AFCrank.Cold.png", width = 8, height = 5, units = "in", res = 600)
print(p_cold)
dev.off()

## 4.2 hot ----
ranks_hot <- hot.snp[order(hot.afc, decreasing = TRUE), .(CHROM, POS)][, rank := .I]

hot.snp.long <- melt(
  hot.snp,
  id.vars = c("CHROM", "POS"),
  measure.vars = c("cold.afc", "hot.afc"),
  variable.name = "Source",
  value.name = "AFC"
)

hot.snp.long <- merge(hot.snp.long, ranks_hot, by = c("CHROM", "POS"))
hot.snp.long[, Source := factor(Source, levels = c("cold.afc", "hot.afc"))]
setorder(hot.snp.long, Source)

h <- ggplot(hot.snp.long, aes(x = rank, y = AFC, color = Source)) +
  geom_point(size = 1.5, alpha = 0.1) +
  scale_color_manual(
    values = c(hot.afc = "maroon", cold.afc = "steelblue"),
    labels = c(cold.afc = "AFC in Cold", hot.afc = "AFC in Hot (focal)")
  ) +
  labs(
    x = "SNP rank by mean hot AFC",
    y = "Mean AFC across replicates",
    color = NULL
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  coord_cartesian(ylim = c(-0.3, 0.6)) +
  theme_bigfig +
  theme(
    legend.position      = c(0.85, 1),
    legend.justification = c("right", "top"),
    panel.grid.major     = element_blank(),
    panel.border         = element_blank(),
    axis.line            = element_blank(),
    plot.background      = element_rect(fill = "white", color = NA),
    plot.margin          = margin(12, 12, 25, 12)
  )

p_hot <- ggMarginal(
  h,
  type = "density",
  margins = "y",
  groupColour = TRUE,
  groupFill = TRUE,
  size = 5
)

png("../Plot/Figure4.AFCrank.Hot.png", width = 8, height = 5, units = "in", res = 600)
print(p_hot)
dev.off()



## 6. Correlation matrix ----
base_cols <- 6:15
cold_cols <- 86:95
hot_cols  <- 166:175

target_all <- unique(rbind(cold.090.target, hot.090.target))
setkey(target_all, CHROM, POS)

base_mat <- as.matrix(xf[target_all, ..base_cols])
cold_mat <- as.matrix(xf[target_all, ..cold_cols])
hot_mat  <- as.matrix(xf[target_all, ..hot_cols])

cold_afc_rep <- cold_mat - base_mat
hot_afc_rep  <- hot_mat - base_mat

colnames(cold_afc_rep) <- paste0("afc.cold.r", 1:ncol(cold_afc_rep))
colnames(hot_afc_rep)  <- paste0("afc.hot.r", 1:ncol(hot_afc_rep))

afc_rep <- cbind(cold_afc_rep, hot_afc_rep)

cor_mat <- cor(afc_rep, use = "pairwise.complete.obs")
cor_dt  <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("Var1", "Var2", "N"))

cor_dt[, Var1 := factor(Var1, levels = colnames(afc_rep))]
cor_dt[, Var2 := factor(Var2, levels = colnames(afc_rep))]
cor_dt[, i := as.integer(Var1)]
cor_dt[, j := as.integer(Var2)]
cor_dt[i > j, N := NA_real_]
cor_dt[, lab := ifelse(is.na(N), NA, sprintf("%.2f", N))]

levels(cor_dt$Var1) <- sub("^afc\\.", "", levels(cor_dt$Var1))
levels(cor_dt$Var2) <- sub("^afc\\.", "", levels(cor_dt$Var2))

k <- ncol(cold_afc_rep)
n <- length(levels(cor_dt$Var2))
x_mid <- k + 0.5
y_mid <- n - k + 0.5

p_cor <- ggplot(cor_dt, aes(x = Var1, y = Var2, fill = N)) +
  geom_tile(color = "white", linewidth = 0.25, na.rm = FALSE) +
  geom_text(aes(label = lab), size = 3, na.rm = TRUE) +
  scale_fill_gradient2(low = "white",mid = "white",high = "red",
    midpoint = 0.2,na.value = "white",limits = c(0, 1)) +
  scale_y_discrete(limits = rev(levels(cor_dt$Var2))) +
  geom_segment(x = x_mid, xend = x_mid, y = 0, yend = y_mid,
    linewidth = 0.8, color = "red", linetype = "dashed") +
  geom_segment(x = 0, xend = x_mid, y = y_mid, yend = y_mid,
    linewidth = 0.8, color = "red", linetype = "dashed") +
  labs(
    title = "Pairwise correlation of AFC between replicate populations",
    x = NULL,y = NULL,fill = NULL) +
  theme_void()+
  theme_bigfig +
  theme(legend.position = "right",
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title  = element_text(size = 18, face = "plain")
  )

png("../Plot/Figure4.CorMatrix.png", width = 9, height = 8, units = "in", res = 600)
print(p_cor)
dev.off()

## 7. Make Figure 4 ----
upper <- plot_grid(
  NULL,NULL,
  p_cold, p_hot,
  labels = c("A", "B", "",""),
  label_size = 24,
  label_fontface = "bold",
  ncol = 2, nrow = 2,
  rel_heights = c(0.1,1)
)
lower <- plot_grid(
  NULL,NULL,
  p_venn, p_cor,
  labels = c("C", "D","",""),
  label_size = 24,
  label_fontface = "bold",
  ncol = 2,  nrow = 2,
  rel_heights = c(0.1,1)
)
mega <- plot_grid(
  upper, lower,
  nrow = 2,
  rel_heights = c(1, 1.1)
)

mega

ggsave(
  "../Plot/Figure4.Mega.png",
  mega,
  width = 16,
  height = 12,
  dpi = 600
)







## S1. Manhattan plot ----
plotMan <- function(res.dt, thres_cmh_genome = 5e-2, alp = 1, top = F, highlightcolor="red") {
  result <- copy(res.dt)
  setnames(result, 1:3, c("CHR","BP","P"))
  chr.name <- data.table(
    oldname = c("2L", "2R", "3L", "3R", "X"),
    newname = c("2",  "3",   "4",  "5", "1")
  )
  setkey(chr.name, "oldname")
  setkey(result, "CHR")
  result[chr.name, CHR := newname]
  result[, CHR := lapply(.SD, as.integer), .SDcols = "CHR"]
  
  don <- result %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP)) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(result, ., by = c("CHR" = "CHR")) %>%
    arrange(CHR, BP) %>%
    mutate(BPcum = BP + tot)
  
  axisdf = don %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2)
  
  resplot <- ggplot(don[P<0.9,], aes(x = BPcum, y = -log10(P))) +
    geom_point(aes(color = as.factor(CHR)), alpha = alp, size = 1) +
    scale_color_manual(values = rep(c("grey60","black"), 4)) +
    scale_x_continuous(name = "Chromosome",label = c("X", "2L", "2R", "3L", "3R"),breaks = axisdf$center,
                       position = ifelse(top,"top", "bottom")) +
    geom_point(data=subset(don, highlight==1), color=highlightcolor,alpha = 0.5, size=1.1) +
    # geom_point(data=subset(don, shared==1), color="red",alpha = 1, size=1.2) +
    geom_hline(yintercept = -log10(thres_cmh_genome),linetype = "dashed",color = "red") +
    theme_minimal() 
  return(resplot)
}
p.cmh.cold <- readRDS("p.cmh.cold.RDS")
p.cmh.hot  <- readRDS("p.cmh.hot.RDS")


p.cmh.cold[p.cmh.hot[p.cmh.hot<0.05, .(CHROM,POS)], highlight :=1]
p.cmh.hot[p.cmh.cold[p.cmh.cold<0.05, .(CHROM,POS)], highlight :=1]
# 
# SNP <- intersect(p.cmh.cold[p.cmh.cold<0.01, .(CHROM, POS)], p.cmh.hot[p.cmh.hot<0.01, .(CHROM, POS)])
# p.cmh.cold[SNP, shared :=1]
# p.cmh.hot[SNP, shared :=1]


library(cowplot)
plot.cold <- plotMan(p.cmh.cold, highlightcolor = "maroon") +
  theme(axis.title.x = element_blank(),
        legend.position = "none") + 
  ylab("-log10(q) Cold ")+
  ylim(c(0, 15))

plot.hot <- plotMan(p.cmh.hot, top = T, highlightcolor = "steelblue") + 
  scale_y_reverse(limits = c(15, 0)) +
  ylab("-log10(q) Hot ")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
title <- ggdraw() + 
  draw_label("Adapted CMH on Cold & Hot evolution F0-F90",fontface = 'bold',x = 0,hjust = 0) +
  theme(    plot.margin = margin(0, 0, 0, 7))
png("../Plot/Suppl.figure.Manhattan.png", width = 8, height = 6, units = "in", res = 600)
cowplot::plot_grid(title, plot.cold,NULL,plot.hot, ncol = 1, rel_heights = c(0.1,1,-0.03,1))
dev.off()


## S2. Ne estimate plot ----

ne_estimates <- readRDS("ne_estimate.RDS")

ne_median <- ne_estimates[,lapply(.SD, FUN = median), by = .(chrom, evo, replicate, start, end)]
ne_median[,`:=`(trial = NULL, ne = round(ne+0.5))]
png("../Plot/Suppl.figure.Ne.png", width = 8, height = 6, units = "in", res = 600)
ggplot(data = ne_median, aes(x = evo, y = ne, color = evo)) +
  geom_hline(yintercept = 1250, linetype = "dashed", colour = "grey20", )+
  geom_boxplot(width = 0.8, size = 1) +
  facet_wrap(~chrom, nrow = 1) +
  scale_color_manual(values = c("steelblue","maroon")) +
  labs(x = "", y = "Ne estimate", color = "Regime", title = "Ne estimate based on temporal AFC F0-F90") 
dev.off()



## 8. Check particular genes' AF ----
### 8.1 Gr28b Gustatory receptor, heat preference ----
## Location 2L:7416508--7423849
af <- readRDS("af.RDS")

gr28b.snps <- af[CHROM=="2L" & POS >= 7411508 & POS <= 7428849, ]

p.cmh.cold <- readRDS("p.cmh.cold.RDS")
p.cmh.cold <- intersect(gr28b.snps[,.(CHROM,POS)], p.cmh.cold[p.cmh.cold<0.3, .(CHROM,POS)])

p.cmh.hot <- readRDS("p.cmh.hot.RDS")
p.cmh.hot <- intersect(gr28b.snps[,.(CHROM,POS)], p.cmh.hot[p.cmh.hot<0.3, .(CHROM,POS)])
targets <- unique(rbind(p.cmh.cold, p.cmh.hot))
gr28b.snps <- gr28b.snps[targets, ]

gr28b.snps[,`:=`(REF = NULL, ALT = NULL, N_ALT = NULL)]
f0_cols      <- grep("^F0_", names(gr28b.snps), value = TRUE)
cold_f0_cols <- paste0("Cold_", f0_cols)
hot_f0_cols  <- paste0("Hot_",  f0_cols)

gr28b.snps[, (cold_f0_cols) := .SD, .SDcols = f0_cols]
gr28b.snps[, (hot_f0_cols)  := .SD, .SDcols = f0_cols]
gr28b.snps[, (f0_cols) := NULL]   # drop original F0_* columns (optional)

gr28b.snps <- melt(gr28b.snps,id.vars = c("CHROM","POS"),variable.name = "Sample",value.name = "AF")

gr28b.snps[, c("Pop","Gen","Rep") := tstrsplit(Sample, "_", fixed = TRUE)]
gr28b.snps[, EvoRep := paste0(Pop, Rep)]
gr28b.snps[, SNP    := paste0(CHROM, POS)]

ggplot(gr28b.snps,aes(x = Gen, y = AF, color = Pop,group = interaction(SNP, EvoRep))) +
  facet_wrap(~ Pop * Rep) +
  geom_line() + 
  scale_color_manual(values = c("steelblue", "maroon"))


## 9. trying other plot than afc.dist----
library(ggpointdensity)
library(viridis)
ggplot(data = hot.snp, aes(y=cold.afc, x=hot.afc))+
  geom_pointdensity(alpha=0.3) +
  scale_color_viridis() +
  xlim(-0.5,0.5)+
  ylim(-0.5,0.5)+
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")+
  geom_abline(slope = -1, intercept = 0, color = "red", linetype = "dashed")+
  theme_minimal()

## ************************************ ----
## ************************************ ----

## 10. instead of CMH, use Chi-square test for all ----

af <- readRDS("af.covfilter.RDS")
xf <- readRDS("xf.covfilter.RDS")
dp <- readRDS("dp.covfilter.RDS")
ne_estimates <- readRDS("ne_estimate.RDS")

ne_median <- ne_estimates[,lapply(.SD, FUN = median), by = .(chrom, evo, replicate, start, end)]
ne_median[,`:=`(trial = NULL, ne = round(ne+0.5))]

library(ACER)

ne_median[chrom == "X", auto:= "X"]
ne_median[chrom != "X", auto:= "auto"]

Ne <- ne_median[  , .(med = median(ne)), by = .(evo,auto, replicate)]
Ne[,med := round(med + 0.5)]

res <- as.data.frame(xf[,.(CHROM,POS)])
setDT(res)
setkey(res,CHROM,POS)

for (evo in c("Hot","Cold")) {
  
  for (repp in 1:10) {
    focal_sample <- c(paste0("F0_r",repp), paste0(evo,"_F90_r",repp))
    af_auto <- as.matrix(xf[CHROM != "X",focal_sample, with = F])
    dp_auto <- as.matrix(dp[CHROM != "X",focal_sample, with = F])
    Ne_auto <- as.integer(Ne[evo == evo & auto == "auto" & replicate == repp, med])
    p.auto <- adapted.chisq.test(freq = af_auto,coverage = dp_auto,gen = c(0,90),mincov = 10,Ne = Ne_auto, poolSize = rep(625, 2),IntGen = F,RetVal = 0)
    res_auto <- as.data.frame(cbind(xf[CHROM != "X",.(CHROM,POS)], p.adjust(p.auto, method = "fdr")))
    colnames(res_auto) <- c("CHROM", "POS", paste0("P.",evo,".r",repp))
  
    af_x <- as.matrix(xf[CHROM == "X",focal_sample, with = F])
    dp_x <- as.matrix(dp[CHROM == "X",focal_sample, with = F])
    Ne_x <- as.integer(Ne[evo == evo & auto == "X" & replicate == repp, med])
    p.x <- adapted.chisq.test(freq = af_x,coverage = dp_x,gen = c(0,90),mincov = 10,Ne = Ne_x, poolSize = rep(625, 2),IntGen = F,RetVal = 0)
    res_x <- as.data.frame(cbind(xf[CHROM == "X",.(CHROM,POS)], p.adjust(p.x, method = "fdr")))
    colnames(res_x) <- c("CHROM", "POS", paste0("P.",evo,".r",repp))
  
    res.temp <- rbind(res_x, res_auto)
    setDT(res.temp)
    setkey(res.temp,CHROM,POS)
    res <- merge(res, res.temp, by = c("CHROM", "POS"),all = T)
    }
}

rm(res_auto, res_x, res.temp, af_auto, af_x, dp_auto, dp_x, Ne, ne_estimates, ne_median, evo, focal_sample, Ne_auto, Ne_x, p.auto, p.x, repp)


hot_cols  <- paste0("P.Hot.r", 1:10)
cold_cols <- paste0("P.Cold.r", 1:10)

res[, hot_sig_n  := rowSums(.SD < 0.001, na.rm = TRUE), .SDcols = hot_cols]
res[, cold_sig_n := rowSums(.SD < 0.001, na.rm = TRUE), .SDcols = cold_cols]

hot.090.target  <- res[hot_sig_n  > 0, .(CHROM, POS)]
cold.090.target <- res[cold_sig_n > 0, .(CHROM, POS)]

saveRDS(res,              "./chi-square.res.RDS")
saveRDS(hot.090.target,   "./chisq.hot.target.RDS")
saveRDS(cold.090.target,  "./chisq.cold.target.RDS")

## 11. Venn Diagram with chisq results ----

library(ggVennDiagram)

hot.090.target <- readRDS("./chisq.hot.target.RDS")
cold.090.target <- readRDS("./chisq.cold.target.RDS")

cold.snp <- cbind(
  cold.090.target,
  cold.afc = apply(xf[cold.090.target, 86:95,  with=F], 1, mean) - apply(xf[cold.090.target, 6:15, with=F], 1, mean),
  hot.afc  = apply(xf[cold.090.target, 166:175,with=F], 1, mean) - apply(xf[cold.090.target, 6:15, with=F], 1, mean)
) |> as.data.table()

hot.snp <- cbind(
  hot.090.target,
  cold.afc = apply(xf[hot.090.target, 86:95,  with=F], 1, mean) - apply(xf[hot.090.target, 6:15, with=F], 1, mean),
  hot.afc  = apply(xf[hot.090.target, 166:175,with=F], 1, mean) - apply(xf[hot.090.target, 6:15, with=F], 1, mean)
) |> as.data.table()

hot_up  <- unique(unlist(hot.snp[hot.afc>0,.(paste0(CHROM,"_",POS))]))
hot_dn  <- unique(unlist(hot.snp[hot.afc<0,.(paste0(CHROM,"_",POS))]))
cold_up  <- unique(unlist(cold.snp[cold.afc>0,.(paste0(CHROM,"_",POS))]))
cold_dn  <- unique(unlist(cold.snp[cold.afc<0,.(paste0(CHROM,"_",POS))]))

venn_obj <- list("Hot Rising"  = hot_up,
                 "Hot Dropping"  = hot_dn,
                 "Cold Rising"  = cold_up,
                 "Cold Dropping"  = cold_dn)

p_venn <- ggVennDiagram(x = venn_obj,label = "both",label_color = "grey5",
                        set_color = c("maroon","maroon","steelblue","steelblue"),
                        label_alpha = 0.15) +
  scale_fill_distiller(palette = "Reds", direction = 1)+
  theme_void(base_family = "sans") +
  coord_fixed(clip = "off") +
  ggtitle("Sharing of targets detected by chi-square tests")
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        legend.position = "none", plot.margin = margin(30, 30, 30, 30))

p_venn
png("../Plot/Chi-sq.Figure4.VennDiagram.png",width = 8, height = 8, units = "in", res = 600)
p_venn
dev.off()

## 12. Correlation matrix with chisq results ----
base_cols <- 6:15
cold_cols <- 86:95
hot_cols  <- 166:175

target_all <- unique(rbind(cold.090.target, hot.090.target))
setkey(target_all, CHROM, POS)

base_mat <- as.matrix(xf[target_all, ..base_cols])
cold_mat <- as.matrix(xf[target_all, ..cold_cols])
hot_mat  <- as.matrix(xf[target_all, ..hot_cols])

cold_afc_rep <- cold_mat - base_mat
hot_afc_rep  <- hot_mat  - base_mat

colnames(cold_afc_rep) <- paste0("afc.cold.r", 1:ncol(cold_afc_rep))
colnames(hot_afc_rep)  <- paste0("afc.hot.r",  1:ncol(hot_afc_rep))

afc_rep <- cbind(cold_afc_rep, hot_afc_rep)

cor_mat <- cor(afc_rep, use = "pairwise.complete.obs")
cor_dt <- as.data.table(as.table(cor_mat))
setnames(cor_dt, c("Var1","Var2","N"))
cor_dt[, Var1 := factor(Var1, levels = colnames(afc_rep))]
cor_dt[, Var2 := factor(Var2, levels = colnames(afc_rep))]

# triangle mask (lower triangle including diagonal)
cor_dt[, i := as.integer(Var1)]
cor_dt[, j := as.integer(Var2)]
cor_dt[i > j, N := NA_real_]
cor_dt[, lab := ifelse(is.na(N), NA, sprintf("%.2f", N))]
levels(cor_dt$Var1) <- sub("^afc\\.", "", levels(cor_dt$Var1))
levels(cor_dt$Var2) <- sub("^afc\\.", "", levels(cor_dt$Var2))
k <- ncol(cold_afc_rep)  # boundary between cold and hot (10 in your case)
n <- length(levels(cor_dt$Var2))
x_mid <- k + 0.5
y_mid <- n - k + 0.5

p_cor <- ggplot(cor_dt, aes(x = Var1, y = Var2, fill = N)) +
  geom_tile(color = "white", linewidth = 0.25, na.rm = FALSE) +
  geom_text(aes(label = lab), size = 2.6, na.rm = TRUE) +
  scale_fill_gradient2(low = "blue",  high = "red",na.value = "white") +
  scale_y_discrete(limits = rev(levels(cor_dt$Var2))) +
  geom_segment(x = x_mid, xend = x_mid, y = 0, yend = y_mid,
               linewidth = 0.6, color = "red", linetype = "dashed") +
  geom_segment(x = 0, xend = x_mid, y = y_mid, yend = y_mid,
               linewidth = 0.6, color = "red", linetype = "dashed") +
  labs(title = "Pairwise cor between replicates on AFC of Chi-sq test targets", x = NULL, y = NULL, fill = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 12))

png("../Plot/Chi-sq.Figure4.CorMatrix.png",width = 9, height = 8, units = "in", res = 600)
p_cor
dev.off()

## 13. AFC rank plot ----
### 13.1 cold ----

cold.snp <- cbind(
  cold.090.target,
  cold.afc = apply(xf[cold.090.target, 86:95,  with=F], 1, mean) - apply(xf[cold.090.target, 6:15, with=F], 1, mean),
  hot.afc  = apply(xf[cold.090.target, 166:175,with=F], 1, mean) - apply(xf[cold.090.target, 6:15, with=F], 1, mean)
) |> as.data.table()
cold.snp[cold.afc < 0, `:=`(cold.afc = abs(cold.afc), hot.afc = -1 * hot.afc)]

hot.snp <- cbind(
  hot.090.target,
  cold.afc = apply(xf[hot.090.target, 86:95,  with=F], 1, mean) - apply(xf[hot.090.target, 6:15, with=F], 1, mean),
  hot.afc  = apply(xf[hot.090.target, 166:175,with=F], 1, mean) - apply(xf[hot.090.target, 6:15, with=F], 1, mean)
) |> as.data.table()
hot.snp[hot.afc < 0, `:=`(hot.afc = abs(hot.afc), cold.afc = -1 * cold.afc)]

ranks_cold <- cold.snp[order(cold.afc, decreasing=T), .(CHROM,POS)][, rank := .I]
cold.snp.long <- melt(cold.snp, id.vars=c("CHROM","POS"), measure.vars=c("cold.afc","hot.afc"),
                      variable.name="Source", value.name="AFC")
cold.snp.long <- merge(cold.snp.long, ranks_cold, by=c("CHROM","POS"))
cold.snp.long$Source <- factor(cold.snp.long$Source, levels=c("hot.afc","cold.afc"))
setorder(cold.snp.long, Source)
 



