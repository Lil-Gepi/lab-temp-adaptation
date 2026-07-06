setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Script/")
library(data.table)
setDTthreads(percent = 80)
getDTthreads()

library(ggplot2)
library(lme4)
library(emmeans)

rf <- fread(
  "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Data/relativefitness.csv",
  quote = FALSE,
  header = TRUE
)

setDT(rf)

rf[, vial := paste(cge, Sample, sep = ".")]

bg <- data.table(
  cge  = c("coldcge", "hotcge"),
  xmin = -Inf, xmax = Inf,
  ymin = -Inf, ymax = Inf
)

###############################################################################
## 1. Raw inspection plots
###############################################################################

## raw cumulative count summary
development_count <- rf[,.(
  mean_Red = mean(Cumulative_Red_Nr, na.rm = TRUE),
  sd_Red   = sd(Cumulative_Red_Nr, na.rm = TRUE)
),by = .(cge, Evo, ColDay)
]

## get final observed count per vial
final_counts <- rf[,
                   .(final_n = max(Cumulative_Red_Nr, na.rm = TRUE)),
                   by = .(cge, Evo, EvoRep, TechRep, Sample, vial)
]

## merge and calculate cumulative proportion per vial
rf_prop <- merge(
  rf,
  final_counts,
  by = c("cge", "Evo", "EvoRep", "TechRep", "Sample", "vial"),
  all.x = TRUE
)

rf_prop[, prop_red := fifelse(final_n > 0, Cumulative_Red_Nr / final_n, NA_real_)]

## summary for proportion plot
development_prop <- rf_prop[,.(
  mean_prop = mean(prop_red, na.rm = TRUE),
  sd_prop   = sd(prop_red, na.rm = TRUE)
),by = .(cge, Evo, ColDay)
]

## optional trimming for display only
development_count_plot <- development_count[
  (cge == "coldcge" & ColDay > 15) |
    (cge == "hotcge"  & ColDay > 7)
]

development_prop_plot <- development_prop[
  (cge == "coldcge" & ColDay > 15) |
    (cge == "hotcge"  & ColDay > 7)
]

png("../Plot/Development_raw_count.png", width = 10, height = 7, units = "in", res = 600)
deve.raw.plot <- ggplot(development_count_plot) +
  facet_wrap(~ cge, nrow = 1, scales = "free_x") +
  geom_ribbon(
    aes(x = ColDay, ymin = mean_Red - sd_Red, ymax = mean_Red + sd_Red, fill = Evo),
    alpha = 0.2
  ) +
  geom_line(
    aes(x = ColDay, y = mean_Red, color = Evo),
    linewidth = 1.5
  ) +
  scale_color_manual(values = c("forestgreen", "steelblue", "maroon")) +
  scale_fill_manual(values = c("forestgreen", "steelblue", "maroon")) +
  labs(
    y = "Cumulative red-eye eclosion",
    x = "Day post egg-laying",
    title = "Raw cumulative counts"
  ) +
  theme_minimal(base_size = 15)
print(deve.raw.plot)
dev.off()
saveRDS(deve.raw.plot, "./deve.raw.plot.RDS")


png("../Plot/Development_raw_proportion.png", width = 10, height = 7, units = "in", res = 600)
ggplot(development_prop_plot) +
  facet_wrap(~ cge, nrow = 1, scales = "free_x") +
  geom_ribbon(
    aes(x = log(ColDay), ymin = pmax(mean_prop - sd_prop, 0), ymax = pmin(mean_prop + sd_prop, 1), fill = Evo),
    alpha = 0.2
  ) +
  geom_line(
    aes(x = log(ColDay), y = mean_prop, color = Evo),
    linewidth = 1.5
  ) +
  scale_color_manual(values = c("forestgreen", "steelblue", "maroon")) +
  scale_fill_manual(values = c("forestgreen", "steelblue", "maroon")) +
  labs(
    y = "Cumulative proportion of final red-eye eclosion",
    x = "Day post egg-laying (log scaled)",
    title = "Raw cumulative proportion"
  ) +
  theme_minimal(base_size = 15)
dev.off()

###############################################################################
## 2. Empirical T10 function
###############################################################################

get_T10 <- function(dt) {
  dt <- copy(dt)
  setorder(dt, ColDay)
  dt <- dt[!is.na(ColDay) & !is.na(Cumulative_Red_Nr)]
  
  final_n <- max(dt$Cumulative_Red_Nr, na.rm = TRUE)
  
  threshold_10 <- 0.10 * final_n
  
  idx_hi <- which(dt$Cumulative_Red_Nr >= threshold_10)[1]
  
  if (idx_hi == 1L) {
    return(data.table(
      final_n = final_n,
      threshold_10 = threshold_10,
      T10 = dt$ColDay[1],
      status = "first_day"
    ))
  }
  
  idx_lo <- idx_hi - 1L
  
  x0 <- dt$ColDay[idx_lo]
  x1 <- dt$ColDay[idx_hi]
  y0 <- dt$Cumulative_Red_Nr[idx_lo]
  y1 <- dt$Cumulative_Red_Nr[idx_hi]
  
  if (y1 == y0) {
    T10 <- x1
    status <- "flat_segment"
  } else {
    T10 <- x0 + (threshold_10 - y0) / (y1 - y0) * (x1 - x0)
    status <- "interpolated"
  }
  
  data.table(
    final_n = final_n,
    threshold_10 = threshold_10,
    T10 = T10,
    status = status
  )
}

###############################################################################
## 3. Compute T10 for each vial
###############################################################################

id_cols <- c("cge", "Evo", "EvoRep", "TechRep", "Sample", "vial")

t10_dt <- rf[ , get_T10(.SD),  by = id_cols]

fwrite(t10_dt, "./development_T10_all.csv", quote = FALSE)

## keep valid values for inference
t10_fit <- t10_dt[is.finite(T10) & T10 > 0]

###############################################################################
## 4. Quick T10 plot
###############################################################################

png("../Plot/Development_T10_boxplot.png", width = 8, height = 6, units = "in", res = 600)
ggplot(t10_fit, aes(x = Evo, y = T10, fill = Evo)) +
  facet_wrap(~ cge, scale = "free_y") +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  ggtitle("Empirical T10 by population") +
  theme_minimal(base_size = 15)
dev.off()

###############################################################################
## 5. Mixed model on T10
###############################################################################

fit_T10 <- lmer(log(T10) ~ Evo * cge + (1 | EvoRep), data = t10_fit)

summary(fit_T10)
plot(fit_T10)
qqnorm(resid(fit_T10))
qqline(resid(fit_T10))

###############################################################################
## 6. Estimated marginal means and pairwise contrasts
###############################################################################

emm.res <- emmeans(
  fit_T10,
  pairwise ~ Evo * cge,
  by = "cge",
  adjust = "BH",
  type = "response"
)

emm.plot <- as.data.table(summary(emm.res, type = "response")$emmeans)
setnames(emm.plot, "response", "emmean")

emm.plot[Evo == "Cold", Evo := "Cold-evolved"]
emm.plot[Evo == "Hot",  Evo := "Hot-evolved"]

pairs <- as.data.table(as.data.frame(emm.res$contrasts))
pairs[, contrast := c(
  "Ancestral - Cold-evolved",
  "Ancestral - Hot-evolved",
  "Cold-evolved - Hot-evolved",
  "Ancestral - Cold-evolved",
  "Ancestral - Hot-evolved",
  "Cold-evolved - Hot-evolved"
)]

pairs[, sig := cut(
  p.value,
  breaks = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  labels = c("****", "***", "**", "*", "n.s."),
  include.lowest = TRUE
)]

lookup_x <- c(
  "Ancestral" = 1,
  "Cold-evolved" = 2,
  "Hot-evolved" = 3
)

panel_span <- emm.plot[
  ,
  .(
    ymax = max(upper.CL, na.rm = TRUE),
    span = diff(range(c(lower.CL, upper.CL), na.rm = TRUE))
  ),
  by = cge
]

pairs[, group1 := sub(" - .*", "", contrast)]
pairs[, group2 := sub(".* - ", "", contrast)]
pairs[, xstart := lookup_x[group1]]
pairs[, xend   := lookup_x[group2]]

pairs[, rank := match(
  contrast,
  rev(c(
    "Ancestral - Cold-evolved",
    "Ancestral - Hot-evolved",
    "Cold-evolved - Hot-evolved"
  ))
)]

ann <- merge(pairs, panel_span, by = "cge", all.x = TRUE)
ann[, y := ymax + (rank + 0.5) * 0.08 * span]
ann[, label_y := fifelse(sig == "n.s.", y + 0.05 * span, y + 0.01 * span)]

###############################################################################
## 7. Final manuscript plot
###############################################################################

png("../Plot/Development_T10.png", width = 8, height = 6, units = "in", res = 600)
deve.plot <- ggplot(emm.plot, aes(x = Evo, y = emmean, colour = Evo)) +
  facet_wrap(~ cge, scales = "free_y", nrow = 1) +
  geom_rect(
    aes(fill = cge),
    data = bg,
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    alpha = 0.25
  ) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.2
  ) +
  geom_point(size = 3) +
  scale_fill_manual(values = c("coldcge" = "lightblue", "hotcge" = "lightpink")) +
  scale_color_manual(values = c(
    "Ancestral"    = "forestgreen",
    "Cold-evolved" = "steelblue",
    "Hot-evolved"  = "maroon"
  )) +
  xlab("") +
  ylab(expression(Estimated ~ T[10])) +
  ggtitle("Egg-to-adult development time (T10)") +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05))) +
  geom_segment(
    data = ann,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = ann,
    aes(x = (xstart + xend) / 2, y = label_y, label = sig),
    inherit.aes = FALSE,
    size = 5
  ) +
  labs(caption = "**** p < 0.0001    *** p < 0.001    n.s. non-significant") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 0.6, vjust = 0.8)
  )
print(deve.plot)
dev.off()

###############################################################################
## 8. Export results
###############################################################################

fwrite(t10_fit, "./development_T10_vial_values.csv", quote = FALSE)
fwrite(emm.plot, "./development_T10_emm.csv", quote = FALSE)
fwrite(as.data.table(as.data.frame(emm.res$contrasts)),
       "./development_T10_contrasts.csv",quote = FALSE)
saveRDS(deve.plot, "./deve.plot.RDS")
