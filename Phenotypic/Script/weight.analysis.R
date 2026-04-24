setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Script/")
rm(list = ls())
gc()

library(data.table)
setDTthreads(percent = 80)
getDTthreads()

library(ggplot2)
library(lme4)
library(emmeans)
library(cowplot)

###############################################################################
## 1. Load dry weight data from respiratory dataset
###############################################################################

resp <- fread(
  "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Data/respiratory.csv",
  quote = FALSE,
  header = TRUE
)

setDT(resp)

resp[, weight.per := NetWeight / FlyNum]

## keep one row per measured vial
resp_w <- unique(
  resp[, .(cge, run, Evo, EvoRep, Sex, TechRep, FlyNum, NetWeight, weight.per)]
)

resp_w[, `:=`(
  cge     = factor(cge, levels = c("coldcge", "hotcge")),
  Evo     = factor(Evo, levels = c("Ancestral", "Cold", "Hot")),
  Sex     = factor(Sex),
  EvoRep  = factor(EvoRep),
  TechRep = factor(TechRep),
  run     = factor(run),
  run_id  = factor(paste(cge, run, sep = "_"))
)]

bg <- data.table(
  cge  = c("coldcge", "hotcge"),
  xmin = -Inf,
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf
)

###############################################################################
## 2. Raw inspection plot
###############################################################################

png("../Plot/resp.Weight.raw.inspection.png", width = 8, height = 6, units = "in", res = 600)
weight_raw_plot <- ggplot(resp_w, aes(Evo, weight.per, color = Evo, shape = Sex)) +
  geom_rect(
    aes(fill = cge),
    data = bg,
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    alpha = 0.25
  ) +
  geom_jitter(width = 0.18, size = 2) +
  facet_wrap(~ cge * Sex, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("coldcge" = "lightblue", "hotcge" = "lightpink")) +
  scale_color_manual(values = c(
    "Ancestral" = "forestgreen",
    "Cold"      = "steelblue",
    "Hot"       = "maroon"
  )) +
  ylab(expression(Weight ~ (mg))) +
  xlab("") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 0.6, vjust = 0.8)
  )
print(weight_raw_plot)
dev.off()

saveRDS(weight_raw_plot, "./weight.raw.plot.RDS")

###############################################################################
## 3. Candidate models for dry weight
##
## true biological replicate = EvoRep
## repeated measurement batch effect = run_id
## optional EvoRep-specific cge response
###############################################################################

## use ML for model comparison
m_w_1 <- lmer(
  (weight.per) ~ Evo * cge * Sex + (1 | EvoRep),
  data = resp_w,
  REML = FALSE
)

m_w_2 <- lmer(
  (weight.per) ~ Evo * cge * Sex + scale(FlyNum) + (1 | EvoRep),
  data = resp_w,
  REML = FALSE
)

m_w_3 <- lmer(
  (weight.per) ~ Evo * cge * Sex + (1 + cge | EvoRep) ,
  data = resp_w,
  REML = FALSE
)

m_w_4 <- lmer(
  (weight.per) ~ Evo * cge * Sex + scale(FlyNum) + (1 + cge | EvoRep) ,
  data = resp_w,
  REML = FALSE
)


m_w_5 <- lm((weight.per) ~ Evo * cge * Sex,data = resp_w)

anova(m_w_1, m_w_2, m_w_3, m_w_4, m_w_5)

###############################################################################
## 4. Final model
## log2(weight.per) ~ Evo * cge * Sex + (1 | EvoRep) 
###############################################################################

fit_weight <- lmer(
  (weight.per) ~ Evo * cge * Sex + (1 | EvoRep),
  data = resp_w,
  REML = TRUE
)

hist(resid(fit_weight), breaks = 20)

pdf("../Plot/resp.Weight.model.check.pdf", width = 10, height = 5)
plot(fit_weight)
qqnorm(resid(fit_weight))
qqline(resid(fit_weight))
dev.off()

###############################################################################
## 5. emmeans and pairwise contrasts
###############################################################################

emm_weight <- emmeans(fit_weight, ~ Evo | cge * Sex, type = "response")
pairs_weight <- contrast(emm_weight, method = "pairwise")

emm_weight_dt <- as.data.table(summary(emm_weight, type = "response"))
setnames(emm_weight_dt, "response", "emmean")
emm_weight_dt[Evo == "Cold", Evo := "Cold-evolved"]
emm_weight_dt[Evo == "Hot",  Evo := "Hot-evolved"]

pairs_weight_dt <- as.data.table(as.data.frame(summary(pairs_weight)))
pairs_weight_dt[, contrast := gsub("Cold", "Cold-evolved", contrast)]
pairs_weight_dt[, contrast := gsub("Hot",  "Hot-evolved",  contrast)]

pairs_weight_dt[, sig := cut(
  p.value,
  breaks = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  labels = c("****", "***", "**", "*", "n.s."),
  include.lowest = TRUE
)]
## keep only significant comparisons
pairs_weight_dt <- pairs_weight_dt[sig != "n.s."]


lookup_x <- c(
  "Ancestral"    = 1,
  "Cold-evolved" = 2,
  "Hot-evolved"  = 3
)

pairs_weight_dt[, group1 := sub(" - .*", "", contrast)]
pairs_weight_dt[, group2 := sub(".* - ", "", contrast)]
pairs_weight_dt[, xstart := lookup_x[group1]]
pairs_weight_dt[, xend   := lookup_x[group2]]

panel_span_weight <- emm_weight_dt[,.(
  ymax = max(upper.CL, na.rm = TRUE), 
  span = diff(range(c(lower.CL, upper.CL), na.rm = TRUE))
  ),by = .(cge, Sex)]

pairs_weight_dt[, rank := match(
  contrast,
  c("Ancestral - Cold-evolved",
    "Ancestral - Hot-evolved",
    "Cold-evolved - Hot-evolved")
)]

ann_weight <- merge(pairs_weight_dt, panel_span_weight, by = c("cge", "Sex"), all.x = TRUE)
ann_weight[, y := ymax + rank * 0.12 * span]
ann_weight[, label_y := ifelse(sig == "n.s.", y + 0.05 * span, y + 0.015 * span)]

###############################################################################
## 6. Final dry weight plot
###############################################################################

png("../Plot/Respiratory.Weight.png", width = 8, height = 6, units = "in", res = 600)
weight_tested_plot <- ggplot(emm_weight_dt, aes(Evo, emmean, color = Evo, shape = Sex)) +
  geom_rect(
    aes(fill = cge),
    data = bg,
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    alpha = 0.25
  ) +
  facet_wrap(~ Sex*cge, scales = "free_y", nrow = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_segment(
    data = ann_weight,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = ann_weight,
    aes(x = (xstart + xend) / 2, y = label_y, label = sig),
    inherit.aes = FALSE,
    size = 4.5
  ) +
  scale_fill_manual(values = c("coldcge" = "lightblue", "hotcge" = "lightpink")) +
  scale_color_manual(values = c(
    "Ancestral"    = "forestgreen",
    "Cold-evolved" = "steelblue",
    "Hot-evolved"  = "maroon"
  )) +
  ylab(expression(Weight ~ (mg))) +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 0.6, vjust = 0.8)
  )
print(weight_tested_plot)
dev.off()

fwrite(emm_weight_dt, "./resp.Weight.emm.csv", quote = FALSE)
fwrite(pairs_weight_dt, "./resp.Weight.contrasts.csv", quote = FALSE)
saveRDS(weight_tested_plot, "./weight.tested.plot.RDS")

###############################################################################
## 7. Combine development and weight plots
###############################################################################
library(grid)

theme_bigfig <- theme(
  plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
  axis.title   = element_text(size = 20),
  axis.text    = element_text(size = 18, colour = "black"),
  strip.text   = element_text(size = 18, face = "plain"),
  legend.title = element_text(size = 18, face = "bold"),
  legend.text  = element_text(size = 16),
  plot.caption = element_text(size = 16),
  legend.key.size = unit(0.7, "cm")
)

dev_raw_plot <- readRDS("./deve.raw.plot.RDS")
dev_raw_plot <- dev_raw_plot +
  ggtitle("Eclosion curve") +
  ylab("Cumulative number of eclosed individual") +
  theme_bigfig +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  ) +
  guides(
    color = guide_legend(nrow = 1),
    fill  = guide_legend(nrow = 1)
  )

dev_t10_plot <- readRDS("./deve.plot.RDS")
dev_t10_plot <- dev_t10_plot +
  ggtitle("Development time") +
  labs(x = NULL, caption = "") +
  theme_bigfig +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )


weight_tested_plot <- readRDS("./weight.tested.plot.RDS")
weight_tested_plot <- weight_tested_plot +
  ggtitle("Dry weight") +
  labs(caption = "**** p < 0.0001    *** p < 0.001    ** p < 0.01    * p < 0.05") +
  theme_bigfig +
  theme(
    plot.caption = element_text(size = 16, hjust = 0.5)
  )

left <- plot_grid(
  dev_t10_plot,
  weight_tested_plot,
  labels = c("A", "C"),
  label_size = 24,
  label_fontface = "bold",
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1.6)
)

right <- plot_grid(
  dev_raw_plot,
  NULL,
  labels = c("B", ""),
  label_size = 24,
  label_fontface = "bold",
  ncol = 1,
  align = "h",
  rel_heights = c(1, 0)
)

final_dev_weight_figure <- plot_grid(
  left,
  right,
  ncol = 2,
  rel_widths = c(1.3, 1)
)

ggsave(
  "../Plot/Figure8.png",
  final_dev_weight_figure,
  width = 15,
  height = 10,
  units = "in",
  dpi = 600
)
