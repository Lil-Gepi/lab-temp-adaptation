setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Script/")
rm(list = ls())
gc()

library(data.table)
setDTthreads(percent = 80)
getDTthreads()

library(ggplot2)
library(lme4)
library(emmeans)

###############################################################################
## 1. Load data and basic setup
###############################################################################

resp <- fread(
  "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Data/respiratory.csv",
  quote = FALSE,  header = TRUE)

setDT(resp)

## derived variables
resp[, `:=`(
  co2.mg     = (co2 * FlyNum) / (NetWeight * (FlyNum - FlyDead)),
  o2.mg      = (o2  * FlyNum) / (NetWeight * (FlyNum - FlyDead)),
  weight.per = NetWeight / FlyNum
)]

## factors
resp[, `:=`(
  cge     = factor(cge, levels = c("coldcge", "hotcge")),
  Evo     = factor(Evo, levels = c("Ancestral", "Cold", "Hot")),
  Sex     = factor(Sex),
  EvoRep  = factor(EvoRep),
  TechRep = factor(TechRep),
  run     = factor(run),
  chamber = factor(chamber)
)]

## IDs
## one biological / technical sample = one vial = one chamber position in one run
resp[, `:=`(
  run_id    = factor(paste(cge, run, sep = "_")),
  sample_id = factor(paste(cge, EvoRep, Sex, TechRep, sep = "_")),
  tech_id   = factor(paste(EvoRep, TechRep, sep = "_"))
)]

bg <- data.table(
  cge  = c("coldcge", "hotcge"),
  xmin = -Inf,
  xmax = Inf,
  ymin = -Inf,
  ymax = Inf
)

###############################################################################
## 2. Raw cycle plots to justify the stability window
###############################################################################

## CO2 across cycles
co2_cycle_sum <- resp[,.(
    mean_co2 = mean(co2.mg, na.rm = TRUE),
    se_co2   = sd(co2.mg, na.rm = TRUE) / sqrt(.N)
  ),by = .(cge, Sex, cycle)
]

png("../Plot/resp.CO2.cycle_summary.png", width = 9, height = 6, units = "in", res = 600)
ggplot(co2_cycle_sum, aes(x = cycle, y = mean_co2, color = Sex, fill = Sex)) +
  facet_wrap(~ cge, nrow = 1, scales = "free_y") +
  geom_ribbon(aes(ymin = mean_co2 - se_co2, ymax = mean_co2 + se_co2), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  xlab("Cycle") +
  ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1}))) +
  ggtitle("CO2 across cycles") +
  theme_minimal(base_size = 15)
dev.off()

ggplot(resp[cge== "coldcge",], aes(x=cycle, y=co2.mg,group = chamber, colour = Sex))+
  facet_wrap(~run)+geom_smooth(method = "loess",linewidth = 1.5)+
  ggtitle("CO2 emission in coldCGE")+
  xlab("Cycles")+ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1})))

## O2 across cycles
o2_cycle_sum <- resp[, .(mean_o2 = mean(o2.mg, na.rm = TRUE),
    se_o2   = sd(o2.mg, na.rm = TRUE) / sqrt(.N)),by = .(cge, Sex, cycle)]

png("../Plot/resp.O2.cycle_summary.png", width = 9, height = 6, units = "in", res = 600)
ggplot(o2_cycle_sum, aes(x = cycle, y = mean_o2, color = Sex, fill = Sex)) +
  facet_wrap(~ cge, nrow = 1, scales = "free_y") +
  geom_ribbon(aes(ymin = mean_o2 - se_o2, ymax = mean_o2 + se_o2), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  xlab("Cycle") +
  ylab(expression(V[O2]~~(mu*L*h^{-1}*mg^{-1}))) +
  ggtitle("O2 across cycles") +
  theme_minimal(base_size = 15)
dev.off()

###############################################################################
## 3. Define resting metabolic rate
##
## Based on the cycle plots:
## - early cycles show settling, especially in males
## - middle cycles look most stable
## - later cycles may drift up again
##
## So use cycles 5:10 and summarize each sample by the median.
###############################################################################

resp_rest <- rbind(resp[run_id == "hotcge_12",],
                   resp[run_id != "hotcge_12" & cycle >= 5 & cycle <= 10,])

resp_rest <- resp_rest[,.(
    co2.rest    = median(co2.mg, na.rm = TRUE),
    o2.rest     = median(o2.mg,  na.rm = TRUE),
    weight.per  = mean(weight.per, na.rm = TRUE),
    FlyNum      = mean(FlyNum, na.rm = TRUE),
    NetWeight   = mean(NetWeight, na.rm = TRUE),
    n_cycles    = .N
  ), by = .(cge, Evo, EvoRep, Sex, TechRep, run_id, sample_id, tech_id)
]

###############################################################################
## 4. Raw inspection plots of resting values
###############################################################################

png("../Plot/resp.CO2.rest.raw.png", width = 8, height = 6, units = "in", res = 600)
ggplot(resp_rest, aes(Evo, co2.rest, color = Evo, shape = Sex)) +
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
  ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1}))) +
  xlab("") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 0.6, vjust = 0.8)
  )
dev.off()

png("../Plot/resp.O2.rest.raw.png", width = 8, height = 6, units = "in", res = 600)
ggplot(resp_rest, aes(Evo, o2.rest, color = Evo, shape = Sex)) +
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
  ylab(expression(V[O2]~~(mu*L*h^{-1}*mg^{-1}))) +
  xlab("") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 0.6, vjust = 0.8)
  )
dev.off()

###############################################################################
## 5. CO2 models
###############################################################################

## Use ML when comparing models
m_co2_1 <- lmer(
  log2(co2.rest) ~ Evo * cge * Sex + (1 | EvoRep),
  data = resp_rest,
  REML = FALSE
)

m_co2_2 <- lmer(
  log2(co2.rest) ~ Evo * cge * Sex + (1 | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

m_co2_3 <- lmer(
  log2(co2.rest) ~ Evo * cge * Sex + (1 + cge | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

m_co2_4 <- lmer(
  log2(co2.rest) ~ Evo * cge * Sex + (1 + Sex | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

m_co2_5 <- lmer(
  log2(co2.rest) ~ Evo * cge * Sex + (1 + cge + Sex | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE)

anova(m_co2_1, m_co2_2, m_co2_3, m_co2_4, m_co2_5)

hist(resid(m_co2_1), breaks=15)

## choose final model
fit_co2 <- lmer(co2.rest ~ Evo * cge * Sex + (1 | EvoRep), data = resp_rest,REML = T)
hist(resid(fit_co2), breaks=20)

pdf("../Plot/resp.CO2.model.check.pdf", width = 10, height = 5)
plot(fit_co2)
qqnorm(resid(fit_co2))
qqline(resid(fit_co2))
dev.off()

emm_co2 <- emmeans(fit_co2, ~ Evo | cge * Sex, type = "response")
pairs_co2 <- contrast(emm_co2, method = "pairwise", adjust = "holm")

emm_co2_dt <- as.data.table(summary(emm_co2, type = "response"))

emm_co2_dt[Evo == "Cold", Evo := "Cold-evolved"]
emm_co2_dt[Evo == "Hot",  Evo := "Hot-evolved"]

pairs_co2_dt <- as.data.table(as.data.frame(summary(pairs_co2)))
pairs_co2_dt[, contrast := gsub("Cold", "Cold-evolved", contrast)]
pairs_co2_dt[, contrast := gsub("Hot",  "Hot-evolved",  contrast)]

pairs_co2_dt[, sig := cut(
  p.value,
  breaks = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  labels = c("****", "***", "**", "*", "n.s."),
  include.lowest = TRUE
)]
pairs_co2_dt <- pairs_co2_dt[sig != "n.s."]

lookup_x <- c(
  "Ancestral"    = 1,
  "Cold-evolved" = 2,
  "Hot-evolved"  = 3
)

pairs_co2_dt[, group1 := sub(" - .*", "", contrast)]
pairs_co2_dt[, group2 := sub(".* - ", "", contrast)]
pairs_co2_dt[, xstart := lookup_x[group1]]
pairs_co2_dt[, xend   := lookup_x[group2]]

panel_span_co2 <- emm_co2_dt[,.(
  ymax = max(upper.CL, na.rm = TRUE),
  span = diff(range(c(lower.CL, upper.CL), na.rm = TRUE))),
  by = .(cge, Sex)]

pairs_co2_dt[, rank := match(
  contrast,
  c("Ancestral - Cold-evolved",
    "Ancestral - Hot-evolved",
    "Cold-evolved - Hot-evolved")
)]

ann_co2 <- merge(pairs_co2_dt, panel_span_co2, by = c("cge", "Sex"), all.x = TRUE)
ann_co2[, y := ymax + rank * 0.12 * span]
ann_co2[, label_y := ifelse(sig == "n.s.", y + 0.05 * span, y + 0.015 * span)]

png("../Plot/Respiratory.CO2.png", width = 8, height = 6, units = "in", res = 600)
co2.plot <- ggplot(emm_co2_dt, aes(Evo, emmean, color = Evo, shape = Sex)) +
  geom_rect(
    aes(fill = cge),
    data = bg,
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    alpha = 0.25
  ) +
  facet_wrap(~ cge * Sex, scales = "free_y", nrow = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_segment(
    data = ann_co2,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = ann_co2,
    aes(x = (xstart + xend) / 2, y = label_y, label = sig),
    inherit.aes = FALSE,
    size = 6
  ) +
  scale_fill_manual(values = c("coldcge" = "lightblue", "hotcge" = "lightpink")) +
  scale_color_manual(values = c(
    "Ancestral"    = "forestgreen",
    "Cold-evolved" = "steelblue",
    "Hot-evolved"  = "maroon"
  )) +
  ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1}))) +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 0.6, vjust = 0.8)
  )
print(co2.plot)
dev.off()

fwrite(emm_co2_dt, "./resp.CO2.emm.csv", quote = FALSE)
fwrite(pairs_co2_dt, "./resp.CO2.contrasts.csv", quote = FALSE)
saveRDS(co2.plot, "./co2.plot.RDS")

###############################################################################
## 6. O2 models
###############################################################################

## Use ML when comparing models
m_o2_1 <- lmer(
  o2.rest ~ Evo * cge * Sex + (1 | EvoRep),
  data = resp_rest,
  REML = FALSE
)

m_o2_2 <- lmer(
  o2.rest ~ Evo * cge * Sex + (1 | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

m_o2_3 <- lmer(
  o2.rest ~ Evo * cge * Sex + (1 + cge | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

m_o2_4 <- lmer(
  o2.rest ~ Evo * cge * Sex + (1 + Sex | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

m_o2_5 <- lmer(
  o2.rest ~ Evo * cge * Sex + (1 + cge + Sex | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

anova(m_o2_1, m_o2_2, m_o2_3, m_o2_4, m_o2_5)

hist(resid(m_o2_2), breaks = 15)

## choose final model
fit_o2 <- lmer(o2.rest ~ Evo * cge * Sex + (1 | EvoRep)+ (1 | run_id), data = resp_rest, REML = TRUE)
hist(resid(fit_o2), breaks = 20)

pdf("../Plot/resp.O2.model.check.pdf", width = 10, height = 5)
plot(fit_o2)
qqnorm(resid(fit_o2))
qqline(resid(fit_o2))
dev.off()

emm_o2 <- emmeans(fit_o2, ~ Evo | cge * Sex, type = "response")
pairs_o2 <- contrast(emm_o2, method = "pairwise", adjust = "holm")

emm_o2_dt <- as.data.table(summary(emm_o2, type = "response"))

emm_o2_dt[Evo == "Cold", Evo := "Cold-evolved"]
emm_o2_dt[Evo == "Hot",  Evo := "Hot-evolved"]

pairs_o2_dt <- as.data.table(as.data.frame(summary(pairs_o2)))
pairs_o2_dt[, contrast := gsub("Cold", "Cold-evolved", contrast)]
pairs_o2_dt[, contrast := gsub("Hot",  "Hot-evolved",  contrast)]

pairs_o2_dt[, sig := cut(
  p.value,
  breaks = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  labels = c("****", "***", "**", "*", "n.s."),
  include.lowest = TRUE
)]
pairs_o2_dt <- pairs_o2_dt[sig != "n.s."]

pairs_o2_dt[, group1 := sub(" - .*", "", contrast)]
pairs_o2_dt[, group2 := sub(".* - ", "", contrast)]
pairs_o2_dt[, xstart := lookup_x[group1]]
pairs_o2_dt[, xend   := lookup_x[group2]]

panel_span_o2 <- emm_o2_dt[,.(
    ymax = max(upper.CL, na.rm = TRUE),
    span = diff(range(c(lower.CL, upper.CL), na.rm = TRUE))
  ),by = .(cge, Sex)]

pairs_o2_dt[, rank := match(
  contrast,
  c("Ancestral - Cold-evolved",
    "Ancestral - Hot-evolved",
    "Cold-evolved - Hot-evolved")
)]

ann_o2 <- merge(pairs_o2_dt, panel_span_o2, by = c("cge", "Sex"), all.x = TRUE)
ann_o2[, y := ymax + rank * 0.12 * span]
ann_o2[, label_y := ifelse(sig == "n.s.", y + 0.05 * span, y + 0.015 * span)]

png("../Plot/Respiratory.O2.png", width = 8, height = 6, units = "in", res = 600)

o2.plot <- ggplot(emm_o2_dt, aes(Evo, emmean, color = Evo, shape = Sex)) +
  geom_rect(
    aes(fill = cge),
    data = bg,
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    alpha = 0.25
  ) +
  facet_wrap(~ cge * Sex, scales = "free_y", nrow = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_segment(
    data = ann_o2,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = ann_o2,
    aes(x = (xstart + xend) / 2, y = label_y, label = sig),
    inherit.aes = FALSE,
    size = 6
  ) +
  scale_fill_manual(values = c("coldcge" = "lightblue", "hotcge" = "lightpink")) +
  scale_color_manual(values = c(
    "Ancestral"    = "forestgreen",
    "Cold-evolved" = "steelblue",
    "Hot-evolved"  = "maroon"
  )) +
  ylab(expression(V[O2]~~(mu*L*h^{-1}*mg^{-1}))) +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 0.6, vjust = 0.8)
  )
print(o2.plot)
dev.off()

fwrite(emm_o2_dt, "./resp.O2.emm.csv", quote = FALSE)
fwrite(pairs_o2_dt, "./resp.O2.contrasts.csv", quote = FALSE)
saveRDS(o2.plot, "./o2.plot.RDS")

###############################################################################
## 7. Respiratory quotient (RQ)
###############################################################################

## resting RQ from resting CO2 and O2
resp_rest[, rq.rest := co2.rest / o2.rest]

png("../Plot/resp.RQ.rest.raw.png", width = 8, height = 6, units = "in", res = 600)
ggplot(resp_rest, aes(Evo, rq.rest, color = Evo, shape = Sex)) +
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
  ylab(expression(RQ == V[CO2] / V[O2])) +
  xlab("") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 0.6, vjust = 0.8)
  )
dev.off()

## candidate models
m_rq_1 <- lmer(
  log2(rq.rest) ~ Evo * cge * Sex + (1 | EvoRep),
  data = resp_rest,
  REML = FALSE
)

m_rq_2 <- lmer(
  log2(rq.rest) ~ Evo * cge * Sex + (1 | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

m_rq_3 <- lmer(
  log2(rq.rest) ~ Evo * cge * Sex + (1 + cge | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

m_rq_4 <- lmer(
  log2(rq.rest) ~ Evo * cge * Sex + (1 + Sex | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

m_rq_5 <- lmer(
  log2(rq.rest) ~ Evo * cge * Sex + (1 + cge + Sex | EvoRep) + (1 | run_id),
  data = resp_rest,
  REML = FALSE
)

anova(m_rq_1, m_rq_2, m_rq_3, m_rq_4, m_rq_5)

hist(resid(m_rq_2), breaks = 20)

## choose final model
fit_rq <- lmer(log2(rq.rest) ~ Evo * cge * Sex + (1 | EvoRep)+ (1 | run_id), data = resp_rest, REML = TRUE)
hist(resid(fit_rq), breaks = 20)

pdf("../Plot/resp.RQ.model.check.pdf", width = 10, height = 5)
plot(fit_rq)
qqnorm(resid(fit_rq))
qqline(resid(fit_rq))
dev.off()

emm_rq <- emmeans(fit_rq, ~ Evo | cge * Sex, type = "response")
pairs_rq <- contrast(emm_rq, method = "pairwise", adjust = "holm")

emm_rq_dt <- as.data.table(summary(emm_rq, type = "response"))
setnames(emm_rq_dt, "response", "emmean")
emm_rq_dt[Evo == "Cold", Evo := "Cold-evolved"]
emm_rq_dt[Evo == "Hot",  Evo := "Hot-evolved"]

pairs_rq_dt <- as.data.table(as.data.frame(summary(pairs_rq)))
pairs_rq_dt[, contrast := gsub("Cold", "Cold-evolved", contrast)]
pairs_rq_dt[, contrast := gsub("Hot",  "Hot-evolved",  contrast)]

pairs_rq_dt[, sig := cut(
  p.value,
  breaks = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  labels = c("****", "***", "**", "*", "n.s."),
  include.lowest = TRUE
)]
pairs_rq_dt <- pairs_rq_dt[sig != "n.s."]

pairs_rq_dt[, group1 := sub(" - .*", "", contrast)]
pairs_rq_dt[, group2 := sub(".* - ", "", contrast)]
pairs_rq_dt[, xstart := lookup_x[group1]]
pairs_rq_dt[, xend   := lookup_x[group2]]

panel_span_rq <- emm_rq_dt[,.(
    ymax = max(upper.CL, na.rm = TRUE),
    span = diff(range(c(lower.CL, upper.CL), na.rm = TRUE))
  ),by = .(cge, Sex)]

pairs_rq_dt[, rank := match(
  contrast,
  c("Ancestral - Cold-evolved",
    "Ancestral - Hot-evolved",
    "Cold-evolved - Hot-evolved")
)]

ann_rq <- merge(pairs_rq_dt, panel_span_rq, by = c("cge", "Sex"), all.x = TRUE)
ann_rq[, y := ymax + rank * 0.12 * span]
ann_rq[, label_y := ifelse(sig == "n.s.", y + 0.05 * span, y + 0.015 * span)]

png("../Plot/Respiratory.RQ.png", width = 8, height = 6, units = "in", res = 600)
rq.plot <- ggplot(emm_rq_dt, aes(Evo, emmean, color = Evo, shape = Sex)) +
  geom_rect(
    aes(fill = cge),
    data = bg,
    inherit.aes = FALSE,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    alpha = 0.25
  ) +
  facet_wrap(~ cge * Sex, nrow = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_segment(
    data = ann_rq,
    aes(x = xstart, xend = xend, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_text(
    data = ann_rq,
    aes(x = (xstart + xend) / 2, y = label_y, label = sig),
    inherit.aes = FALSE,
    size = 6
  ) +
  scale_fill_manual(values = c("coldcge" = "lightblue", "hotcge" = "lightpink")) +
  scale_color_manual(values = c(
    "Ancestral"    = "forestgreen",
    "Cold-evolved" = "steelblue",
    "Hot-evolved"  = "maroon"
  )) +
  ylab(expression(RQ == V[CO2] / V[O2])) +
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 15, hjust = 0.6, vjust = 0.8)
  )+
  labs(caption = "**** p < 0.0001    *** p < 0.001    ** p < 0.01    * p < 0.05    n.s. non-significant")
print(rq.plot)
dev.off()

fwrite(emm_rq_dt, "./resp.RQ.emm.csv", quote = FALSE)
fwrite(pairs_rq_dt, "./resp.RQ.contrasts.csv", quote = FALSE)
saveRDS(rq.plot, "./rq.plot.RDS")

###############################################################################
## 8. Combine Figure 7 
###############################################################################

library(cowplot)
library(grid)

theme_bigfig <- theme_minimal(base_size = 18, base_family = "sans") +
  theme(
    plot.title      = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title      = element_text(size = 22),
    axis.text       = element_text(size = 18, colour = "black"),
    strip.text      = element_text(size = 18, face = "plain"),
    legend.title    = element_text(size = 18, face = "bold"),
    legend.text     = element_text(size = 18),
    plot.caption    = element_text(size = 18, hjust = 0.5),
    legend.key.size = unit(0.7, "cm"),
    panel.spacing   = unit(0.6, "lines")
  )

theme_no_strips <- theme(
  strip.text = element_blank(),
  strip.background = element_blank()
)

theme_no_x <- theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank()
)

evo_cols <- c(
  "Ancestral"    = "forestgreen",
  "Cold-evolved" = "steelblue",
  "Hot-evolved"  = "maroon"
)

sex_shapes <- c(
  "Female" = 16,
  "Male"   = 17
)

## ---------------------------------------------------------------------------
## Panel A: CO2
## ---------------------------------------------------------------------------
co2_panel <- co2.plot +
  facet_wrap(~ cge * Sex, scales = "free_y", nrow = 1) +
  scale_color_manual(values = evo_cols) +
  scale_shape_manual(values = sex_shapes) +
  labs(caption = NULL) +
  theme_bigfig +
  theme_no_x +
  theme(
    legend.position = "none"
  )

## ---------------------------------------------------------------------------
## Panel B: O2
## ---------------------------------------------------------------------------
o2_panel <- o2.plot +
  facet_wrap(~ cge * Sex, scales = "free_y", nrow = 1) +
  scale_color_manual(values = evo_cols) +
  scale_shape_manual(values = sex_shapes) +
  labs(caption = NULL) +
  theme_bigfig +
  theme_no_strips +
  theme_no_x +
  theme(
    legend.position = "none"
  )

## ---------------------------------------------------------------------------
## Panel C: RQ
## ---------------------------------------------------------------------------
rq_panel <- rq.plot +
  facet_wrap(~ cge * Sex, nrow = 1) +
  scale_color_manual(values = evo_cols) +
  scale_shape_manual(values = sex_shapes) +
  labs(caption = NULL) +
  theme_bigfig +
  theme_no_strips +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 18, hjust = 1)
  )

## ---------------------------------------------------------------------------
## Shared legend
## ---------------------------------------------------------------------------
legend_plot <- ggplot(
  emm_co2_dt,
  aes(x = Evo, y = emmean, color = Evo, shape = Sex)
) +
  geom_point(size = 5) +
  scale_color_manual(values = evo_cols, name = "Evo") +
  scale_shape_manual(values = sex_shapes, name = "Sex") +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(shape = 16, size = 5)
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(color = "grey20", size = 5)
    )
  ) +
  theme_void(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  )

shared_legend <- get_legend(legend_plot)

## significance footnote below legend
sig_note <- ggdraw() +
  draw_label(
    "**** p < 0.0001    *** p < 0.001    ** p < 0.01    * p < 0.05",
    x = 0.5, y = 0.5,
    hjust = 0.5, vjust = 0.5,
    size = 16
  )

legend_block <- plot_grid(
  shared_legend,
  sig_note,
  ncol = 1,
  rel_heights = c(1, 0.35)
)

## ---------------------------------------------------------------------------
## Stack panels vertically
## ---------------------------------------------------------------------------
panels_only <- plot_grid(
  co2_panel,
  o2_panel,
  rq_panel,
  labels = c("A", "B", "C"),
  label_size = 24,
  label_fontface = "bold",
  ncol = 1,
  align = "v",
  rel_heights = c(1.25, 1, 1.3)
)

final_resp_figure <- plot_grid(
  panels_only,
  legend_block,
  ncol = 1,
  rel_heights = c(1, 0.11)
)

ggsave(
  "../Plot/Figure7.png",
  final_resp_figure,
  width = 16,
  height = 12.5,
  units = "in",
  dpi = 600
)
