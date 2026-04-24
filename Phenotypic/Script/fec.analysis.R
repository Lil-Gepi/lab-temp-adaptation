setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Script/")
rm(list = ls());gc()
library(data.table)
setDTthreads(percent = 80)
getDTthreads()
library(stringr)
library(ggplot2)
library(lme4)
library(emmeans)
library(dplyr)

fec <- fread(file = "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Data/fecundity.csv",quote = F,header = T)
fec[cge == "Cold", cge:= "coldcge"][cge == "Hot", cge:= "hotcge"]
## 1. modeling  ----
fit <- lmer(log2(cum.fec) ~ cge*Evo + scale(flynr) + (1+cge|EvoRep), data = fec, REML = T)
# fit0 <- lmer(log2(cum.fec) ~ cge+Lab +cge:Lab + cge:Lab:Temp + scale(flynr) + (1|EvoRep), data = fec, REML = T)
# anova(fit, fit0)
# fit1 <- lmer(log2(cum.fec) ~ cge+Lab +cge:Lab + cge:Lab:Temp  + (1+cge|EvoRep), data = fec, REML = T)
# anova(fit, fit1) 
# fit2 <- lmer(log2(cum.fec) ~ cge+Lab +cge:Lab + cge:Lab:Temp + cge*Lab:Temp*scale(flynr) + (1+cge|EvoRep), data = fec, REML = T)
# anova(fit, fit2)
# fit3 <- lmer(log2(cum.fec) ~ cge+Lab +cge:Lab  + scale(flynr) + (1+cge|EvoRep), data = fec, REML = T)
# anova(fit, fit3)
# fit4 <- lmer(log2(cum.fec) ~ cge+Lab + cge:Lab:Temp + scale(flynr) + (1+cge|EvoRep), data = fec, REML = T)
# anova(fit, fit4)
# fit5 <- lmer(log2(cum.fec) ~ cge + cge:Lab + cge:Lab:Temp + scale(flynr) + (1+cge|EvoRep), data = fec, REML = T)
# anova(fit, fit5)
# fit6 <- lmer(log2(cum.fec) ~ cge+Lab +cge:Lab +Lab:Temp+ cge:Lab:Temp + scale(flynr) + (1+cge|EvoRep), data = fec, REML = T)
# anova(fit, fit6)
### all these alternative models are either significantly worse or not significantly better with additional parameters than `fit`
summary(fit)
coef(fit)
hist(resid(fit), breaks = 30)
plot(fit)
qqnorm(resid(fit)); qqline(resid(fit))

emmip(fit, ~ Evo | cge)

emm.res <- emmeans(fit, ~ Evo | cge)   
emm.res
contrast(emm.res,by = "cge", method = "pairwise")
# Get the effect size estimated in SD unit----
fit <- lmer(log2(cum.fec) ~ Evo*cge +  scale(flynr) + (1+cge|EvoRep), data = fec, REML = T)
emm.res <- emmeans(fit, ~ Evo*cge,by = "cge", adjust = "none",type = "response")
emm.res
contrast(emm.res,method = "pairwise")

eff_evo <- eff_size(emm.res,
                    sigma = sigma(fit),
                    edf   = df.residual(fit))         #cohen's d
(eff_evo <- as.data.frame(eff_evo))
fwrite(eff_evo, file = "./fecundity.eff.csv", quote = F, sep = ",", col.names = T)

## 2. plotting  ----
emm.plot <- summary(emm.res, type = "response")
setDT(emm.plot)
colnames(emm.plot)
setnames(emm.plot,"response", "emmean")
emm.plot[Evo == "Cold", Evo:="Cold-evolved"][Evo == "Hot", Evo:="Hot-evolved"]

pairs <- contrast(emm.res, by = "cge", method = "pairwise") %>%
  as.data.frame()

pairs <- pairs %>%
  mutate(sig = cut(p.value,
                   breaks = c(0, .0001, .001, .01, .05, Inf),
                   labels = c("****","***","**","*","n.s.")))
setDT(pairs)
pairs$contrast <- c("Ancestral - Hot-evolved","Ancestral - Cold-evolved","Hot-evolved - Cold-evolved",
                    "Ancestral - Hot-evolved","Ancestral - Cold-evolved","Hot-evolved - Cold-evolved")

lookup_x <- c(Ancestral = 1,
              `Cold-evolved` = 2,
              `Hot-evolved`  = 3)

ann <- pairs %>%
  transmute(cge,
            xstart = lookup_x[sub(" - .*",  "", contrast)],
            xend   = lookup_x[sub(".* - ",  "", contrast)],
            y      = 1.05 * tapply(emm.plot$emmean, emm.plot$cge, max)[cge] +
              as.numeric(factor(contrast, levels = rev(unique(contrast))))*3,
            label  = sig)

bg <- data.frame(
  cge  = c("coldcge", "hotcge"),      
  xmin = -Inf, xmax =  Inf,           
  ymin = -Inf, ymax =  Inf)

fec.plot <- ggplot(emm.plot, aes(Evo, emmean, color = Evo, fill = cge)) +
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  facet_wrap(~ cge) +
  ylim(0,80)+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),position = position_dodge(0.7),width = .2) +
  geom_point(size = 3,position = position_dodge(0.7)) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                               "Cold-evolved" = "steelblue",
                               "Hot-evolved"="maroon")) +
  scale_fill_manual(values = c("coldcge" = "lightblue",
                    "hotcge" = "lightpink"))+
  ylab("Eggs laid per Female") +
  ggtitle("Fecundity")+
  xlab("") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 0.6, vjust = 0.8),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0),
        legend.position = "none") +
  geom_segment(data = ann,
             aes(x = xstart, xend = xend, y = y*1.035,  yend = y*1.035),
             inherit.aes = FALSE, linewidth = .6) +
  geom_text(data = ann,
            aes(x = (xstart + xend)/2, y = ifelse(label == "n.s.", y*1.04+2, y*1.04+0.5), label = label),
            inherit.aes = FALSE, size = 5)+
   labs(caption = " ") 
# labs(caption = " **** p < 0.0001    *** p < 0.001    ** p < 0.01    * p < 0.05    n.s. non-significant") removed for plot aggregating
png("../Plot/Fecundity.png", width = 6, height = 6, units = "in", res = 400)
print(fec.plot)
dev.off()

## 3. export the emmeans res for the lab:temp proportion barplot ----
fwrite(emm.plot, file = "./fec.emm.csv", quote = F, sep = ",", col.names = T)
saveRDS(object = fec.plot, file = "./fec.plot.RDS")


## 4. effect size plot----
fec.eff <- fread("./fecundity.eff.csv")
fec.eff[,phenotype:="fecundity"]
setnames(fec.eff, 3, "effect.size")
fec.eff <- fec.eff[,.(phenotype, cge,contrast, effect.size)]

fec.eff.wide <- dcast(fec.eff, phenotype + cge ~ contrast,value.var = "effect.size")

fec.eff.wide[,lab.effect := min(abs(`Ancestral - Cold`),abs(`Ancestral - Hot`)) * sign(`Ancestral - Hot`*`Ancestral - Cold`), by = .(phenotype, cge)]
fec.eff.wide[,temp.effect := abs(`Cold - Hot`)]
fec.eff.wide <- fec.eff.wide[,.(phenotype,cge,lab.effect,temp.effect)]

fec.eff <- melt(fec.eff.wide, id.vars = c("phenotype", "cge"), variable.name = "effect.type", value.name = "effect.size")
#plotting
fec.eff.plot <- ggplot(data = fec.eff, aes(x = cge, y = effect.size, fill = effect.type))+
  geom_col(position = position_dodge(width = 0.6)) +
  labs(x = NULL,y="Cohen's d (effect size in unit of SDs)", title = "Fecundity") +
  scale_fill_manual(values = c("lab.effect"  = "purple",
                               "temp.effect" = "gold"),
                    name   = "",
                    labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 1, vjust = 1),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0))

png(paste0("../Plot/fecundity.effect.size.lab.temp.png"), width = 4, height = 6, units = "in", res = 600)
print(fec.eff.plot)
dev.off()

saveRDS(fec.eff.plot, "./fec.eff.plot.RDS")

