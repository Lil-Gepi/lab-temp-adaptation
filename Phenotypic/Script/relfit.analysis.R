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
library(minpack.lm)
rf <- fread(file = "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Data/relativefitness.csv",quote = F,header = T)
bg <- data.frame(cge  = c("coldcge", "hotcge"), xmin = -Inf, xmax =  Inf, ymin = -Inf, ymax =  Inf)
rf[, vial:= paste(cge, Sample, sep = ".")]
## raw data inspection ----
### red eclosion by Day
pdf("../Plot/relfit.raw.inspection.pdf")
ggplot(data=rf, aes(x=ColDay, y=Red_Nr, group=Sample, col = Evo)) +
  facet_wrap(~Evo*cge,ncol = 2, scale = "free_x")+
  ggtitle("Red eye fly eclosion by day")+
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  geom_line()+ scale_color_manual(values = c("forestgreen", "steelblue","maroon"))+theme_minimal()
### white eclosion by Day
ggplot(data=rf, aes(x=ColDay, y=White_Nr, group=Sample, col = Evo)) +
  facet_wrap(~Evo*cge,ncol = 2, scale = "free_x")+
  ggtitle("White eye fly eclosion by day")+
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  geom_line()+ scale_color_manual(values = c("forestgreen", "steelblue","maroon"))+theme_minimal()
### Cumulative red eclosion by Day
ggplot(rf, aes(x=ColDay, y=Cumulative_Red_Nr, group=Sample, col = Evo))+
  facet_wrap(~Evo*cge, ncol=2, scale = "free_x")+
  ggtitle("Cumulative red eye fly eclosion by day")+
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  geom_line()+scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("forestgreen", "steelblue","maroon"))+theme_minimal()
### Cumulative white eclosion by Day
ggplot(rf, aes(x=ColDay, y=Cumulative_White_Nr, group=Sample, col = Evo))+
  facet_wrap(~Evo*cge, ncol=2, scale = "free_x")+
  ggtitle("Cumulative white eye fly eclosion by day")+
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  geom_line()+
  #geom_smooth(method = "loess",linewidth = 1.5)+
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("forestgreen", "steelblue","maroon"))+theme_minimal()

ggplot(rf, aes(x=ColDay, y=ratio, group=Sample, col = Evo))+
  facet_wrap(~Evo*cge, ncol=2, scale = "free_x")+
  ggtitle("Red:White ratio")+
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  geom_line()+
  #geom_smooth(method = "loess",linewidth = 1.5)+
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("forestgreen", "steelblue","maroon"))+theme_minimal()
dev.off()


## 2. test the relative fitness (ratio) ----
png(file = "../Plot/RelativeFitness.raw.png", width = 7, height = 5, units = "in", res = 600)
ggplot(data = rf[(cge=="hotcge"&ColDay == 19)|(cge=="coldcge"&ColDay == 47),], aes( x =Evo, y = ratio,  color = Evo))+
  facet_wrap(~ cge, scale = "free_x")+
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  geom_boxplot()+geom_jitter(width = 0.1)+
  ggtitle("Relative fitness in ColdCGE")+
  ylab("red-eye %")+
  xlab("")+
  scale_color_manual(values = c("forestgreen","steelblue","maroon"))+
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )
dev.off()


rf.plot <- rf[!is.na(ratio) & ratio > 0]
rf.plot[, ColDay := factor(ColDay)]

ggplot(rf.plot, aes(ColDay, ratio, fill = Evo, colour = Evo)) +
  facet_wrap(~ cge, nrow = 2, scales = "free_x") +
  geom_rect(aes(fill = cge), data = bg, inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .25) +
  geom_boxplot(position = position_dodge(0.8),width = .7, outlier.shape = NA) +
  geom_jitter(position = position_dodge(0.8),  size = 1.4, alpha = .8) +
  scale_fill_manual(values = c( coldcge = "lightblue", hotcge = "lightpink")) +
  scale_colour_manual(values = c(Ancestral = "forestgreen", Cold = "steelblue",
                                 Hot = "maroon")) +
  labs(title = "Relative fitness by day",
       x = "Day post egg-laying",
       y = "Red-eye %") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

rf_endpoint <- rf[, total := Red_Nr+White_Nr][(cge=="hotcge"&ColDay == 20)|(cge=="coldcge"&ColDay == 47),.(cge, Evo, EvoRep, TechRep, Sample,vial, ratio, Cumulative_Red_Nr,Cumulative_White_Nr)]
# the end points were chosen so that we are sure there are no second generation flies among the eclosed
# rf_endpoint[,ratio := ratio /100]
colnames(rf_endpoint)
fit <- glmer(cbind(Cumulative_Red_Nr,Cumulative_White_Nr ) ~ Evo * cge + (0+cge|EvoRep) +(1|Sample) ,family = binomial, data = rf_endpoint)
fit1 <- glmer(cbind(Cumulative_Red_Nr,Cumulative_White_Nr ) ~ Evo * cge + (1|EvoRep) +(1|Sample) ,family = binomial, data = rf_endpoint)
fit2 <- glmer(cbind(Cumulative_Red_Nr,Cumulative_White_Nr ) ~ Evo * cge + (1+cge|EvoRep)+(1|Sample),family = binomial, data = rf_endpoint)
anova(fit, fit1, fit2)

hist(resid(fit), breaks = 50)
summary(fit)
coef(fit)
plot(fit)
qqnorm(resid(fit)); qqline(resid(fit))

emm.res<-emmeans(fit, pairwise ~ Evo*cge, by = "cge", adjust = "none",type = "response")#transform back  odd-ratio are multiplicative
emmip(fit, ~ Evo | cge, CIs = T)
emm.res

emm_link <- emmeans(fit, ~ Evo | cge, type = "link")
ct <- contrast(emm_link, method = "trt.vs.ctrlk", ref = "Ancestral", adjust = "none")
k <- sqrt(3) / pi
ct_df <- as.data.frame(ct)
ct_df <- ct_df %>%
  mutate(
    d        = estimate * k,
    SE_d     = SE * k,
    lower_d  = (estimate - 1.96 * SE) * k,
    upper_d  = (estimate + 1.96 * SE) * k
  )


pd <- position_dodge(width = 0.4)

ggplot(ct_df, aes(x = cge, y = d, color = contrast, group = contrast)) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(aes(ymin = lower_d, ymax = upper_d), width = 0.5, linewidth = 1,position = pd) +
  labs(x = NULL,y = "Effect size (Cohen's d, SD units)",
    title = "Evolved vs. Ancestral (within each CGE)") +
  theme_classic()

eff_evo <- eff_size(emm.res,
                    sigma = sigma(fit),
                    edf   = df.residual(fit))          # "d" = Cohen; "hedges_g" available too
(eff_evo <- as.data.frame(eff_evo))
fwrite(eff_evo, file = "./relfit.eff.csv", quote = F, sep = ",", col.names = T)


emm.plot <- summary(emm.res, type = "response")
emm.plot <- emm.plot$emmeans
setDT(emm.plot)


### plotting ----
setnames(emm.plot,"prob", "emmean")
setnames(emm.plot,"asymp.LCL", "lower.CL")
setnames(emm.plot,"asymp.UCL", "upper.CL")

emm.plot[Evo == "Cold", Evo:="Cold-evolved"][Evo == "Hot", Evo:="Hot-evolved"]
bg <- data.frame(
  cge  = c("coldcge", "hotcge"),      
  xmin = -Inf, xmax =  Inf,           
  ymin = -Inf, ymax =  Inf)

pairs <- as.data.frame(emm.res$contrasts)
pairs <- pairs %>%
  mutate(sig = cut(p.value,
                   breaks = c(0, .0001, .001, .01, .05, Inf),
                   labels = c("****","***","**","*","n.s.")))
setDT(pairs)
pairs$contrast <- c("Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved")

lookup_x <- c(Ancestral = 1,`Cold-evolved` = 2,`Hot-evolved`  = 3)
panel_span <- emm.plot |>
  group_by(cge) |>
  summarise(ymax  = max(upper.CL), span  = diff(range(c(lower.CL, upper.CL))), .groups = "drop")

ann <- pairs |>
  mutate(xstart = lookup_x[str_remove(contrast, " - .*")],
         xend   = lookup_x[str_remove(contrast, ".* - ")],
         rank   = as.numeric(factor(contrast,levels = rev(unique(contrast))))) |>
  left_join(panel_span, by = "cge") |>
  mutate(y     = ymax + (rank+0.5)*0.08*span,label = sig)


rf.plot <- ggplot(emm.plot, aes( x= Evo, emmean, color = Evo)) +
  facet_wrap(~ cge, nrow = 1) +
  geom_rect(aes(fill = cge), data = bg, inherit.aes = FALSE,
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .25) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),position = position_dodge(0.7),width = .2) +
  geom_point(size = 3,position = position_dodge(0.7)) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                                "Cold-evolved" = "steelblue",
                                "Hot-evolved"="maroon"))+
  xlab("") +ylab("Red-eye offspring proportion")+ggtitle("Relative fitness")+
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05))) +
  geom_segment(data = ann,
               aes(x = xstart, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = .6) +
  geom_text(data = ann,
            aes(x = (xstart + xend)/2, y = ifelse(label == "n.s.", y + 0.05*span,y + 0.01*span), label = label),
            inherit.aes = FALSE, size = 5) +
  labs(caption = " **** p < 0.0001    *** p < 0.001    ** p < 0.01    * p < 0.05    n.s. non-significant")+
  theme_minimal(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x        = element_text(angle = 15, hjust = 0.6, vjust = 0.8),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18)
  )

png("../Plot/RelativeFitness.png", width = 6, height = 6, units = "in", res = 600)
print(rf.plot)
dev.off()

## 2.4. export the emmeans res for the lab:temp proportion barplot ----
fwrite(emm.plot, file = "./relfit.emm.csv", quote = F, sep = ",", col.names = T)
saveRDS(object = rf.plot, file = "./rf.plot.RDS")

## 4. effect size plot----
rf.eff <- fread("./relfit.eff.csv")
rf.eff[,phenotype:="relfit"]
rf.eff[contrast=="(Ancestral - Cold)", contrast:="Ancestral - Cold"][
  contrast=="(Ancestral - Hot)", contrast:="Ancestral - Hot"][
    contrast=="(Cold - Hot)", contrast:="Cold - Hot"]

setnames(rf.eff, 3, "effect.size")
rf.eff <- rf.eff[,.(phenotype, cge,contrast, effect.size)]

rf.eff.wide <- dcast(rf.eff, phenotype + cge ~ contrast,value.var = "effect.size")

rf.eff.wide[,lab.effect := min(abs(`Ancestral - Cold`),abs(`Ancestral - Hot`)) * sign(`Ancestral - Hot`*`Ancestral - Cold`), by = .(phenotype, cge)]
rf.eff.wide[,temp.effect := abs(`Cold - Hot`)]
rf.eff.wide <- rf.eff.wide[,.(phenotype,cge,lab.effect,temp.effect)]

rf.eff <- melt(rf.eff.wide, id.vars = c("phenotype", "cge"), variable.name = "effect.type", value.name = "effect.size")
#plotting
rf.eff.plot <- ggplot(data = rf.eff, aes(x = cge, y = effect.size, fill = effect.type))+
  geom_col(position = position_dodge(width = 0.6)) +
  labs(x = NULL,y="Cohen's d (effect size in unit of SDs)", title = "Relative Fitness") +
  scale_fill_manual(values = c("lab.effect"  = "purple",
                               "temp.effect" = "gold"),
                    name   = "",
                    labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 1, vjust = 1),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0))

png(paste0("../Plot/relative.fitness.effect.size.lab.temp.png"), width = 4, height = 6, units = "in", res = 600)
print(rf.eff.plot)
dev.off()

saveRDS(rf.eff.plot, "./rf.eff.plot.RDS")

