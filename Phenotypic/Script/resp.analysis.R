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

resp <- fread(file = "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Data/respiratory.csv",quote = F,header = T)
resp[,`:=`(co2.mg = (co2*FlyNum)/(NetWeight*(FlyNum-FlyDead)), o2.mg = (o2*FlyNum)/(NetWeight*(FlyNum-FlyDead)))]

bg <- data.frame(cge  = c("coldcge", "hotcge"), xmin = -Inf, xmax =  Inf, ymin = -Inf, ymax =  Inf)
## CO2 first ----

##1. inspect raw data----
pdf("../Plot/resp.CO2.raw.inspection.pdf")
ggplot(resp, aes(Evo, co2.mg, color = Evo, shape = Sex)) +
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  facet_wrap(~ Sex*cge, scales = "free_y") +
  geom_jitter(size = 1, width = 0.2) +
  # geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
  #               width = .25) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                                "Cold" = "steelblue",
                                "Hot"="maroon"))+
  ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1})))+
  xlab("") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 0.6, vjust = 0.8),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0),
        legend.position = "none")

ggplot(resp[cge== "hotcge",], aes(x=cycle, y=co2.mg,group = chamber, colour = Sex))+
  facet_wrap(~run)+geom_smooth(method = "loess",inewidth = 1.5)+
  ggtitle("CO2 emission in hotCGE")+
  xlab("Cycles")+ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1})))
## from the plot above we see that the male flies spent quite some time to finally calm down,
##  whereas the females are quite stable throughout the experiment
## male vs female pattern confirmed, also we do see that the two evolved tends to be above the ancestral 
ggplot(resp[cge== "coldcge",], aes(x=cycle, y=co2.mg,group = chamber, colour = Sex))+
  facet_wrap(~run)+geom_smooth(method = "loess",linewidth = 1.5)+
  ggtitle("CO2 emission in coldCGE")+
  xlab("Cycles")+ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1})))
## in coldcge, the male vs. female pattern is not so obvious anymore
ggplot(resp, aes(x=cycle, y=co2.mg,group = interaction(run,chamber), color = Evo))+
  facet_grid(~cge*Sex)+geom_smooth(method = "loess",linewidth = 1.5, alpha = 0.1)+
  ggtitle("CO2 emission by cycles")+xlab("Cycles")+ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1})))
dev.off()
## but once again, we see that ancestral population is below the two evolved populations.
## note that here this graph shows crazy noise, but it is just because that the flies
## don't do much so their co2 emission is 1/4 as in hot. Therefore more visually susceptible to 
## noise or technical variation

ggplot(resp, aes(x=cge, y=co2.mg,group = interaction(run,chamber), color = EvoRep))+
  facet_grid(~Sex)+geom_jitter()+
  ggtitle("CO2 emission in coldCGE")+xlab("Cycles")+ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1})))
## in summary, we have to model the cge (fixed), sex (fixed), Evolution (fixed, of interest)
## the cooling pattern, seen only in Male samples from teh hotcge will be included in the random effect. 
## it would have its random slop and intercept.
## each run could also have a different intercept, so we should model that by random effect.
## probably each EvoRep has its own intercept, but probably different among cge, I need to test it.
## I am not sure about the effect of FlyNum, so I will test whether I can drop it.

## 2. fit model, remove extreme outlier if necessary ----
fit <- lmer(data = resp, formula = log2(co2.mg) ~ cge*Evo*Sex + scale(FlyNum)+ (1|run) +(1|cycle:run:Sex:cge) + (1+cge|EvoRep))

pdf("../Plot/resp.CO2.model.fit.raw.pdf")
hist(resid(fit), breaks = 300)
summary(fit)
coef(fit)
plot(fit)
qqnorm(resid(fit)); qqline(resid(fit))
dev.off()
resp[, residuals := residuals(fit)]
threshold <- 2 * sd(resp$residuals)
resp <- resp[abs(residuals) < threshold,]

fit<-lmer(data = resp, formula = log2(co2.mg) ~ cge*Evo*Sex + scale(FlyNum)+ (1|run) +(1|cycle:run:Sex:cge) + (1+cge|EvoRep))
pdf("../Plot/resp.CO2.model.fit.cleaned.pdf")
hist(resid(fit), breaks = 100)
summary(fit)
coef(fit)
plot(fit)
qqnorm(resid(fit)); qqline(resid(fit))
dev.off()

# emmip(fit, ~ Sex:Evo | cge, CIs = T)
emm.res <- emmeans(fit, ~ Evo | cge*Sex)
emm.res
contrast(emm.res, method = "pairwise")

# Get the effect size estimated in SD unit----
eff_evo <- eff_size(emm.res,
                    sigma = sigma(fit),
                    edf   = df.residual(fit),method = )          # "d" = Cohen; "hedges_g" available too
eff_evo <- as.data.frame(eff_evo)
fwrite(eff_evo, file = "./resp.CO2.eff.csv", quote = F, sep = ",", col.names = T)


## 3. plotting  ----
emm.plot <- summary(emm.res, type = "response")
setDT(emm.plot)
setnames(emm.plot,"response", "emmean")
emm.plot[Evo == "Cold", Evo:="Cold-evolved"][Evo == "Hot", Evo:="Hot-evolved"]
bg <- data.frame(
  cge  = c("coldcge", "hotcge"),      
  xmin = -Inf, xmax =  Inf,           
  ymin = -Inf, ymax =  Inf)

pairs <- contrast(emm.res, method = "pairwise") %>%
  as.data.frame()

pairs <- pairs %>%
  mutate(sig = cut(p.value,
                   breaks = c(0, .0001, .001, .01, .05, Inf),
                   labels = c("****","***","**","*","n.s.")))
setDT(pairs)
pairs$contrast <- c("Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved")

lookup_x <- c(Ancestral = 1,`Cold-evolved` = 2,`Hot-evolved`  = 3)
panel_span <- emm.plot |>
  group_by(cge, Sex) |>
  summarise(ymax  = max(upper.CL),                # highest point we draw
    span  = diff(range(c(lower.CL, upper.CL))),  # height of the panel
    .groups = "drop")

ann <- pairs |>
  mutate(xstart = lookup_x[str_remove(contrast, " - .*")],
    xend   = lookup_x[str_remove(contrast, ".* - ")],
    rank   = as.numeric(factor(contrast,
                               levels = rev(unique(contrast))))) |>
  left_join(panel_span, by = c("cge","Sex")) |>
  mutate(y     = ymax + (rank+0.5)*0.2*span,
    label = sig
  )


png("../Plot/Respiratory.CO2.png", width = 8, height = 6, units = "in", res = 400)
co2.plot <- ggplot(emm.plot, aes(Evo, emmean, color = Evo, shape = Sex)) +
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  facet_wrap(~ cge*Sex, scales = "free_y", nrow= 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = .2) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                                "Cold-evolved" = "steelblue",
                                "Hot-evolved"="maroon"))+
  ylab(expression(V[CO2]~~(mu*L*h^{-1}*mg^{-1})))+
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  geom_segment(data = ann,
               aes(x = xstart, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = .6) +
  geom_text(data = ann,
            aes(x = (xstart + xend)/2, y = ifelse(label == "n.s.", y + 0.11*span,y + 0.01*span), label = label),
            inherit.aes = FALSE, size = 5) +
  labs(caption = " ")+
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 0.6, vjust = 0.8),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0),
        legend.position = "none")
print(co2.plot)
dev.off()
## 4. export the emmeans res for the lab:temp proportion barplot ----
fwrite(emm.plot, file = "./resp.CO2.emm.csv", quote = F, sep = ",", col.names = T)
saveRDS(co2.plot, file = "./co2.plot.RDS")

## O2 then ----
rm(list = ls());gc()
resp <- fread(file = "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Data/respiratory.csv",quote = F,header = T)
resp[,`:=`(co2.mg = (co2*FlyNum)/(NetWeight*(FlyNum-FlyDead)), o2.mg = (o2*FlyNum)/(NetWeight*(FlyNum-FlyDead)))]
bg <- data.frame(cge  = c("coldcge", "hotcge"), xmin = -Inf, xmax =  Inf, ymin = -Inf, ymax =  Inf)

##1. inspect raw data----
## the pattern would be most likely the same as CO2
pdf("../Plot/resp.O2.raw.inspection.pdf")
ggplot(resp, aes(Evo, o2.mg, color = Evo, shape = Sex)) +
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  facet_wrap(~ Sex*cge, scales = "free_y") +
  geom_jitter(size = 1, width = 0.2) +
  # geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
  #               width = .25) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                                "Cold" = "steelblue",
                                "Hot"="maroon"))+
  ylab(expression(V[O2]~~(mu*L*h^{-1}*mg^{-1})))+
  xlab("") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 0.6, vjust = 0.8),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0),
        legend.position = "none")

ggplot(resp[cge== "hotcge",], aes(x=cycle, y=o2.mg,group = chamber, colour = Sex))+
  facet_wrap(~run)+geom_smooth(method = "loess",inewidth = 1.5)+
  ggtitle("O2 emission in hotCGE")+
  xlab("Cycles")+ylab(expression(V[O2]~~(mu*L*h^{-1}*mg^{-1})))
## from the plot above we see that the male flies spent quite some time to finally calm down,
##  whereas the females are quite stable throughout the experiment
## male vs female pattern confirmed, also we do see that the two evolved tends to be above the ancestral 
ggplot(resp[cge== "coldcge",], aes(x=cycle, y=o2.mg,group = chamber, colour = Sex))+
  facet_wrap(~run)+geom_smooth(method = "loess",linewidth = 1.5)+
  ggtitle("O2 emission in coldCGE")+
  xlab("Cycles")+ylab(expression(V[O2]~~(mu*L*h^{-1}*mg^{-1})))
## in coldcge, the male vs. female pattern is not so obvious anymore
ggplot(resp, aes(x=cycle, y=o2.mg,group = interaction(run,chamber), color = Evo))+
  facet_grid(~cge*Sex)+geom_smooth(method = "loess",linewidth = 1.5, alpha = 0.1)+
  ggtitle("O2 emission by cycles")+xlab("Cycles")+ylab(expression(V[O2]~~(mu*L*h^{-1}*mg^{-1})))
dev.off()


## 2. fit model, remove extreme outlier if necessary ----
fit<-lmer(data = resp, formula = log2(o2.mg) ~ cge*Evo*Sex + scale(FlyNum)+ (1|run) +(1|cycle:run:Sex:cge) + (1+cge|EvoRep))
pdf("../Plot/resp.O2.model.fit.raw.pdf")
hist(resid(fit), breaks = 300)
summary(fit)
coef(fit)
plot(fit)
qqnorm(resid(fit)); qqline(resid(fit))
dev.off()
resp[, residuals := residuals(fit)]
threshold <- 2 * sd(resp$residuals)
resp<- resp[abs(residuals) < threshold,]

fit<-lmer(data = resp, formula = log2(o2.mg) ~ cge*Evo*Sex + scale(FlyNum)+ (1|run) +(1|cycle:run:Sex:cge) + (1+cge|EvoRep))
pdf("../Plot/resp.O2.model.fit.cleaned.pdf")
hist(resid(fit), breaks = 100)
summary(fit)
coef(fit)
plot(fit)
qqnorm(resid(fit)); qqline(resid(fit))
dev.off()

# emmip(fit, ~ Sex:Evo | cge, CIs = T)
emm.res <- emmeans(fit, ~ Evo | cge*Sex)
emm.res
contrast(emm.res, method = "pairwise")

# Get the effect size estimated in SD unit----
eff_evo <- eff_size(emm.res,
                    sigma = sigma(fit),
                    edf   = df.residual(fit),method = )          # "d" = Cohen; "hedges_g" available too
eff_evo <- as.data.frame(eff_evo)
fwrite(eff_evo, file = "./resp.O2.eff.csv", quote = F, sep = ",", col.names = T)


## 3. plotting  ----
emm.plot <- summary(emm.res, type = "response")
setDT(emm.plot)
setnames(emm.plot,"response", "emmean")
emm.plot[Evo == "Cold", Evo:="Cold-evolved"][Evo == "Hot", Evo:="Hot-evolved"]
bg <- data.frame(
  cge  = c("coldcge", "hotcge"),      
  xmin = -Inf, xmax =  Inf,           
  ymin = -Inf, ymax =  Inf)

pairs <- contrast(emm.res, method = "pairwise") %>%
  as.data.frame()

pairs <- pairs %>%
  mutate(sig = cut(p.value,
                   breaks = c(0, .0001, .001, .01, .05, Inf),
                   labels = c("****","***","**","*","n.s.")))
setDT(pairs)
pairs$contrast <- c("Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved")

lookup_x <- c(Ancestral = 1,
              `Cold-evolved` = 2,
              `Hot-evolved`  = 3)
panel_span <- emm.plot |>
  group_by(cge, Sex) |>
  summarise(
    ymax  = max(upper.CL),                # highest point we draw
    span  = diff(range(c(lower.CL, upper.CL))),  # height of the panel
    .groups = "drop"
  )

ann <- pairs |>
  mutate(
    xstart = lookup_x[str_remove(contrast, " - .*")],
    xend   = lookup_x[str_remove(contrast, ".* - ")],
    rank   = as.numeric(factor(contrast,
                               levels = rev(unique(contrast))))
  ) |>
  left_join(panel_span, by = c("cge","Sex")) |>
  mutate(
    y     = ymax + (rank+0.5)*0.2*span,  # 0.06 = 6 % – tweak to taste
    label = sig
  )


png("../Plot/Respiratory.O2.png", width = 8, height = 6, units = "in", res = 400)
o2.plot <- ggplot(emm.plot, aes(Evo, emmean, color = Evo, shape = Sex)) +
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  facet_wrap(~ cge*Sex, scales = "free_y", nrow= 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = .2) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                                "Cold-evolved" = "steelblue",
                                "Hot-evolved"="maroon"))+
  ylab(expression(V[O2]~~(mu*L*h^{-1}*mg^{-1})))+
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  geom_segment(data = ann,
               aes(x = xstart, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = .6) +
  geom_text(data = ann,
            aes(x = (xstart + xend)/2, y = ifelse(label == "n.s.", y + 0.11*span,y + 0.01*span), label = label),
            inherit.aes = FALSE, size = 5) +
  labs(caption = " ")+
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 0.6, vjust = 0.8),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0),
        legend.position = "none")
print(o2.plot)
dev.off()
## 4. export the emmeans res for the lab:temp proportion barplot ----
fwrite(emm.plot, file = "./resp.O2.emm.csv", quote = F, sep = ",", col.names = T)
saveRDS(o2.plot, "./o2.plot.RDS")

## RER CO2/O2 ratio ----
rm(list = ls());gc()
resp <- fread(file = "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Data/respiratory.csv",quote = F,header = T)
resp[,`:=`(co2.mg = (co2*FlyNum)/(NetWeight*(FlyNum-FlyDead)), o2.mg = (o2*FlyNum)/(NetWeight*(FlyNum-FlyDead)))]
bg <- data.frame(cge  = c("coldcge", "hotcge"), xmin = -Inf, xmax =  Inf, ymin = -Inf, ymax =  Inf)

##1. inspect raw data----

pdf("../Plot/resp.RER.raw.inspection.pdf")
ggplot(resp, aes(Evo, RER, color = Evo, shape = Sex)) +
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  facet_wrap(~ Sex*cge, scales = "free_y") +
  geom_jitter(size = 1, width = 0.2) +
  # geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
  #               width = .25) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                                "Cold" = "steelblue",
                                "Hot"="maroon"))+
  ylab(expression(V[CO2] / V[O2]))+
  xlab("") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 0.6, vjust = 0.8),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0),
        legend.position = "none")

ggplot(resp[cge== "hotcge",], aes(x=cycle, y=RER,group = chamber, colour = Sex))+
  facet_wrap(~run)+geom_smooth(method = "loess",inewidth = 1.5)+
  ggtitle("RER in hotCGE")+
  xlab("Cycles")+ylab(expression(V[CO2] / V[O2]))
## from the plot above we see that the male flies spent quite some time to finally calm down,
##  whereas the females are quite stable throughout the experiment
## male vs female pattern confirmed, also we do see that the two evolved tends to be above the ancestral 
ggplot(resp[cge== "coldcge",], aes(x=cycle, y=RER,group = chamber, colour = Sex))+
  facet_wrap(~run)+geom_smooth(method = "loess",linewidth = 1.5)+
  ggtitle("RER in coldCGE")+
  xlab("Cycles")+ylab(expression(V[CO2] / V[O2]))
## in coldcge, the male vs. female pattern is not so obvious anymore
ggplot(resp, aes(x=cycle, y=RER,group = interaction(run,chamber), color = Evo))+
  facet_grid(~cge*Sex)+geom_smooth(method = "loess",linewidth = 1.5, alpha = 0.1)+
  ylim(0.6, 1.2)+
  ggtitle("RER by cycles")+xlab("Cycles")+ylab(expression(RER = V[CO2] / V[O2]))
dev.off()


## 2. fit model, remove extreme outlier if necessary ----
fit<-lmer(data = resp, formula = log2(RER) ~ cge*Evo*Sex + (1|run) +(1|cycle:run:cge) + (1+cge|EvoRep))
pdf("../Plot/resp.RER.model.fit.cleaned.pdf")
hist(resid(fit), breaks = 100)
summary(fit)
coef(fit)
plot(fit)
qqnorm(resid(fit)); qqline(resid(fit))
dev.off()

# emmip(fit, ~ Sex:Evo | cge, CIs = T)
emm.res <- emmeans(fit, ~ Evo | cge*Sex)
emm.res
contrast(emm.res, method = "pairwise")


## 3. plotting  ----
emm.plot <- summary(emm.res, type = "response")
setDT(emm.plot)
setnames(emm.plot,"response", "emmean")
emm.plot[Evo == "Cold", Evo:="Cold-evolved"][Evo == "Hot", Evo:="Hot-evolved"]
bg <- data.frame(
  cge  = c("coldcge", "hotcge"),      
  xmin = -Inf, xmax =  Inf,           
  ymin = -Inf, ymax =  Inf)

png("../Plot/Respiratory.RER.png", width = 8, height = 6, units = "in", res = 400)
rer.plot <- ggplot(emm.plot, aes(Evo, emmean, color = Evo, shape = Sex)) +
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  facet_wrap(~ cge*Sex, nrow =1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = .2) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                                "Cold-evolved" = "steelblue",
                                "Hot-evolved"="maroon"))+
  ylab(expression(RER = V[CO2] / V[O2]))+
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.20))) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 0.6, vjust = 0.8),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0),
        legend.position = "none")+
  labs(caption = "**** p < 0.0001    *** p < 0.001    ** p < 0.01    * p < 0.05    n.s. non-significant")

print(rer.plot)
dev.off()
## 4. export the emmeans res for the lab:temp proportion barplot ----
fwrite(emm.plot, file = "./resp.RER.emm.csv", quote = F, sep = ",", col.names = T)
saveRDS(rer.plot, "./rer.plot.RDS")





## Lastly, dry weight ----
rm(list = ls());gc()
resp <- fread(file = "~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Data/respiratory.csv",quote = F,header = T)
resp[,`:=`(weight.per = NetWeight / FlyNum)]
resp <- unique(resp[,.(cge, run, Evo, EvoRep, Sex, TechRep, FlyNum, NetWeight, weight.per)])
bg <- data.frame(cge  = c("coldcge", "hotcge"), xmin = -Inf, xmax =  Inf, ymin = -Inf, ymax =  Inf)

##1. inspect raw data----
pdf("../Plot/resp.Weight.raw.inspection.pdf")
ggplot(resp, aes(Evo, weight.per, color = Evo, shape = Sex)) +
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  facet_wrap(~ cge*Sex, scales = "free_y") +
  geom_jitter(size = 3, width = 0.2) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                                "Cold" = "steelblue",
                                "Hot"="maroon"))+
  ylab(expression(Weight~(mg)))+
  xlab("") +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 15, hjust = 0.6, vjust = 0.8),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0),
        legend.position = "none")

dev.off()


## 2. fit model, remove extreme outlier if necessary ----
fit<-lmer(data = resp, formula = log2(weight.per) ~ cge*Evo*Sex + (1|run)+ (0+cge|EvoRep))
fit1<-lmer(data = resp, formula = log2(weight.per) ~ cge*Evo*Sex + (0+cge|EvoRep))
anova(fit1, fit)

fit2<-lmer(data = resp, formula = log2(weight.per) ~ cge*Evo*Sex + (1|run) )
anova(fit2, fit)
fit3<-lmer(data = resp, formula = log2(weight.per) ~ cge*Evo*Sex + (1|run)+ (1|EvoRep))
anova(fit3, fit)

fit<-lmer(data = resp, formula = log2(weight.per) ~ cge*Evo*Sex + (1|run)+ (0+cge|EvoRep))

pdf("../Plot/resp.Weight.model.fit.raw.pdf")
hist(resid(fit), breaks = 50)
summary(fit)
coef(fit)
plot(fit)
qqnorm(resid(fit)); qqline(resid(fit))
dev.off()

# emmip(fit, ~ Sex:Evo | cge, CIs = T)
emm.res <- emmeans(fit, ~ Evo | cge*Sex)
emm.res
contrast(emm.res, method = "pairwise")

# Get the effect size estimated in SD unit----
eff_evo <- eff_size(emm.res,
                    sigma = sigma(fit),
                    edf   = df.residual(fit),method = )          # "d" = Cohen; "hedges_g" available too
(eff_evo <- as.data.frame(eff_evo))
fwrite(eff_evo, file = "./resp.weight.eff.csv", quote = F, sep = ",", col.names = T)



## 3. plotting  ----
emm.plot <- summary(emm.res, type = "response")
setDT(emm.plot)
setnames(emm.plot,"response", "emmean")
emm.plot[Evo == "Cold", Evo:="Cold-evolved"][Evo == "Hot", Evo:="Hot-evolved"]
bg <- data.frame(
  cge  = c("coldcge", "hotcge"),      
  xmin = -Inf, xmax =  Inf,           
  ymin = -Inf, ymax =  Inf)

pairs <- contrast(emm.res, method = "pairwise") %>%
  as.data.frame()

pairs <- pairs %>%
  mutate(sig = cut(p.value,
                   breaks = c(0, .05, 0.10,1),
                   labels = c("*","•","n.s.")))
setDT(pairs)
pairs$contrast <- c("Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved",
                    "Ancestral - Cold-evolved","Ancestral - Hot-evolved","Cold-evolved - Hot-evolved")

panel_span <- emm.plot |>
  group_by(cge, Sex) |>
  summarise(
    ymax  = max(upper.CL),                # highest point we draw
    span  = diff(range(c(lower.CL, upper.CL))),  # height of the panel
    .groups = "drop"
  )

y_span <- diff(range(c(emm.plot$lower.CL, emm.plot$upper.CL))) # total range

lookup_x <- c(Ancestral = 1,
              `Cold-evolved` = 2,
              `Hot-evolved`  = 3)

ann <- pairs |>
  mutate(
    xstart = lookup_x[str_remove(contrast, " - .*")],
    xend   = lookup_x[str_remove(contrast, ".* - ")],
    rank   = as.numeric(factor(contrast,
                               levels = rev(unique(contrast))))
  ) |>
  left_join(panel_span, by = c("cge","Sex")) |>
  mutate(
    y     = ymax + (rank+0.5)*0.1*y_span,  # 0.06 = 6 % – tweak to taste
    label = sig
  )


png("../Plot/Respiratory.Weight.png", width = 8, height = 6, units = "in", res = 600)
ggplot(emm.plot, aes(Evo, emmean, color = Evo, shape = Sex)) +
  geom_rect(aes(fill = cge),data = bg, inherit.aes = FALSE,
            xmin=-Inf, xmax=Inf, ymin = -Inf, ymax=Inf, alpha = 0.3) +
  facet_wrap(~ cge*Sex,nrow =1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = .2) +
  scale_fill_manual(values = c("coldcge" = "lightblue","hotcge" = "lightpink")) +
  scale_color_manual(values = c("Ancestral" ="forestgreen",
                                "Cold-evolved" = "steelblue",
                                "Hot-evolved"="maroon"))+
  ylab(expression(Weight~(mg)))+
  xlab("") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  geom_segment(data = ann,
               aes(x = xstart, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = .6) +
  geom_text(data = ann,
            aes(x = (xstart + xend)/2, y = ifelse(label == "n.s.", y + 0.07*y_span, y + 0.05*y_span), label = label),
            inherit.aes = FALSE, size = 5) +
  labs(caption = " * p < 0.05    • p < 0.10    n.s. non-significant")+
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0),
        legend.position = "none")
dev.off()
 ## 4. export the emmeans res for the lab:temp proportion barplot ----
fwrite(emm.plot, file = "./resp.Weight.emm.csv", quote = F, sep = ",", col.names = T)




## plotting Figure 7, CO2 O2 RER----
co2.plot <- readRDS("co2.plot.RDS")
o2.plot <- readRDS("o2.plot.RDS")
rer.plot <- readRDS("rer.plot.RDS")

library(cowplot)
top_row <- plot_grid(co2.plot, o2.plot,  labels = c("A", "B"),
                     label_size = 20,label_fontface = "bold",
                     ncol = 2,  align = "hv",  rel_widths = c(1, 1))

bottom_row <- plot_grid(  rer.plot,  labels = c("C"), label_size = 20,label_fontface = "bold", ncol = 1)

final_resp_figure <- plot_grid(  top_row, bottom_row,  ncol = 1,  rel_heights = c(1, 0.7))
print(final_resp_figure)
ggsave("../Plot/Figure7.png", final_resp_figure, width = 14, height = 11, dpi = 600)

