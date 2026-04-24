setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Script/")
rm(list = ls());gc()
library(data.table)
setDTthreads(percent = 80)
getDTthreads()
library(ggplot2)
development <- fread("./development.eff.csv"); development[,phenotype:="development"]; setnames(development, "estimate","effect.size" )
fecundity <- fread("./fecundity.eff.csv");fecundity[,phenotype:="fecundity"];setnames(fecundity, 3, "effect.size")
relfit <- fread("./relfit.eff.csv");relfit[,phenotype:="relfit"];setnames(relfit, c("estimate","asymp.LCL","asymp.UCL"),c("effect.size","lower.CL","upper.CL"))
resp.co2 <- fread("./resp.CO2.eff.csv")
resp.o2 <- fread("./resp.O2.eff.csv")
weight <- fread("./resp.Weight.eff.csv")
resp.co2.male <- resp.co2[Sex == "Male",];resp.co2.male[,Sex := NULL];resp.co2.male[,phenotype:="resp.co2.male"]
resp.co2.female <- resp.co2[Sex == "Female",];resp.co2.female[,Sex := NULL];resp.co2.female[,phenotype:="resp.co2.female"]
resp.o2.male <- resp.o2[Sex == "Male",];resp.o2.male[,Sex := NULL];resp.o2.male[,phenotype:="resp.o2.male"]
resp.o2.female <- resp.o2[Sex == "Female",];resp.o2.female[,Sex := NULL];resp.o2.female[,phenotype:="resp.o2.female"]
weight.male <- weight[Sex == "Male",];weight.male[,Sex := NULL];weight.male[,phenotype:="weight.male"]
weight.female <- weight[Sex == "Female",];weight.female[,Sex := NULL];weight.female[,phenotype:="weight.female"]
rm(resp.co2, resp.o2, weight)
phenotype <- rbind(relfit, fecundity, resp.co2.female, resp.co2.male, resp.o2.female, resp.o2.male, weight.female, weight.male, development)
rm(relfit, fecundity, resp.co2.female, resp.co2.male, resp.o2.female, resp.o2.male, weight.female, weight.male, development)
# phenotype[phenotype == "development", effect.size := -1L*effect.size] ## so that the adaptive response, aka, shorter development is positive value now.
phenotype[contrast=="(Ancestral - Cold)", contrast:="Ancestral - Cold"][
  contrast=="(Ancestral - Hot)", contrast:="Ancestral - Hot"][
    contrast=="(Cold - Hot)", contrast:="Cold - Hot"]

phenotype <- phenotype[,.(phenotype, cge,contrast, effect.size)]

wide_phenotype <- dcast(phenotype, phenotype + cge ~ contrast,value.var = "effect.size")

wide_phenotype[,lab.effect := min(abs(`Ancestral - Cold`),abs(`Ancestral - Hot`)) * sign(`Ancestral - Hot`*`Ancestral - Cold`), by = .(phenotype, cge)]
wide_phenotype[,temp.effect := abs(`Cold - Hot`)]

plot.dt <- wide_phenotype[,.(phenotype, cge, lab.effect, temp.effect)]
plot.dt <- melt(plot.dt, id.vars = c("phenotype","cge"), variable.name = "type", value.name = "effect.size")
plot.dt[,pheno.cge := paste(phenotype, cge, sep = ".")]
plot.dt[, type:= factor( type, levels = c("temp.effect", "lab.effect"))]
plot.dt[phenotype == "relfit", phenotype := "relative.fitness"]
order_vec <- c("fecundity","relative.fitness","resp.co2.male","resp.co2.female",
               "resp.o2.male", "resp.o2.female", "weight.male", "weight.female", "development")
plot.dt[ , phenotype := factor(phenotype, levels = order_vec) ]

png("../Plot/Phenotype.effect.size.lab.temp.png", width = 8, height = 6, units = "in", res = 600)
ggplot(data = plot.dt, aes(x = phenotype, y = effect.size, fill = type))+
  geom_col(position = position_dodge(width = 0.6)) +
  facet_wrap(~cge)+
  labs(x = NULL,y="Effect size (in unit of SDs)",
       title = "Temp- vs. Lab-driven response",
       fill = "Type") +
  scale_fill_manual(values = c("lab.effect"  = "purple",
                               "temp.effect" = "gold"),
                    name   = "Type",
                    labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0))
dev.off()

for (p in order_vec) {
  print(paste0("dealing with ", p, "..."))
  temp.plot <- plot.dt[phenotype == p, ]
  plot <- ggplot(data = temp.plot, aes(x = cge, y = effect.size, fill = type))+
    geom_col(position = position_dodge(width = 0.6)) +
    labs(x = NULL,y="Effect size (in unit of SDs)", title = p,
         fill = "Type") +
    scale_fill_manual(values = c("lab.effect"  = "purple",
                                 "temp.effect" = "gold"),
                      name   = "",
                      labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
    theme_minimal(base_size = 16) +
    theme(panel.grid.minor.x = element_line(colour = "grey60"),
          axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.caption = element_text(size = 10, vjust = 5, hjust = 0))
  png(paste0("../Plot/",p,".effect.size.lab.temp.png"), width = 4, height = 6, units = "in", res = 600)
  print(plot)
  dev.off()
}


# 
# 
# 
# wide_phenotype[,lab.effect := (min(`Cold-evolved`, `Hot-evolved`) - Ancestral), by = .(phenotype, cge)]
# wide_phenotype[,temp.effect := (max(`Cold-evolved`, `Hot-evolved`) - Ancestral), by = .(phenotype, cge)]
# wide_phenotype[, temp.prop := 1 - (lab.effect/temp.effect)][temp.prop >1, temp.prop := 1]
# wide_phenotype[, lab.prop:= 1 - temp.prop]
# 
# plot.dt <- wide_phenotype[,.(phenotype, cge, lab.prop, temp.prop)]
# setorder(plot.dt, lab.prop)
# plot.dt <- melt(plot.dt, id.vars = c("phenotype","cge"), variable.name = "type", value.name = "proportion")
# plot.dt[,pheno.cge := paste(phenotype, cge, sep = ".")]
# plot.dt[, type:= factor( type, levels = c("temp.prop", "lab.prop"))]
# plot.dt[phenotype == "relfit", phenotype := "relative.fitness"]
# order_vec <- unique(plot.dt[type == "lab.prop"][order(-proportion), phenotype])
# plot.dt[ , phenotype := factor(phenotype, levels = order_vec) ]
# 
# 
# png("../Plot/Phenotype.proportion.lab.temp.png", width = 8, height = 6, units = "in", res = 600)
# ggplot(data = plot.dt, aes(x = phenotype, y = proportion, fill = type))+
#   geom_col() +
#   facet_wrap(~cge)+
#   scale_y_continuous(labels = percent_format(accuracy = 1)) +
#   labs(x = NULL,
#        title = "Temp- vs. Lab-driven response",
#        fill = "Type") +
#   scale_fill_manual(values = c("lab.prop"  = "purple",
#                                "temp.prop" = "gold"),
#                     name   = "Type",
#                     labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
#   theme_minimal(base_size = 16) +
#   theme(panel.grid.minor.x = element_line(colour = "grey60"),
#         axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
#         plot.caption = element_text(size = 10, vjust = 5, hjust = 0))
# dev.off()


################################################################################
################################################################################
########### Change the proportion calculation to FC proportion #################
################################################################################
################################################################################

wide_phenotype <- dcast(phenotype, phenotype + cge ~ Evo,value.var = "emmean")
wide_phenotype[,lab.effect := log(min(`Cold-evolved`, `Hot-evolved`) / Ancestral), by = .(phenotype, cge)]
wide_phenotype[,temp.effect := log(max(`Cold-evolved`, `Hot-evolved`) / Ancestral), by = .(phenotype, cge)]
wide_phenotype[, temp.prop := (temp.effect - lab.effect) / temp.effect][temp.prop >1, temp.prop := 1]
wide_phenotype[, lab.prop:= 1 - temp.prop]

plot.dt <- wide_phenotype[,.(phenotype, cge, lab.prop, temp.prop)]
setorder(plot.dt, lab.prop)
plot.dt <- melt(plot.dt, id.vars = c("phenotype","cge"), variable.name = "type", value.name = "proportion")
plot.dt[,pheno.cge := paste(phenotype, cge, sep = ".")]
plot.dt[, type:= factor( type, levels = c("temp.prop", "lab.prop"))]
plot.dt[phenotype == "relfit", phenotype := "relative.fitness"]
order_vec <- unique(plot.dt[type == "lab.prop"][order(-proportion), phenotype])
plot.dt[ , phenotype := factor(phenotype, levels = order_vec) ]


png("../Plot/Phenotype.proportion.lab.temp.FC.png", width = 8, height = 6, units = "in", res = 600)
ggplot(data = plot.dt, aes(x = phenotype, y = proportion, fill = type))+
  geom_col() +
  facet_wrap(~cge)+
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL,
       title = "Temp- vs. Lab-driven response",
       fill = "Type") +
  scale_fill_manual(values = c("lab.prop"  = "purple",
                               "temp.prop" = "gold"),
                    name   = "Type",
                    labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0))
dev.off()


#############################################################################
#############################################################################
############### Old script for plotting the proportion ########################
#############################################################################
#############################################################################
setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Script/")
rm(list = ls());gc()
library(data.table)
setDTthreads(percent = 80)
getDTthreads()
library(ggplot2)
development <- fread("./development.emm.csv"); development[,phenotype:="development"]
fecundity <- fread("./fec.emm.csv");fecundity[,Lab := NULL][,Temp:=NULL];fecundity[,phenotype:="fecundity"]
relfit <- fread("./relfit.emm.csv");relfit[,phenotype:="relfit"]
resp.co2 <- fread("./resp.CO2.emm.csv")
resp.o2 <- fread("./resp.O2.emm.csv")
weight <- fread("./resp.Weight.emm.csv")
resp.co2.male <- resp.co2[Sex == "Male",];resp.co2.male[,Sex := NULL];resp.co2.male[,phenotype:="resp.co2.male"]
resp.co2.female <- resp.co2[Sex == "Female",];resp.co2.female[,Sex := NULL];resp.co2.female[,phenotype:="resp.co2.female"]
resp.o2.male <- resp.o2[Sex == "Male",];resp.o2.male[,Sex := NULL];resp.o2.male[,phenotype:="resp.o2.male"]
resp.o2.female <- resp.o2[Sex == "Female",];resp.o2.female[,Sex := NULL];resp.o2.female[,phenotype:="resp.o2.female"]
weight.male <- weight[Sex == "Male",];weight.male[,Sex := NULL];weight.male[,phenotype:="weight.male"]
weight.female <- weight[Sex == "Female",];weight.female[,Sex := NULL];weight.female[,phenotype:="weight.female"]
rm(resp.co2, resp.o2, weight)
phenotype <- rbind(relfit, fecundity, resp.co2.female, resp.co2.male, resp.o2.female, resp.o2.male, weight.female, weight.male, development)
rm(relfit, fecundity, resp.co2.female, resp.co2.male, resp.o2.female, resp.o2.male, weight.female, weight.male, development)
phenotype[phenotype == "development", emmean := -1L*emmean] ## so that the adaptive response, aka, shorter development is positive value now.


phenotype <- phenotype[,.(phenotype, cge, Evo,emmean)]
wide_phenotype <- dcast(phenotype, phenotype + cge ~ Evo,value.var = "emmean")
wide_phenotype[,lab.effect := (min(`Cold-evolved`, `Hot-evolved`) - Ancestral), by = .(phenotype, cge)]
wide_phenotype[,temp.effect := (max(`Cold-evolved`, `Hot-evolved`) - Ancestral), by = .(phenotype, cge)]
wide_phenotype[, temp.prop := 1 - (lab.effect/temp.effect)][temp.prop >1, temp.prop := 1]
wide_phenotype[, lab.prop:= 1 - temp.prop]


plot.dt <- wide_phenotype[,.(phenotype, cge, lab.prop, temp.prop)]
setorder(plot.dt, lab.prop)
plot.dt <- melt(plot.dt, id.vars = c("phenotype","cge"), variable.name = "type", value.name = "proportion")
plot.dt[,pheno.cge := paste(phenotype, cge, sep = ".")]
plot.dt[, type:= factor( type, levels = c("temp.prop", "lab.prop"))]
plot.dt[phenotype == "relfit", phenotype := "relative.fitness"]
order_vec <- unique(plot.dt[type == "lab.prop"][order(-proportion), phenotype])
plot.dt[ , phenotype := factor(phenotype, levels = order_vec) ]




png("../Plot/Phenotype.proportion.lab.temp.png", width = 8, height = 6, units = "in", res = 600)
ggplot(data = plot.dt, aes(x = phenotype, y = proportion, fill = type))+
  geom_col() +
  facet_wrap(~cge)+
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL,
       title = "Temp- vs. Lab-driven response",
       fill = "Type") +
  scale_fill_manual(values = c("lab.prop"  = "purple",
                               "temp.prop" = "gold"),
                    name   = "Type",
                    labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0))
dev.off()




################################################################################
################################################################################
########### Change the proportion calculation to FC proportion #################
################################################################################
################################################################################


wide_phenotype <- dcast(phenotype, phenotype + cge ~ Evo,value.var = "emmean")
wide_phenotype[,lab.effect := log(min(`Cold-evolved`, `Hot-evolved`) / Ancestral), by = .(phenotype, cge)]
wide_phenotype[,temp.effect := log(max(`Cold-evolved`, `Hot-evolved`) / Ancestral), by = .(phenotype, cge)]
wide_phenotype[, temp.prop := (temp.effect - lab.effect) / temp.effect][temp.prop >1, temp.prop := 1]
wide_phenotype[, lab.prop:= 1 - temp.prop]


plot.dt <- wide_phenotype[,.(phenotype, cge, lab.prop, temp.prop)]
setorder(plot.dt, lab.prop)
plot.dt <- melt(plot.dt, id.vars = c("phenotype","cge"), variable.name = "type", value.name = "proportion")
plot.dt[,pheno.cge := paste(phenotype, cge, sep = ".")]
plot.dt[, type:= factor( type, levels = c("temp.prop", "lab.prop"))]
plot.dt[phenotype == "relfit", phenotype := "relative.fitness"]
order_vec <- unique(plot.dt[type == "lab.prop"][order(-proportion), phenotype])
plot.dt[ , phenotype := factor(phenotype, levels = order_vec) ]




png("../Plot/Phenotype.proportion.lab.temp.FC.png", width = 8, height = 6, units = "in", res = 600)
ggplot(data = plot.dt, aes(x = phenotype, y = proportion, fill = type))+
  geom_col() +
  facet_wrap(~cge)+
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL,
       title = "Temp- vs. Lab-driven response",
       fill = "Type") +
  scale_fill_manual(values = c("lab.prop"  = "purple",
                               "temp.prop" = "gold"),
                    name   = "Type",
                    labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0))
dev.off()


png("../Plot/Fitness.proportion.lab.temp.FC.png", width = 5, height = 6, units = "in", res = 600)
ggplot(data = plot.dt[phenotype %in% c("relative.fitness", "fecundity")], aes(x = phenotype, y = proportion, fill = type))+
  geom_col() +
  facet_wrap(~cge)+
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Fitness Proxies",
       fill = "Type") +
  scale_fill_manual(values = c("lab.prop"  = "purple",
                               "temp.prop" = "gold"),
                    name   = "Type",
                    labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),
        axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0))
dev.off()


png("../Plot/Pheno.noFitness.proportion.lab.temp.FC.png", width = 7, height = 6, units = "in", res = 600)
ggplot(data = plot.dt[!phenotype %in% c("relative.fitness", "fecundity")], aes(x = phenotype, y = proportion, fill = type))+
  geom_col() +
  facet_wrap(~cge)+
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = NULL,
       title = "Phenotypes, exlcuding fitness proxies",
       fill = "Type") +
  scale_fill_manual(values = c("lab.prop"  = "purple",
                               "temp.prop" = "gold"),
                    name   = "Type",
                    labels = c("lab.prop"  ="Lab-driven","temp.prop" =  "Temp-driven")) +
  theme_minimal(base_size = 16) +
  theme(panel.grid.minor.x = element_line(colour = "grey60"),legend.position = "none",
        axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.caption = element_text(size = 10, vjust = 5, hjust = 0))
dev.off()


