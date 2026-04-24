## This script is purely for generating Figrue 2
setwd("~/Dropbox/Manuscript/pleiotropy_yiwen/Analyses/Phenotypic/Script/")
rm(list = ls());gc()
library(data.table)
setDTthreads(percent = 80)
getDTthreads()
library(ggplot2)
library(cowplot)
fec.plot <- readRDS("./fec.plot.RDS")
fec.eff.plot <- readRDS("./fec.eff.plot.RDS")
rf.plot <- readRDS("./rf.plot.RDS")
rf.eff.plot <- readRDS("./rf.eff.plot.RDS")

## optional small theme tweaks for panel consistency
fec.plot <- fec.plot +
  theme(plot.title = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 13))

rf.plot <- rf.plot +
  theme(plot.title = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 13))

fec.eff.plot <- fec.eff.plot +
  theme(plot.title = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 12),
        legend.position = "right")

rf.eff.plot <- rf.eff.plot +
  theme(plot.title = element_text(size = 18, face = "plain"),
        axis.title = element_text(size = 16),
        axis.text  = element_text(size = 12),
        legend.position = "right")

left_col <- plot_grid(fec.plot,rf.plot,labels = c("A", "C"),
  label_size = 20,label_fontface = "bold",ncol = 1,align = "v",rel_heights = c(1, 1)
)

right_col <- plot_grid(fec.eff.plot, rf.eff.plot, labels = c("B", "D"),
  label_size = 20,label_fontface = "bold",ncol = 1,align = "v",rel_heights = c(1, 1)
)

final_fig <- plot_grid(left_col,right_col,
  ncol = 2,rel_widths = c(2.2, 1.15),align = "h")

ggsave(filename = "../Plot/Figure2.png",
  plot = final_fig,width = 14,height = 10,dpi = 600)
