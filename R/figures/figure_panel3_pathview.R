## ============================================================================
## Figure panel 3 - KEGG Pathview maps: Oxidative phosphorylation (hsa00190)
## and Citrate cycle (hsa00020)
##
## Requires the pathview PNGs written by R/pt1/01_WT_vs_MUT_DESeq2_enrichment.Rmd
## ============================================================================

source("R/00_config.R")

library(ggplot2)
library(ggpubr)
library(png)

img_oxphos <- readPNG(file.path(figures_dir, "hsa00190.pathview.png"))
img_tca    <- readPNG(file.path(figures_dir, "hsa00020.pathview.png"))

panel3 <- ggarrange(ggplot() + background_image(img_oxphos),
                    ggplot() + background_image(img_tca),
                    ncol = 1, nrow = 2, labels = c("A", "B"))

ggsave(file.path(figures_dir, "panel3.png"), panel3,
       width = 15, height = 20, units = "cm", dpi = 2000)
