## ============================================================================
## Figure panel 2 - xCell immune infiltration (WT vs MUT)
## Common lymphoid progenitor (CLP) and Hematopoietic stem cell are excluded;
## Naive B cell is always kept.
##
## Requires: R/pt1/02_xCell_immune_infiltration.R (xcell_objects.RData)
## ============================================================================

source("R/00_config.R")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

load(file.path(results_dir, "xcell_objects.RData"))  # xcell_long_grouped, xcell_df, celltype_map

exclude_cells <- c("Common lymphoid progenitor (CLP)", "Hematopoietic stem cell")

## ---- A. Stacked bar plot -----------------------------------------------------
threshold <- 0.015   # keep cell types with mean proportion >= 1.5%

xcell_long_grouped_mod <- xcell_long_grouped %>%
  filter(!CellType %in% exclude_cells) %>%
  group_by(CellType) %>%
  mutate(mean_prop = mean(Proportion)) %>%
  ungroup() %>%
  filter(mean_prop >= threshold | CellType == "Naive B cell") %>%
  group_by(Sample) %>%                         # rescale to 1 after removal
  mutate(Proportion = Proportion / sum(Proportion)) %>%
  ungroup()

## samples ordered by macrophage proportion within each group
sample_order <- xcell_long_grouped_mod %>%
  filter(CellType == "Macrophage") %>%
  arrange(Group, desc(Proportion)) %>%
  pull(Sample)
xcell_long_grouped_mod$Sample <- factor(xcell_long_grouped_mod$Sample, levels = sample_order)

other_celltypes <- setdiff(unique(as.character(xcell_long_grouped_mod$CellType)), "Macrophage")
xcell_long_grouped_mod$CellType <- factor(xcell_long_grouped_mod$CellType,
                                          levels = c(sort(other_celltypes), "Macrophage"))

p_stack <- ggplot(xcell_long_grouped_mod, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  labs(x = "Samples", y = "Proportion", fill = "Cell Type") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2))

## ---- B. Z-score heatmap of cell types with Wilcoxon p < 0.05 ------------------
inv_map      <- setNames(names(celltype_map), celltype_map)
exclude_orig <- na.omit(as.character(inv_map[exclude_cells]))
cell_types_filt <- setdiff(colnames(xcell_df), c("Sample", "comp_3", exclude_orig))

pvalues <- lapply(cell_types_filt, function(cell) {
  data.frame(CellType = cell, p.value = wilcox.test(xcell_df[[cell]] ~ xcell_df$comp_3)$p.value)
}) %>%
  bind_rows() %>%
  arrange(p.value) %>%
  mutate(FDR = p.adjust(p.value, method = "fdr"))
print(pvalues)   # note: selection below is on nominal p, FDR reported for transparency

sig_cells <- pvalues %>% filter(p.value < 0.05) %>% pull(CellType)

heat_long <- xcell_df[, c("Sample", "comp_3", sig_cells)] %>%
  pivot_longer(all_of(sig_cells), names_to = "CellType", values_to = "Score") %>%
  mutate(CellType = ifelse(CellType %in% names(celltype_map), celltype_map[CellType], CellType)) %>%
  group_by(CellType) %>%
  mutate(Zscore = as.numeric(scale(Score))) %>%
  ungroup()
heat_long$Sample   <- factor(heat_long$Sample, levels = xcell_df$Sample[order(xcell_df$comp_3)])
heat_long$CellType <- factor(heat_long$CellType,
                             levels = heat_long %>% group_by(CellType) %>%
                               summarise(m = mean(Zscore)) %>% arrange(desc(m)) %>% pull(CellType))

p_heat <- ggplot(heat_long, aes(x = Sample, y = CellType, fill = Zscore)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Z-score") +
  facet_grid(. ~ comp_3, scales = "free_x", space = "free_x") +
  labs(x = "Samples", y = "Immune Cell Type") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 9), panel.grid = element_blank())

## ---- Assemble and save -------------------------------------------------------
panel2 <- ggarrange(p_stack, p_heat, ncol = 1, nrow = 2, labels = c("A", "B"))
ggsave(file.path(figures_dir, "panel2_noCLP_noHSC.png"), panel2,
       width = 30, height = 25, units = "cm", dpi = 1000)
