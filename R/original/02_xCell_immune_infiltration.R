## ============================================================================
## 02_xCell_immune_infiltration.R
##
## Cell-type enrichment scores with xCell (Aran et al., Genome Biol 2017,
## https://github.com/dviraran/xCell) on the filtered count matrix, compared
## between WT and MUT (comp_3). Produces the objects used by
## R/figures/figure_panel2_xCell.R.
##
## Requires: R/original/01_WT_vs_MUT_DESeq2_enrichment.Rmd (saves `filtered`, `ids2`)
## ============================================================================

source("R/00_config.R")

library(xCell)          # devtools::install_github("dviraran/xCell")
library(SummarizedExperiment)
library(tidyverse)
library(ggpubr)
library(pheatmap)
library(Polychrome)

load(original_rdata)    # filtered, ids2, ...

## ---- 1. xCell scores ---------------------------------------------------------
exprMatrix <- assay(filtered, "counts")
rownames(exprMatrix) <- ids2
xcell <- xCellAnalysis(exprMatrix)
write.csv(xcell, file.path(results_dir, "xcell_scores.csv"))

## ---- 2. Immune / stromal cell types shown in the paper -----------------------
## xCell name -> label used in the figures
celltype_map <- c(
  "B-cells" = "B cell", "Memory B-cells" = "Memory B cell", "naive B-cells" = "Naive B cell",
  "Plasma cells" = "Plasma B cell", "Class-switched memory B-cells" = "Class-switched memory B cell",
  "CLP" = "Common lymphoid progenitor (CLP)", "CMP" = "Common myeloid progenitor",
  "Endothelial cells" = "Endothelial cell", "Eosinophils" = "Eosinophil",
  "GMP" = "Granulocyte-monocyte progenitor", "HSC" = "Hematopoietic stem cell",
  "Macrophages" = "Macrophage", "Macrophages M1" = "M1 macrophage", "Macrophages M2" = "M2 macrophage",
  "Mast cells" = "Mast cell", "Monocytes" = "Monocyte", "aDC" = "Activated myeloid dendritic cell",
  "iDC" = "Myeloid dendritic cell", "Neutrophils" = "Neutrophil", "NK cells" = "NK cell",
  "pDC" = "Plasmacytoid dendritic cell", "CD4+ memory T-cells" = "Memory CD4+ T cell",
  "CD4+ naive T-cells" = "Naive CD4+ T cell", "CD8+ naive T-cells" = "Naive CD8+ T cell",
  "CD8+ T-cells" = "CD8+ T cell", "CD4+ T-cells" = "Non-regulatory CD4+ T cell",
  "CD4+ Tcm" = "Central memory CD4+ T cell", "CD4+ Tem" = "Effector memory CD4+ T cell",
  "CD8+ Tcm" = "Central memory CD8+ T cell", "CD8+ Tem" = "Effector memory CD8+ T cell",
  "Tgd cells" = "Gamma delta T cell", "NKT" = "NK T cell", "Tregs" = "Regulatory T cell",
  "Th1 cells" = "CD4+ Th1", "Th2 cells" = "CD4+ Th2"
)

matched_cells  <- names(celltype_map)[names(celltype_map) %in% rownames(xcell)]
xcell_filtered <- xcell[matched_cells, ]
rownames(xcell_filtered) <- celltype_map[matched_cells]

## proportions per sample (scores rescaled to sum to 1) for the stacked bar plot
xcell_norm <- as.data.frame(t(xcell_filtered))
xcell_norm <- xcell_norm / rowSums(xcell_norm)
xcell_long <- xcell_norm %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "Proportion")

## sample metadata: sample name + comp_3 group
metainfo_df <- data.frame(Sample = filtered$`Sample Name`, Group = filtered$comp_3)
xcell_long_grouped <- merge(xcell_long, metainfo_df, by = "Sample")

ggplot(xcell_long_grouped, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = unname(palette36.colors(36))) +
  labs(title = "Immune cell infiltration in WT vs MUT", x = "Samples", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2))

## mean proportion per cell type and group
celltype_means_wide <- xcell_long_grouped %>%
  group_by(CellType, Group) %>%
  summarise(mean_prop = mean(Proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_prop, values_fill = 0, names_prefix = "Mean_")
write.csv(celltype_means_wide, file.path(results_dir, "xcell_mean_proportion_by_group.csv"), row.names = FALSE)

## ---- 3. WT vs MUT comparison on raw xCell scores -----------------------------
xcell_df <- as.data.frame(t(xcell[matched_cells, ]))
xcell_df$Sample <- rownames(xcell_df)
xcell_df <- left_join(xcell_df, data.frame(Sample = metainfo_df$Sample, comp_3 = metainfo_df$Group),
                      by = "Sample")
cell_types <- setdiff(colnames(xcell_df), c("Sample", "comp_3"))

## standardised mean difference (positive = higher in WT)
stats_wide <- xcell_df %>%
  pivot_longer(all_of(cell_types), names_to = "CellType", values_to = "score") %>%
  group_by(CellType, comp_3) %>%
  summarise(mean = mean(score), sd = sd(score), .groups = "drop") %>%
  pivot_wider(names_from = comp_3, values_from = c(mean, sd)) %>%
  mutate(Z_score = (mean_WT - mean_MUT) / sqrt((sd_WT^2 + sd_MUT^2) / 2))

ggplot(stats_wide, aes(x = reorder(CellType, Z_score), y = Z_score)) +
  geom_bar(stat = "identity") + coord_flip() +
  labs(title = "Standardised difference WT vs MUT", x = "Cell type", y = "Z-score") +
  theme_minimal()

## Wilcoxon test per cell type, FDR-adjusted
safe_wilcox <- function(cell) {
  values <- xcell_df[[cell]]
  group  <- xcell_df$comp_3
  vals   <- split(values, group)
  if (length(vals) < 2 || any(sapply(vals, function(x) all(is.na(x)) || length(unique(x)) <= 1)))
    return(data.frame(CellType = cell, p.value = NA))
  data.frame(CellType = cell, p.value = wilcox.test(values ~ group)$p.value)
}

pvalues <- lapply(cell_types, safe_wilcox) %>%
  bind_rows() %>%
  filter(!is.na(p.value)) %>%
  arrange(p.value) %>%
  mutate(FDR = p.adjust(p.value, method = "fdr"))
write.csv(pvalues, file.path(results_dir, "xcell_wilcoxon_WT_vs_MUT.csv"), row.names = FALSE)

sig_cells <- pvalues %>% filter(p.value < 0.05) %>% pull(CellType)
for (cell in sig_cells) {
  print(
    ggplot(xcell_df, aes(x = comp_3, y = .data[[cell]], fill = comp_3)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.2, size = 1) +
      stat_compare_means(method = "wilcox.test") +
      labs(title = cell, x = "", y = "xCell score") +
      theme_minimal()
  )
}

## ---- 4. Objects for the figure script ----------------------------------------
save(xcell, xcell_long_grouped, xcell_df, celltype_map,
     file = file.path(results_dir, "xcell_objects.RData"))
