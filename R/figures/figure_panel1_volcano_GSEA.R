## ============================================================================
## Figure panel 1 - volcano plot (MUT vs WT) and KEGG GSEA bubble plot
##
## Requires: R/original/01_WT_vs_MUT_DESeq2_enrichment.Rmd (res, res_df, kk2)
##
## The DESeq2 contrast is WT_vs_MUT; here every sign is flipped so that the
## figure reads MUT vs WT (positive log2FC / "upregulated" = higher in MUT).
## ============================================================================

source("R/00_config.R")

library(ggplot2)
library(ggrepel)
library(dplyr)
library(forcats)
library(ggpubr)
library(org.Hs.eg.db)
library(clusterProfiler)

load(original_rdata)

## ---- A. Volcano plot ---------------------------------------------------------
res_df$log2FoldChange <- -res_df$log2FoldChange
res_df$diffexpressed  <- recode(res_df$diffexpressed, "Down" = "Up", "Up" = "Down")
table(res_df$diffexpressed)

## DEG table, MUT vs WT
geni_de <- na.omit(as.data.frame(res)[res$padj < 0.05, ])
geni_de$log2FoldChange <- -geni_de$log2FoldChange
geni_de$type   <- ifelse(geni_de$log2FoldChange > 0, "upregulated", "downregulated")
geni_de$symbol <- mapIds(org.Hs.eg.db, rownames(geni_de), "SYMBOL", "ENSEMBL")
write.csv(geni_de, file.path(results_dir, "DEG_MUT_vs_WT_padj005.csv"))

## genes labelled in the volcano plot
ensg_ids <- c(
  "ENSG00000077264", "ENSG00000077274", "ENSG00000103449", "ENSG00000111249",
  "ENSG00000112541", "ENSG00000128710", "ENSG00000145808", "ENSG00000148734",
  "ENSG00000164120", "ENSG00000166265", "ENSG00000173809", "ENSG00000182389",
  "ENSG00000186281", "ENSG00000203727", "ENSG00000231246", "ENSG00000255910",
  "ENSG00000257017", "ENSG00000259417", "ENSG00000268621", "ENSG00000273079",
  "ENSG00000276476", "ENSG00000091483", "ENSG00000187837", "ENSG00000101639"
)
volcano_plot_gene <- mapIds(org.Hs.eg.db, ensg_ids, "SYMBOL", "ENSEMBL")
res_df$label     <- volcano_plot_gene[rownames(res_df)]
res_df$label_ens <- ifelse(rownames(res_df) %in% ensg_ids, rownames(res_df), NA)
res_df$diffexpressed[res_df$label == "CEP192"] <- "Up"
res_df$label[res_df$label_ens == "ENSG00000255910"] <- "ENSG00000255910"  # no symbol

volcano <- function(label_col) {
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) +
    scale_color_manual(values = c("#4F9BF2", "#BEBEBE", "#BF011D")) +
    labs(col = "") + xlim(-5, 5) + ylim(0, 4) +
    geom_point(size = 2.5) +
    geom_text_repel(aes(label = .data[[label_col]]), color = "black",
                    box.padding = unit(0.5, "lines"), na.rm = TRUE) +
    theme_classic() +
    theme(axis.title = element_text(size = 15), axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 12), legend.text = element_text(size = 15))
}
b_symbol <- volcano("label")
b_ensg   <- volcano("label_ens")

## ---- B. KEGG GSEA bubble plot ------------------------------------------------
k <- kk2@result
k$type <- ifelse(k$NES < 0, "upregulated", "downregulated")   # flipped: MUT vs WT
k$Count     <- sapply(strsplit(k$core_enrichment, "/"), length)
k$GeneRatio <- k$Count / k$setSize
k$Description <- sapply(k$Description, function(d) paste(strwrap(d, 30), collapse = "\n"))

kegg_final <- read.csv(kegg_categories, sep = "")
k2 <- merge(k, kegg_final, by = "ID", all.x = TRUE)
k2 <- k2[k2$p.adjust < 0.35, ]
k2$type <- factor(k2$type, levels = c("upregulated", "downregulated"))

## leading-edge genes as symbols (for the supplementary table)
k2$core_enrichment_symbol <- sapply(k2$core_enrichment, function(g) {
  paste(bitr(strsplit(g, "/")[[1]], fromType = "ENTREZID", toType = "SYMBOL",
             OrgDb = org.Hs.eg.db)$SYMBOL, collapse = "/")
})
write.csv(k2, file.path(results_dir, "GSEA_KEGG_MUT_vs_WT_padj035.csv"), row.names = FALSE)

c_allpath_padj <- ggplot(k2, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 7.75)) +
  scale_colour_gradient(limits = c(0.001, 0.50), low = "#2600D1FF", high = "#EE3F3FFF") +
  ylab(NULL) + xlab("GeneRatio") +
  facet_grid(. ~ type)

c_all_path_type <- ggplot(k2, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = Main)) +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 6.75)) +
  ylab(NULL) + xlab("GeneRatio") +
  facet_grid(. ~ type)

## ---- Assemble and save -------------------------------------------------------
save_panel <- function(volc, gsea, file) {
  panel <- ggarrange(volc, gsea, ncol = 1, nrow = 2, labels = c("A", "B"), heights = c(1, 2))
  ggsave(file.path(figures_dir, file), panel, width = 35, height = 34, units = "cm", dpi = 1000)
}
save_panel(b_symbol, c_allpath_padj,  "panel1_symbol_padj.png")
save_panel(b_symbol, c_all_path_type, "panel1_symbol_type.png")
save_panel(b_ensg,   c_allpath_padj,  "panel1_ens_padj.png")
save_panel(b_ensg,   c_all_path_type, "panel1_ens_type.png")

## Legend (draft)
## A Volcano plot of MUT vs WT. Each dot is a gene: grey = not significant, red = up-regulated,
##   blue = down-regulated in MUT (padj < 0.05, |log2FC| > 0.6).
## B KEGG GSEA (genes ranked by the DESeq2 Wald statistic). Colour = adjusted p-value
##   (BH), size = number of leading-edge genes, GeneRatio = leading-edge genes / set size.
