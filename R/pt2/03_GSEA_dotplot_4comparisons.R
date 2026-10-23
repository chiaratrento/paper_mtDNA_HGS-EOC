# ------------------------------------------------------------------
# GSEA (KEGG) dot plot: LOW vs WT, MID vs WT, HIGH vs LOW, HIGH vs MID
#
# Two hues (direction of enrichment), shaded by adj. p - no NES gradient:
#   red  = NES > 0  -> enriched in the FIRST-named group (LOW, MID, HIGH, HIGH)
#   blue = NES < 0  -> enriched in the SECOND-named group (WT, WT, LOW, MID)
#   size = -log10(p.adjust), floored at 1e-4
#   colour intensity = -log10(p.adjust) (darker = more significant), red/blue = direction
#
# Input: the GSEA KEGG csv files written by section 8 of
# R/pt2/01_DESeq2_stratified_heteroplasmy.R (gseKEGG, ranked by DESeq2 `stat`, BH-adjusted).
# HIGH vs LOW / HIGH vs MID are obtained from LOW_vs_HIGH / MID_vs_HIGH
# by flipping the NES sign (pvalue and p.adjust unchanged).
# ------------------------------------------------------------------
source("R/00_config.R")
suppressPackageStartupMessages({
  library(ggplot2)
})

in_dir  <- file.path(results_dir, "pt2")   # GSEA_*.csv from 01_DESeq2_stratified_heteroplasmy.R
out_dir <- in_dir

# ---- 1. Comparisons, files and sign ------------------------------
comparisons <- data.frame(
  comparison = c("LOW vs WT", "MID vs WT", "HIGH vs LOW", "HIGH vs MID"),
  file       = c("GSEA_LOW_vs_WT.csv", "GSEA_MID_vs_WT.csv",
                 "GSEA_LOW_vs_HIGH.csv", "GSEA_MID_vs_HIGH.csv"),
  sign       = c(1, 1, -1, -1),
  x          = c(0, 1, 2.7, 3.7),   # gap between the WT-referenced and HIGH-referenced blocks
  stringsAsFactors = FALSE
)

# ---- 2. Pathways (KEGG IDs), theme and label ---------------------
# Selection made by Chiara (38 pathways, adj. p < 0.05 in >= 1 of the 4 contrasts,
# redundant disease-named KEGG sets dropped). Selected by ID, not by name.
mk <- function(theme, ids, labels) data.frame(theme = theme, ID = ids, label = labels,
                                              stringsAsFactors = FALSE)
paths <- rbind(
  mk("Mitochondrial\nmetabolism",
     c("hsa00190", "hsa00020", "hsa04714"),
     c("Oxidative phosphorylation", "Citrate cycle (TCA cycle)", "Thermogenesis")),
  mk("Proliferation &\nDNA repair",
     c("hsa04110", "hsa03030", "hsa03440", "hsa03460", "hsa03420", "hsa03410", "hsa03430"),
     c("Cell cycle", "DNA replication", "Homologous recombination", "Fanconi anemia pathway",
       "Nucleotide excision repair", "Base excision repair", "Mismatch repair")),
  mk("RNA & protein\nprocessing",
     c("hsa03010", "hsa03008", "hsa03040", "hsa03015", "hsa03013", "hsa03082", "hsa03050"),
     c("Ribosome", "Ribosome biogenesis", "Spliceosome", "mRNA surveillance pathway",
       "Nucleocytoplasmic transport", "Chromatin remodeling (ATP-dep.)", "Proteasome")),
  mk("Immune &\ninflammation",
     c("hsa04640", "hsa04060", "hsa04658", "hsa04659", "hsa04650", "hsa04662",
       "hsa04660", "hsa04064", "hsa05235", "hsa04062", "hsa04610"),
     c("Hematopoietic cell lineage", "Cytokine–cytokine receptor interaction",
       "Th1 and Th2 cell differentiation", "Th17 cell differentiation",
       "NK cell mediated cytotoxicity", "B cell receptor signaling pathway",
       "T cell receptor signaling pathway", "NF-kappa B signaling pathway",
       "PD-L1 / PD-1 checkpoint pathway", "Chemokine signaling pathway",
       "Complement and coagulation cascades")),
  mk("Adhesion &\nsignalling",
     c("hsa04510", "hsa04512", "hsa04518", "hsa04514", "hsa04151",
       "hsa04014", "hsa04010", "hsa04015", "hsa04020", "hsa05200"),
     c("Focal adhesion", "ECM-receptor interaction", "Integrin signaling",
       "Cell adhesion molecules (CAM)", "PI3K-Akt signaling pathway", "Ras signaling pathway",
       "MAPK signaling pathway", "Rap1 signaling pathway", "Calcium signaling pathway",
       "Pathways in cancer"))
)
stopifnot(!anyDuplicated(paths$ID), !anyDuplicated(paths$label))

# ---- 3. Read, flip sign, combine ---------------------------------
res <- do.call(rbind, lapply(seq_len(nrow(comparisons)), function(i) {
  d <- read.csv(file.path(in_dir, comparisons$file[i]), stringsAsFactors = FALSE)
  missing_ids <- setdiff(paths$ID, d$ID)
  if (length(missing_ids)) stop("IDs not found in ", comparisons$file[i], ": ",
                                paste(missing_ids, collapse = ", "))
  d <- d[d$ID %in% paths$ID, c("ID", "NES", "pvalue", "p.adjust")]
  d$NES        <- d$NES * comparisons$sign[i]
  d$comparison <- comparisons$comparison[i]
  d$x          <- comparisons$x[i]
  d
}))
res <- merge(res, paths, by = "ID")
stopifnot(nrow(res) == nrow(paths) * nrow(comparisons))

# ---- 4. Derived variables ----------------------------------------
lev_up   <- "NES > 0 (higher in first group)"
lev_down <- "NES < 0 (higher in second group)"
res$direction <- factor(ifelse(res$NES > 0, lev_up, lev_down), levels = c(lev_up, lev_down))

# size: -log10(adj. p), adj. p floored at 1e-4
res$logq <- pmin(-log10(pmax(res$p.adjust, 1e-4)), 4)
# colour: direction (sign of NES) x strength of evidence (-log10 adj. p)
res$signed_logq <- sign(res$NES) * res$logq

lev_sig <- "adj. p < 0.05"
lev_nom <- "nominal p < 0.05 (adj. p ≥ 0.05)"
lev_ns  <- "p ≥ 0.05"
res$tier <- factor(ifelse(res$p.adjust < 0.05, lev_sig,
                          ifelse(res$pvalue < 0.05, lev_nom, lev_ns)),
                   levels = c(lev_sig, lev_nom, lev_ns))

res$label <- factor(res$label, levels = rev(paths$label))   # first pathway of each theme on top
res$theme <- factor(res$theme, levels = unique(paths$theme))

col_up   <- "#C0392B"
col_down <- "#2166AC"

# ---- 5. Plot ------------------------------------------------------
p <- ggplot(res, aes(x = x, y = label)) +
  geom_vline(xintercept = 1.85, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_point(aes(fill = signed_logq, size = logq), shape = 21, colour = "grey30", stroke = 0.25) +
  scale_fill_gradientn(
    colours = c(col_down, "#D1E5F0", "#FDDBC7", col_up),
    values  = scales::rescale(c(-4, -0.01, 0.01, 4)),
    limits  = c(-4, 4),
    breaks  = c(-4, -3, -2, -1.3, 1.3, 2, 3, 4),
    labels  = c("\u22640.0001", "0.001", "0.01", "0.05", "0.05", "0.01", "0.001", "\u22640.0001"),
    name    = "Adj. p-value\n(red = up in first group,\nblue = up in second)") +
  scale_size_continuous(range = c(1, 7), limits = c(0, 4),
                        breaks = c(1.3, 2, 3, 4),
                        labels = c("0.05", "0.01", "0.001", "≤0.0001"),
                        name = "Adj. p-value") +
  scale_x_continuous(breaks = comparisons$x, labels = comparisons$comparison,
                     position = "top", limits = c(-0.6, 4.3), expand = c(0, 0)) +
  facet_grid(theme ~ ., scales = "free_y", space = "free_y", switch = "y") +
  guides(
    fill = guide_colourbar(order = 1, barheight = unit(150, "pt")),
    size = guide_legend(order = 2, override.aes = list(fill = "grey60"))) +
  labs(x = NULL, y = NULL,
       caption = paste0(
         "KEGG GSEA (gseKEGG on genes ranked by the DESeq2 Wald statistic; p-values BH-adjusted).\n",
         "Red: pathway enriched in the first-named group of each comparison; blue: enriched in the second.\n",
         "Groups by mtDNA heteroplasmy (VAF): WT (none detected), LOW (0–50%), MID (50–90%), HIGH (≥90%).\n",
         "Pathways: adj. p < 0.05 in at least one comparison; redundant disease-named sets omitted.")) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3),
    panel.background   = element_rect(fill = NA, colour = "grey88", linewidth = 0.3),
    panel.spacing.y    = unit(7, "pt"),
    strip.placement    = "outside",
    strip.background   = element_blank(),
    strip.text.y.left  = element_text(angle = 90, face = "bold", size = 8.5, margin = margin(r = 8)),
    axis.text.x.top    = element_text(angle = 45, hjust = 0, vjust = 0, size = 9, colour = "black"),
    axis.text.y        = element_text(size = 8.5, colour = "black"),
    axis.ticks         = element_blank(),
    plot.caption       = element_text(hjust = 0, size = 7, colour = "grey20"),
    plot.caption.position = "plot",
    legend.position    = "right",
    legend.key.height  = unit(14, "pt"))

# ---- 5b. Vertical bar next to each theme label --------------------
# Adds a grey line on the right edge of every left strip (theme label),
# spanning the rows of that theme. Uses only ggplot2 + grid.
g <- ggplotGrob(p)
for (i in grep("^strip-l", g$layout$name)) {
  g$grobs[[i]] <- grid::grobTree(
    g$grobs[[i]],
    grid::segmentsGrob(x0 = unit(1, "npc") - unit(2, "pt"), x1 = unit(1, "npc") - unit(2, "pt"),
                       y0 = unit(0.02, "npc"), y1 = unit(0.98, "npc"),
                       gp = grid::gpar(col = "grey35", lwd = 1.6, lineend = "butt")))
}
g$layout$clip[grep("^strip-l", g$layout$name)] <- "off"
p_final <- g   # ggsave accepts a gtable

# ---- 6. Save ------------------------------------------------------
ggsave(file.path(out_dir, "GSEA_dotplot_4comparisons_updown.png"), p_final,
       width = 8.5, height = 12.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "GSEA_dotplot_4comparisons_updown.pdf"), p_final,
       width = 8.5, height = 12.5, bg = "white")

# values used, to check against Chiara's dotplot_values_4comparisons.csv
write.csv(res[order(res$theme, match(res$label, levels(res$label)), res$x),
              c("theme", "ID", "label", "comparison", "NES", "pvalue", "p.adjust", "tier")],
          file.path(out_dir, "dotplot_values_from_R.csv"), row.names = FALSE)
