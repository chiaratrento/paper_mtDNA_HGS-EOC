## ============================================================================
## Stratified RNA-seq re-analysis by mtDNA heteroplasmy level
##
## Re-runs the WT vs MUT DESeq2 comparison as a multi-level heteroplasmy
## design (WT / LOW 0-50% / MID 50-90% / HIGH >=90%) and, as a sensitivity
## analysis, with heteroplasmy percentage as a continuous covariate.
##
## Purpose: test whether the transcriptional / KEGG pathway signature
## currently attributed to near-homoplasmic (HIGH) variants is already
## present at intermediate (MID) heteroplasmy, and report enrichment
## results with explicit multiple-testing status.
## ============================================================================

library(clusterProfiler)   # load before tidyverse (rlang conflicts)
library(DESeq2)
library(tximport)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ReactomePA)
library(org.Hs.eg.db)
library(EDASeq)
library(edgeR)
library(readxl)

source("R/00_config.R")

## ----------------------------------------------------------------------
## 0. Heteroplasmy metadata
##
## Normalized (tumour-area corrected) VAF per sample for the 62-sample
## cohort (BEVA-BS 1-30, CTR-BS 31-62).
## ----------------------------------------------------------------------

heteroplasmy_meta <- read.csv(heteroplasmy_csv, stringsAsFactors = FALSE)
# expected columns: sample_n, vaf (numeric, %), complex, cohort (BEVA/CTR)

heteroplasmy_meta <- heteroplasmy_meta %>%
  mutate(
    heteroplasmy_group = case_when(
      vaf == 0            ~ "WT",
      vaf > 0  & vaf < 50  ~ "LOW",
      vaf >= 50 & vaf < 90 ~ "MID",
      vaf >= 90            ~ "HIGH"
    ),
    heteroplasmy_group = factor(heteroplasmy_group, levels = c("WT", "LOW", "MID", "HIGH")),
    heteroplasmy_percent = vaf,
    # Reconstruct the sample ID as used in the Salmon quant / colData
    # (B_1..B_30 for BEVA, C_31..C_62 for CTR), to allow a direct join
    # below without relying on numeric sample_n alone.
    sample_id = ifelse(cohort == "BEVA", paste0("B_", sample_n), paste0("C_", sample_n))
  )

## QC exclusion: C_34 removed for low sequencing quality (documented here
## for Methods transparency, per Reviewer 1 #4). Heteroplasmy for C_34 was
## 63.4% (MID group) - flag this explicitly in the rebuttal/Methods, since
## excluding it shrinks an already small MID group (n=8 -> n=7).
heteroplasmy_meta <- heteroplasmy_meta %>% filter(sample_id != "C_34")

table(heteroplasmy_meta$heteroplasmy_group)

## Output directory
out_dir <- file.path(results_dir, "pt2/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ----------------------------------------------------------------------
## 1. Import counts (tximport / Salmon)
## ----------------------------------------------------------------------

samples <- list.files(path = salmon_dir, full.names = TRUE, pattern = "B_")
files_B <- file.path(samples, "quant.sf")
names(files_B) <- gsub(".*\\/B_(\\d+)_.*", "B_\\1", files_B)

samples <- list.files(path = salmon_dir, full.names = TRUE, pattern = "C_")
files_C <- file.path(samples, "quant.sf")
files_C <- files_C[-c(1:24)]   # first 24 C_ folders belong to a different project
names(files_C) <- gsub(".*\\/C_(\\d+)_.*", "C_\\1", files_C)

files <- c(files_B, files_C)
tx2gene <- read.delim(tx2gene_path)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c(1, 3)])
txi_counts <- round(txi$counts)

## Merge sample-level metadata (must match column names / order used in colnames(txi_counts))
meta_data <- read_excel(sample_list, sheet = "Foglio2")
meta_data <- meta_data[-39, ]   # C_39 excluded (low sequencing quality)
meta_data <- as.data.frame(meta_data)
meta_data <- meta_data[order(meta_data$`Sample Name`), ]

## Safety check: ensure sample_id/heteroplasmy_percent exist even if this
## object was reloaded from a stale CSV / previous session.
if (!"sample_id" %in% colnames(heteroplasmy_meta)) {
  heteroplasmy_meta <- heteroplasmy_meta %>%
    mutate(
      heteroplasmy_percent = vaf,
      sample_id = ifelse(cohort == "BEVA", paste0("B_", sample_n), paste0("C_", sample_n))
    )
}

meta_data <- left_join(meta_data, heteroplasmy_meta, by = c("Sample Name" = "sample_id"))

## Flag any sample present in the clinical metadata but missing a VAF match
## (or vice versa) - a silent mismatch here would drop samples unnoticed.
if (any(is.na(meta_data$heteroplasmy_group))) {
  message("Samples in RNA-seq metadata with NO matching heteroplasmy data (will be dropped): ",
          paste(meta_data$`Sample Name`[is.na(meta_data$heteroplasmy_group)], collapse = ", "))
}
missing_from_rnaseq <- setdiff(heteroplasmy_meta$sample_id, meta_data$`Sample Name`)
if (length(missing_from_rnaseq) > 0) {
  message("Samples with heteroplasmy data but NO RNA-seq metadata entry: ",
          paste(missing_from_rnaseq, collapse = ", "))
}


meta_data <- meta_data %>% filter(!is.na(heteroplasmy_group))
txi_counts <- txi_counts[, meta_data$`Sample Name`]
txi$counts <- txi$counts[, meta_data$`Sample Name`]
txi$abundance <- txi$abundance[, meta_data$`Sample Name`]
txi$length <- txi$length[, meta_data$`Sample Name`]

stopifnot(all(colnames(txi_counts) == meta_data$`Sample Name`))

## ----------------------------------------------------------------------
## 2. Build DESeqDataSet - categorical multi-level design
## ----------------------------------------------------------------------

dds <- DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ heteroplasmy_group)

## Low-count filtering (documented threshold for Methods: mean count >= 10)
keep <- rowMeans(counts(dds)) >= 10
dds  <- dds[keep, ]

dds <- DESeq(dds)
resultsNames(dds)


## Save for reuse in downstream scripts (02_PCA_BRCA_annotation.R)
## without re-running tximport/DESeq from scratch.
saveRDS(dds, paste0(out_dir, "dds_heteroplasmy_group.rds"))

## ----------------------------------------------------------------------
## 3. QC: PCA / RLE annotated by heteroplasmy group (and cohort, if relevant)
##
## Addresses reviewer request for PCA/clustering to assess sample outliers,
## and to check whether cohort (BEVA vs CTR) structure explains variance
## independently of heteroplasmy status.
## ----------------------------------------------------------------------

vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, paste0(out_dir, "vsd_heteroplasmy_group.rds"))

pcaData <- plotPCA(vsd, intgroup = c("heteroplasmy_group", "cohort"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, colour = heteroplasmy_group, shape = cohort)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - vst, annotated by heteroplasmy group and cohort")

ggsave(paste0(out_dir, "PCA_heteroplasmy_group.png"), width = 7, height = 6, dpi = 600)

pal <- brewer.pal(4, "Set2")
plotRLE(counts(dds), col = pal[as.factor(dds$heteroplasmy_group)], outline = FALSE, las = 2)

## If BRCA1/2, HRD or CCNE1 status become available, add them to intgroup
## above and/or as additional colour/shape aesthetics to check confounding.

## ----------------------------------------------------------------------
## 4. Differential expression - pairwise contrasts
##
## DEG thresholds (consistent with original manuscript, now stated
## explicitly for the Methods): padj < 0.05 and |log2FC| > 0.6
## ----------------------------------------------------------------------

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 0.6

get_deg <- function(dds, contrast, variable = "heteroplasmy_group") {
  res <- results(dds, contrast = c(variable, contrast[1], contrast[2]), alpha = PADJ_CUTOFF)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    filter(!is.na(padj))
  res_df$symbol <- mapIds(org.Hs.eg.db, keys = sub("\\..*", "", res_df$gene_id),
                          column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  res_df$sig <- with(res_df, padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF)
  res_df
}

res_HIGH_vs_WT <- get_deg(dds, c("HIGH", "WT"))   # original comparison
res_MID_vs_WT  <- get_deg(dds, c("MID",  "WT"))   # new
res_LOW_vs_WT  <- get_deg(dds, c("LOW",  "WT"))   # new
res_MID_vs_HIGH <- get_deg(dds, c("MID", "HIGH")) # key contrast: is MID already like HIGH?
res_LOW_vs_HIGH <- get_deg(dds, c("LOW", "HIGH")) # completes all 6 pairwise contrasts
res_LOW_vs_MID  <- get_deg(dds, c("LOW", "MID"))  # completes all 6 pairwise contrasts

n_deg <- c(
  HIGH_vs_WT  = sum(res_HIGH_vs_WT$sig),
  MID_vs_WT   = sum(res_MID_vs_WT$sig),
  LOW_vs_WT   = sum(res_LOW_vs_WT$sig),
  MID_vs_HIGH = sum(res_MID_vs_HIGH$sig),
  LOW_vs_HIGH = sum(res_LOW_vs_HIGH$sig),
  LOW_vs_MID  = sum(res_LOW_vs_MID$sig)
)
print(n_deg)

write.csv(res_HIGH_vs_WT,  paste0(out_dir, "DEG_HIGH_vs_WT_symb.csv"),  row.names = FALSE)
write.csv(res_MID_vs_WT,   paste0(out_dir, "DEG_MID_vs_WT_symb.csv"),   row.names = FALSE)
write.csv(res_LOW_vs_WT,   paste0(out_dir, "DEG_LOW_vs_WT_symb.csv"),   row.names = FALSE)
write.csv(res_MID_vs_HIGH, paste0(out_dir, "DEG_MID_vs_HIGH_symb.csv"), row.names = FALSE)
write.csv(res_LOW_vs_HIGH, paste0(out_dir, "DEG_LOW_vs_HIGH_symb.csv"), row.names = FALSE)
write.csv(res_LOW_vs_MID,  paste0(out_dir, "DEG_LOW_vs_MID_symb.csv"),  row.names = FALSE)

## ----------------------------------------------------------------------
## 4b. Sensitivity groupings, added to directly match reviewer wording
##     and test robustness of the WT/LOW boundary choice.
## ----------------------------------------------------------------------

## (a) Reviewer 3's literal 3-tier scheme: undetectable / intermediate (10-90%) / high (>90%)
heteroplasmy_meta <- heteroplasmy_meta %>%
  mutate(
    reviewer3_group = case_when(
      vaf < 10             ~ "undetectable",
      vaf >= 10 & vaf < 90  ~ "intermediate",
      vaf >= 90             ~ "high"
    ),
    reviewer3_group = factor(reviewer3_group, levels = c("undetectable", "intermediate", "high")),
    ## (b) Sensitivity merge: WT + LOW collapsed into a single reference group,
    ## on the basis that VAF < ~40-50% is not expected to be functionally
    ## relevant (Mahmood et al. 2024)
    heteroplasmy_group_merged = case_when(
      heteroplasmy_group %in% c("WT", "LOW") ~ "WT_LOW",
      TRUE ~ as.character(heteroplasmy_group)
    ),
    heteroplasmy_group_merged = factor(heteroplasmy_group_merged, levels = c("WT_LOW", "MID", "HIGH"))
  )

## Re-attach to meta_data (already joined above) - only needed if this block
## runs after the join; otherwise these columns arrive automatically.
if (!"reviewer3_group" %in% colnames(meta_data)) {
  meta_data <- left_join(meta_data,
                         heteroplasmy_meta %>% dplyr::select(sample_id, reviewer3_group, heteroplasmy_group_merged),
                         by = c("Sample Name" = "sample_id"))
}

## Reviewer-3-literal contrast: high vs undetectable, high vs intermediate
dds_r3 <- DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ reviewer3_group)
dds_r3 <- dds_r3[keep, ]
dds_r3 <- DESeq(dds_r3)
res_r3_high_vs_undetectable <- get_deg(dds_r3, c("high", "undetectable"), variable = "reviewer3_group")
res_r3_high_vs_intermediate <- get_deg(dds_r3, c("high", "intermediate"), variable = "reviewer3_group")
cat("Reviewer 3 literal grouping - DEG counts:\n")
cat("high vs undetectable:", sum(res_r3_high_vs_undetectable$sig), "\n")
cat("high vs intermediate:", sum(res_r3_high_vs_intermediate$sig), "\n")

## WT+LOW merged sensitivity contrast
dds_merged <- DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ heteroplasmy_group_merged)
dds_merged <- dds_merged[keep, ]
dds_merged <- DESeq(dds_merged)
res_HIGH_vs_WTLOW <- get_deg(dds_merged, c("HIGH", "WT_LOW"), variable = "heteroplasmy_group_merged")
res_MID_vs_WTLOW  <- get_deg(dds_merged, c("MID",  "WT_LOW"), variable = "heteroplasmy_group_merged")
cat("WT+LOW merged reference - DEG counts:\n")
cat("HIGH vs WT_LOW:", sum(res_HIGH_vs_WTLOW$sig), "\n")
cat("MID vs WT_LOW:",  sum(res_MID_vs_WTLOW$sig),  "\n")

write.csv(res_r3_high_vs_undetectable, paste0(out_dir, "DEG_reviewer3_high_vs_undetectable.csv"), row.names = FALSE)
write.csv(res_r3_high_vs_intermediate, paste0(out_dir, "DEG_reviewer3_high_vs_intermediate.csv"), row.names = FALSE)
write.csv(res_HIGH_vs_WTLOW, paste0(out_dir, "DEG_HIGH_vs_WTLOW_merged.csv"), row.names = FALSE)
write.csv(res_MID_vs_WTLOW,  paste0(out_dir, "DEG_MID_vs_WTLOW_merged.csv"),  row.names = FALSE)

## ----------------------------------------------------------------------
## 5. Sensitivity analysis - heteroplasmy percentage as continuous covariate
##
## Avoids relying on an arbitrary categorical cutoff; tests a dose-response
## relationship between % heteroplasmy and gene expression via LRT.
## ----------------------------------------------------------------------

dds_cont <- DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ heteroplasmy_percent)
dds_cont <- dds_cont[keep, ]
dds_cont <- DESeq(dds_cont, test = "LRT", reduced = ~ 1)

res_cont <- results(dds_cont, alpha = PADJ_CUTOFF)
res_cont_df <- as.data.frame(res_cont) %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj))
res_cont_df$symbol <- mapIds(org.Hs.eg.db, keys = sub("\\..*", "", res_cont_df$gene_id),
                             column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
sum(res_cont_df$padj < PADJ_CUTOFF)

write.csv(res_cont_df, paste0(out_dir, "DEG_continuous_heteroplasmy_LRT.csv"), row.names = FALSE)

## ----------------------------------------------------------------------
## 6. Pathway enrichment per contrast (ORA + GSEA), with explicit
## multiple-testing reporting (addresses Reviewer 3 #7)
## ----------------------------------------------------------------------

run_kegg_ora <- function(res_df, universe_ids) {
  sig_ids <- mapIds(org.Hs.eg.db, res_df$gene_id[res_df$sig], "ENTREZID", "ENSEMBL")
  enrichKEGG(gene = na.omit(sig_ids), universe = universe_ids, organism = "hsa", pvalueCutoff = 1) # keep all, filter after
}

universe_ids <- mapIds(org.Hs.eg.db, rownames(dds), "ENTREZID", "ENSEMBL")

kegg_HIGH_vs_WT  <- run_kegg_ora(res_HIGH_vs_WT,  universe_ids)
kegg_MID_vs_WT   <- run_kegg_ora(res_MID_vs_WT,   universe_ids)
kegg_MID_vs_HIGH <- run_kegg_ora(res_MID_vs_HIGH, universe_ids)
kegg_LOW_vs_WT   <- run_kegg_ora(res_LOW_vs_WT,   universe_ids)
kegg_LOW_vs_HIGH <- run_kegg_ora(res_LOW_vs_HIGH, universe_ids)
kegg_LOW_vs_MID  <- run_kegg_ora(res_LOW_vs_MID,  universe_ids)

## Explicitly report how many pathways survive BH correction at 0.05,
## vs how many are only nominal (p < 0.05, padj >= 0.05) - this is the
## table/sentence Reviewer 3 #7 is asking for.
report_significance <- function(kegg_res, label) {
  df <- as.data.frame(kegg_res)
  data.frame(
    contrast              = label,
    n_nominal_p05          = sum(df$pvalue < 0.05, na.rm = TRUE),
    n_significant_padj05   = sum(df$p.adjust < 0.05, na.rm = TRUE),
    ## 0.15 = threshold actually used for Fig 2B in the manuscript
    n_significant_padj15   = sum(df$p.adjust < 0.15, na.rm = TRUE)
  )
}

kegg_significance_summary <- bind_rows(
  report_significance(kegg_HIGH_vs_WT,   "HIGH vs WT"),
  report_significance(kegg_MID_vs_WT,    "MID vs WT"),
  report_significance(kegg_MID_vs_HIGH,  "MID vs HIGH"),
  report_significance(kegg_LOW_vs_WT,    "LOW vs WT"),
  report_significance(kegg_LOW_vs_HIGH,  "LOW vs HIGH"),
  report_significance(kegg_LOW_vs_MID,   "LOW vs MID")
)
print(kegg_significance_summary)
write.csv(kegg_significance_summary, paste0(out_dir, "KEGG_significance_summary.csv"), row.names = FALSE)

## ----------------------------------------------------------------------
## 7. Pathway overlap: is the HIGH-group KEGG signature unique, or
## already present at MID? (directly addresses Reviewer 2, lines 245-247)
## ----------------------------------------------------------------------

pathways_HIGH <- as.data.frame(kegg_HIGH_vs_WT) %>% filter(p.adjust < 0.15) %>% pull(Description)
pathways_MID  <- as.data.frame(kegg_MID_vs_WT)  %>% filter(p.adjust < 0.15) %>% pull(Description)

shared_pathways <- intersect(pathways_HIGH, pathways_MID)
unique_to_HIGH  <- setdiff(pathways_HIGH, pathways_MID)

cat("Pathways shared between MID and HIGH:", length(shared_pathways), "\n")
cat("Pathways unique to HIGH (>=90%):", length(unique_to_HIGH), "\n")

# Optional: UpSet-style visualisation for the revised Figure 2B
library(UpSetR)
upset_input <- fromList(list(HIGH = pathways_HIGH, MID = pathways_MID))
upset(upset_input, order.by = "freq")

## ----------------------------------------------------------------------
## 8. GSEA (rank-based) per contrast - addresses Reviewer 3 #7 with more
## power than ORA on small DEG lists: uses the full ranked gene list
## (DESeq2 `stat` column) instead of a hard significance cutoff.
##
## BASELINE REFERENCE (important for reading the CSVs/collaborators):
## DESeq2 contrast c("heteroplasmy_group", "A", "B") tests A vs B, with
## B as the denominator/reference level. A positive NES / log2FC means
## "higher in A relative to B" - NOT necessarily relative to WT.
##   HIGH_vs_WT   : baseline = WT
##   MID_vs_WT    : baseline = WT
##   LOW_vs_WT    : baseline = WT
##   MID_vs_HIGH  : baseline = HIGH   <- NOT WT
##   LOW_vs_HIGH  : baseline = HIGH   <- NOT WT
##   LOW_vs_MID   : baseline = MID    <- NOT WT
## All 6 pairwise contrasts among the 4 groups (WT/LOW/MID/HIGH) are
## produced below.
## ----------------------------------------------------------------------

## Full (unfiltered) results per contrast, needed to build the ranked list
res_HIGH_vs_WT_full  <- results(dds, contrast = c("heteroplasmy_group", "HIGH", "WT"))
res_MID_vs_WT_full   <- results(dds, contrast = c("heteroplasmy_group", "MID",  "WT"))
res_MID_vs_HIGH_full <- results(dds, contrast = c("heteroplasmy_group", "MID",  "HIGH"))

build_ranked_list <- function(res_full) {
  res_df <- as.data.frame(res_full) %>%
    rownames_to_column("gene_id") %>%
    filter(!is.na(stat))
  ids <- mapIds(org.Hs.eg.db, res_df$gene_id, "ENTREZID", "ENSEMBL")
  ranked <- res_df$stat
  names(ranked) <- ids
  ranked <- ranked[!is.na(names(ranked))]
  ranked <- ranked[!duplicated(names(ranked))]
  sort(ranked, decreasing = TRUE)
}

rank_HIGH_vs_WT  <- build_ranked_list(res_HIGH_vs_WT_full)
rank_MID_vs_WT   <- build_ranked_list(res_MID_vs_WT_full)
rank_MID_vs_HIGH <- build_ranked_list(res_MID_vs_HIGH_full)

## Remaining pairwise contrasts, added for a complete 4-group picture
## (all 6 pairs among WT/LOW/MID/HIGH)
res_LOW_vs_WT_full   <- results(dds, contrast = c("heteroplasmy_group", "LOW", "WT"))
res_LOW_vs_HIGH_full <- results(dds, contrast = c("heteroplasmy_group", "LOW", "HIGH"))
res_LOW_vs_MID_full  <- results(dds, contrast = c("heteroplasmy_group", "LOW", "MID"))

rank_LOW_vs_WT   <- build_ranked_list(res_LOW_vs_WT_full)
rank_LOW_vs_HIGH <- build_ranked_list(res_LOW_vs_HIGH_full)
rank_LOW_vs_MID  <- build_ranked_list(res_LOW_vs_MID_full)

set.seed(1)
gsea_HIGH_vs_WT <- gseKEGG(geneList = rank_HIGH_vs_WT, organism = "hsa",
                           minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1, seed = TRUE)
gsea_MID_vs_WT <- gseKEGG(geneList = rank_MID_vs_WT, organism = "hsa",
                          minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1, seed = TRUE)
gsea_MID_vs_HIGH <- gseKEGG(geneList = rank_MID_vs_HIGH, organism = "hsa",
                            minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1, seed = TRUE)
gsea_LOW_vs_WT <- gseKEGG(geneList = rank_LOW_vs_WT, organism = "hsa",
                          minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1, seed = TRUE)
gsea_LOW_vs_HIGH <- gseKEGG(geneList = rank_LOW_vs_HIGH, organism = "hsa",
                            minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1, seed = TRUE)
gsea_LOW_vs_MID <- gseKEGG(geneList = rank_LOW_vs_MID, organism = "hsa",
                           minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1, seed = TRUE)

## Same explicit nominal-vs-corrected reporting as the ORA summary above
gsea_significance_summary <- bind_rows(
  report_significance(gsea_HIGH_vs_WT,   "GSEA HIGH vs WT"),
  report_significance(gsea_MID_vs_WT,    "GSEA MID vs WT"),
  report_significance(gsea_MID_vs_HIGH,  "GSEA MID vs HIGH"),
  report_significance(gsea_LOW_vs_WT,    "GSEA LOW vs WT"),
  report_significance(gsea_LOW_vs_HIGH,  "GSEA LOW vs HIGH"),
  report_significance(gsea_LOW_vs_MID,   "GSEA LOW vs MID")
)
print(gsea_significance_summary)
write.csv(gsea_significance_summary,
          paste0(out_dir, "GSEA_KEGG_significance_summary.csv"), row.names = FALSE)

write.csv(as.data.frame(gsea_HIGH_vs_WT),   paste0(out_dir, "GSEA_HIGH_vs_WT.csv"),   row.names = FALSE)
write.csv(as.data.frame(gsea_MID_vs_WT),    paste0(out_dir, "GSEA_MID_vs_WT.csv"),    row.names = FALSE)
write.csv(as.data.frame(gsea_MID_vs_HIGH),  paste0(out_dir, "GSEA_MID_vs_HIGH.csv"),  row.names = FALSE)
write.csv(as.data.frame(gsea_LOW_vs_WT),    paste0(out_dir, "GSEA_LOW_vs_WT.csv"),    row.names = FALSE)
write.csv(as.data.frame(gsea_LOW_vs_HIGH),  paste0(out_dir, "GSEA_LOW_vs_HIGH.csv"),  row.names = FALSE)
write.csv(as.data.frame(gsea_LOW_vs_MID),   paste0(out_dir, "GSEA_LOW_vs_MID.csv"),   row.names = FALSE)

## Quick dotplots per contrast, split by up/down (NES sign) - same style as
## the original comp_3.Rmd GSE KEGG section
for (obj_name in c("gsea_HIGH_vs_WT", "gsea_MID_vs_WT", "gsea_MID_vs_HIGH",
                   "gsea_LOW_vs_WT", "gsea_LOW_vs_HIGH", "gsea_LOW_vs_MID")) {
  obj <- get(obj_name)
  if (nrow(as.data.frame(obj)) > 0) {
    print(dotplot(obj, showCategory = 15, split = ".sign") + facet_grid(. ~ .sign) +
            ggtitle(obj_name))
  } else {
    message(obj_name, ": no enriched terms returned (empty result).")
  }
}

## ----------------------------------------------------------------------
## 8b. GSEA (rank-based) on Reactome pathways, for all 6 contrasts
##
## Same ranked gene lists as the KEGG GSEA above (section 8), same
## rationale, using ReactomePA::gsePathway - matches the Reactome GSEA
## block already present in the original comp_3.Rmd (`pw2 <- gsePathway(...)`).
## pvalueCutoff = 0.35 here keeps everything through to reporting, exactly
## as in the original script, so BH-corrected vs nominal can be compared
## explicitly afterwards (as for KEGG above) rather than pre-filtered away.
## ----------------------------------------------------------------------

set.seed(123)
run_gsea_reactome <- function(ranked_list) {
  gsePathway(geneList = ranked_list,
             organism = "human",
             exponent = 1,
             #nPermSimple = 10000,
             minGSSize = 10,
             maxGSSize = 1000,
             pvalueCutoff = 0.35,
             pAdjustMethod = "BH",
             verbose = TRUE)
}

reactome_HIGH_vs_WT  <- run_gsea_reactome(rank_HIGH_vs_WT)
reactome_MID_vs_WT   <- run_gsea_reactome(rank_MID_vs_WT)
reactome_MID_vs_HIGH <- run_gsea_reactome(rank_MID_vs_HIGH)
reactome_LOW_vs_WT   <- run_gsea_reactome(rank_LOW_vs_WT)
reactome_LOW_vs_HIGH <- run_gsea_reactome(rank_LOW_vs_HIGH)
reactome_LOW_vs_MID  <- run_gsea_reactome(rank_LOW_vs_MID)

## Same explicit nominal-vs-corrected reporting as the KEGG summaries above
reactome_significance_summary <- bind_rows(
  report_significance(reactome_HIGH_vs_WT,  "Reactome HIGH vs WT"),
  report_significance(reactome_MID_vs_WT,   "Reactome MID vs WT"),
  report_significance(reactome_MID_vs_HIGH, "Reactome MID vs HIGH"),
  report_significance(reactome_LOW_vs_WT,   "Reactome LOW vs WT"),
  report_significance(reactome_LOW_vs_HIGH, "Reactome LOW vs HIGH"),
  report_significance(reactome_LOW_vs_MID,  "Reactome LOW vs MID")
)
print(reactome_significance_summary)
write.csv(reactome_significance_summary,
          paste0(out_dir, "GSEA_Reactome_significance_summary.csv"), row.names = FALSE)

write.csv(as.data.frame(reactome_HIGH_vs_WT),  paste0(out_dir, "GSEA_Reactome_HIGH_vs_WT.csv"),  row.names = FALSE)
write.csv(as.data.frame(reactome_MID_vs_WT),   paste0(out_dir, "GSEA_Reactome_MID_vs_WT.csv"),   row.names = FALSE)
write.csv(as.data.frame(reactome_MID_vs_HIGH), paste0(out_dir, "GSEA_Reactome_MID_vs_HIGH.csv"), row.names = FALSE)
write.csv(as.data.frame(reactome_LOW_vs_WT),   paste0(out_dir, "GSEA_Reactome_LOW_vs_WT.csv"),   row.names = FALSE)
write.csv(as.data.frame(reactome_LOW_vs_HIGH), paste0(out_dir, "GSEA_Reactome_LOW_vs_HIGH.csv"), row.names = FALSE)
write.csv(as.data.frame(reactome_LOW_vs_MID),  paste0(out_dir, "GSEA_Reactome_LOW_vs_MID.csv"),  row.names = FALSE)

for (obj_name in c("reactome_HIGH_vs_WT", "reactome_MID_vs_WT", "reactome_MID_vs_HIGH",
                   "reactome_LOW_vs_WT", "reactome_LOW_vs_HIGH", "reactome_LOW_vs_MID")) {
  obj <- get(obj_name)
  if (nrow(as.data.frame(obj)) > 0) {
    print(dotplot(obj, showCategory = 15, split = ".sign") + facet_grid(. ~ .sign) +
            ggtitle(obj_name))
  } else {
    message(obj_name, ": no enriched Reactome terms returned (empty result).")
  }
}



## ============================================================================
## Notes:
## - Run from the repository root; all paths are set in R/00_config.R.
## - Outputs go to results/pt2/.
## - Comparisons involving the small MID group (n = 7) showed poorly calibrated
##   p-value distributions on diagnostic inspection: they are reported as
##   hypothesis-generating only (see README).
## ============================================================================
