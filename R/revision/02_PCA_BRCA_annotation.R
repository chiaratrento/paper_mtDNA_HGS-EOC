## ============================================================================
## PCA annotated by BRCA1/2 status
##
## Addresses Reviewer 3 #6: check whether major oncogenic features (starting
## with BRCA1/2, the only one with data currently available - HRD and CCNE1
## are not in the current clinical database) explain transcriptional
## variance independently of mtDNA heteroplasmy group.
##
## Loads the dds/vsd objects saved by 01_DESeq2_stratified_heteroplasmy.R -
## does NOT re-run tximport/DESeq. Run 03 first if the RDS files below
## don't exist yet.
## ============================================================================

source("R/00_config.R")

library(tidyverse)
library(readxl)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)

out_dir <- file.path(results_dir, "revision/")

## ----------------------------------------------------------------------
## 1. Load dds/vsd from script 01 (no re-import, no re-running DESeq)
## ----------------------------------------------------------------------

dds <- readRDS(paste0(out_dir, "dds_heteroplasmy_group.rds"))
vsd <- readRDS(paste0(out_dir, "vsd_heteroplasmy_group.rds"))

## ----------------------------------------------------------------------
## 2. Load clinical database and extract/simplify BRCA status
## ----------------------------------------------------------------------

clinical_path <- clinical_brca

read_clinical_sheet <- function(sheet_name, prefix) {
  df <- read_excel(clinical_path, sheet = sheet_name)
  names(df)[1:2] <- c("codice_paziente", "sample_n")
  df <- df %>%
    filter(!is.na(sample_n)) %>%
    mutate(
      sample_id = paste0(prefix, sample_n),
      BRCA_raw = .data[["BRCA"]]
    )
  df
}

clin_beva <- read_clinical_sheet("CASI-BEVA", "B_")
clin_ctr  <- read_clinical_sheet("CONTROLLI NO-BEVA", "C_")

clinical <- bind_rows(
  clin_beva %>% dplyr::select(sample_id, BRCA_raw),
  clin_ctr  %>% dplyr::select(sample_id, BRCA_raw)
)

## Simplify free-text BRCA status (mixes germline "g"/somatic "s" testing,
## multiple lines per patient, WT/mutated/VUS/not-determined) into a single
## categorical variable per patient:
##   BRCA_mutated  - at least one line indicates a pathogenic mutation ("m")
##   BRCA_VUS      - no mutation found, but at least one VUS reported
##   BRCA_WT       - all reported results are wild-type
##   not_tested    - "nd" / missing
classify_brca <- function(raw) {
  if (is.na(raw) || trimws(raw) == "") return("not_tested")
  raw_low <- tolower(raw)
  if (raw_low == "nd") return("not_tested")
  lines <- unlist(strsplit(raw, "\n"))
  lines <- trimws(lines)
  has_mut <- any(grepl("\\bm\\b", lines, ignore.case = TRUE))
  has_vus <- any(grepl("vus", lines, ignore.case = TRUE))
  if (has_mut) return("BRCA_mutated")
  if (has_vus) return("BRCA_VUS")
  return("BRCA_WT")
}

clinical <- clinical %>%
  mutate(BRCA_status = map_chr(BRCA_raw, classify_brca),
         BRCA_status = factor(BRCA_status,
                              levels = c("BRCA_WT", "BRCA_VUS", "BRCA_mutated", "not_tested")))

table(clinical$BRCA_status)
## Sanity check: inspect raw -> classified mapping before trusting it
clinical %>% distinct(BRCA_raw, BRCA_status) %>% print(n = 30)

write.csv(clinical, paste0(out_dir, "BRCA_status_by_sample.csv"), row.names = FALSE)

## ----------------------------------------------------------------------
## 3. Join BRCA status to the sample metadata and plot PCA
## ----------------------------------------------------------------------

colData(vsd)$BRCA_status <- clinical$BRCA_status[
  match(colData(vsd)$`Sample Name`, clinical$sample_id)
]

if (any(is.na(colData(vsd)$BRCA_status))) {
  missing_ids <- colData(vsd)$`Sample Name`[is.na(colData(vsd)$BRCA_status)]
  message("Samples in RNA-seq data with NO BRCA status match: ", paste(missing_ids, collapse = ", "))
}

pcaData <- plotPCA(vsd, intgroup = c("heteroplasmy_group", "BRCA_status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

## (a) PCA coloured by BRCA status alone - the primary check requested by
## Reviewer 3: does BRCA status explain the transcriptional variance?
ggplot(pcaData, aes(PC1, PC2, colour = BRCA_status)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - vst, annotated by BRCA1/2 status")

ggsave(paste0(out_dir, "PCA_BRCA_status.png"), width = 7, height = 6, dpi = 600)

## (b) PCA with heteroplasmy group (colour) and BRCA status (shape) together -
## checks for confounding between the two variables in the same plot.
ggplot(pcaData, aes(PC1, PC2, colour = heteroplasmy_group, shape = BRCA_status)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - vst, heteroplasmy group x BRCA1/2 status")

ggsave(paste0(out_dir, "PCA_heteroplasmy_x_BRCA.png"), width = 8, height = 6, dpi = 600)

## ----------------------------------------------------------------------
## 4. Quick cross-tabulation: is BRCA status balanced across heteroplasmy
## groups, or is it a potential confounder that needs to be called out
## explicitly in the rebuttal regardless of what the PCA shows visually?
## ----------------------------------------------------------------------

crosstab <- table(colData(vsd)$heteroplasmy_group, colData(vsd)$BRCA_status)
print(crosstab)
write.csv(as.data.frame.matrix(crosstab), paste0(out_dir, "crosstab_heteroplasmy_BRCA.csv"))

## ----------------------------------------------------------------------
## 5. Fisher's exact test: is BRCA-mutated proportion imbalanced across
## heteroplasmy groups? (small-cell-count safe, unlike chi-squared here)
## ----------------------------------------------------------------------

## Collapse to BRCA_mutated vs everything-else (WT + VUS + not_tested),
## since the question that matters for confounding is specifically whether
## BRCA-mutated status tracks heteroplasmy group.
brca_mut_tab <- table(
  heteroplasmy_group = colData(vsd)$heteroplasmy_group,
  BRCA_mutated = ifelse(colData(vsd)$BRCA_status == "BRCA_mutated", "mutated", "other")
)
print(brca_mut_tab)

fisher_overall <- fisher.test(brca_mut_tab, simulate.p.value = TRUE, B = 100000)
cat("\nFisher's exact test (BRCA-mutated vs other, across all 4 heteroplasmy groups):\n")
print(fisher_overall)

## Pairwise HIGH vs WT specifically (the comparison flagged visually above:
## 31% BRCA-mutated in HIGH vs 9% in WT)
fisher_high_vs_wt <- fisher.test(brca_mut_tab[c("HIGH", "WT"), ])
cat("\nFisher's exact test (BRCA-mutated proportion, HIGH vs WT only):\n")
print(fisher_high_vs_wt)

fisher_summary <- data.frame(
  comparison = c("All 4 groups (BRCA-mutated vs other)", "HIGH vs WT (BRCA-mutated vs other)"),
  p_value    = c(fisher_overall$p.value, fisher_high_vs_wt$p.value)
)
print(fisher_summary)
write.csv(fisher_summary, paste0(out_dir, "fisher_BRCA_heteroplasmy_summary.csv"), row.names = FALSE)



## ============================================================================
## Notes:
## - Requires R/revision/01_DESeq2_stratified_heteroplasmy.R to have been run
##   first (dds_heteroplasmy_group.rds / vsd_heteroplasmy_group.rds).
## - classify_brca() collapses germline/somatic testing into one status.
## - HRD and CCNE1 status (also requested by Reviewer 3 #6) are not available
##   in the clinical database.
## ============================================================================
