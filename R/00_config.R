## ============================================================================
## 00_config.R - paths shared by all scripts
##
## Every script in R/ starts with source("R/00_config.R") and is meant to be
## run from the repository root. Raw data and clinical metadata are NOT
## distributed with this repository (patient privacy, see README): place them
## under data/ following the layout below, or edit the paths here.
## ============================================================================

## ---- Inputs (not in the repository) ----------------------------------------
salmon_dir      <- "data/salmon"              # one folder per sample, each with quant.sf
tx2gene_path    <- "data/tx2gene.txt"         # GENCODE GRCh38.p13 transcript -> gene table
sample_list     <- "data/clinical/Lista RNA pazienti OC.xlsx"         # sheet "Foglio2"
heteroplasmy_csv <- "data/clinical/heteroplasmy_metadata.csv"         # sample_n, vaf, complex, cohort
clinical_brca   <- "data/clinical/Data-Base_Indraccolo_con BRCA e p53_Dic2024_clean.xlsx"

## ---- Annotation tables (pathway -> category, used only for figures) --------
kegg_categories     <- "annotation/kegg_final.txt"      # columns: ID, Main, Second, ...
reactome_categories <- "annotation/reactome_final.txt"  # columns: ID, Group, ...

## ---- Outputs ---------------------------------------------------------------
results_dir  <- "results"                         # tables / RDS objects
figures_dir  <- file.path(results_dir, "figures") # figures
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

## RData produced by R/original/01_WT_vs_MUT_DESeq2_enrichment.Rmd and reused
## by the xCell and figure scripts (replaces the old comp_3.RData workspace).
original_rdata <- file.path(results_dir, "WT_vs_MUT_objects.RData")
