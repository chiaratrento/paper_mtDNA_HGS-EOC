# Homoplasmic mitochondrial DNA variants associate with favorable prognosis in high-grade serous ovarian cancer

This repository contains the analysis code supporting Trento et al., *"Homoplasmic mitochondrial DNA variants associate with favorable prognosis in high-grade serous ovarian cancer"*.


## Repository structure

```
.
├── R/
│   ├── 00_config.R                                   # all input/output paths (edit here)
│   ├── 01_preprocessing_notes.md                     # FastQC / MultiQC / cutadapt / Salmon (shell, not R)
│   ├── pt1/                                     # analyses of the submitted manuscript
│   │   ├── 01_WT_vs_MUT_DESeq2_enrichment.Rmd        # tximport, QC, DESeq2 (~ comp_3), ORA, GSEA, Pathview
│   │   └── 02_xCell_immune_infiltration.R            # xCell cell-type scores, WT vs MUT
│   └── pt2/                                     # analyses added in response to the reviewers
│       ├── 01_DESeq2_stratified_heteroplasmy.R       # WT / LOW / MID / HIGH design, sensitivity groupings,
│       │                                             #   continuous VAF (LRT), ORA + GSEA (KEGG, Reactome)
│       ├── 02_PCA_BRCA_annotation.R                  # PCA by BRCA1/2 status, Fisher's exact test
│       └── 03_GSEA_dotplot_4comparisons.R            # KEGG GSEA dot plot across heteroplasmy groups
├── annotation/                                       # KEGG / Reactome pathway categories used in figures
├── metadata/
│   └── sample_metadata_description.md                # variables used (no patient-identifiable data)
├── results/                                          # created by the scripts (not tracked)
├── sessionInfo.txt
└── README.md
```

## How to run

All scripts are run from the repository root. Paths to raw data and clinical
metadata are set once in `R/00_config.R`; these files are not distributed
(see *Data availability*) and are expected under `data/` (git-ignored).

| Order | Script | Main outputs |
|---|---|---|
| 1 | `R/pt1/01_WT_vs_MUT_DESeq2_enrichment.Rmd` | DEG and GSEA tables, Pathview maps, `WT_vs_MUT_objects.RData` |
| 2 | `R/pt1/02_xCell_immune_infiltration.R` | xCell scores, Wilcoxon tests, `xcell_objects.RData` |
| 3 | `R/figures/figure_panel*.R` | figure panels |
| 4 | `R/pt2/01_DESeq2_stratified_heteroplasmy.R` | DEG, ORA/GSEA tables, significance summaries, RDS objects |
| 5 | `R/pt2/02_PCA_BRCA_annotation.R` | PCA by BRCA status, Fisher's test |
| 6 | `R/pt2/03_GSEA_dotplot_4comparisons.R` | GSEA dot plot (LOW vs WT, MID vs WT, HIGH vs LOW, HIGH vs MID) |

## Analysis pipeline

1. **Sequencing & preprocessing**: 150 bp paired-end reads, Illumina platform. Quality control with FastQC and MultiQC; adapter/quality trimming with cutadapt.
2. **Quantification**: Salmon against the GENCODE GRCh38.p13 human transcriptome; gene-level counts with `tximport`.
3. **Sample QC**: samples C_34 and C_39 excluded for low sequencing quality. Library sizes, RLE plots, PCA (rlog/vst) and comparison of normalisation methods (upper-quartile, full-quantile, TMM, RLE).
4. **Differential expression**: `DESeq2`, genes with mean count ≥ 10. DEGs: `padj < 0.05` (Benjamini–Hochberg) and `|log2FC| > 0.6`.
   - *pt1*: near-homoplasmic variants vs all other samples (`~ comp_3`).
   - *pt2*: four heteroplasmy groups by tumour-corrected VAF, WT (no variant), LOW (0–50%), MID (50–90%), HIGH (≥ 90%), all six pairwise contrasts; sensitivity analyses with the three-tier grouping undetectable (< 10%) / intermediate (10–90%) / high (≥ 90%), with WT and LOW merged, and with VAF as a continuous covariate (likelihood-ratio test).
5. **Pathway enrichment**: over-representation analysis (`enrichGO`, `enrichKEGG`, `enrichPathway`) and GSEA on genes ranked by the DESeq2 Wald statistic (`gseGO`, `gseKEGG`, `gsePathway`). All p-values are BH-adjusted; the number of pathways significant at nominal p < 0.05, padj < 0.05 and padj < 0.15 is reported for every contrast. Leading-edge overlap between pathways is assessed with the Jaccard index.
6. **Immune infiltration**: xCell scores compared between groups with the Wilcoxon rank-sum test (FDR-adjusted).


## Software requirements

R (≥ 4.5) with: `DESeq2`, `tximport`, `edgeR`, `EDASeq`, `clusterProfiler`, `ReactomePA`, `org.Hs.eg.db`, `pathview`, `xCell`, `tidyverse`, `readxl`, `ggplot2`, `ggrepel`, `ggpubr`, `pheatmap`, `UpSetR`, `Polychrome`, `png`. Full versions in [`sessionInfo.txt`](sessionInfo.txt).

## Data availability

- **Raw RNA-seq data**: deposited in the NCBI Sequence Read Archive (SRA) under BioProject accession number `[PRJNA1401716]`, publicly accessible upon publication.
- **Targeted mtDNA sequencing data**: not currently required for public deposition per journal policy; available from the corresponding author upon reasonable request.
- **Processed count matrices / sample-level clinical metadata**: not included in this repository due to patient privacy restrictions; available from the corresponding author upon reasonable request, subject to institutional data sharing agreements.

## Citation

Trento et al. Homoplasmic mitochondrial DNA variants associate with favorable prognosis in high-grade serous ovarian cancer.

## Contact

Corresponding author: Stefano Indraccolo — stefano.indraccolo@unipd.it

## Session Info

<details>
<summary>
Click here for Session Info
</summary>

``` r

sessionInfo()
# sessionInfo()
# R version 4.6.1 (2026-06-24)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 24.04.4 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Rome
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils    
# [6] datasets  methods   base     
# 
# other attached packages:
#   [1] UpSetR_1.4.1               
# [2] readxl_1.5.0               
# [3] edgeR_4.10.3               
# [4] limma_3.68.5               
# [5] EDASeq_2.46.0              
# [6] ShortRead_1.70.0           
# [7] GenomicAlignments_1.48.0   
# [8] Rsamtools_2.28.0           
# [9] Biostrings_2.80.1          
# [10] XVector_0.52.0             
# [11] BiocParallel_1.46.0        
# [12] org.Hs.eg.db_3.23.1        
# [13] AnnotationDbi_1.74.0       
# [14] ReactomePA_1.56.0          
# [15] clusterProfiler_4.20.0     
# [16] RColorBrewer_1.1-3         
# [17] pheatmap_1.0.13            
# [18] lubridate_1.9.5            
# [19] forcats_1.0.1              
# [20] stringr_1.6.0              
# [21] dplyr_1.2.1                
# [22] purrr_1.2.2                
# [23] readr_2.2.0                
# [24] tidyr_1.3.2                
# [25] tibble_3.3.1               
# [26] ggplot2_4.0.3              
# [27] tidyverse_2.0.0            
# [28] tximport_1.40.0            
# [29] DESeq2_1.52.0              
# [30] SummarizedExperiment_1.42.0
# [31] Biobase_2.72.0             
# [32] MatrixGenerics_1.24.0      
# [33] matrixStats_1.5.0          
# [34] GenomicRanges_1.64.0       
# [35] Seqinfo_1.2.0              
# [36] IRanges_2.46.0             
# [37] S4Vectors_0.50.1           
# [38] BiocGenerics_0.58.1        
# [39] generics_0.1.4             
# 
# loaded via a namespace (and not attached):
#   [1] R.methodsS3_1.8.2          
# [2] dichromat_2.0-1            
# [3] vroom_1.7.1                
# [4] progress_1.2.3             
# [5] goftest_1.2-3              
# [6] vctrs_0.7.3                
# [7] ggtangle_0.1.2             
# [8] spatstat.random_3.5-1      
# [9] digest_0.6.39              
# [10] png_0.1-9                  
# [11] ggrepel_0.9.8              
# [12] deldir_2.0-4               
# [13] parallelly_1.48.0          
# [14] MASS_7.3-66                
# [15] fontLiberation_0.1.0       
# [16] reshape2_1.4.5             
# [17] httpuv_1.6.17              
# [18] qvalue_2.44.0              
# [19] withr_3.0.3                
# [20] xfun_0.60                  
# [21] ggfun_0.2.1                
# [22] survival_3.8-11            
# [23] memoise_2.0.1              
# [24] gson_0.2.1                 
# [25] systemfonts_1.3.2          
# [26] ragg_1.5.2                 
# [27] tidytree_0.4.8             
# [28] zoo_1.9-0                  
# [29] pbapply_1.7-4              
# [30] R.oo_1.27.1                
# [31] prettyunits_1.2.0          
# [32] KEGGREST_1.52.2            
# [33] promises_1.5.0             
# [34] otel_0.2.0                 
# [35] httr_1.4.8                 
# [36] restfulr_0.0.17            
# [37] globals_0.19.1             
# [38] fitdistrplus_1.2-6         
# [39] ps_1.9.3                   
# [40] rstudioapi_0.19.0          
# [41] miniUI_0.1.2               
# [42] DOSE_4.6.0                 
# [43] reactome.db_1.96.0         
# [44] processx_3.9.0             
# [45] curl_7.1.0                 
# [46] ScaledMatrix_1.20.0        
# [47] ggraph_2.2.2               
# [48] polyclip_1.10-7            
# [49] SparseArray_1.12.2         
# [50] xtable_1.8-8               
# [51] evaluate_1.0.5             
# [52] S4Arrays_1.12.0            
# [53] BiocFileCache_3.2.0        
# [54] hms_1.1.4                  
# [55] irlba_2.3.7                
# [56] filelock_1.0.3             
# [57] ROCR_1.0-12                
# [58] reticulate_1.46.0          
# [59] spatstat.data_3.1-9        
# [60] magrittr_2.0.5             
# [61] lmtest_0.9-40              
# [62] later_1.4.8                
# [63] viridis_0.6.5              
# [64] ggtree_4.2.0               
# [65] lattice_0.23-1             
# [66] spatstat.geom_3.8-2        
# [67] future.apply_1.20.2        
# [68] XML_3.99-0.24              
# [69] scattermore_1.2            
# [70] scuttle_1.22.0             
# [71] cowplot_1.2.0              
# [72] RcppAnnoy_0.0.23           
# [73] pillar_1.11.1              
# [74] nlme_3.1-170               
# [75] pwalign_1.8.0              
# [76] compiler_4.6.1             
# [77] beachmat_2.28.0            
# [78] RSpectra_0.16-2            
# [79] stringi_1.8.9              
# [80] tensor_1.5.1               
# [81] plyr_1.8.9                 
# [82] BiocIO_1.22.0              
# [83] crayon_1.5.3               
# [84] abind_1.4-8                
# [85] gridGraphics_0.5-1         
# [86] locfit_1.5-9.12            
# [87] sp_2.2-3                   
# [88] graphlayouts_1.2.5         
# [89] bit_4.6.0                  
# [90] textshaping_1.0.5          
# [91] codetools_0.2-20           
# [92] BiocSingular_1.28.0        
# [93] plotly_4.12.1              
# [94] mime_0.13                  
# [95] splines_4.6.1              
# [96] Rcpp_1.1.2                 
# [97] fastDummies_1.7.6          
# [98] tidydr_0.0.6               
# [99] dbplyr_2.6.0               
# [100] scrapper_1.6.3             
# [101] cellranger_1.1.0           
# [102] interp_1.1-6               
# [103] knitr_1.51                 
# [104] blob_1.3.0                 
# [105] fs_2.1.0                   
# [106] listenv_1.0.0              
# [107] ggplotify_0.1.3            
# [108] Matrix_1.7-6               
# [109] callr_3.8.0                
# [110] statmod_1.5.2              
# [111] tzdb_0.5.0                 
# [112] tweenr_2.0.3               
# [113] pkgconfig_2.0.3            
# [114] tools_4.6.1                
# [115] cachem_1.1.0               
# [116] cigarillo_1.2.1            
# [117] RSQLite_3.53.3             
# [118] viridisLite_0.4.3          
# [119] DBI_1.3.0                  
# [120] graphite_1.58.0            
# [121] fastmap_1.2.0              
# [122] rmarkdown_2.31             
# [123] scales_1.4.0               
# [124] grid_4.6.1                 
# [125] ica_1.0-3                  
# [126] Seurat_5.5.1               
# [127] patchwork_1.3.2            
# [128] dotCall64_1.2              
# [129] graph_1.90.0               
# [130] RANN_2.6.2                 
# [131] farver_2.1.2               
# [132] tidygraph_1.3.1            
# [133] scatterpie_0.2.6           
# [134] yaml_2.3.12                
# [135] latticeExtra_0.6-31        
# [136] rtracklayer_1.72.0         
# [137] cli_3.6.6                  
# [138] lifecycle_1.0.5            
# [139] uwot_0.2.4                 
# [140] bluster_1.22.0             
# [141] rjson_0.2.23               
# [142] timechange_0.4.0           
# [143] gtable_0.3.6               
# [144] ggridges_0.5.7             
# [145] progressr_1.0.0            
# [146] parallel_4.6.1             
# [147] ape_5.8-1                  
# [148] enrichit_0.2.1             
# [149] jsonlite_2.0.0             
# [150] RcppHNSW_0.7.0             
# [151] bitops_1.1-0               
# [152] bit64_4.8.4                
# [153] Rtsne_0.17                 
# [154] yulab.utils_0.2.4          
# [155] spatstat.utils_3.2-4       
# [156] BiocNeighbors_2.6.0        
# [157] SeuratObject_5.4.0         
# [158] aisdk_1.4.12               
# [159] metapod_1.20.0             
# [160] GOSemSim_2.38.3            
# [161] dqrng_0.4.1                
# [162] spatstat.univar_3.2-0      
# [163] R.utils_2.13.0             
# [164] lazyeval_0.2.3             
# [165] shiny_1.14.0               
# [166] htmltools_0.5.9            
# [167] enrichplot_1.32.0          
# [168] GO.db_3.23.1               
# [169] sctransform_0.4.3          
# [170] rappdirs_0.3.4             
# [171] glue_1.8.1                 
# [172] spam_2.11-4                
# [173] httr2_1.3.0                
# [174] RCurl_1.98-1.20            
# [175] gdtools_0.5.1              
# [176] treeio_1.36.1              
# [177] scran_1.40.0               
# [178] jpeg_0.1-11                
# [179] gridExtra_2.3.1            
# [180] igraph_2.3.3               
# [181] R6_2.6.1                   
# [182] SingleCellExperiment_1.34.0
# [183] ggiraph_0.9.6              
# [184] labeling_0.4.3             
# [185] GenomicFeatures_1.64.0     
# [186] cluster_2.1.8.3            
# [187] aplot_0.3.1                
# [188] DelayedArray_0.38.2        
# [189] tidyselect_1.2.1           
# [190] ggforce_0.5.0              
# [191] fontBitstreamVera_0.1.1    
# [192] future_1.75.0              
# [193] rsvd_1.0.5                 
# [194] KernSmooth_2.23-27         
# [195] S7_0.2.2                   
# [196] fontquiver_0.2.1           
# [197] data.table_1.18.4          
# [198] htmlwidgets_1.6.4          
# [199] aroma.light_3.42.0         
# [200] hwriter_1.3.2.1            
# [201] biomaRt_2.68.0             
# [202] rlang_1.3.0                
# [203] spatstat.sparse_3.2-0      
# [204] spatstat.explore_3.8-2     
# [205] ggnewscale_0.5.2
```

</details>
