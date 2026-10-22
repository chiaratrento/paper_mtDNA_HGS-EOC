# Sample metadata (description only)

No patient-level data are distributed. The scripts expect these files under `data/clinical/` (paths in `R/00_config.R`).

| File | Variables used | Used by |
|---|---|---|
| `Lista RNA pazienti OC.xlsx`, sheet `Foglio2` | `Sample Name` (B_1–B_30 = BEVA cohort, C_31–C_62 = CTR cohort), `comp_3` (WT / MUT: near-homoplasmic mtDNA variant vs other) | all RNA-seq scripts |
| `heteroplasmy_metadata.csv` | `sample_n`, `vaf` (tumour-area corrected variant allele frequency, %), `complex`, `cohort` (BEVA / CTR) | `R/revision/01_*` |
| `Data-Base_Indraccolo_con BRCA e p53_Dic2024_clean.xlsx`, sheets `CASI-BEVA`, `CONTROLLI NO-BEVA` | `BRCA` (free-text germline/somatic test results, classified as WT / VUS / mutated / not tested) | `R/revision/02_*` |

Heteroplasmy groups: WT = no variant detected; LOW = 0 < VAF < 50%; MID = 50% ≤ VAF < 90%; HIGH = VAF ≥ 90%.

Excluded samples (low sequencing quality): C_34, C_39.
