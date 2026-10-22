# Preprocessing (shell)

Steps run before the R analysis. Tool versions marked `<version>` still need to be filled in.

1. **Quality control**: FastQC (`1.14`) on raw FASTQ files, summarised with MultiQC (`1.14`).
2. **Trimming**: cutadapt (`1.18`), adapter and quality trimming of 150 bp paired-end reads.
3. **Quantification**: Salmon (`1.4.0`) against the GENCODE GRCh38.p13 transcriptome, one output folder per sample (`quant.sf`).
4. **Gene-level counts**: `tximport` in R with a transcript-to-gene table built from the GENCODE GTF (`data/tx2gene.txt`, see `R/00_config.R`).
