# Output Reference

This page describes both tracked repository outputs and larger files generated
when the workflow is rerun.

## Tracked Outputs

| Path | Description |
|---|---|
| `figures/pca.png` | PCA of variance-stabilized counts by condition |
| `figures/volcano_plot.png` | Log2 fold change versus adjusted significance |
| `figures/heatmap.png` | Top 15 genes ranked by adjusted p-value |
| `figures/ma_plot.png` | Mean normalized abundance versus log2 fold change |
| `results/differential_expression/DESeq2_results.csv` | Complete DESeq2 result table |
| `results/differential_expression/SaltStress_Annotated_DEGs.csv` | Significant DEGs with annotation |
| `data/metadata/sample_metadata.csv` | Sample-to-condition and accession mapping |

## Pipeline-Generated Directories

| Directory | Expected contents |
|---|---|
| `raw/` | SRA downloads and converted `.fastq` files |
| `trimmed/` | Single-end reads after `TRAILING:10` trimming |
| `fastqc_raw/` | Raw FastQC reports and a MultiQC summary |
| `fastqc_trimmed/` | Post-trimming FastQC reports and a MultiQC summary |
| `reference/` | Genome FASTA, annotation, and HISAT2 index files |
| `alignment/` | Coordinate-sorted BAM files and BAM indexes |
| `counts/` | featureCounts matrix and summary file |
| `logs/` | Reserved for command logs |

## R-Generated Files Mentioned in `DEG.R`

The R script also writes the following files when executed:

| File | Description |
|---|---|
| `Significant_DEGs.csv` | Unannotated DEGs passing both thresholds |
| `Upregulated_genes.csv` | Unannotated genes with `padj < 0.05` and log2FC > 1 |
| `Downregulated_genes.csv` | Unannotated genes with `padj < 0.05` and log2FC < -1 |
| `upregulated_annotated.csv` | Annotated upregulated genes |
| `downregulated_annotated.csv` | Annotated downregulated genes |

These supplementary tables were not among the supplied project files and are
therefore not tracked in this repository snapshot.

