# Differential-Expression Results

## `DESeq2_results.csv`

Contains all tested features.

| Column | Description |
|---|---|
| First unnamed column | Feature or locus identifier |
| `baseMean` | Mean normalized count across all samples |
| `log2FoldChange` | Salt relative to control; positive is higher under salt |
| `lfcSE` | Standard error of the log2 fold change |
| `stat` | Wald test statistic |
| `pvalue` | Unadjusted Wald test p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |

## `SaltStress_Annotated_DEGs.csv`

Contains features satisfying:

```text
padj < 0.05 and abs(log2FoldChange) > 1
```

It adds:

| Column | Description |
|---|---|
| `locus_tag` | Reference annotation locus tag |
| `gene` | Gene symbol, when available |
| `product` | Annotated gene product, when available |

Counts in the supplied files:

- 2,376 tested features
- 2,372 features with non-missing adjusted p-values
- 357 significant features
- 127 upregulated features
- 230 downregulated features

