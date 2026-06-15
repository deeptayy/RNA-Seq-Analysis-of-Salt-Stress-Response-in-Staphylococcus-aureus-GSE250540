# GSE250540: S. aureus salt stress
# setwd("~/Staph_rnase")

library(DESeq2)
library(tidyverse)
library(rtracklayer)

counts <- read.delim('gene_counts.txt', comment.char='#')

count_matrix <- counts[,7:ncol(counts)]
rownames(count_matrix) <- counts$Geneid
head(count_matrix)

metadata<-data.frame(
  row.names = colnames(count_matrix) <- c(
    "Salt_1",
    "Salt_2",
    "Salt_3",
    "Control_1",
    "Control_2",
    "Control_3"
  ),
  
  condition = c(
    "Salt",
    "Salt",
    "Salt",
    "Control",
    "Control",
    "Control"
  )
)

dds <-DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~condition
)

dds <- dds[rowSums(counts(dds)) > 10,]

dds <- DESeq(dds)

res <- results(dds)

sig <- subset(
  res,
  padj <0.05 &
  abs(log2FoldChange) > 1
)

write.csv(as.data.frame(res),'DESeq2_results.csv')

write.csv(as.data.frame(sig),'Significant_DEGs.csv')

vsd <- vst(dds)

plotPCA(vsd,intgroup = 'condition')

topgenes <-head (rownames(sig[order(sig$padj),]),15)

library(EnhancedVolcano)

EnhancedVolcano(
  res,
  lab = ifelse(rownames(res) %in% topgenes, rownames(res),''),
  x='log2FoldChange',
  y='padj'
)

library(pheatmap)

pheatmap(assay(vsd)[topgenes,], fontsize_row =8 ,fontsize_col = 10)

up<- subset(res,padj<0.05 & log2FoldChange>1)
down<- subset(res, padj<0.05 & log2FoldChange< -1)

write.csv(as.data.frame(up),'Upregulated_genes.csv')

write.csv(as.data.frame(down),'Downregulated_genes.csv')

plotMA(res)

head(res[order(res$padj),c('log2FoldChange','padj')],50)
write.csv(as.data.frame(res), 'DESeq2_results.csv')

# Import GFF file
gff <- import("GCF_000013425.1_ASM1342v1_genomic.gff")

# Convert to data frame
cds_features <- as.data.frame(gff)

# Keep only CDS features
cds_features <- cds_features %>%
  filter(type == "CDS")

# Check available columns
colnames(cds_features)

# Create annotation table
annotation <- cds_features %>%
  select(locus_tag,gene,product) %>%
  distinct()

# Check annotation quality
sum(is.na(annotation$gene))
sum(is.na(annotation$product))
nrow(annotation)

# Convert DESeq2 results to data frame
res_df <- as.data.frame(res)

# Move rownames into a column
res_df$locus_tag <- rownames(res_df)

# Merge DESeq2 results with annotation
annotated_res <- merge(res_df,annotation,by = "locus_tag",all.x = TRUE)

# Check merged result
head(annotated_res)

sum(is.na(annotation$product))

top20<- annotated_res[order(annotated_res$padj),]

top20<- top20[,c('locus_tag','gene','product','log2FoldChange','padj')]
head(top20,20)

annotated_res[annotated_res$locus_tag %in% c('SAOUHSC_00094','SAOUHSC_00187','SAOUHSC_00188','SAOUHSC_00356','SAOUHSC_00544',
                                             'SAOUHSC_00561','SAOUHSC_00634','SAOUHSC_00636','SAOUHSC_00660','SAOUHSC_00898','SAOUHSC_00899',
                                             'SAOUHSC_00923','SAOUHSC_00924','SAOUHSC_00925','SAOUHSC_00927','SAOUHSC_01991','SAOUHSC_02119', 
                                             'SAOUHSC_02333','SAOUHSC_02990','SAOUHSC_01053' ),c('locus_tag','gene','product','log2FoldChange','padj')]

sig <- subset(annotated_res, padj<0.05 & abs(log2FoldChange)>1)

sig <- sig[order(sig$padj),]

write.csv(sig, 'SaltStress_Annotated_DEGs.csv', row.names =FALSE)

up_annotated <- subset(annotated_res,padj<0.05 & log2FoldChange >1)

down_annotated <- subset(annotated_res, padj <0.05 & log2FoldChange < -1)

write.csv(up_annotated,'upregulated_annotated.csv', row.names = FALSE)

write.csv(down_annotated,'downregulated_annotated.csv', row.names = FALSE)

head(up_annotated[order(-up_annotated$log2FoldChange),],20)

head(down_annotated[order(down_annotated$log2FoldChange),],20)

clusterProfiler
