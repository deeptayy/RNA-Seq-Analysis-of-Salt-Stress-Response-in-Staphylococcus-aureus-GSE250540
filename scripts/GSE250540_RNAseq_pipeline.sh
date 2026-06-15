#!/bin/bash

set -euo pipefail

SECONDS=0

# change working directory
cd ~/Desktop/GSE250540_RNASeq/

mkdir -p raw trimmed fastqc_raw fastqc_trimmed alignment counts reference logs scripts

SAMPLES=(SRR27291863 SRR27291864 SRR27291865 SRR27291866 SRR27291867 SRR27291868)

#Download SRA data (single-end)

echo "Downloding SRA"

for sample in "${SAMPLES[@]}"
do
  prefetch $sample -O raw/
  fasterq-dump raw/$sample -O raw/ --threads 4
done
echo "SRA download finished"

# Run FastQC (raw)

echo "Running FastQC (raw)"

for sample in "${SAMPLES[@]}"
do
  fastqc raw/${sample}.fastq -o fastqc_raw/
done
echo "FastQC (raw) finished"

# MultiQC summary

echo "Running MultiQC(raw)"

multiqc fastqc_raw/ -o fastqc_raw/

echo "MultiQC(raw) finished"

# Trimmomatic

echo "Running Trimmomatic"

for sample in "${SAMPLES[@]}"
do
  java -jar ~/Desktop/demo/tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 \
    raw/${sample}.fastq trimmed/${sample}_trimmed.fastq TRAILING:10 -phred33
done
echo "Trimmomatic finished"

# FastQC (trimmed)

echo "Running FastQC (trimmed)"
for sample in "${SAMPLES[@]}"
do
  fastqc trimmed/${sample}_trimmed.fastq -o fastqc_trimmed/
done
echo "FastQC (trimmed) finished"

# MultiQC summary

echo "Running MultiQC (trimmed)"

multiqc fastqc_trimmed/ -o fastqc_trimmed/

echo "MultiQC(trimmed) finished"

# Download Staphylococcus aureus genome & GTF

echo "Downloding Genome & GTF"

# Genome FASTA (NCBI GCF_000013425.1)
wget -O reference/GCF_000013425.1.fna.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz"
gunzip reference/GCF_000013425.1.fna.gz

# GTF annotation
wget -O reference/GCF_000013425.1.gtf.gz \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gtf.gz"
gunzip reference/GCF_000013425.1.gtf.gz

echo " Genome and GTF download finished"

# Build HISAT2 index
hisat2-build reference/GCF_000013425.1.fna reference/GCF_000013425.1_index
echo " HISAT2 index build finished"

# HISAT2 alignment

for sample in "${SAMPLES[@]}"
do
  hisat2 -q --threads 4 -x reference/GCF_000013425.1_index \
    -U trimmed/${sample}_trimmed.fastq \
    | samtools sort -@ 4 -o alignment/${sample}.bam
  samtools index alignment/${sample}.bam
done
echo "alignment finished"

# featureCounts 
featureCounts -T 4 -a reference/GCF_000013425.1.gtf \
  -o counts/featurecounts.txt \
  alignment/*.bam
echo "featureCounts finished"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
