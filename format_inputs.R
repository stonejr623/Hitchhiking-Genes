#!/usr/bin/R
#Script to add genomic locations, strand, and biotype info from the gtf to the DESeq2 results data frame

#build library
library(readr)
library(tidyverse)
library(dplyr)
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(optparse)

#parse command line arguments
option_list <- list(
  make_option(c("--deseq"), type="character", help="Path to DESeq2 CSV output"),
  make_option(c("--gene_gtf"), type="character", help="Path to gene GTF file"),
  make_option(c("--te_gtf"), type="character", help="Path to TE GTF file"),
  make_option(c("--id_col"), type="character", default="gene", help="Column name for gene/TE ID"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
gene_gtf <- rtracklayer::import(opt$gene_gtf)
te_gtf <- rtracklayer::import(opt$te_gtf)
id_col   <- opt$id_col
deseq <- read.csv(opt$deseq, row.names = 1)

#clean deseq file
deseq <- deseq[!duplicated(rownames(deseq)), ]
rownames(deseq) <- sub(":.*", "", rownames(deseq))
deseq$`gene/TE` <- rownames(deseq)

#load gtf files separately
gene_gtf <- as.data.frame(gene_gtf)
genes_to_keep <- gene_gtf %>%
  filter(gene_name %in% rownames(deseq)) %>%
  select(seqnames, start, end, strand, gene_name, gene_biotype)
#keep longest annotation (full seauence not just one exon)
genes_to_keep <- genes_to_keep %>%
  group_by(gene_name) %>%
  mutate(length = end - start) %>%
  slice_max(order_by = length, n = 1, with_ties = FALSE) %>%
  ungroup()
colnames(genes_to_keep)[5] <- "gene/TE"

#repeat for TE gtf
te_gtf <- as.data.frame(te_gtf)
colnames(te_gtf)[13]="gene_biotype"
TEs_to_keep <- te_gtf %>%
  filter(transcript_id %in% rownames(deseq)) %>%
  select (seqnames, start, end, strand, transcript_id, gene_biotype)
TEs_to_keep <- TEs_to_keep %>%
  group_by(transcript_id) %>%
  mutate(length = end - start) %>%
  slice_max(order_by = length, n = 1, with_ties = FALSE) %>%
  ungroup()
colnames(TEs_to_keep)[5] <- "gene/TE"

#merge desired columns from gtf to deseq
genes_to_keep <- column_to_rownames(genes_to_keep, var="gene/TE")
TEs_to_keep <- column_to_rownames(TEs_to_keep, var="gene/TE")
big_df <- rbind(TEs_to_keep, genes_to_keep)
big_df$`gene/TE` <- rownames(big_df)
big_df <- merge(deseq, big_df, by = "gene/TE", all.x = TRUE)
big_df <- big_df %>%
  select(`gene/TE`, log2FoldChange, pvalue, chr = seqnames, start, end, strand, gene_biotype)

#separate data frame into coding and noncoding
coding <- subset(big_df, big_df$gene_biotype == "protein_coding")
noncoding <- subset(big_df, big_df$gene_biotype != "protein_coding")

#save as csv files to import into the next step
write.csv(coding, file= "coding_genes.csv", row.names = FALSE)
write.csv(noncoding, file= "noncoding_genes.csv", row.names = FALSE)
