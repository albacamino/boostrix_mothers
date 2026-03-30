rm(list=ls())
library("apeglm")
library("tidyverse")
library("DESeq2")
library("RColorBrewer")
library("pheatmap")

setwd("/Users/alba/repos/boostrix_mothers")

# Load data (counts and metadata)
counts <- read.csv("data/counts.csv", header = TRUE, row.names = 1)
metadata <- read.csv("data/metadata.csv", header = TRUE)

# Factorize the variables
metadata$Ext_ID <- factor(metadata$Ext_ID)
metadata$Condition <- factor(metadata$Condition)
metadata$Timepoint <- factor(metadata$Timepoint)

rownames(metadata) <- metadata$NIM_ID
names(counts) <- sub("^X", "", names(counts))

# Count the number of timepoints in each individual per condition
condition <- "Vaccinated"
metadata_vac <- metadata[metadata$Condition == condition, ]
table_tp <- table(metadata_vac$Ext_ID)

# Keep individuals with 2 or more samples
valid_ids <- names(table_tp[table_tp > 1])

# Filter metadata and counts
metadata_filter<- metadata_vac[metadata_vac$Ext_ID %in% valid_ids, ]
counts_filter <- counts[, rownames(metadata_vac)]

# Ensure row names of metadata match column names of counts
all(rownames(metadata_filter) %in% colnames(counts_filter))

# Check the order of the columns in the counts and the metadata
counts_filter <- counts_filter[, rownames(metadata_filter)]
all(rownames(metadata_filter) == colnames(counts_filter))

# Create the DESeq2 object
design <- ~Ext_ID + Timepoint
dds <- DESeqDataSetFromMatrix(countData = counts_filter,
                              colData = metadata_filter,
                              design = design)

# Pre-filtering: Remove low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# PCA plot
vsd <- vst(dds, blind = TRUE)  # Variance-stabilizing transformation
plotPCA(vsd, intgroup = "Timepoint")

# Exportar matriz VST
vsd_mat <- assay(vsd)
write.csv(vsd_mat, paste0("data/vsd_", condition,"_matrix.csv"))

# Differential Expression Analysis
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Timepoint_3_vs_1")

# Shrinkage of effect size (LFC estimates)
resLFC <- lfcShrink(dds, coef="Timepoint_3_vs_1", type="apeglm")
resLFC

# p-values and p-adjusted
resOrdered <- res[order(res$pvalue),]
summary(res)

sum(res$padj < 0.05, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

sum(res05$padj < 0.05 & abs(res05$log2FoldChange)>1.5, na.rm=TRUE)

# MA plot
plotMA(resLFC, ylim=c(-2,2))

# Heatmap
sig <- rownames(res)[which(res$padj < 0.01)]
annotation_col <- data.frame(Timepoint = metadata_filter$Timepoint)
rownames(annotation_col) <- rownames(metadata_filter)  # fila = muestra
mat <- assay(vst(dds))[sig, ]
all(rownames(annotation_col) == colnames(mat)) 

pheatmap(mat,
         annotation_col = annotation_col,
         scale = "row")  # opcional para escalar genes
