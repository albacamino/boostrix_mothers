rm(list=ls())
library("apeglm")
library("tidyverse")
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("limma")
library("ggrepel")

setwd("/Users/alba/repos/boostrix_mothers")

plot_volcano <- function(res, titulo) {
  as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(
      significance = case_when(
        padj < 0.05 & log2FoldChange >  1.5 ~ "Upregulated",
        padj < 0.05 & log2FoldChange < -1.5 ~ "Downregulated",
        TRUE                               ~ "NS"
      ),
      label = ifelse(padj < 0.05 & abs(log2FoldChange) > 1.5, gene, NA)
    ) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj),
               color = significance, label = label)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = c("Upregulated"   = "#534AB7",
                                  "Downregulated" = "#0F6E56",
                                  "NS"            = "grey70")) +
    labs(x = "log2 Fold Change",
         y = "-log10(adjusted p-value)", color = NULL) +
    theme_bw() +
    theme_gray(base_size = 15) +
    theme(legend.position = "none") +
    guides(fill = FALSE)
}

plot_pca <- function(res, var1 = "Condition", var2 = "Timepoint", 
                     opciones1 = NULL, opciones2 = NULL) {
  
  # Si no se especifican opciones, usar todos los niveles disponibles
  if (is.null(opciones1)) opciones1 <- unique(metadata_filter[[var1]])
  if (is.null(opciones2)) opciones2 <- unique(metadata_filter[[var2]])
  
  idx <- metadata_filter[[var1]] %in% opciones1 & 
    metadata_filter[[var2]] %in% opciones2
  
  # Verificar que hay suficientes muestras
  if (sum(idx) < 3) stop("Menos de 3 muestras seleccionadas, revisa los filtros")
  
  pca_res <- prcomp(t(res[, idx]), scale. = FALSE)
  pca_df  <- as.data.frame(pca_res$x[, 1:2])
  
  pca_df$Condition  <- metadata_filter$Condition[idx]
  pca_df$Timepoint  <- paste0("T", metadata_filter$Timepoint[idx])
  pca_df$Individual <- metadata_filter$Ext_ID[idx]
  
  var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  
  ggplot(pca_df, aes(x = PC1, y = PC2,
                     color = Timepoint,
                     shape = Condition,
                     label = Individual)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("T1" = "#E8593C",
                                  "T2" = "#534AB7", 
                                  "T3" = "#0F6E56")) +
    labs(x     = paste0("PC1 (", var_exp[1], "%)"),
         y     = paste0("PC2 (", var_exp[2], "%)")) +
    theme_bw() +
    theme_gray(base_size = 15)
}

plot_heatmap <- function(res, metadata, opciones1 = NULL, opciones2 = NULL) {
  
  res_df <- as.data.frame(res)
  
  sig <- rownames(res_df)[which(!is.na(res_df$padj) & 
                                  res_df$padj < 0.05 & 
                                  abs(res_df$log2FoldChange) > 1.5)]
  
  if (length(sig) == 0) stop("No hay genes significativos con estos criterios")
  cat("Genes significativos:", length(sig), "\n")
  
  # Usar vst_corrected en lugar de recalcular vst
  mat <- vst_corrected[sig, ]
  
  # Filtrar muestras
  if (is.null(opciones1) | is.null(opciones2)) {
    samples_to_plot <- rownames(metadata)
  } else {
    samples_to_plot <- rownames(metadata)[metadata$Timepoint %in% opciones1 & metadata$Condition %in% opciones2]
  }
  
  mat_filtered <- mat[, samples_to_plot]
  
  annotation_col <- data.frame(
    Condition = metadata[samples_to_plot, "Condition"],
    Timepoint = metadata[samples_to_plot, "Timepoint"],
    row.names = samples_to_plot
  )
  
  pheatmap(mat_filtered,
           annotation_col = annotation_col,
           scale          = "row",
           show_colnames  = FALSE,
           fontsize_row   = 6)
}

# Load data (counts and metadata)
counts <- read.csv("data/counts.csv", header = TRUE, row.names = 1)
metadata <- read.csv("data/metadata.csv", header = TRUE)

# Factorize the variables
metadata$Condition <- factor(metadata$Condition)
metadata$Timepoint <- factor(metadata$Timepoint)

rownames(metadata) <- metadata$NIM_ID
names(counts) <- sub("^X", "", names(counts))

# Count the number of timepoints in each individual per condition
table_tp <- table(metadata$Ext_ID)

# Keep individuals with 2 or more samples
valid_ids <- names(table_tp[table_tp == 3])

# Filter metadata and counts
metadata_filter <- metadata[metadata$Ext_ID %in% valid_ids, ]

# Remove sample
metadata_filter <- metadata_filter[which(metadata_filter$Ext_ID!="A4316"),]

# Ensure counts have the same samples as metadata
counts_filter <- counts[, rownames(metadata_filter)]

nrow(metadata_filter)
ncol(counts_filter)

all(table(metadata_filter$Ext_ID) == 3)

# Ensure row names of metadata match column names of counts
all(rownames(metadata_filter) %in% colnames(counts_filter))

# Check the order of the columns in the counts and the metadata
counts_filter <- counts_filter[, rownames(metadata_filter)]
all(rownames(metadata_filter) == colnames(counts_filter))

# Create the DESeq2 object
design <- ~ Condition + Timepoint + Condition:Timepoint
dds <- DESeqDataSetFromMatrix(countData = counts_filter,
                              colData = metadata_filter,
                              design = design)

# Pre-filtering: Remove low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Differential Expression Analysis
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

# Exportar matriz VST
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation
vst_mat <- assay(vsd)

# Remove Batch Effect 
vst_corrected <- removeBatchEffect(
  vst_mat,
  batch  = metadata_filter$Ext_ID,
  design = model.matrix(~ Timepoint, data = metadata_filter)
)

# PCA plot
plot_pca(vst_corrected, opciones1 = c("Placebo","Vaccinated"), opciones2 = c("1","2","3"))

# Ver qué niveles tiene cada factor (el primero es la referencia)
levels(metadata_filter$Condition)
levels(metadata_filter$Timepoint)

resultsNames(dds)

# 1. T1 vs T3 en Placebo
#    (efecto directo de Timepoint en la referencia = Placebo)
res_placebo_T1vsT3 <- results(dds, 
                              name = "Timepoint_3_vs_1")

# 2. T1 vs T3 en Vacunados
#    (efecto de Timepoint en Placebo + término de interacción)
res_vac_T1vsT3 <- results(dds,
                          contrast = list(c("Timepoint_3_vs_1",
                                            "ConditionVaccinated.Timepoint3")))

# 3. Vacunados vs Placebo en T1
#    (efecto principal de Condition, en T1 = referencia)
res_T1_vacVsPlacebo <- results(dds,
                               name = "Condition_Vaccinated_vs_Placebo")

# 4. Vacunados vs Placebo en T3
#    (efecto de Condition en T1 + término de interacción en T3)
res_T3_vacVsPlacebo <- results(dds,
                               contrast = list(c("Condition_Vaccinated_vs_Placebo",
                                                 "ConditionVaccinated.Timepoint3")))

# Volcano Plot para cada contraste
plot_volcano(res_T3_vacVsPlacebo,  "Vaccinated vs Placebo - T3")
plot_volcano(res_T1_vacVsPlacebo,  "Vaccinated vs Placebo - T1")
plot_volcano(res_vac_T1vsT3,       "Vaccinated: T1 vs T3")
plot_volcano(res_placebo_T1vsT3,   "Placebo: T1 vs T3")

# p-values and p-adjusted
resOrdered <- res_T3_vacVsPlacebo[order(res_T3_vacVsPlacebo$pvalue),]
summary(res_T3_vacVsPlacebo)

sum(res_T3_vacVsPlacebo$padj < 0.05, na.rm=TRUE)
sum(res_T3_vacVsPlacebo$padj < 0.05 & 
      abs(res_T3_vacVsPlacebo$log2FoldChange) > 2, na.rm = TRUE)

plot_heatmap(res_T3_vacVsPlacebo, metadata_filter, opciones1 = c("3"), opciones2 = c("Vaccinated", "Placebo"))

  
